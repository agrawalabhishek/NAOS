/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <vector>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/misc.hpp"
#include "NAOS/gslIntegratorOrbiterAroundUREPointMassGravity.hpp"

namespace naos
{

//! equations of motion
/*!
 * first order differential equations describing the motion of a particle or spacecraft around a
 * central body modeled as a point mass gravity. Look at GSL odeiv2 examples to see the format
 * of how the equations are written for a c++ application.
 */
int equationsOfMotion( double currentTime,
                       const double stateVector[ ],
                       double stateTimeDerivative[ ],
                       void *parameters )
{
    ( void )( currentTime ); // to avoid warning for an unused parameter value

    // type-casting from void to asteroidParameters. This may give some warnings on compilation.
    asteroidParameters param = *( asteroidParameters *) ( parameters );
    const double Ux = param.xGravAcceleration;
    const double Uy = param.yGravAcceleration;
    const double Uz = param.zGravAcceleration;
    const double Wz = param.zRotation;

    // write down the first order equations
    stateTimeDerivative[ xPositionIndex ] = stateVector[ xVelocityIndex ];
    stateTimeDerivative[ yPositionIndex ] = stateVector[ yVelocityIndex ];
    stateTimeDerivative[ zPositionIndex ] = stateVector[ zVelocityIndex ];

    stateTimeDerivative[ xVelocityIndex ] = Wz * Wz * stateVector[ xPositionIndex ]
                                            + Ux
                                            + 2.0 * Wz * stateVector[ yVelocityIndex ];

    stateTimeDerivative[ yVelocityIndex ] = Wz * Wz * stateVector[ yPositionIndex ]
                                            + Uy
                                            - 2.0 * Wz * stateVector[ xVelocityIndex ];

    stateTimeDerivative[ zVelocityIndex ] = Uz;

    // return success flag
    return GSL_SUCCESS;
}

//! Get the jacobian for the integrator
int jacobianForIntegrator( double currentTime,
                           const double stateVector[ ],
                           double dfdy[ ],
                           double dfdt[ ],
                           void *parameters )
{
    ( void )( currentTime ); // to avoid warning for unused parameter

    asteroidParameters param = *( asteroidParameters *) ( parameters );
    const double Ux = param.xGravAcceleration;
    const double Uy = param.yGravAcceleration;
    const double Uz = param.zGravAcceleration;
    const double Wz = param.zRotation;

    // form the dfdy matrix
    gsl_matrix_view dfdyMatrix = gsl_matrix_view_array( dfdy, 6, 6 );
    gsl_matrix * m = &dfdyMatrix.matrix;

    gsl_matrix_set_all( m, 0.0 );

    gsl_matrix_set( m, 0, 3, 1.0 );
    gsl_matrix_set( m, 1, 4, 1.0 );
    gsl_matrix_set( m, 2, 5, 1.0 );
    gsl_matrix_set( m, 3, 0, Wz * Wz );
    gsl_matrix_set( m, 3, 4, 2.0 * Wz );
    gsl_matrix_set( m, 4, 1, Wz * Wz );
    gsl_matrix_set( m, 4, 3, -2.0 * Wz );

    // form the dfdt array
    dfdt[ 0 ] = 0.0;
    dfdt[ 1 ] = 0.0;
    dfdt[ 2 ] = 0.0;
    dfdt[ 3 ] = 0.0;
    dfdt[ 4 ] = 0.0;
    dfdt[ 5 ] = 0.0;

    // return success
    return GSL_SUCCESS;
}

//! Calculate the jacobian
/*!
 * calculate the jacobian for the point mass gravity potential case.
 */
double calculateJacobianPointMassGravity( const double stateVector[ ],
                                          const double angularVelocity,
                                          const double gravParameter )
{
    double xVelocitySquare = stateVector[ xVelocityIndex ] * stateVector[ xVelocityIndex ];

    double yVelocitySquare = stateVector[ yVelocityIndex ] * stateVector[ yVelocityIndex ];

    double zVelocitySquare = stateVector[ zVelocityIndex ] * stateVector[ zVelocityIndex ];

    double xPositionSquare = stateVector[ xPositionIndex ] * stateVector[ xPositionIndex ];

    double yPositionSquare = stateVector[ yPositionIndex ] * stateVector[ yPositionIndex ];

    double zPositionSquare = stateVector[ zPositionIndex ] * stateVector[ zPositionIndex ];

    double radialDistance = std::sqrt( xPositionSquare + yPositionSquare + zPositionSquare );

    double angularVelocitySquare = angularVelocity * angularVelocity;

    double gravPotential = gravParameter / radialDistance;

    double jacobian = 0.5 * ( xVelocitySquare + yVelocitySquare + zVelocitySquare )
                        - 0.5 * angularVelocitySquare * ( xPositionSquare + yPositionSquare )
                        - gravPotential;

    return jacobian;
}

//! Execute orbiter around an asteroid (URE) modeled as a point mass
/*!
 * integrate the equations of motion for a particle around an asteroid using the GSL odeiv2 library
 * for a given set of initial conditions given in the form of classical orbital elements. The
 * cartesian state vector is store in a csv file.
 */
void gslIntegratorOrbiterAroundUREPointMassGravity( const double gravParameter,
                                                    std::vector< double > &asteroidRotationVector,
                                                    std::vector< double > &initialOrbitalElements,
                                                    const double initialStepSize,
                                                    const double startTime,
                                                    const double endTime,
                                                    std::ostringstream &outputFilePath )
{
    //! open the output csv file to save data. Declare file headers.
    std::ofstream outputFile;
    outputFile.open( outputFilePath.str( ) );
    outputFile << "x" << ",";
    outputFile << "y" << ",";
    outputFile << "z" << ",";
    outputFile << "vx" << ",";
    outputFile << "vy" << ",";
    outputFile << "vz" << ",";
    outputFile << "t" << ",";
    outputFile << "jacobian" << ",";
    outputFile << "stepSize" << std::endl;

    //! convert the initial orbital elements to cartesian state (expressed in the body frame)
    std::vector< double > initialStateInertialFrame( 6, 0.0 );
    initialStateInertialFrame = convertKeplerianElementsToCartesianCoordinates(
                                                            initialOrbitalElements,
                                                            gravParameter );

    //! get the intial state in body frame coordinates
    std::vector< double > initialStateBodyFrame( 6, 0.0 );

    // account for a non zero start time value
    double rotationAngle = asteroidRotationVector[ zPositionIndex ] * startTime;

    initialStateBodyFrame[ xPositionIndex ]
                    = initialStateInertialFrame[ xPositionIndex ] * std::cos( rotationAngle )
                    + initialStateInertialFrame[ yPositionIndex ] * std::sin( rotationAngle );

    initialStateBodyFrame[ yPositionIndex ]
                    = -1.0 * initialStateInertialFrame[ xPositionIndex ] * std::sin( rotationAngle )
                    + initialStateInertialFrame[ yPositionIndex ] * std::cos( rotationAngle );

    initialStateBodyFrame[ zPositionIndex ] = initialStateInertialFrame[ zPositionIndex ];

    // get the omega cross position vector term of the transport theorem
    std::vector< double > omegaCrossPosition( 3, 0.0 );
    std::vector< double > bodyFramePositionVector = { initialStateBodyFrame[ xPositionIndex ],
                                                      initialStateBodyFrame[ yPositionIndex ],
                                                      initialStateBodyFrame[ zPositionIndex ] };

    omegaCrossPosition = crossProduct( asteroidRotationVector, bodyFramePositionVector );

    // get the time derivative of the body frame position vector in inertial frame by using the
    // principal Z rotation matrix
    double xBodyFrameVelocityInertialCoordinates
                    = initialStateInertialFrame[ xVelocityIndex ] * std::cos( rotationAngle )
                    + initialStateInertialFrame[ yVelocityIndex ] * std::sin( rotationAngle );

    double yBodyFrameVelocityInertialCoordinates
                    = -1.0 * initialStateInertialFrame[ xVelocityIndex ] * std::sin( rotationAngle )
                    + initialStateInertialFrame[ yVelocityIndex ] * std::cos( rotationAngle );

    double zBodyFrameVelocityInertialCoordinates = initialStateInertialFrame[ zVelocityIndex ];

    // get the velocity in body frame using the transport theorem and previous result
    initialStateBodyFrame[ xVelocityIndex ]
        = xBodyFrameVelocityInertialCoordinates - omegaCrossPosition[ xPositionIndex ];

    initialStateBodyFrame[ yVelocityIndex ]
        = yBodyFrameVelocityInertialCoordinates - omegaCrossPosition[ yPositionIndex ];

    initialStateBodyFrame[ zVelocityIndex ]
        = zBodyFrameVelocityInertialCoordinates - omegaCrossPosition[ zPositionIndex ];

    //! initialize the step size and current time variables
    double stepSize = initialStepSize;
    double currentTime = startTime;

    //! set the current state vector
    double currentStateVector[ 6 ] = { initialStateBodyFrame[ xPositionIndex ],
                                       initialStateBodyFrame[ yPositionIndex ],
                                       initialStateBodyFrame[ zPositionIndex ],
                                       initialStateBodyFrame[ xVelocityIndex ],
                                       initialStateBodyFrame[ yVelocityIndex ],
                                       initialStateBodyFrame[ zVelocityIndex ] };

    //! initial jacobian
    double jacobian = calculateJacobianPointMassGravity( currentStateVector,
                                                         asteroidRotationVector[ zPositionIndex ],
                                                         gravParameter );

    //! save the initial body frame cartesian state vector in the csv file
    outputFile << initialStateBodyFrame[ xPositionIndex ] << ",";
    outputFile << initialStateBodyFrame[ yPositionIndex ] << ",";
    outputFile << initialStateBodyFrame[ zPositionIndex ] << ",";
    outputFile << initialStateBodyFrame[ xVelocityIndex ] << ",";
    outputFile << initialStateBodyFrame[ yVelocityIndex ] << ",";
    outputFile << initialStateBodyFrame[ zVelocityIndex ] << ",";
    outputFile << currentTime << ",";
    outputFile << jacobian << ",";
    outputFile << initialStepSize << std::endl;

    //! GSL integrator
    // set the numerical integrator algorithm
    const gsl_odeiv2_step_type * stepType = gsl_odeiv2_step_msadams;
    gsl_odeiv2_step * stepper = gsl_odeiv2_step_alloc( stepType, 6 );

    // set the tolerances
    const double absoluteTolerance = 10.0e-12;
    const double relativeTolerance = 10.0e-12;

    // set the step size control method
    gsl_odeiv2_control * stepSizeControl = gsl_odeiv2_control_standard_new( absoluteTolerance,
                                                                            relativeTolerance,
                                                                            1.0,
                                                                            1.0 );

    // set the step evolution function
    gsl_odeiv2_evolve * stepEvolution = gsl_odeiv2_evolve_alloc( 6 );

    // assign the gravitational acceleration values to an object of asteroid paramter struct
    asteroidParameters * parameters = ( asteroidParameters *)std::malloc( sizeof( asteroidParameters ) );
    parameters->xGravAcceleration = 0.0;
    parameters->yGravAcceleration = 0.0;
    parameters->zGravAcceleration = 0.0;
    parameters->zRotation = asteroidRotationVector[ zPositionIndex ];

    // set up the integration system
    gsl_odeiv2_system stepSystem = { equationsOfMotion,
                                     jacobianForIntegrator,
                                     6,
                                     parameters };

    // set up the driver wrapper for the integrator
    gsl_odeiv2_driver * driver =  gsl_odeiv2_driver_alloc_standard_new( &stepSystem,
                                                                        stepType,
                                                                        stepSize,
                                                                        absoluteTolerance,
                                                                        relativeTolerance,
                                                                        1.0,
                                                                        1.0 );

    int stepSetDriverStatus = gsl_odeiv2_step_set_driver( stepper, driver );

    // start the integration outer loop
    while( currentTime != endTime )
    {
        // calculate the gravitational acceleration values at current state values
        std::vector< double > positionVector = { currentStateVector[ xPositionIndex ],
                                                 currentStateVector[ yPositionIndex ],
                                                 currentStateVector[ zPositionIndex ] };

        double radialDistance = vectorNorm( positionVector );
        double radialDistanceCube = radialDistance * radialDistance * radialDistance;

        double xAcceleration = -1.0 * gravParameter * positionVector[ xPositionIndex ] / radialDistanceCube;

        double yAcceleration = -1.0 * gravParameter * positionVector[ yPositionIndex ] / radialDistanceCube;

        double zAcceleration = -1.0 * gravParameter * positionVector[ zPositionIndex ] / radialDistanceCube;

        // update the asteroid parameters object
        parameters->xGravAcceleration = xAcceleration;
        parameters->yGravAcceleration = yAcceleration;
        parameters->zGravAcceleration = zAcceleration;
        parameters->zRotation = asteroidRotationVector[ zPositionIndex ];

        // reset the step evolve to account for the new parametric values in the equations of motion
        int resetStatus = gsl_odeiv2_step_reset( stepper );

        // perform one step integration (currentTime and currentStateVector are auto-updated)
        int stepStatus = gsl_odeiv2_evolve_apply( stepEvolution,
                                                  stepSizeControl,
                                                  stepper,
                                                  &stepSystem,
                                                  &currentTime,
                                                  endTime,
                                                  &stepSize,
                                                  currentStateVector );

        // throw exception if integration failed
        if( stepStatus != GSL_SUCCESS )
        {
            std::ostringstream errorMessage;
            errorMessage << std::endl;
            errorMessage << "ERROR: gsl integration routine failed" << std::endl;
            errorMessage << std::endl;
            throw std::runtime_error( errorMessage.str( ) );
        } // end of if loop for runtime exception

        // calculate the current jacobian
        jacobian = calculateJacobianPointMassGravity( currentStateVector,
                                                      asteroidRotationVector[ zPositionIndex ],
                                                      gravParameter );
        // store the step result
        outputFile << currentStateVector[ xPositionIndex ] << ",";
        outputFile << currentStateVector[ yPositionIndex ] << ",";
        outputFile << currentStateVector[ zPositionIndex ] << ",";
        outputFile << currentStateVector[ xVelocityIndex ] << ",";
        outputFile << currentStateVector[ yVelocityIndex ] << ",";
        outputFile << currentStateVector[ zVelocityIndex ] << ",";
        outputFile << currentTime << ",";
        outputFile << jacobian << ",";
        outputFile << stepSize << std::endl;

    } // end of outer while loop for integration

    //! end of everything
    gsl_odeiv2_evolve_free( stepEvolution );
    gsl_odeiv2_control_free( stepSizeControl );
    gsl_odeiv2_step_free( stepper );
    gsl_odeiv2_driver_free( driver );
    outputFile.close( );
}

} // namespace naos
