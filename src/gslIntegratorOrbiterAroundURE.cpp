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

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

#include "NAOS/constants.hpp"
#include "NAOS/ellipsoidGravitationalAcceleration.hpp"
#include "NAOS/orbiterEquationsOfMotion.hpp"
#include "NAOS/rk4.hpp"
#include "NAOS/rk54.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/ellipsoidPotential.hpp"
#include "NAOS/misc.hpp"
#include "NAOS/gslIntegratorOrbiterAroundURE.hpp"

namespace naos
{

//! Equations of motion for spacecraft or particle around a uniformly rotating ellipsoid (URE)
/*!
 * function defining the equation of motion which will be used in the integrator routine of GSL
 *
 */

 int equationsOfMotion( double t, const double X[], double dXdt[], void *parameters )
 {
    ( void )( t );

    asteroidParameters param = *(asteroidParameters *) ( parameters );
    const double Ux = param.xGravAcceleration;
    const double Uy = param.yGravAcceleration;
    const double Uz = param.zGravAcceleration;
    const double Wz = param.zRotation;

    dXdt[ naos::xPositionIndex ] = X[ naos::xVelocityIndex ];
    dXdt[ naos::yPositionIndex ] = X[ naos::yVelocityIndex ];
    dXdt[ naos::zPositionIndex ] = X[ naos::zVelocityIndex ];
    dXdt[ naos::xVelocityIndex ] = Ux + 2.0 * Wz * X[ naos::yVelocityIndex ]
                                    + Wz * Wz * X[ naos::xPositionIndex ];
    dXdt[ naos::yVelocityIndex ] = Uy - 2.0 * Wz * X[ naos::xVelocityIndex ]
                                    + Wz * Wz * X[ naos::yPositionIndex ];
    dXdt[ naos::zVelocityIndex ] = Uz;
    return GSL_SUCCESS;
 }

//! Execute routine to compute orbiter/particle motion around a uniformly rotating ellipsoid (URE)
/*!
 * This routine will compute the motion of the orbiter or particle around an asteroid modelled as
 * a uniformly rotating triaxial ellipsoid. The routine is called in main for execution.
 * @sa ellipsoidGravitationalAcceleration.cpp, orbiterEquationsOfMotion.hpp, rk4.hpp
 *
 * @param[in] alpha                         largest semi axis of the ellipsoid
 * @param[in] beta                          intermediate semi axis of the ellipsoid
 * @param[in] gamma                         smallest largest semi axis of the ellipsoid
 * @param[in] gravParameter                 gravitational parameter of the ellipsoid
 * @param[in] density                       density of the ellipsoid
 * @param[in] Wvector                       angular velocity vector for the ellipsoid
 * @param[in] Wmagnitude                    magnitude of the angular velocity
 * @param[in] initialVectorIsCartesian      boolean flag to indicate type of initial vector
 *                                          value is false for kepler elements, and true for
 *                                          cartesian elements.
 * @param[in] initialVector                 vector containing the initial state of the orbiter
 * @param[in] integrationStepSize           step size to be used for integration
 * @param[in] startTime                     integration start time
 * @param[in] endTime                       integration end time
 * @param[in] filePath                      path to csv file where all the results are stored
 */
 void gslIntegratorOrbiterAroundURE( const double alpha,
                                     const double beta,
                                     const double gamma,
                                     const double gravParameter,
                                     const double density,
                                     std::vector< double > Wvector,
                                     const double Wmagnitude,
                                     const bool initialVectorIsCartesian,
                                     const std::vector< double > &initialVector,
                                     const double integrationStepSize,
                                     const double startTime,
                                     const double endTime,
                                     std::ostringstream &filePath )
 {
    // Use a numerical integrator to evaluate the equations of motion for an orbiter around a uniformly
    // rotating triaxial ellipsoid (URE). Store the solutions for the entire time of integration in
    // a CSV file.
    std::ofstream eomOrbiterUREFile;
    eomOrbiterUREFile.open( filePath.str( ) );
    eomOrbiterUREFile << "x" << "," << "y" << "," << "z" << ",";
    eomOrbiterUREFile << "vx" << "," << "vy" << "," << "vz" << "," << "t" << ",";
    eomOrbiterUREFile << "jacobian" << "," << "semiMajor" << "," << "eccentricity" << ",";
    eomOrbiterUREFile << "inclination" << "," << "RAAN" << "," << "AOP" << "," << "TA" << ",";
    eomOrbiterUREFile << "stepSize" << std::endl;

    // get the initial state from the input arguments
    Vector6 initialStateVector( 6 );
    Vector6 orbitalElements( 6 );
    if( initialVectorIsCartesian == false )
    {
        for( int i = 0; i < 6; i++ )
        {
            orbitalElements[ i ] = initialVector[ i ];
        }

        // Specify initial values for the state vector of the orbiter in the inertial frame.
        Vector6 inertialStateVector( 6 );
        inertialStateVector = convertKeplerianElementsToCartesianCoordinates< Vector6 >(
                                orbitalElements,
                                gravParameter );

        // convert position in inertial frame to position in body frame at the starting time
        double phiAngle = Wvector[ zPositionIndex ] * startTime; // angle is in radians

        double xPositionBodyFrame = inertialStateVector[ xPositionIndex ] * std::cos( phiAngle )
                                    + inertialStateVector[ yPositionIndex ] * std::sin( phiAngle );

        double yPositionBodyFrame = -inertialStateVector[ xPositionIndex ] * std::sin( phiAngle )
                                    + inertialStateVector[ yPositionIndex ] * std::cos( phiAngle );

        double zPositionBodyFrame = inertialStateVector[ zPositionIndex ];

        initialStateVector[ xPositionIndex ] = xPositionBodyFrame;
        initialStateVector[ yPositionIndex ] = yPositionBodyFrame;
        initialStateVector[ zPositionIndex ] = zPositionBodyFrame;

        // get velocity in the body frame from inertial frame for the starting time/initial epoch
        Vector3 positionVector { initialStateVector[ xPositionIndex ],
                                 initialStateVector[ yPositionIndex ],
                                 initialStateVector[ zPositionIndex ] };

        Vector3 rotationTerm = crossProduct< Vector3 >( Wvector, positionVector );

        double xBodyFrameVelocityInInertialFrame =
                                    inertialStateVector[ xVelocityIndex ] * std::cos( phiAngle )
                                    + inertialStateVector[ yVelocityIndex ] * std::sin( phiAngle );

        double yBodyFrameVelocityInInertialFrame =
                                    -inertialStateVector[ xVelocityIndex ] * std::sin( phiAngle )
                                    + inertialStateVector[ yVelocityIndex ] * std::cos( phiAngle );

        double zBodyFrameVelocityInInertialFrame = inertialStateVector[ zVelocityIndex ];

        double xVelocityBodyFrame = xBodyFrameVelocityInInertialFrame - rotationTerm[ 0 ];
        double yVelocityBodyFrame = yBodyFrameVelocityInInertialFrame - rotationTerm[ 1 ];
        double zVelocityBodyFrame = zBodyFrameVelocityInInertialFrame - rotationTerm[ 2 ];

        initialStateVector[ xVelocityIndex ] = xVelocityBodyFrame;
        initialStateVector[ yVelocityIndex ] = yVelocityBodyFrame;
        initialStateVector[ zVelocityIndex ] = zVelocityBodyFrame;
    }
    else
    {
        // assuming initialVector has body frame coordinates
        initialStateVector = initialVector;

        // convert the initial state vector in body frame to inertial frame and obtain orbital elements
        double phiAngle = Wvector[ zPositionIndex ] * startTime; // the angle is in radians

        double xInertial = initialStateVector[ xPositionIndex ] * std::cos( phiAngle )
                            - initialStateVector[ yPositionIndex ] * std::sin( phiAngle );
        double yInertial = initialStateVector[ xPositionIndex ] * std::sin( phiAngle )
                            + initialStateVector[ yPositionIndex ] * std::cos( phiAngle );
        double zInertial = initialStateVector[ zPositionIndex ];

        // body frame position vector
        Vector3 positionVector = { initialStateVector[ xPositionIndex ],
                                   initialStateVector[ yPositionIndex ],
                                   initialStateVector[ zPositionIndex ] };

        Vector3 wCrossR( 3 );
        wCrossR = crossProduct( Wvector, positionVector );
        double xBodyFrameVelocityInInertialFrame = initialStateVector[ xVelocityIndex ] + wCrossR[ 0 ];
        double yBodyFrameVelocityInInertialFrame = initialStateVector[ yVelocityIndex ] + wCrossR[ 1 ];
        double zBodyFrameVelocityInInertialFrame = initialStateVector[ zVelocityIndex ] + wCrossR[ 2 ];

        double vxInertial = xBodyFrameVelocityInInertialFrame * std::cos( phiAngle )
                            - yBodyFrameVelocityInInertialFrame * std::sin( phiAngle );
        double vyInertial = xBodyFrameVelocityInInertialFrame * std::sin( phiAngle )
                            + yBodyFrameVelocityInInertialFrame * std::cos( phiAngle );
        double vzInertial = zBodyFrameVelocityInInertialFrame;

        Vector6 inertialCartesianElements = { xInertial,
                                              yInertial,
                                              zInertial,
                                              vxInertial,
                                              vyInertial,
                                              vzInertial };

        convertCartesianCoordinatesToKeplerianElements( inertialCartesianElements,
                                                        gravParameter,
                                                        orbitalElements );
    }

    // Specify the step size value [s]
    static double stepSize = integrationStepSize;

    // Specify the integration time limits [s]
    const double tStart = startTime;
    const double tEnd = endTime;
    double tCurrent = tStart;

    // Calculate the initial jacobian
    double xVelocitySquare = initialStateVector[ xVelocityIndex ]
                                * initialStateVector[ xVelocityIndex ];
    double yVelocitySquare = initialStateVector[ yVelocityIndex ]
                                * initialStateVector[ yVelocityIndex ];
    double zVelocitySquare = initialStateVector[ zVelocityIndex ]
                                * initialStateVector[ zVelocityIndex ];

    double xPositionSquare = initialStateVector[ xPositionIndex ]
                                * initialStateVector[ xPositionIndex ];
    double yPositionSquare = initialStateVector[ yPositionIndex ]
                                * initialStateVector[ yPositionIndex ];

    double wSquare = Wmagnitude * Wmagnitude;

    double gravPotential = 0.0;
    naos::computeEllipsoidGravitationalPotential< double >( alpha, beta, gamma,
                                                  gravParameter,
                                                  initialStateVector[ xPositionIndex ],
                                                  initialStateVector[ yPositionIndex ],
                                                  initialStateVector[ zPositionIndex ],
                                                  gravPotential );

    double jacobian = 0.5 * ( xVelocitySquare + yVelocitySquare + zVelocitySquare )
                        - 0.5 * wSquare * ( xPositionSquare + yPositionSquare )
                        - gravPotential;

    // Save initial values in the CSV file
    eomOrbiterUREFile << initialStateVector[ xPositionIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ yPositionIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ zPositionIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ xVelocityIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ yVelocityIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ zVelocityIndex ] << ",";
    eomOrbiterUREFile << tCurrent << ",";
    eomOrbiterUREFile << jacobian << ",";
    eomOrbiterUREFile << orbitalElements[ 0 ] << ",";
    eomOrbiterUREFile << orbitalElements[ 1 ] << ",";
    eomOrbiterUREFile << orbitalElements[ 2 ] << ",";
    eomOrbiterUREFile << orbitalElements[ 3 ] << ",";
    eomOrbiterUREFile << orbitalElements[ 4 ] << ",";
    eomOrbiterUREFile << orbitalElements[ 5 ] << ",";
    eomOrbiterUREFile << stepSize << std::endl;

    // Define a vector to store latest/current state values in the integration loop
    Vector6 currentStateVector = initialStateVector;

    // Define a vector to store the integrated state vector values
    Vector6 nextStateVector = initialStateVector;

    // SET UP THE GSL INTEGRATOR
    // use RK8(9) dormand prince integrator (eight order rk integrator)
    const gsl_odeiv2_step_type * stepType = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_step * stepper = gsl_odeiv2_step_alloc( stepType, 6 );

    // set up the adaptive step size controller
    const double absoluteTolerance = 10.0e-15;
    const double relativeTolerance = 10.0e-15;

    gsl_odeiv2_control * stepControl = gsl_odeiv2_control_standard_new( absoluteTolerance,
                                                                        relativeTolerance,
                                                                        1.0,
                                                                        1.0 );

    // set up the step evolution function
    gsl_odeiv2_evolve * stepEvolution = gsl_odeiv2_evolve_alloc( 6 );

    // Start the integration outer loop
    while( tCurrent != tEnd )
    {
        // Evaluate gravitational acceleration values at the current state values
        Vector3 currentGravAcceleration( 3 );
        computeEllipsoidGravitationalAcceleration(
            alpha,
            beta,
            gamma,
            gravParameter,
            currentStateVector[ xPositionIndex ],
            currentStateVector[ yPositionIndex ],
            currentStateVector[ zPositionIndex ],
            currentGravAcceleration );

        asteroidParameters * parameters = ( asteroidParameters * )std::malloc( sizeof( asteroidParameters ) );
        parameters->xGravAcceleration = currentGravAcceleration[ xPositionIndex ];
        parameters->yGravAcceleration = currentGravAcceleration[ yPositionIndex ];
        parameters->zGravAcceleration = currentGravAcceleration[ zPositionIndex ];
        parameters->zRotation = Wvector[ zPositionIndex ];

        gsl_odeiv2_system stepSystem = { equationsOfMotion,
                                         NULL,
                                         6,
                                         parameters };

        // int stepperResetStatus = gsl_odeiv2_step_reset( stepper );
        // int evolveResetStatus = gsl_odeiv2_evolve_reset( stepEvolution );
        // if( evolveResetStatus != GSL_SUCCESS )
        // {
        //     std::ostringstream errorMessage;
        //     errorMessage << std::endl;
        //     errorMessage << "ERROR: GSL evolve reset routine failed" << std::endl;
        //     errorMessage << std::endl;
        //     throw std::runtime_error( errorMessage.str( ) );
        // }

        // perform the gsl integration step here
        double stateHolder[ 6 ] = { currentStateVector[ xPositionIndex ],
                                    currentStateVector[ yPositionIndex ],
                                    currentStateVector[ zPositionIndex ],
                                    currentStateVector[ xVelocityIndex ],
                                    currentStateVector[ yVelocityIndex ],
                                    currentStateVector[ zVelocityIndex ] };
        double timeHolder = tCurrent;

        int gslStatus = gsl_odeiv2_evolve_apply( stepEvolution,
                                                 stepControl,
                                                 stepper,
                                                 &stepSystem,
                                                 &timeHolder,
                                                 tEnd,
                                                 &stepSize,
                                                 stateHolder );

        if( gslStatus != GSL_SUCCESS )
        {
            std::ostringstream errorMessage;
            errorMessage << std::endl;
            errorMessage << "ERROR: gsl integration routine failed" << std::endl;
            errorMessage << std::endl;
            throw std::runtime_error( errorMessage.str( ) );
        }

        for( int i = 0; i < 6; i++ )
        {
            nextStateVector[ i ] = stateHolder[ i ];
        }
        double tNext = timeHolder;

        // convert the next state vector in body frame to inertial frame
        double phiAngle = Wvector[ zPositionIndex ] * tNext; // the angle is in radians

        double xInertial = nextStateVector[ xPositionIndex ] * std::cos( phiAngle )
                            - nextStateVector[ yPositionIndex ] * std::sin( phiAngle );
        double yInertial = nextStateVector[ xPositionIndex ] * std::sin( phiAngle )
                            + nextStateVector[ yPositionIndex ] * std::cos( phiAngle );
        double zInertial = nextStateVector[ zPositionIndex ];

        // body frame position vector
        Vector3 positionVector = { nextStateVector[ xPositionIndex ],
                                   nextStateVector[ yPositionIndex ],
                                   nextStateVector[ zPositionIndex ] };

        Vector3 wCrossR( 3 );
        wCrossR = crossProduct( Wvector, positionVector );
        double xBodyFrameVelocityInInertialFrame = nextStateVector[ xVelocityIndex ] + wCrossR[ 0 ];
        double yBodyFrameVelocityInInertialFrame = nextStateVector[ yVelocityIndex ] + wCrossR[ 1 ];
        double zBodyFrameVelocityInInertialFrame = nextStateVector[ zVelocityIndex ] + wCrossR[ 2 ];

        double vxInertial = xBodyFrameVelocityInInertialFrame * std::cos( phiAngle )
                            - yBodyFrameVelocityInInertialFrame * std::sin( phiAngle );
        double vyInertial = xBodyFrameVelocityInInertialFrame * std::sin( phiAngle )
                            + yBodyFrameVelocityInInertialFrame * std::cos( phiAngle );
        double vzInertial = zBodyFrameVelocityInInertialFrame;

        // check whether length of the vectors are the same between the two frames
        // frame conversions should not affect the vector lengths (magnitudes)
        Vector3 inertialPositionVector = { xInertial,
                                           yInertial,
                                           zInertial };
        double inertialFramePositionNorm = vectorNorm( inertialPositionVector );

        Vector3 bodyFramePositionVector = { nextStateVector[ xPositionIndex ],
                                            nextStateVector[ yPositionIndex ],
                                            nextStateVector[ zPositionIndex ] };
        double bodyFramePositionNorm = vectorNorm( bodyFramePositionVector );

        double precision = 10.0e-10; // in metre
        if( std::fabs( bodyFramePositionNorm - inertialFramePositionNorm ) > precision )
        {
            std::cout.precision( 20 );
            std::cout << "integration time value: " << tNext << std::endl;
            std::cout << "difference in position norms: ";
            std::cout << std::fabs( bodyFramePositionNorm - inertialFramePositionNorm );
            std::cout << std::endl;

            std::ostringstream errorMessage;
            errorMessage << std::endl;
            errorMessage << "ERROR!: the body and inertial position vector norms do not match";
            errorMessage << std::endl;
            throw std::runtime_error( errorMessage.str( ) );
        }

        Vector6 inertialCartesianElements = { xInertial,
                                              yInertial,
                                              zInertial,
                                              vxInertial,
                                              vyInertial,
                                              vzInertial };

        convertCartesianCoordinatesToKeplerianElements( inertialCartesianElements,
                                                        gravParameter,
                                                        orbitalElements );

        // Check if the new state is valid or not i.e. the particle/orbiter is not
        // inside the surface of the asteroid. If it is then terminate the outer loop. Evaluate the
        // function phi(x,y,z;0) and check if it has a positive sign (i.e. point is outside)
        double xCoordinateSquare = nextStateVector[ xPositionIndex ]
                                    * nextStateVector[ xPositionIndex ];
        double yCoordinateSquare = nextStateVector[ yPositionIndex ]
                                    * nextStateVector[ yPositionIndex ];
        double zCoordinateSquare = nextStateVector[ zPositionIndex ]
                                    * nextStateVector[ zPositionIndex ];
        double phiCheck = xCoordinateSquare / ( alpha * alpha )
                            + yCoordinateSquare / ( beta * beta )
                            + zCoordinateSquare / ( gamma * gamma )
                            - 1.0;

        // calculate the jacobian after each integration step
        xVelocitySquare = nextStateVector[ xVelocityIndex ]
                            * nextStateVector[ xVelocityIndex ];
        yVelocitySquare = nextStateVector[ yVelocityIndex ]
                            * nextStateVector[ yVelocityIndex ];
        zVelocitySquare = nextStateVector[ zVelocityIndex ]
                            * nextStateVector[ zVelocityIndex ];

        xPositionSquare = nextStateVector[ xPositionIndex ]
                            * nextStateVector[ xPositionIndex ];
        yPositionSquare = nextStateVector[ yPositionIndex ]
                            * nextStateVector[ yPositionIndex ];

        wSquare = Wmagnitude * Wmagnitude;

        gravPotential = 0.0;
        naos::computeEllipsoidGravitationalPotential< double >( alpha, beta, gamma,
                                                      gravParameter,
                                                      nextStateVector[ xPositionIndex ],
                                                      nextStateVector[ yPositionIndex ],
                                                      nextStateVector[ zPositionIndex ],
                                                      gravPotential );

        jacobian = 0.5 * ( xVelocitySquare + yVelocitySquare + zVelocitySquare )
                            - 0.5 * wSquare * ( xPositionSquare + yPositionSquare )
                            - gravPotential;

        if( phiCheck <= 0.0 )
        {
            // point is either inside or on the surface of the ellipsoid, so terminate the
            // simulator and exit the outer loop after saving the data.
            eomOrbiterUREFile << nextStateVector[ xPositionIndex ] << ",";
            eomOrbiterUREFile << nextStateVector[ yPositionIndex ] << ",";
            eomOrbiterUREFile << nextStateVector[ zPositionIndex ] << ",";
            eomOrbiterUREFile << nextStateVector[ xVelocityIndex ] << ",";
            eomOrbiterUREFile << nextStateVector[ yVelocityIndex ] << ",";
            eomOrbiterUREFile << nextStateVector[ zVelocityIndex ] << ",";
            eomOrbiterUREFile << tNext << ",";
            eomOrbiterUREFile << jacobian << ",";
            eomOrbiterUREFile << orbitalElements[ 0 ] << ",";
            eomOrbiterUREFile << orbitalElements[ 1 ] << ",";
            eomOrbiterUREFile << orbitalElements[ 2 ] << ",";
            eomOrbiterUREFile << orbitalElements[ 3 ] << ",";
            eomOrbiterUREFile << orbitalElements[ 4 ] << ",";
            eomOrbiterUREFile << orbitalElements[ 5 ] << ",";
            eomOrbiterUREFile << stepSize << std::endl;
            std::cout << "Houston, we've got a problem!" << std::endl;
            std::cout << "Particle at or inside the ellipsoidal surface" << std::endl;
            std::cout << "Event occurence at: " << tNext << " seconds" << std::endl;
            break;
        }

        // Store the integrated state values in the CSV file
        eomOrbiterUREFile << nextStateVector[ xPositionIndex ] << ",";
        eomOrbiterUREFile << nextStateVector[ yPositionIndex ] << ",";
        eomOrbiterUREFile << nextStateVector[ zPositionIndex ] << ",";
        eomOrbiterUREFile << nextStateVector[ xVelocityIndex ] << ",";
        eomOrbiterUREFile << nextStateVector[ yVelocityIndex ] << ",";
        eomOrbiterUREFile << nextStateVector[ zVelocityIndex ] << ",";
        eomOrbiterUREFile << tNext << ",";
        eomOrbiterUREFile << jacobian << ",";
        eomOrbiterUREFile << orbitalElements[ 0 ] << ",";
        eomOrbiterUREFile << orbitalElements[ 1 ] << ",";
        eomOrbiterUREFile << orbitalElements[ 2 ] << ",";
        eomOrbiterUREFile << orbitalElements[ 3 ] << ",";
        eomOrbiterUREFile << orbitalElements[ 4 ] << ",";
        eomOrbiterUREFile << orbitalElements[ 5 ] << ",";
        eomOrbiterUREFile << stepSize << std::endl;

        // save the new state values in the vector of current state values. these will be used in
        // the next loop iteration
        tCurrent = tNext;
        currentStateVector = nextStateVector;
    }
 }

} // namespace naos
