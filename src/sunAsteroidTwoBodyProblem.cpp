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
#include <stdexcept>

#include <boost/numeric/odeint.hpp>

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/misc.hpp"

namespace naos
{

//! equations of motion (sun-asteroid restricted two body problem)
/*!
 * first order differential equations describing the motion of sun around
 * asteroid modeled as a point mass gravity. EOMs are defined in an inertial frame centred at the
 * asteroid.
 */
class sunAsteroidEquationsOfMotion
{
    // declare parameters (gravitational parameter of the Sun)
    const double gravParameter;

public:
    // Default constructor with member initializer list, get the gravitational parameter
    sunAsteroidEquationsOfMotion( const double aGravParameter )
                                : gravParameter( aGravParameter )
    { }
    void operator() ( const std::vector< double > &stateVector,
                      std::vector< double > &dXdt,
                      const double currentTime )
    {
        // calculate the gravitational accelerations first
        std::vector< double > positionVector = { stateVector[ xPositionIndex ],
                                                 stateVector[ yPositionIndex ],
                                                 stateVector[ zPositionIndex ] };

        double radialDistance = vectorNorm( positionVector );
        double radialDistanceCube = radialDistance * radialDistance * radialDistance;

        double Ux = -1.0 * gravParameter * positionVector[ xPositionIndex ] / radialDistanceCube;

        double Uy = -1.0 * gravParameter * positionVector[ yPositionIndex ] / radialDistanceCube;

        double Uz = -1.0 * gravParameter * positionVector[ zPositionIndex ] / radialDistanceCube;

        // now calculate the derivatives
        dXdt[ xPositionIndex ] = stateVector[ xVelocityIndex ];
        dXdt[ yPositionIndex ] = stateVector[ yVelocityIndex ];
        dXdt[ zPositionIndex ] = stateVector[ zVelocityIndex ];
        dXdt[ xVelocityIndex ] = Ux;
        dXdt[ yVelocityIndex ] = Uy;
        dXdt[ zVelocityIndex ] = Uz;
    }
};

//! Store intermediate state values and time( if needed )
/*!
 * This structure contains members that will save all intermediate state values and times when
 * an object of this structure is passed as an argument to the integrator function.
 */
struct pushBackStateAndTime
{
    // declare containers to store state and time
    std::vector< std::vector< double > > &stateContainer;
    std::vector< double > &timeContainer;

    //member initializer list
    pushBackStateAndTime( std::vector< std::vector< double > > &aState,
                          std::vector< double > &aTime )
                : stateContainer( aState ),
                  timeContainer( aTime )
    { }

    void operator() ( const std::vector< double > &singleStateVector, const double singleTime )
    {
        // store the intermediate state and time values in the containers
        stateContainer.push_back( singleStateVector );
        timeContainer.push_back( singleTime );
    }
};

//! convert state vector from inertial frame to body frame
/*!
 * function converts state vector from inertial frame centred at the asteroid to a rotating frame
 * which is alligned with the asteroids principle axes.
 */
void convertInertialFrameVectorToBodyFrame( const std::vector< double > &asteroidRotationVector,
                                            const std::vector< double > &inertialStateVector,
                                            const double currentTime,
                                            std::vector< double > &bodyFrameStateVector )
{
    // calculate the angle
    double rotationAngle = asteroidRotationVector[ zPositionIndex ] * currentTime;

    // convert the inertial position vector to body frame position vector
    bodyFrameStateVector[ xPositionIndex ]
        = inertialStateVector[ xPositionIndex ] * std::cos( rotationAngle )
        + inertialStateVector[ yPositionIndex ] * std::sin( rotationAngle );

    bodyFrameStateVector[ yPositionIndex ]
        = -1.0 * inertialStateVector[ xPositionIndex ] * std::sin( rotationAngle )
        + inertialStateVector[ yPositionIndex ] * std::cos( rotationAngle );

    bodyFrameStateVector[ zPositionIndex ] = inertialStateVector[ zPositionIndex ];

    // apply transport theorem - get body frame velocity in inertial coordinates
    std::vector< double > inertialPositionVector = { inertialStateVector[ xPositionIndex ],
                                                     inertialStateVector[ yPositionIndex ],
                                                     inertialStateVector[ zPositionIndex ] };

    std::vector< double > omegaCrossPosition( 3, 0.0 );
    omegaCrossPosition = crossProduct( asteroidRotationVector, inertialPositionVector );

    double xBodyFrameVelocityInertialCoordinates
        = inertialStateVector[ xVelocityIndex ] - omegaCrossPosition[ 0 ];

    double yBodyFrameVelocityInertialCoordinates
        = inertialStateVector[ yVelocityIndex ] - omegaCrossPosition[ 1 ];

    double zBodyFrameVelocityInertialCoordinates
        = inertialStateVector[ zVelocityIndex ] - omegaCrossPosition[ 2 ];

    // use rotatin matrix to get the body frame velocity in body coordinates
    bodyFrameStateVector[ xVelocityIndex ]
        = xBodyFrameVelocityInertialCoordinates * std::cos( rotationAngle )
        + yBodyFrameVelocityInertialCoordinates * std::sin( rotationAngle );

    bodyFrameStateVector[ yVelocityIndex ]
        = -1.0 * xBodyFrameVelocityInertialCoordinates * std::sin( rotationAngle )
        + yBodyFrameVelocityInertialCoordinates * std::cos( rotationAngle );

    bodyFrameStateVector[ zVelocityIndex ] = zBodyFrameVelocityInertialCoordinates;
}

//! compute jacobi constant
/*!
 * compute the jacobi integral for the two body problem. It should be conserved since eom's have no
 * explicit time dependance.
 */
double computeSunAsteroidJacobiConstant( const std::vector< double > &bodyFrameStateVector,
                                         const std::vector< double > &asteroidRotationVector,
                                         const double gravParameter )
{
    // square of angular velocity magnitude
    const double omegaSquare = asteroidRotationVector[ 2 ] * asteroidRotationVector[ 2 ];

    // velocity and position squares
    const double xVelocitySquare
        = bodyFrameStateVector[ xVelocityIndex ] * bodyFrameStateVector[ xVelocityIndex ];

    const double yVelocitySquare
        = bodyFrameStateVector[ yVelocityIndex ] * bodyFrameStateVector[ yVelocityIndex ];

    const double zVelocitySquare
        = bodyFrameStateVector[ zVelocityIndex ] * bodyFrameStateVector[ zVelocityIndex ];

    const double xPositionSquare
        = bodyFrameStateVector[ xPositionIndex ] * bodyFrameStateVector[ xPositionIndex ];

    const double yPositionSquare
        = bodyFrameStateVector[ yPositionIndex ] * bodyFrameStateVector[ yPositionIndex ];

    const double zPositionSquare
        = bodyFrameStateVector[ zPositionIndex ] * bodyFrameStateVector[ zPositionIndex ];

    // point-mass grav potential
    double radialDistance = std::sqrt( xPositionSquare + yPositionSquare + zPositionSquare );

    double gravPotential = gravParameter / radialDistance;

    // jacobi constant
    double jacobiIntegral
        = 0.5 * ( xVelocitySquare + yVelocitySquare + zVelocitySquare )
        - 0.5 * omegaSquare * ( xPositionSquare + yPositionSquare )
        - gravPotential;

    // return final value
    return jacobiIntegral;
}

//! compute energy
/*!
 * computes the energy for sun orbiting the asteroid
 */
double computeSunAsteroidEnergy( const std::vector< double > &inertialStateVector,
                                 const double gravParameter )
{
    // velocity and position squares
    const double xVelocitySquare
        = inertialStateVector[ xVelocityIndex ] * inertialStateVector[ xVelocityIndex ];

    const double yVelocitySquare
        = inertialStateVector[ yVelocityIndex ] * inertialStateVector[ yVelocityIndex ];

    const double zVelocitySquare
        = inertialStateVector[ zVelocityIndex ] * inertialStateVector[ zVelocityIndex ];

    const double xPositionSquare
        = inertialStateVector[ xPositionIndex ] * inertialStateVector[ xPositionIndex ];

    const double yPositionSquare
        = inertialStateVector[ yPositionIndex ] * inertialStateVector[ yPositionIndex ];

    const double zPositionSquare
        = inertialStateVector[ zPositionIndex ] * inertialStateVector[ zPositionIndex ];

    // total energy
    const double velocity = std::sqrt( xVelocitySquare + yVelocitySquare + zVelocitySquare );
    const double distance = std::sqrt( xPositionSquare + yPositionSquare + zPositionSquare );

    const double velocitySquare = velocity * velocity;

    double energy
        = ( velocitySquare / 2.0 ) - ( gravParameter / distance );

    return energy;
}

//! sun-asteroid Restricted two body problem integration
/*!
 * integrate the equations of motion for the sun around the asteroid.
 */
void executeSunAsteroidTwoBodyProblem( const double gravParameter,
                                       const std::vector< double > &asteroidRotationVector,
                                       std::vector< double > &initialOrbitalElements,
                                       const double initialStepSize,
                                       const double startTime,
                                       const double endTime,
                                       std::ostringstream &outputFilePath,
                                       const int dataSaveIntervals )
{
    //! open the output csv file to save data. Declare file headers.
    std::ofstream outputFile;
    outputFile.open( outputFilePath.str( ) );
    outputFile << "x_body_frame" << ",";
    outputFile << "y_body_frame" << ",";
    outputFile << "z_body_frame" << ",";
    outputFile << "vx_body_frame" << ",";
    outputFile << "vy_body_frame" << ",";
    outputFile << "vz_body_frame" << ",";
    outputFile << "t" << ",";

    outputFile << "x_inertial_frame" << ",";
    outputFile << "y_inertial_frame" << ",";
    outputFile << "z_inertial_frame" << ",";
    outputFile << "vx_inertial_frame" << ",";
    outputFile << "vy_inertial_frame" << ",";
    outputFile << "vz_inertial_frame" << ",";

    outputFile << "sma" << ",";
    outputFile << "eccentricity" << ",";
    outputFile << "inclination" << ",";
    outputFile << "raan" << ",";
    outputFile << "aop" << ",";
    outputFile << "ta" << ",";

    outputFile << "jacobi" << ",";
    outputFile << "energy" << std::endl;

    // outputFile.precision( 16 );

    //! convert the initial orbital elements to cartesian state
    std::vector< double > initialState( 6, 0.0 );
    initialState = convertKeplerianElementsToCartesianCoordinates( initialOrbitalElements,
                                                                   gravParameter );

    // set up boost odeint
    const double absoluteTolerance = 1.0e-15;
    const double relativeTolerance = 1.0e-15;
    typedef boost::numeric::odeint::runge_kutta_fehlberg78< std::vector< double > > stepperType;

    // state step size guess (at each step this initial guess will be used)
    double stepSizeGuess = initialStepSize;

    // initialize the ode system
    sunAsteroidEquationsOfMotion sunAsteroidTwoBodyProblem( gravParameter );

    // initialize current state vector and time
    std::vector< double > currentStateVector = initialState;
    double currentTime = startTime;
    double intermediateEndTime = currentTime + dataSaveIntervals;

    // save the initial state in body frame
    std::vector< double > bodyFrameStateVector( 6, 0.0 );
    convertInertialFrameVectorToBodyFrame( asteroidRotationVector,
                                           currentStateVector,
                                           currentTime,
                                           bodyFrameStateVector );

    outputFile << bodyFrameStateVector[ xPositionIndex ] << ",";
    outputFile << bodyFrameStateVector[ yPositionIndex ] << ",";
    outputFile << bodyFrameStateVector[ zPositionIndex ] << ",";
    outputFile << bodyFrameStateVector[ xVelocityIndex ] << ",";
    outputFile << bodyFrameStateVector[ yVelocityIndex ] << ",";
    outputFile << bodyFrameStateVector[ zVelocityIndex ] << ",";
    outputFile << currentTime << ",";

    // save the initial state vector (inertial frame)
    outputFile << currentStateVector[ xPositionIndex ] << ",";
    outputFile << currentStateVector[ yPositionIndex ] << ",";
    outputFile << currentStateVector[ zPositionIndex ] << ",";
    outputFile << currentStateVector[ xVelocityIndex ] << ",";
    outputFile << currentStateVector[ yVelocityIndex ] << ",";
    outputFile << currentStateVector[ zVelocityIndex ] << ",";

    // save the initial orbital elements, jacobi and energy
    std::vector< double > orbitalElements( 6, 0.0 );
    convertCartesianCoordinatesToKeplerianElements( currentStateVector,
                                                    gravParameter,
                                                    orbitalElements );

    outputFile << orbitalElements[ 0 ] << ",";
    outputFile << orbitalElements[ 1 ] << ",";
    outputFile << orbitalElements[ 2 ] << ",";
    outputFile << orbitalElements[ 3 ] << ",";
    outputFile << orbitalElements[ 4 ] << ",";
    outputFile << orbitalElements[ 5 ] << ",";

    double jacobi = computeSunAsteroidJacobiConstant( bodyFrameStateVector,
                                                      asteroidRotationVector,
                                                      gravParameter );

    double energy = computeSunAsteroidEnergy( currentStateVector,
                                              gravParameter );

    outputFile << jacobi << ",";
    outputFile << energy << std::endl;

    // start the integration outer loop
    while( intermediateEndTime <= endTime )
    {
        // perform integration, integrated result stored in currentStateVector
        size_t steps = boost::numeric::odeint::integrate_adaptive(
                            make_controlled( absoluteTolerance, relativeTolerance, stepperType( ) ),
                            sunAsteroidTwoBodyProblem,
                            currentStateVector,
                            currentTime,
                            intermediateEndTime,
                            stepSizeGuess );

        // update the time variables
        currentTime = intermediateEndTime;
        intermediateEndTime = currentTime + dataSaveIntervals;

        // save data
        convertInertialFrameVectorToBodyFrame( asteroidRotationVector,
                                               currentStateVector,
                                               currentTime,
                                               bodyFrameStateVector );

        outputFile << bodyFrameStateVector[ xPositionIndex ] << ",";
        outputFile << bodyFrameStateVector[ yPositionIndex ] << ",";
        outputFile << bodyFrameStateVector[ zPositionIndex ] << ",";
        outputFile << bodyFrameStateVector[ xVelocityIndex ] << ",";
        outputFile << bodyFrameStateVector[ yVelocityIndex ] << ",";
        outputFile << bodyFrameStateVector[ zVelocityIndex ] << ",";
        outputFile << currentTime << ",";

        outputFile << currentStateVector[ xPositionIndex ] << ",";
        outputFile << currentStateVector[ yPositionIndex ] << ",";
        outputFile << currentStateVector[ zPositionIndex ] << ",";
        outputFile << currentStateVector[ xVelocityIndex ] << ",";
        outputFile << currentStateVector[ yVelocityIndex ] << ",";
        outputFile << currentStateVector[ zVelocityIndex ] << ",";

        convertCartesianCoordinatesToKeplerianElements( currentStateVector,
                                                        gravParameter,
                                                        orbitalElements );

        outputFile << orbitalElements[ 0 ] << ",";
        outputFile << orbitalElements[ 1 ] << ",";
        outputFile << orbitalElements[ 2 ] << ",";
        outputFile << orbitalElements[ 3 ] << ",";
        outputFile << orbitalElements[ 4 ] << ",";
        outputFile << orbitalElements[ 5 ] << ",";

        jacobi = computeSunAsteroidJacobiConstant( bodyFrameStateVector,
                                                   asteroidRotationVector,
                                                   gravParameter );

        energy = computeSunAsteroidEnergy( currentStateVector,
                                           gravParameter );

        outputFile << jacobi << ",";
        outputFile << energy << std::endl;
    } // end of outer while loop for integration
    outputFile.close( );
}

} // namespace naos
