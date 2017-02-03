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
#include <utility>

#include <boost/math/tools/roots.hpp>

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/misc.hpp"

namespace naos
{

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

//! define the equations for the kepler problem
/*!
 * the structure returns, in the form of a st::pair, the values for the eccentric-mean anomaly
 * equation and its first derivative. These are used in the function to convert mean anomaly
 * to eccentric anomaly. Provide the angle values in rads for input. The output will be in rads as
 * well.
 */
struct keplerProblem
{
    // declare parameters
    double eccentricity = 0.0;
    double meanAnomaly = 0.0; // radians

public:
    // default constructor with member initializer list
    keplerProblem(
        const double anEccentricityValue,
        const double aMeanAnomalyValue )
        : eccentricity( anEccentricityValue ),
          meanAnomaly( aMeanAnomalyValue )
    { }
    std::pair< double, double > operator() ( double const& eccentricAnomaly )
    {
        // NOTE- all angles are in radians
        // compute f(E)
        double eccentricAnomalyFunction
            = eccentricAnomaly - ( eccentricity * std::sin( eccentricAnomaly ) ) - meanAnomaly;

        // compute the first derivative
        double firstDerivative
            = 1.0 - ( eccentricity * std::cos( eccentricAnomaly ) );

        // return the values
        return std::make_pair( eccentricAnomalyFunction, firstDerivative );
    }
};

struct keplerProblemSecondOrder
{
    // declare parameters
    double eccentricity = 0.0;
    double meanAnomaly = 0.0; // radians

public:
    // default constructor with member initializer list
    keplerProblemSecondOrder(
        const double anEccentricityValue,
        const double aMeanAnomalyValue )
        : eccentricity( anEccentricityValue ),
          meanAnomaly( aMeanAnomalyValue )
    { }
    std::tuple< double, double, double > operator() ( double const& eccentricAnomaly )
    {
        // NOTE- all angles are in radians
        // compute f(E)
        double eccentricAnomalyFunction
            = eccentricAnomaly - ( eccentricity * std::sin( eccentricAnomaly ) ) - meanAnomaly;

        // compute the first derivative
        double firstDerivative
            = 1.0 - ( eccentricity * std::cos( eccentricAnomaly ) );

        // compute second derivative
        double secondDerivative
            = eccentricity * std::sin( eccentricAnomaly );

        // return the values
        return std::make_tuple( eccentricAnomalyFunction, firstDerivative, secondDerivative );
    }
};

//! convert mean anomaly to eccentric anomaly
/*!
 * This function solves the kepler problem by using newton-rhapson method to solve for eccentric
 * anomaly from the mean anomaly. Only the case for eccentric and circular orbits is covered here.
 * the mean anomaly input has to be in radians here.
 */
double convertMeanAnomalyToEccentricAnomaly( const double eccentricity,
                                             const double meanAnomaly )
{
    //NOTE- input mean anomaly to be in radian
    // initialize for newton-rhapson method
    double initialGuess = meanAnomaly;
    double minEccentricAnomaly = -1.0 * naos::PI;
    double maxEccentricAnomaly = naos::PI;

    // set the desired accuracy
    const int digits = std::numeric_limits< double >::digits;
    // Accuracy doubles with each step, so stop when we have
    // just over half the digits correct. (BOOST lib)
    // int binaryAccuracy = static_cast< int >( 0.6 * digits );
    int binaryAccuracy = digits;

    // set max iterations
    const boost::uintmax_t maxIterationsAllowed = 50;
    // set variable to store actual iterations undertaken in newton-rhapson
    boost::uintmax_t iterationCounter = maxIterationsAllowed;

    // form an object of the struct keplerProblem
    keplerProblem keplerProblemObject( eccentricity, meanAnomaly );
    // keplerProblemSecondOrder keplerProblemObject( eccentricity, meanAnomaly );

    // perform the newton-rhapson iteration
    double eccentricAnomaly = boost::math::tools::newton_raphson_iterate(
                                    keplerProblemObject,
                                    initialGuess,
                                    minEccentricAnomaly,
                                    maxEccentricAnomaly,
                                    binaryAccuracy,
                                    iterationCounter );

    // perform the halley's iteration
    // double eccentricAnomaly = boost::math::tools::halley_iterate(
    //                                 keplerProblemObject,
    //                                 initialGuess,
    //                                 minEccentricAnomaly,
    //                                 maxEccentricAnomaly,
    //                                 binaryAccuracy,
    //                                 iterationCounter );

    // perform the schroder's iteration method
    // double eccentricAnomaly = boost::math::tools::schroder_iterate(
    //                                 keplerProblemObject,
    //                                 initialGuess,
    //                                 minEccentricAnomaly,
    //                                 maxEccentricAnomaly,
    //                                 binaryAccuracy,
    //                                 iterationCounter );

    return eccentricAnomaly;
}

//! sun-asteroid Restricted two body problem integration
/*!
 * integrate the equations of motion for the sun around the asteroid.
 */
void executeSunAsteroidKeplerSolver( const double gravParameter,
                                     const std::vector< double > &asteroidRotationVector,
                                     std::vector< double > &initialOrbitalElements,
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
    outputFile << "eccentric_anomaly" << ",";
    outputFile << "mean_anomaly" << ",";

    outputFile << "jacobi" << ",";
    outputFile << "energy" << std::endl;

    outputFile.precision( 20 );

    //! convert the initial orbital elements to cartesian state
    std::vector< double > initialState( 6, 0.0 );
    initialState = convertKeplerianElementsToCartesianCoordinates( initialOrbitalElements,
                                                                   gravParameter );

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
    outputFile << initialOrbitalElements[ 0 ] << ",";
    outputFile << initialOrbitalElements[ 1 ] << ",";
    outputFile << initialOrbitalElements[ 2 ] << ",";
    outputFile << initialOrbitalElements[ 3 ] << ",";
    outputFile << initialOrbitalElements[ 4 ] << ",";
    outputFile << initialOrbitalElements[ 5 ] << ",";

    // get the initial eccentric anomaly
    double trueAnomalyRadian = naos::convertDegreeToRadians( initialOrbitalElements[ 5 ] );
    double eccentricity = initialOrbitalElements[ 1 ];
    double eccentricitySquare = eccentricity * eccentricity;

    double sineOfEccentricAnomaly
        = ( std::sqrt( 1.0 - eccentricitySquare ) * std::sin( trueAnomalyRadian ) )
        / ( 1.0 + eccentricity * std::cos( trueAnomalyRadian ) );

    double cosineOfEccentricAnomaly
        = ( eccentricity + std::cos( trueAnomalyRadian ) )
        / ( 1.0 + eccentricity * std::cos( trueAnomalyRadian ) );

    // value returned in the range of -pi to +pi radians
    double initialEccentricAnomalyRadian
        = std::atan2( sineOfEccentricAnomaly, cosineOfEccentricAnomaly );

    // get the initial mean anomaly
    double initialMeanAnomalyRadian
        = initialEccentricAnomalyRadian - eccentricity * std::sin( initialEccentricAnomalyRadian );

    outputFile << naos::convertRadiansToDegree( initialEccentricAnomalyRadian ) << ",";
    outputFile << naos::convertRadiansToDegree( initialMeanAnomalyRadian ) << ",";

    // get the mean motion
    double semiMajorAxis = initialOrbitalElements[ 0 ];
    double semiMajorAxisCube = semiMajorAxis * semiMajorAxis * semiMajorAxis;
    double meanMotion = std::sqrt( gravParameter / semiMajorAxisCube );

    // get the jacobi and the energy
    double jacobi = computeSunAsteroidJacobiConstant( bodyFrameStateVector,
                                                      asteroidRotationVector,
                                                      gravParameter );

    double energy = computeSunAsteroidEnergy( currentStateVector,
                                              gravParameter );

    outputFile << jacobi << ",";
    outputFile << energy << std::endl;

    // start the propagator outer loop
    while( intermediateEndTime <= endTime )
    {
        // get the mean anomaly for the next time value
        double timeDifference = intermediateEndTime - startTime;
        double meanAnomalyRadian = meanMotion * timeDifference + initialMeanAnomalyRadian;

        // get the corresponding eccentric anomaly
        double eccentricAnomalyRadian
            = convertMeanAnomalyToEccentricAnomaly( eccentricity,
                                                    meanAnomalyRadian );

        // get the corresponding true anomaly
        double sineOfTrueAnomaly
            = std::sqrt( 1.0 - eccentricitySquare ) * std::sin( eccentricAnomalyRadian )
            / ( 1.0 - eccentricity * std::cos( eccentricAnomalyRadian ) );

        double cosineOfTrueAnomaly
            = ( std::cos( eccentricAnomalyRadian ) - eccentricity )
            / ( 1.0 - eccentricity * std::cos( eccentricAnomalyRadian ) );

        trueAnomalyRadian
            = std::atan2( sineOfTrueAnomaly, cosineOfTrueAnomaly );

        double trueAnomaly = naos::convertRadiansToDegree( trueAnomalyRadian );

        if( trueAnomaly < 0.0 )
        {
            trueAnomaly = trueAnomaly + 360.0;
        }

        // form the orbital elements vector
        // basically update the orbital elements true anomaly,
        // everything else remains the same
        std::vector< double > orbitalElements( 6, 0.0 );
        orbitalElements[ 0 ] = initialOrbitalElements[ 0 ];
        orbitalElements[ 1 ] = initialOrbitalElements[ 1 ];
        orbitalElements[ 2 ] = initialOrbitalElements[ 2 ];
        orbitalElements[ 3 ] = initialOrbitalElements[ 3 ];
        orbitalElements[ 4 ] = initialOrbitalElements[ 4 ];
        orbitalElements[ 5 ] = trueAnomaly;

        // get the inertial state vector from the updated orbital elements vector
        std::vector< double > currentInertialStateVector( 6, 0.0 );
        currentInertialStateVector
            = convertKeplerianElementsToCartesianCoordinates( orbitalElements,
                                                              gravParameter );

        // update the time variables
        currentTime = intermediateEndTime;
        intermediateEndTime = currentTime + dataSaveIntervals;

        // save data
        convertInertialFrameVectorToBodyFrame( asteroidRotationVector,
                                               currentInertialStateVector,
                                               currentTime,
                                               bodyFrameStateVector );

        outputFile << bodyFrameStateVector[ xPositionIndex ] << ",";
        outputFile << bodyFrameStateVector[ yPositionIndex ] << ",";
        outputFile << bodyFrameStateVector[ zPositionIndex ] << ",";
        outputFile << bodyFrameStateVector[ xVelocityIndex ] << ",";
        outputFile << bodyFrameStateVector[ yVelocityIndex ] << ",";
        outputFile << bodyFrameStateVector[ zVelocityIndex ] << ",";
        outputFile << currentTime << ",";

        outputFile << currentInertialStateVector[ xPositionIndex ] << ",";
        outputFile << currentInertialStateVector[ yPositionIndex ] << ",";
        outputFile << currentInertialStateVector[ zPositionIndex ] << ",";
        outputFile << currentInertialStateVector[ xVelocityIndex ] << ",";
        outputFile << currentInertialStateVector[ yVelocityIndex ] << ",";
        outputFile << currentInertialStateVector[ zVelocityIndex ] << ",";

        outputFile << orbitalElements[ 0 ] << ",";
        outputFile << orbitalElements[ 1 ] << ",";
        outputFile << orbitalElements[ 2 ] << ",";
        outputFile << orbitalElements[ 3 ] << ",";
        outputFile << orbitalElements[ 4 ] << ",";
        outputFile << orbitalElements[ 5 ] << ",";

        outputFile << naos::convertRadiansToDegree( eccentricAnomalyRadian ) << ",";
        outputFile << naos::convertRadiansToDegree( meanAnomalyRadian ) << ",";

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
