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

#include "NAOS/basicAstro.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/constants.hpp"
#include "NAOS/misc.hpp"
#include "NAOS/sunAsteroidKeplerProblemSolver.hpp"

namespace naos
{

//! open sun ephemeris and extract data for a specified time value
/*!
 * This file opens the ephemeris data for the sun around the asteroid at specified time values in
 * multiples of ten. If a finer ephemeris is required then the sunAsteroidTwoBodyProblem must be
 * re-run prior to using this mode, with a finer data save interval value.
 */
std::vector< double > extractSunEphemeris( const double timeValue,
                                           std::ostringstream &ephemerisFilePath )
{
    // open the file
    std::ifstream data;
    data.open( ephemerisFilePath.str( ) );

    // initialize the variables to read the data
    std::string line;
    std::getline( data, line ); // reads the header of the input csv file

    const int timeIndex = 6;

    std::vector< double > bodyFrameStateVectorAndTime( 7, 0.0 );

    bool valueFoundFlag = false;

    // read data line by line and convert to orbital elements
    while( std::getline( data, line ) )
    {
        std::stringstream linestream( line );
        std::string cell;
        std::vector< double > dataVector( 21, 0.0 );
        int index = 0;

        while( std::getline( linestream, cell, ',' ) )
        {
            std::stringstream stringData( cell );
            double floatingPointData = 1.0;
            stringData >> floatingPointData;
            dataVector[ index ] = floatingPointData;
            index++;
        }

        if( std::fabs( dataVector[ timeIndex ] - timeValue ) <= 1.0e-10 )
        {
            valueFoundFlag = true;
            for( int i = 0; i < 7; i++ )
            {
                bodyFrameStateVectorAndTime[ i ] = dataVector[ i ];
            }
            // break out of outer most while loop
            break;
        }
    }

    if( !valueFoundFlag )
    {
        std::ostringstream errorMessage;
        errorMessage << std::endl;
        errorMessage << "Error in calculating perturbing accelerations cpp file" << std::endl;
        errorMessage << "Sun's ephemeris could not be obtained because of illegal time value" << std::endl;
        errorMessage << "Time value = " << timeValue << std::endl;
        errorMessage << std::endl;
        throw std::runtime_error( errorMessage.str( ) );
    }
    // close the csv file
    data.close( );

    // return the data
    return bodyFrameStateVectorAndTime;
}

//! Third-body effect
/*!
 * Compute the perturbing acceleration from the thrid body effect of the sun.
 */
std::vector< double > computeSunThirdBodyEffectAcceleration( const std::vector< double > &regolithPositionVector,
                                                             std::vector< double > &asteroidRotationVector,
                                                             const double initialTime,
                                                             const double initialSunMeanAnomalyRadian,
                                                             const std::vector< double > &initialSunOrbitalElements,
                                                             const double timeValue,
                                                             const double sunMeanMotion )
{
    // extract body frame state vector for the sun's motion around the asteroid for the given time
    // value
    // std::ostringstream sunAsteroidFilePath;
    // sunAsteroidFilePath << "../../data/sun_asteroid_2BP/sunAsteroid2BP.csv";

    // std::vector< double > bodyFrameStateVector( 7, 0.0 );
    // bodyFrameStateVector = extractSunEphemeris( timeValue, sunAsteroidFilePath );

    // accessed 3 jan 2016 from:
    // http://ssd.jpl.nasa.gov/?constants
    const double sunGravParameter = 1.32712440018 * 1.0e+20;

    // solve the kepler problem to get the sun's position for the required time value
    // get the mean anomaly for the time value
    // double timeDifference = std::fabs( timeValue - initialTime );
    double timeDifference = std::fabs( timeValue );
    double meanAnomalyRadian = sunMeanMotion * timeDifference + initialSunMeanAnomalyRadian;

    double eccentricity = initialSunOrbitalElements[ 1 ];
    double eccentricitySquare = eccentricity * eccentricity;

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

    double trueAnomalyRadian
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
    orbitalElements[ 0 ] = initialSunOrbitalElements[ 0 ];
    orbitalElements[ 1 ] = initialSunOrbitalElements[ 1 ];
    orbitalElements[ 2 ] = initialSunOrbitalElements[ 2 ];
    orbitalElements[ 3 ] = initialSunOrbitalElements[ 3 ];
    orbitalElements[ 4 ] = initialSunOrbitalElements[ 4 ];
    orbitalElements[ 5 ] = trueAnomaly;

    // get the inertial state vector from the updated orbital elements vector
    std::vector< double > currentInertialStateVector( 6, 0.0 );
    currentInertialStateVector
        = convertKeplerianElementsToCartesianCoordinates( orbitalElements,
                                                          sunGravParameter );

    // get the asteroid-body frame state vector for the Sun
    std::vector< double > sunStateVectorInRotatingFrame( 6, 0.0 );
    convertInertialFrameVectorToBodyFrame( asteroidRotationVector,
                                           currentInertialStateVector,
                                           timeValue,
                                           sunStateVectorInRotatingFrame );

    // form the body frame position vector of the sun
    std::vector< double > sunPositionVector = { sunStateVectorInRotatingFrame[ xPositionIndex ],
                                                sunStateVectorInRotatingFrame[ yPositionIndex ],
                                                sunStateVectorInRotatingFrame[ zPositionIndex ] };

    double sunPositionVectorMagnitude = vectorNorm( sunPositionVector );

    double sunPositionVectorMagnitudeCube
        = sunPositionVectorMagnitude * sunPositionVectorMagnitude * sunPositionVectorMagnitude;

    // compute the perturbing acceleration
    std::vector< double > regolithPositionFromSun( 3, 0.0 );
    regolithPositionFromSun[ xPositionIndex ]
        = regolithPositionVector[ xPositionIndex ] - sunPositionVector[ xPositionIndex ];

    regolithPositionFromSun[ yPositionIndex ]
        = regolithPositionVector[ yPositionIndex ] - sunPositionVector[ yPositionIndex ];

    regolithPositionFromSun[ zPositionIndex ]
        = regolithPositionVector[ zPositionIndex ] - sunPositionVector[ zPositionIndex ];

    double regolithPositionFromSunMagnitude = vectorNorm( regolithPositionFromSun );

    double regolithPositionFromSunMagnitudeCube
        = regolithPositionFromSunMagnitude * regolithPositionFromSunMagnitude * regolithPositionFromSunMagnitude;

    std::vector< double > accelerationTermOne( 3, 0.0 );
    std::vector< double > accelerationTermTwo( 3, 0.0 );

    accelerationTermOne[ xPositionIndex ]
        = ( -1.0 * sunGravParameter / regolithPositionFromSunMagnitudeCube ) * regolithPositionFromSun[ xPositionIndex ];

    accelerationTermOne[ yPositionIndex ]
        = ( -1.0 * sunGravParameter / regolithPositionFromSunMagnitudeCube ) * regolithPositionFromSun[ yPositionIndex ];

    accelerationTermOne[ zPositionIndex ]
        = ( -1.0 * sunGravParameter / regolithPositionFromSunMagnitudeCube ) * regolithPositionFromSun[ zPositionIndex ];

    accelerationTermTwo[ xPositionIndex ]
        = ( -1.0 * sunGravParameter / sunPositionVectorMagnitudeCube ) * sunPositionVector[ xPositionIndex ];

    accelerationTermTwo[ yPositionIndex ]
        = ( -1.0 * sunGravParameter / sunPositionVectorMagnitudeCube ) * sunPositionVector[ yPositionIndex ];

    accelerationTermTwo[ zPositionIndex ]
        = ( -1.0 * sunGravParameter / sunPositionVectorMagnitudeCube ) * sunPositionVector[ zPositionIndex ];

    std::vector< double > thirdBodyPerturbingAcceleration( 3, 0.0 );
    thirdBodyPerturbingAcceleration[ xPositionIndex ]
        = accelerationTermOne[ xPositionIndex ] + accelerationTermTwo[ xPositionIndex ];

    thirdBodyPerturbingAcceleration[ yPositionIndex ]
        = accelerationTermOne[ yPositionIndex ] + accelerationTermTwo[ yPositionIndex ];

    thirdBodyPerturbingAcceleration[ zPositionIndex ]
        = accelerationTermOne[ zPositionIndex ] + accelerationTermTwo[ zPositionIndex ];

    return thirdBodyPerturbingAcceleration;
}

//! Solar radiation pressure
/*!
 * comupte the perturbing acceleration from solar radition pressure.
 */
std::vector< double > computeSolarRadiationPressureAcceleration( const std::vector< double > &regolithPositionVector,
                                                                 std::vector< double > &asteroidRotationVector,
                                                                 const double initialTime,
                                                                 const double initialSunMeanAnomalyRadian,
                                                                 const std::vector< double > &initialSunOrbitalElements,
                                                                 const double timeValue,
                                                                 const double sunMeanMotion )
{
    // accessed 3 jan 2016 from:
    // http://ssd.jpl.nasa.gov/?constants
    const double sunGravParameter = 1.32712440018 * 1.0e+20;

    // solve the kepler problem to get the sun's position for the required time value
    // get the mean anomaly for the time value
    // double timeDifference = std::fabs( timeValue - initialTime );
    double timeDifference = std::fabs( timeValue );
    double meanAnomalyRadian = sunMeanMotion * timeDifference + initialSunMeanAnomalyRadian;

    double eccentricity = initialSunOrbitalElements[ 1 ];
    double eccentricitySquare = eccentricity * eccentricity;

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

    double trueAnomalyRadian
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
    orbitalElements[ 0 ] = initialSunOrbitalElements[ 0 ];
    orbitalElements[ 1 ] = initialSunOrbitalElements[ 1 ];
    orbitalElements[ 2 ] = initialSunOrbitalElements[ 2 ];
    orbitalElements[ 3 ] = initialSunOrbitalElements[ 3 ];
    orbitalElements[ 4 ] = initialSunOrbitalElements[ 4 ];
    orbitalElements[ 5 ] = trueAnomaly;

    // get the inertial state vector from the updated orbital elements vector
    std::vector< double > currentInertialStateVector( 6, 0.0 );
    currentInertialStateVector
        = convertKeplerianElementsToCartesianCoordinates( orbitalElements,
                                                          sunGravParameter );

    // get the asteroid-body frame state vector for the Sun
    std::vector< double > sunStateVectorInRotatingFrame( 6, 0.0 );
    convertInertialFrameVectorToBodyFrame( asteroidRotationVector,
                                           currentInertialStateVector,
                                           timeValue,
                                           sunStateVectorInRotatingFrame );

    // form the body frame position vector of the sun
    std::vector< double > sunPositionVector = { sunStateVectorInRotatingFrame[ xPositionIndex ],
                                                sunStateVectorInRotatingFrame[ yPositionIndex ],
                                                sunStateVectorInRotatingFrame[ zPositionIndex ] };

    // regolith position from sun
    std::vector< double > regolithPositionFromSun( 3, 0.0 );
    regolithPositionFromSun[ xPositionIndex ]
        = -1.0 * regolithPositionVector[ xPositionIndex ] + sunPositionVector[ xPositionIndex ];

    regolithPositionFromSun[ yPositionIndex ]
        = -1.0 * regolithPositionVector[ yPositionIndex ] + sunPositionVector[ yPositionIndex ];

    regolithPositionFromSun[ zPositionIndex ]
        = -1.0 * regolithPositionVector[ zPositionIndex ] + sunPositionVector[ zPositionIndex ];

    double regolithPositionFromSunMagnitude = vectorNorm( regolithPositionFromSun );

    double regolithPositionFromSunMagnitudeCube
        = regolithPositionFromSunMagnitude * regolithPositionFromSunMagnitude * regolithPositionFromSunMagnitude;

    // specify the albedo value, for max acceleration from SRP, set rho=1
    const double rho = 1.0;

    // specify the solar constant value in SI units [kg m/s^2]
    const double solarConstant = 1.0e17;

    // specify the area to mass ratio for the regolith
    const double regolithGrainDensity = 3.2 * 1.0e3;
    const double regolithGrainRadius = 1.0 * 1.0e-2;
    const double areaToMassRatio = 3.0 / ( regolithGrainRadius * regolithGrainDensity );

    // compute the acceleration
    const double multiplicationConstant = -1.0 * ( 1.0 + rho ) * solarConstant * areaToMassRatio;

    std::vector< double > solarRadiationPerturbingAcceleration( 3, 0.0 );

    solarRadiationPerturbingAcceleration[ xPositionIndex ]
        = multiplicationConstant * regolithPositionFromSun[ xPositionIndex ] / regolithPositionFromSunMagnitudeCube;

    solarRadiationPerturbingAcceleration[ yPositionIndex ]
        = multiplicationConstant * regolithPositionFromSun[ yPositionIndex ] / regolithPositionFromSunMagnitudeCube;

    solarRadiationPerturbingAcceleration[ zPositionIndex ]
        = multiplicationConstant * regolithPositionFromSun[ zPositionIndex ] / regolithPositionFromSunMagnitudeCube;

    return solarRadiationPerturbingAcceleration;
}

} // namespace naos
