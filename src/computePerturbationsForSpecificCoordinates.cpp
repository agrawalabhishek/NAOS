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

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/misc.hpp"
#include "NAOS/ellipsoidGravitationalAcceleration.hpp"
#include "NAOS/perturbingAccelerations.hpp"

namespace naos
{

//! compute the gravitational and perturbing accelerations
/*!
 * this function computes and stores the gravitational and perturbing accelerations on the surface
 * an asteroid. The data can be plotted to analyze the relative strengths of the perturbing
 * accelerations.
 */
void computeGravitationalAndPerturbingAccelerationsForSpecificCoordinates( const double alpha,
                                                                           const double beta,
                                                                           const double gamma,
                                                                           const double gravParameter,
                                                                           std::vector< double > &asteroidRotationVector,
                                                                           std::ostringstream &filePath )
{
    // open the csv file for writing data
    std::ofstream outputData;
    outputData.open( filePath.str( ) );

    // declare data headers
    outputData << "xPosition" << ",";
    outputData << "yPosition" << ",";
    outputData << "zPosition" << ",";

    outputData << "xGravAcceleration" << ",";
    outputData << "yGravAcceleration" << ",";
    outputData << "zGravAcceleration" << ",";
    outputData << "gravAcceleration" << ",";

    outputData << "xThirdBodyAcceleration" << ",";
    outputData << "yThirdBodyAcceleration" << ",";
    outputData << "zThirdBodyAcceleration" << ",";
    outputData << "thirdBodyAcceleration" << ",";

    outputData << "xSRPAcceleration" << ",";
    outputData << "ySRPAcceleration" << ",";
    outputData << "zSRPAcceleration" << ",";
    outputData << "srpAcceleration" << std::endl;

    int xPoint = 119474;
    int yPoint = -158118;
    int zPoint = -26946;

    // start position iterator and compute accelerations for those points
    for( int xPosition = xPoint - 500; xPosition <= xPoint + 500; xPosition = xPosition + 10 )
    {
        for( int yPosition = yPoint - 500; yPosition <= yPoint + 500; yPosition = yPosition + 10 )
        {
            for( int zPosition = zPoint - 500; zPosition <= zPoint + 500; zPosition = zPosition + 10 )
            {
                // compute the surface coordinate
                double xCoordinate = xPosition;
                double yCoordinate = yPosition;
                double zCoordinate = zPosition;

                // form a vector for the position
                const std::vector< double > regolithPositionVector = { xCoordinate,
                                                                       yCoordinate,
                                                                       zCoordinate };

                // save the coordinates and the lat long values in the output file
                outputData << xCoordinate << ",";
                outputData << yCoordinate << ",";
                outputData << zCoordinate << ",";

                // compute gravitational acceleration
                std::vector< double > gravAcceleration( 3, 0.0 );
                computeEllipsoidGravitationalAcceleration( alpha,
                                                           beta,
                                                           gamma,
                                                           gravParameter,
                                                           xCoordinate,
                                                           yCoordinate,
                                                           zCoordinate,
                                                           gravAcceleration );

                double gravAccelerationMagnitude = vectorNorm( gravAcceleration );

                // save the gravitational acceleration in output file
                outputData << gravAcceleration[ 0 ] << ",";
                outputData << gravAcceleration[ 1 ] << ",";
                outputData << gravAcceleration[ 2 ] << ",";
                outputData << gravAccelerationMagnitude << ",";

                // set initial conditions for the Sun (to compute perturbations)
                // NOTE: initial true anomaly for sun's position is independant of the starting time for regolith
                // simulation
                const double oneAstronomicalUnit = 149597870700.0;
                const double initialTimeForSun = 0.0;
                std::vector< double > initialSunOrbitalElements = { 1.457945652635353 * oneAstronomicalUnit,
                                                                    0.2225680937603629,
                                                                    10.82771477612614,
                                                                    0.0,
                                                                    0.0,
                                                                    0.0 };

                // accessed 3 jan 2016 from:
                // http://ssd.jpl.nasa.gov/?constants
                const double sunGravParameter = 1.32712440018 * 10.0e+20;

                // get the initial eccentric anomaly
                double trueAnomalyRadian = naos::convertDegreeToRadians( initialSunOrbitalElements[ 5 ] );
                double eccentricity = initialSunOrbitalElements[ 1 ];
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
                double initialSunMeanAnomalyRadian
                    = initialEccentricAnomalyRadian - eccentricity * std::sin( initialEccentricAnomalyRadian );

                // get the mean motion
                double semiMajorAxis = initialSunOrbitalElements[ 0 ];
                double semiMajorAxisCube = semiMajorAxis * semiMajorAxis * semiMajorAxis;
                double sunMeanMotion = std::sqrt( sunGravParameter / semiMajorAxisCube );

                // compute the perturbing accelerations due to third body effect
                double perturbationTimeValue = 165420.0;
                std::vector< double > thirdBodyEffect( 3, 0.0 );
                thirdBodyEffect = computeSunThirdBodyEffectAcceleration(
                                                regolithPositionVector,
                                                asteroidRotationVector,
                                                initialTimeForSun,
                                                initialSunMeanAnomalyRadian,
                                                initialSunOrbitalElements,
                                                perturbationTimeValue,
                                                sunMeanMotion );

                double thirdBodyEffectMagnitude = vectorNorm( thirdBodyEffect );

                outputData << thirdBodyEffect[ 0 ] << ",";
                outputData << thirdBodyEffect[ 1 ] << ",";
                outputData << thirdBodyEffect[ 2 ] << ",";
                outputData << thirdBodyEffectMagnitude << ",";

                // compute the perturbing acceleration due to SRP effect
                perturbationTimeValue = 165420.0;
                std::vector< double > srpEffect( 3, 0.0 );
                srpEffect = computeSolarRadiationPressureAcceleration(
                                                regolithPositionVector,
                                                asteroidRotationVector,
                                                initialTimeForSun,
                                                initialSunMeanAnomalyRadian,
                                                initialSunOrbitalElements,
                                                perturbationTimeValue,
                                                sunMeanMotion );

                double srpEffectMagnitude = vectorNorm( srpEffect );

                outputData << srpEffect[ 0 ] << ",";
                outputData << srpEffect[ 1 ] << ",";
                outputData << srpEffect[ 2 ] << ",";
                outputData << srpEffectMagnitude << std::endl;
            } // end zPosition loop
        } // end yPosition loop
    } // end xPosition loop
    outputData.close( );
}

} // namespace naos
