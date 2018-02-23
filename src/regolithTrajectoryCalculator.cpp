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
#include <vector>
#include <limits>
#include <stdexcept>

#include <SQLiteCpp/SQLiteCpp.h>
#include <sqlite3.h>

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/misc.hpp"
#include "NAOS/particleAroundUniformlyRotatingEllipsoid.hpp"
#include "NAOS/regolithTrajectoryCalculator.hpp"
#include "NAOS/ellipsoidPotential.hpp"

namespace naos
{

void computeRegolithVelocityVector2( std::vector< double > regolithPositionVector,
                                     const double velocityMagnitude,
                                     const double coneAngleAzimuth,
                                     const double coneAngleDeclination,
                                     std::vector< double > &unitNormalVector,
                                     std::vector< double > &regolithVelocityVector,
                                     std::vector< double > &regolithDirectionUnitVector )
{
    // get the z basis vector of the surface frame
    std::vector< double > zUnitVector = unitNormalVector;

    std::vector< double > xUnitVector( 3 );
    std::vector< double > yUnitVector( 3 );

    // body frame principal axis Z
    std::vector< double > zPrincipalAxisBodyFrame { 0.0, 0.0, 1.0 };
    std::vector< double > zNegativePrincipalAxisBodyFrame { 0.0, 0.0, -1.0 };

    // check if the position vector is along the poles
    const double positionDotPrincipalZ
            = dotProduct( normalize( regolithPositionVector ), zPrincipalAxisBodyFrame );
    const double positionDotNegativePrincipalZ
            = dotProduct( normalize( regolithPositionVector ), zNegativePrincipalAxisBodyFrame );

    if( positionDotPrincipalZ == 1.0 )
    {
        // the position vector is pointing to the poles, hence x basis vector pointing to the north
        // direction wouldn't work
        xUnitVector = { 1.0, 0.0, 0.0 };
        yUnitVector = { 0.0, 1.0, 0.0 };
    }
    else if( positionDotNegativePrincipalZ == 1.0 )
    {
        xUnitVector = { -1.0, 0.0, 0.0 };
        yUnitVector = { 0.0, 1.0, 0.0 };
    }
    else
    {
        // do the regular thing where the x basis is pointing to north
        // get the intermediate RTN frame at the surface point
        std::vector< double > unitR = normalize( regolithPositionVector );

        std::vector< double > unitT = normalize( crossProduct( unitR, zPrincipalAxisBodyFrame ) );

        // get the x basis vector, pointing to north
        xUnitVector = normalize( crossProduct( unitT, zUnitVector ) );

        // get the y basis vector
        yUnitVector = normalize( crossProduct( zUnitVector, xUnitVector ) );
    }

    const double cosDelta = std::cos( coneAngleDeclination );
    const double sinDelta = std::sin( coneAngleDeclination );
    const double cosGamma = std::cos( coneAngleAzimuth );
    const double sinGamma = std::sin( coneAngleAzimuth );

    regolithVelocityVector[ 0 ] = velocityMagnitude * ( cosDelta * zUnitVector[ 0 ]
                                                        + sinDelta * cosGamma * xUnitVector[ 0 ]
                                                        + sinDelta * sinGamma * yUnitVector[ 0 ] );

    regolithVelocityVector[ 1 ] = velocityMagnitude * ( cosDelta * zUnitVector[ 1 ]
                                                        + sinDelta * cosGamma * xUnitVector[ 1 ]
                                                        + sinDelta * sinGamma * yUnitVector[ 1 ] );

    regolithVelocityVector[ 2 ] = velocityMagnitude * ( cosDelta * zUnitVector[ 2 ]
                                                        + sinDelta * cosGamma * xUnitVector[ 2 ]
                                                        + sinDelta * sinGamma * yUnitVector[ 2 ] );

    regolithDirectionUnitVector[ 0 ] = 1.0 * ( cosDelta * zUnitVector[ 0 ]
                                               + sinDelta * cosGamma * xUnitVector[ 0 ]
                                               + sinDelta * sinGamma * yUnitVector[ 0 ] );

    regolithDirectionUnitVector[ 1 ] = 1.0 * ( cosDelta * zUnitVector[ 1 ]
                                               + sinDelta * cosGamma * xUnitVector[ 1 ]
                                               + sinDelta * sinGamma * yUnitVector[ 1 ] );

    regolithDirectionUnitVector[ 2 ] = 1.0 * ( cosDelta * zUnitVector[ 2 ]
                                               + sinDelta * cosGamma * xUnitVector[ 2 ]
                                               + sinDelta * sinGamma * yUnitVector[ 2 ] );

    regolithDirectionUnitVector = normalize( regolithDirectionUnitVector );
}

//! Compute trajectory for regolith ejected from the surface of an asteroid (data saved in SQL db)
/*!
 * This routine computes the trajectory for a single regolith ejected from the surface of an
 * asteroid by first computing the appropriate initial conditions (IC) and then integrating the
 * equations of motion using those ICs.
 *
 */
void executeRegolithTrajectoryCalculation( const double alpha,
                                           const double beta,
                                           const double gamma,
                                           const double gravitationalParameter,
                                           const std::vector< double > W,
                                           const double Wmagnitude,
                                           const double aXValue,
                                           const double aYValue,
                                           const double aZValue,
                                           const double coneAngleAzimuth,
                                           const double coneAngleDeclination,
                                           const double velocityMagnitudeFactor,
                                           const double integrationStepSize,
                                           const double startTime,
                                           const double endTime,
                                           const double dataSaveIntervals,
                                           SQLite::Statement &databaseQuery,
                                           std::ofstream &escapeSpeedFileHandle )
{
    // get the initial position coordinates for the particle on the surface
    std::vector< double > aVector { aXValue, aYValue, aZValue };
    // std::vector< double > positionUnitVector = naos::normalize( aVector );
    std::vector< double > positionUnitVector = aVector;

    double xTempTerm
        = ( positionUnitVector[ xPositionIndex ] * positionUnitVector[ xPositionIndex ] )
            / ( alpha * alpha );
    double yTempTerm
        = ( positionUnitVector[ yPositionIndex ] * positionUnitVector[ yPositionIndex ] )
            / ( beta * beta );
    double zTempTerm
        = ( positionUnitVector[ zPositionIndex ] * positionUnitVector[ zPositionIndex ] )
            / ( gamma * gamma );
    double positionMagnitude = std::sqrt( 1.0 / ( xTempTerm + yTempTerm + zTempTerm ) );

    std::vector< double > regolithPositionVector( 3 );
    regolithPositionVector[ xPositionIndex ]
                            = positionMagnitude * positionUnitVector[ xPositionIndex ];
    regolithPositionVector[ yPositionIndex ]
                            = positionMagnitude * positionUnitVector[ yPositionIndex ];
    regolithPositionVector[ zPositionIndex ]
                            = positionMagnitude * positionUnitVector[ zPositionIndex ];

    // get the normal vector for the point on the surface from where the particle is launched
    // this normal vector is not in the same direction as the position vector unless the central
    // body is a sphere.
    std::vector< double > regolithNormalVector
                               = { 2.0 * regolithPositionVector[ xPositionIndex ] / ( alpha * alpha ),
                                   2.0 * regolithPositionVector[ yPositionIndex ] / ( beta * beta ),
                                   2.0 * regolithPositionVector[ zPositionIndex ] / ( gamma * gamma ) };

    // get the unit normal vector
    std::vector< double > regolithNormalUnitVector = naos::normalize( regolithNormalVector );

    // get the velocity vector
    // calculate local normal escape speed (ref: DJ Scheeres, orbits close to asteroid 4769 castalia)
    double velocityMagnitude = 0.0;

    std::vector< double > omegaCrossPosition( 3, 0.0 );
    omegaCrossPosition = crossProduct( W, regolithPositionVector );
    double omegaCrossPositionSquare = dotProduct( omegaCrossPosition, omegaCrossPosition );

    double normalDotOmegaCrossPosition = dotProduct( regolithNormalUnitVector, omegaCrossPosition );
    double normalDotOmegaCrossPositionSquare = normalDotOmegaCrossPosition * normalDotOmegaCrossPosition;

    positionMagnitude = vectorNorm( regolithPositionVector );
    double pointMassGravPotential = gravitationalParameter / positionMagnitude;

    double ellipsoidGravPotential;
    computeEllipsoidGravitationalPotential( alpha,
                                            beta,
                                            gamma,
                                            gravitationalParameter,
                                            regolithPositionVector[ xPositionIndex ],
                                            regolithPositionVector[ yPositionIndex ],
                                            regolithPositionVector[ zPositionIndex ],
                                            ellipsoidGravPotential );
    double maxPotentialValue = 0.0;
    if( ellipsoidGravPotential > pointMassGravPotential )
    {
        maxPotentialValue = ellipsoidGravPotential;
        // std::cout << "maxPotentialValue = ellipsoid model = " << maxPotentialValue << std::endl;
    }
    else
    {
        maxPotentialValue = pointMassGravPotential;
        // std::cout << "maxPotentialValue = point mass model" << maxPotentialValue << std::endl;
    }

    double localNormalEscapeSpeed
        = -1.0 * normalDotOmegaCrossPosition
            + std::sqrt( normalDotOmegaCrossPositionSquare + 2.0 * maxPotentialValue - omegaCrossPositionSquare );

    velocityMagnitude = velocityMagnitudeFactor * localNormalEscapeSpeed;

    // Sun in equatorial plane of asteroid, circular orbit, varying only the longitude of Sun
    for( int sunLongitudeIterator = 225; sunLongitudeIterator <= 225; sunLongitudeIterator = sunLongitudeIterator + 90 )
    {
        double sunLongitude = 1.0 * sunLongitudeIterator;
        std::vector< double > initialSunOrbitalElements = { 1.0 * naos::oneAstronomicalUnit,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            sunLongitude };

        for( int velocityIterator = 9; velocityIterator <= 9; velocityIterator = velocityIterator + 1 )
        {
            // all particles, irrespective of surface location will have the same launch velocity as computed
            // below.
            velocityMagnitude = 1.0 * velocityIterator;

            std::vector< double > regolithVelocityVector( 3 );

            // get the regolith launch direction unit vector
            std::vector< double > regolithDirectionUnitVector( 3, 0.0 );

            computeRegolithVelocityVector2( regolithPositionVector,
                                            velocityMagnitude,
                                            coneAngleAzimuth,
                                            coneAngleDeclination,
                                            regolithNormalUnitVector,
                                            regolithVelocityVector,
                                            regolithDirectionUnitVector );

            // compute local escape velocity magnitude in the direction of the launched regolith
            // i.e. not in the normal direction!
            double regolithDirectionDotOmegaCrossPosition = dotProduct( regolithDirectionUnitVector,
                                                                        omegaCrossPosition );
            double regolithDirectionDotOmegaCrossPositionSquare
                        = regolithDirectionDotOmegaCrossPosition * regolithDirectionDotOmegaCrossPosition;

            // conservative inertial escape speed calculation
            double conservativeInertialEscapeSpeedSquared = 2.0 * maxPotentialValue;

            // calculate the non conservative inertial escape speed magnitude and store it in
            // a seperate .csv file instead of the database
            double angularRate = W[ 2 ];

            std::vector< double > initialStateVectorInertialFrame( 6, 0.0 );
            std::vector< double > initialStateVectorRotatingFrame( 6, 0.0 );
            initialStateVectorRotatingFrame[ 0 ] = regolithPositionVector[ xPositionIndex ];
            initialStateVectorRotatingFrame[ 1 ] = regolithPositionVector[ yPositionIndex ];
            initialStateVectorRotatingFrame[ 2 ] = regolithPositionVector[ zPositionIndex ];
            initialStateVectorRotatingFrame[ 3 ] = regolithVelocityVector[ 0 ];
            initialStateVectorRotatingFrame[ 4 ] = regolithVelocityVector[ 1 ];
            initialStateVectorRotatingFrame[ 5 ] = regolithVelocityVector[ 2 ];
            convertBodyFrameVectorToInertialFrame( W,
                                                   initialStateVectorRotatingFrame,
                                                   0.0,
                                                   initialStateVectorInertialFrame );

            std::vector< double > initialPositionVectorInertialFrame( 3, 0.0 );
            std::vector< double > initialVelocityVectorInertialFrame( 3, 0.0 );

            initialPositionVectorInertialFrame[ 0 ] = initialStateVectorInertialFrame[ xPositionIndex ];
            initialPositionVectorInertialFrame[ 1 ] = initialStateVectorInertialFrame[ yPositionIndex ];
            initialPositionVectorInertialFrame[ 2 ] = initialStateVectorInertialFrame[ zPositionIndex ];

            initialVelocityVectorInertialFrame[ 0 ] = initialStateVectorInertialFrame[ xVelocityIndex ];
            initialVelocityVectorInertialFrame[ 1 ] = initialStateVectorInertialFrame[ yVelocityIndex ];
            initialVelocityVectorInertialFrame[ 2 ] = initialStateVectorInertialFrame[ zVelocityIndex ];

            double omegaDot_Hcap = dotProduct( W,
                                        normalize(
                                            crossProduct( initialPositionVectorInertialFrame,
                                                          initialVelocityVectorInertialFrame ) ) );

            double tempTerm1 = positionMagnitude * ( std::sin( coneAngleDeclination ) ) * omegaDot_Hcap;
            double tempTerm1Squared = tempTerm1 * tempTerm1;

            int qInfinityUpperLimit = static_cast< int >( alpha );

            double qInfinity = 0.3 * alpha;
            double qInfinityTerm
                = 2.0 * angularRate * std::sqrt( 2.0 * gravitationalParameter * qInfinity );

            double inertialEscapeSpeedWithPlus
                = tempTerm1 + std::sqrt( tempTerm1Squared - qInfinityTerm + 2.0 * ellipsoidGravPotential );

            double inertialEscapeSpeedWithPlusSquared
                = inertialEscapeSpeedWithPlus * inertialEscapeSpeedWithPlus;

            double inertialEscapeSpeedWithMinus
                = tempTerm1 - std::sqrt( tempTerm1Squared - qInfinityTerm + 2.0 * ellipsoidGravPotential );

            double bodyFrameDirectionalEscapeSpeedWithMinus = 0.0;
            if( inertialEscapeSpeedWithMinus >= 0.0 )
            {
                 double inertialEscapeSpeedWithMinusSquared
                      = inertialEscapeSpeedWithMinus * inertialEscapeSpeedWithMinus;

                 bodyFrameDirectionalEscapeSpeedWithMinus
                     = -1.0 * regolithDirectionDotOmegaCrossPosition
                     + std::sqrt( regolithDirectionDotOmegaCrossPositionSquare
                     + inertialEscapeSpeedWithMinusSquared - omegaCrossPositionSquare );
             }
             else
             {
                bodyFrameDirectionalEscapeSpeedWithMinus = 0.0;
             }

             double bodyFrameDirectionalEscapeSpeedWithPlus
                 = -1.0 * regolithDirectionDotOmegaCrossPosition
                 + std::sqrt( regolithDirectionDotOmegaCrossPositionSquare
                 + inertialEscapeSpeedWithPlusSquared - omegaCrossPositionSquare );

             // escapeSpeedFileHandle << convertRadiansToDegree( coneAngleDeclination ) << ",";
             // escapeSpeedFileHandle << qInfinity << ",";
             // escapeSpeedFileHandle << inertialEscapeSpeedWithPlus << ",";
             // escapeSpeedFileHandle << inertialEscapeSpeedWithMinus << ",";
             // escapeSpeedFileHandle << bodyFrameDirectionalEscapeSpeedWithPlus << ",";
             // escapeSpeedFileHandle << bodyFrameDirectionalEscapeSpeedWithMinus << std::endl;

            const double localRegolithDirectionEscapeSpeed
                                = -1.0 * regolithDirectionDotOmegaCrossPosition
                                + std::sqrt( regolithDirectionDotOmegaCrossPositionSquare
                                + conservativeInertialEscapeSpeedSquared - omegaCrossPositionSquare );

            // get the inertial directional escape velocity magnitude.
            // Note - accounting for non-zero start time
            std::vector< double > inertialDirectionalEscapeVelocityInBodyComponents( 3, 0.0 );
            std::vector< double > inertialDirectionalEscapeVelocity( 3, 0.0 );

            inertialDirectionalEscapeVelocityInBodyComponents[ 0 ]
                = localRegolithDirectionEscapeSpeed * regolithDirectionUnitVector[ 0 ] + omegaCrossPosition[ 0 ];

            inertialDirectionalEscapeVelocityInBodyComponents[ 1 ]
                = localRegolithDirectionEscapeSpeed * regolithDirectionUnitVector[ 1 ] + omegaCrossPosition[ 1 ];

            inertialDirectionalEscapeVelocityInBodyComponents[ 2 ]
                = localRegolithDirectionEscapeSpeed * regolithDirectionUnitVector[ 2 ] + omegaCrossPosition[ 2 ];

            // use the rotation matrix to get the escape velocity in the inertial frame
            double rotationAngle = W[ 2 ] * startTime;

            inertialDirectionalEscapeVelocity[ 0 ]
                        = inertialDirectionalEscapeVelocityInBodyComponents[ 0 ] * std::cos( rotationAngle )
                        - inertialDirectionalEscapeVelocityInBodyComponents[ 1 ] * std::sin( rotationAngle );

            inertialDirectionalEscapeVelocity[ 1 ]
                        = inertialDirectionalEscapeVelocityInBodyComponents[ 0 ] * std::sin( rotationAngle )
                        + inertialDirectionalEscapeVelocityInBodyComponents[ 1 ] * std::cos( rotationAngle );

            inertialDirectionalEscapeVelocity[ 2 ] = inertialDirectionalEscapeVelocityInBodyComponents[ 2 ];

            const double inertialDirectionalEscapeSpeed = vectorNorm( inertialDirectionalEscapeVelocity );

            // form the initial state vector defined in body frame
            naos::Vector6 initialVector { regolithPositionVector[ xPositionIndex ],
                                          regolithPositionVector[ yPositionIndex ],
                                          regolithPositionVector[ zPositionIndex ],
                                          regolithVelocityVector[ 0 ],
                                          regolithVelocityVector[ 1 ],
                                          regolithVelocityVector[ 2 ] };

            // calculate the trajectory of the regolith (uses SQL db to save data)
            executeSingleRegolithTrajectoryCalculation( alpha,
                                                        beta,
                                                        gamma,
                                                        gravitationalParameter,
                                                        W,
                                                        initialVector,
                                                        coneAngleAzimuth,
                                                        coneAngleDeclination,
                                                        localRegolithDirectionEscapeSpeed,
                                                        inertialDirectionalEscapeSpeed,
                                                        integrationStepSize,
                                                        startTime,
                                                        endTime,
                                                        initialSunOrbitalElements,
                                                        databaseQuery,
                                                        dataSaveIntervals );
        } // end of for loop iterating over different magnitudes for the regolith velocity
    } // end of solar phase angle iterator loop
}

} // namespace naos
