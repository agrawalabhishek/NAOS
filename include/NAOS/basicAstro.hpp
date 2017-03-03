/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef BASIC_ASTRO_HPP
#define BASIC_ASTRO_HPP

#include <cmath>
#include <limits>
#include <iostream>

#include "NAOS/basicMath.hpp"
#include "NAOS/constants.hpp"

namespace naos
{

//! Convert orbital elements to cartesian elements
/*!
 * Given a set of classical orbital elements, convert them to equivalent cartesian elements. The
 * final cartesian elements are defined in the body centred inertial frame. The formulas were taken
 * from DA Vallado's book - Fundamentals of astrodynamics and applications
 *
 * @tparam Vector6                  Template paramter for a 6 element vector
 * @param[in] keplerianElements     A 6 element vector containing keplerian elements in this order:
 *                                  semi major axis [m] (or semilatus rectum for a parabolic orbit),
 *                                  eccentricity [-],
 *                                  inclination [deg],
 *                                  RAAN (right ascension of ascending node) [deg],
 *                                  AOP (argument of periapsis) [deg],
 *                                  TA (true anomaly) [deg]
 * @param[in] gravParameter         Gravitational parameter of the central body in SI units
 * @return    result                A 6 element vector containing the converted cartesian elements
 */
template< typename Vector6 >
Vector6 convertKeplerianElementsToCartesianCoordinates(
    const Vector6 &keplerianElements,
    const double gravParameter )
{
    // Initialize the results vector
    Vector6 cartesianElements = keplerianElements;

    // Extract the keplerian elements
    const double semiMajorAxis  = keplerianElements[ 0 ];
    const double eccentricity   = keplerianElements[ 1 ];
    const double inclination    = naos::convertDegreeToRadians< double >( keplerianElements[ 2 ] );
    const double RAAN           = naos::convertDegreeToRadians< double >( keplerianElements[ 3 ] );
    const double AOP            = naos::convertDegreeToRadians< double >( keplerianElements[ 4 ] );
    const double TA             = naos::convertDegreeToRadians< double >( keplerianElements[ 5 ] );

    // Compute semi-latus rectum
    double semiLatus = 0.0;
    // semiLatus for non parabolic orbits (eccentricity != 1)
    const double tolerance = 10.0 * std::numeric_limits< double >::epsilon( );
    if( std::fabs( eccentricity - 1.0 ) > tolerance )
    {
        semiLatus = semiMajorAxis * ( 1.0 - eccentricity * eccentricity );
    }
    // for parabolic orbits
    else
    {
        semiLatus = keplerianElements[ 0 ];
    }

    // compute position and velocity in the perifocal coordinate system
    const double radius = semiLatus / ( 1.0 + eccentricity * std::cos( TA ) );
    const double xPositionPerifocal = radius * std::cos( TA );
    const double yPositionPerifocal = radius * std::sin( TA );
    const double xVelocityPerifocal = -1.0 * std::sqrt( gravParameter / semiLatus ) * std::sin( TA );
    const double yVelocityPerifocal = std::sqrt( gravParameter / semiLatus )
                                        * ( eccentricity + std::cos( TA ) );

    // Calculate components of the rotation matrix
    const double rotation11 = std::cos( RAAN ) * std::cos( AOP )
                                - std::sin( RAAN ) * std::sin( AOP ) * std::cos( inclination );
    const double rotation12 = -1.0 * std::cos( RAAN ) * std::sin( AOP )
                                - std::sin( RAAN ) * std::cos( AOP ) * std::cos( inclination );
    const double rotation21 = std::sin( RAAN ) * std::cos( AOP )
                                + std::cos( RAAN ) * std::sin( AOP ) * std::cos( inclination );
    const double rotation22 = -1.0 * std::sin( RAAN ) * std::sin( AOP )
                                + std::cos( RAAN ) * std::cos( AOP ) * std::cos( inclination );
    const double rotation31 = std::sin( AOP ) * std::sin( inclination );
    const double rotation32 = std::cos( AOP ) * std::sin( inclination );

    // Compute the cartesian position and velocity in the inertial frame
    cartesianElements[ naos::xPositionIndex ] = rotation11 * xPositionPerifocal
                                                    + rotation12 * yPositionPerifocal;
    cartesianElements[ naos::yPositionIndex ] = rotation21 * xPositionPerifocal
                                                    + rotation22 * yPositionPerifocal;
    cartesianElements[ naos::zPositionIndex ] = rotation31 * xPositionPerifocal
                                                    + rotation32 * yPositionPerifocal;
    cartesianElements[ naos::xVelocityIndex ] = rotation11 * xVelocityPerifocal
                                                    + rotation12 * yVelocityPerifocal;
    cartesianElements[ naos::yVelocityIndex ] = rotation21 * xVelocityPerifocal
                                                    + rotation22 * yVelocityPerifocal;
    cartesianElements[ naos::zVelocityIndex ] = rotation31 * xVelocityPerifocal
                                                    + rotation32 * yVelocityPerifocal;

    return cartesianElements;
}

//! Convert Cartesian coordinates to orbital elements
/*!
 * Convert a set of cartesian position and velocity (defined in the inertial frame) to orbital
 * elements. The equations were obtained from DA Vallodo's book - fundamentals of astrodynamics and
 * applications.
 *
 * @param[in] cartesianElements     6 element vector containing position and velocity of the orbiter
 * @param[in] gravParameter         the gravitational parameter of the central body
 * @param[in/out] orbitalElements   output vector containing the final orbital elements in the
 *                                  following order:
 *                                  semi-major axis / semi-latus rectum (semiMajor/ semiLatus)
 *                                  eccentricity
 *                                  inclination
 *                                  Right ascension of ascending node (RAAN)
 *                                  argument of periapsis (AOP)
 *                                  true anomaly (TA)
 */
template< typename Vector6 >
void convertCartesianCoordinatesToKeplerianElements(
    const Vector6 &cartesianElements,
    const double gravParameter,
    Vector6 &orbitalElements )
{
    // get position and velocity vectors
    std::vector< double > position( 3 );
    position[ 0 ] = cartesianElements[ xPositionIndex ];
    position[ 1 ] = cartesianElements[ yPositionIndex ];
    position[ 2 ] = cartesianElements[ zPositionIndex ];
    double positionMagnitude = vectorNorm( position );
    std::vector < double > positionUnitVector( 3 );
    for( int i = 0; i <= 2; i++ )
    {
        positionUnitVector[ i ] = position[ i ] / positionMagnitude;
    }

    std::vector< double > velocity( 3 );
    velocity[ 0 ] = cartesianElements[ xVelocityIndex ];
    velocity[ 1 ] = cartesianElements[ yVelocityIndex ];
    velocity[ 2 ] = cartesianElements[ zVelocityIndex ];
    double velocityMagnitude = vectorNorm( velocity );
    std::vector < double > velocityUnitVector( 3 );
    for( int i = 0; i <= 2; i++ )
    {
        velocityUnitVector[ i ] = velocity[ i ] / velocityMagnitude;
    }

    // get the angular momentum vector
    std::vector< double > angularMomentum( 3 );
    angularMomentum = crossProduct( position, velocity );

    // get the unit vector for the ascending node
    std::vector< double > node( 3 );
    std::vector< double > zUnitVector = { 0.0, 0.0, 1.0 };
    node = normalize( crossProduct( zUnitVector, normalize( angularMomentum ) ) );

    // get the eccentricity vector
    std::vector< double > eccentricityVector( 3 );
    std::vector< double > tempVector = crossProduct( velocity, angularMomentum );
    for( int i = 0; i <= 2; i++ )
    {
        eccentricityVector[ i ] = tempVector[ i ] / gravParameter - positionUnitVector[ i ];
    }

    // get the specific mechanical energy
    double specificMechanicalEnergy = velocityMagnitude * velocityMagnitude / 2.0
                                        - gravParameter / positionMagnitude;

    // calculate eccentricity from the eccentricity vector
    const double eccentricity = vectorNorm( eccentricityVector );
    orbitalElements[ 1 ] = eccentricity;

    // calculate the semiAxis or semiLatus depending on the eccentricity value
    const double tolerance = 10.0 * std::numeric_limits< double >::epsilon( );
    if( std::fabs( eccentricity - 1.0 ) < tolerance )
    {
        // case of a parabolic orbit, so semiLatus is calculated
        double semiLatus = vectorNorm( angularMomentum ) * vectorNorm( angularMomentum )
                             / gravParameter;
        orbitalElements[ 0 ] = semiLatus;
    }
    else
    {
        // double semiAxis = -1.0 * gravParameter / ( 2 * specificMechanicalEnergy );
        double semiLatus = vectorNorm( angularMomentum ) * vectorNorm( angularMomentum )
                             / gravParameter;
        double semiAxis = semiLatus / ( 1.0 - eccentricity * eccentricity );
        orbitalElements[ 0 ] = semiAxis;
    }

    // calculate the inclination
    double inclinationRadians = std::acos( angularMomentum[ 2 ] / vectorNorm( angularMomentum ) );
    double inclination = convertRadiansToDegree( inclinationRadians );
    orbitalElements[ 2 ] = inclination;

    // calculate the right ascenion of ascending node
    if( inclination < tolerance )
    {
        // orbit is circular hence RAAN = 0
        double RAAN = 0.0;
        orbitalElements[ 3 ] = RAAN;
    }
    else
    {
        double RAAN = std::acos( node[ 0 ] );
        RAAN = convertRadiansToDegree( RAAN );
        // quadrant check
        if( node[ 1 ] < 0.0 )
        {
            RAAN = 360.0 - RAAN;
        }
        orbitalElements[ 3 ] = RAAN;
    }

    // calculate the argument of periapsis
    if( eccentricity < tolerance )
    {
        // orbit is circular hence AOP = 0
        double AOP = 0.0;
        orbitalElements[ 4 ] = AOP;
    }
    else if( inclination < tolerance )
    {
        // elliptical equatorial special case orbit
        double AOP = 0.0;
        AOP = std::acos( eccentricityVector[ 0 ] / eccentricity );
        AOP = convertRadiansToDegree( AOP );
        if( eccentricityVector[ 1 ] < 0 )
        {
            AOP = 360.0 - AOP;
        }
        orbitalElements[ 4 ] = AOP;
    }
    else
    {
        double AOP = std::acos( dotProduct( node, eccentricityVector )
                        / ( vectorNorm( node ) * vectorNorm( eccentricityVector ) ) );
        AOP = convertRadiansToDegree( AOP );
        // quadrant check
        if( eccentricityVector[ 2 ] < 0 )
        {
            AOP = 360.0 - AOP;
        }
        orbitalElements[ 4 ] = AOP;
    }

    // calculate the true anomaly
    double TA = std::acos( dotProduct( eccentricityVector, position )
                    / ( vectorNorm( eccentricityVector ) * vectorNorm( position ) ) );
    TA = convertRadiansToDegree( TA );
    // quadrant check
    if( dotProduct( position, velocity ) < 0 )
    {
        TA = 360.0 - TA;
    }
    orbitalElements[ 5 ] = TA;
}

} // namespace naos

#endif // BASIC_ASTRO_HPP
