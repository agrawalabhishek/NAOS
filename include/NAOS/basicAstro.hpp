/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef BASIC_ASTRO_HPP
#define BASIC_ASTRO_HPP

#include <cmath>
#include <limits>

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

} // namespace naos

#endif // BASIC_ASTRO_HPP
