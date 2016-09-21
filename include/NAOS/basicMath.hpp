/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef BASIC_MATH_HPP
#define BASIC_MATH_HPP

#include <cmath>

#include "NAOS/constants.hpp"

namespace naos
{

//! Compute vector cross product
/*!
 * Compute cross product of two vectors, each with 3 elements.
 *
 * @tparam Vector3      Template paramter for a 3 element vector
 * @param[in] vector1   A 3 element vector
 * @param[in] vector2   A 3 element vector
 * @return result       The resulting vector from a cross product
 */
template< typename Vector3 >
Vector3 crossProduct( const Vector3 &vector1, const Vector3 &vector2 )
{
    Vector3 result = vector1;
    // Compute components of resulting 3-vector.
    result[ 0 ] = vector1[ 1 ] * vector2[ 2 ] - vector1[ 2 ] * vector2[ 1 ];
    result[ 1 ] = vector1[ 2 ] * vector2[ 0 ] - vector1[ 0 ] * vector2[ 2 ];
    result[ 2 ] = vector1[ 0 ] * vector2[ 1 ] - vector1[ 1 ] * vector2[ 0 ];

    return result;
}

//! Compute vector dot product
/*!
 * Compute dot product of two vectors, each with 3 elements.
 *
 * @tparam Vector3      Template paramter for a 3 element vector
 * @param[in] vector1   A 3 element vector
 * @param[in] vector2   A 3 element vector
 * @return result       The resulting vector from a dot product
 */
template< typename Vector3 >
double dotProduct( const Vector3 &vector1, const Vector3 &vector2 )
{
    double result = 0.0;
    // Compute dot product.
    result = vector1[ 0 ] * vector2[ 0 ]
                + vector1[ 1 ] * vector2[ 1 ]
                + vector1[ 2 ] * vector2[ 2 ];

    return result;
}

//! Compute vector norm
/*!
 * Compute norm of a vector of 3 elements.
 *
 * @tparam Vector3      Template paramter for a 3 element vector
 * @param[in] vector1   A 3 element vector
 * @return result       The resulting norm of the vector
 */
template< typename Vector3 >
double vectorNorm( const Vector3 &vector1 )
{
    double result = 0.0;
    // Compute Norm.
    double xSquare = vector1[ 0 ] * vector1[ 0 ];
    double ySquare = vector1[ 1 ] * vector1[ 1 ];
    double zSquare = vector1[ 2 ] * vector1[ 2 ];
    result = std::sqrt( xSquare + ySquare + zSquare );

    return result;
}

//! Convert degrees to radians
/*!
 * convert degrees to radians
 *
 * @param[in] degreeAngle   angle in degrees
 * @return    radAngle      angle in radians
 */
 template< typename Real >
 Real convertDegreeToRadians( const Real degreeAngle )
 {
    Real result = degreeAngle * naos::PI / 180.0;
    return result;
 }

//! Convert radians to degrees
/*!
 * convert radians to degrees
 *
 * @param[in] radAngle      angle in radians
 * @return    degreeAngle   angle in degrees
 */
 template< typename Real >
 Real convertRadiansToDegree( const Real radAngle )
 {
    Real result = radAngle * 180.0 / naos::PI;
    return result;
 }

} // namespace naos

#endif // BASIC_MATH_HPP
