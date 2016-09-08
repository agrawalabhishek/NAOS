/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef BASIC_MATH_HPP
#define BASIC_MATH_HPP

#include <stdexcept>
#include <cmath>

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

} // namespace naos

#endif // BASIC_MATH_HPP
