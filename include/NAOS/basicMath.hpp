/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef BASIC_MATH_HPP
#define BASIC_MATH_HPP

#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <vector>

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

//! Compute unit vector
/*!
 * normalize a vector of 3 elements.
 *
 * @tparam Vector3      Template paramter for a 3 element vector
 * @param[in] vector    A 3 element vector
 * @return result       The resulting unit vector
 */
template< typename Vector3 >
Vector3 normalize( const Vector3 &vector )
{
    Vector3 result = vector;
    // normalize
    const double magnitude = vectorNorm< Vector3 >( vector );
    for( int i = 0; i < 3; i++ )
    {
        result[ i ] = vector[ i ] / magnitude;
    }

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

//! calculate mod 360 of an angle
/*!
 * calculate the mod 360 of an angle
 *
 * @param[in] degreeAngle   angle in degrees
 * @return    result        mod 360 of the angle
 */
 template< typename Real >
 Real mod360( const Real degreeAngle )
 {
    Real result = std::fmod( degreeAngle, 360.0 );
    return result;
 }

//! Multiply two matrices
/*!
 * multiply two matrices and return the final resultant matrix. The matrices, both input and output
 * are in the form of c++11 std::vectors. for the sake of simplicity, all the matrices should be
 * given in a 2D format, even if they are just row or column vectors.
 *
 * @param[in] firstMatrix       The matrix on the left hand side of the multiplication symbol
 * @param[in] secondMatrix      The matrix on the right hand side of the multiplication symbol
 * @param[in/out] result        The final product matrix
 */
 template< typename Vector2D >
 void matrixMultiplication( Vector2D &firstMatrix,
                            Vector2D &secondMatrix,
                            Vector2D &result,
                            const int firstMatrixRowSize,
                            const int firstMatrixColumnSize,
                            const int secondMatrixRowSize,
                            const int secondMatrixColumnSize )
 {
    // check if the input matrices are of proper dimensions
    if( firstMatrixColumnSize != secondMatrixRowSize )
    {
        std::ostringstream errorMessage;
        errorMessage << std::endl;
        errorMessage << "ERROR in matrix multiplication, dimension mismatch" << std::endl;
        errorMessage << std::endl;
        throw std::runtime_error( errorMessage.str( ) );
    }

    // perform the multiplication
    for( int i = 0; i < firstMatrixRowSize; i++ )
    {
        for( int j = 0; j < secondMatrixColumnSize; j++ )
        {
            for( int k = 0; k < firstMatrixColumnSize; k++ )
            {
                result[ i ][ j ] += firstMatrix[ i ][ k ] * secondMatrix[ k ][ j ];
            }
        }
    }
 }

//! transpose of a matrix
/*!
 * This routine computes the transpose of a NxN matrix
 *
 * @param[in] matrix        The input matrix whose transpose has to be computed
 * @param[in/out] result    The transposed matrix as output
 */
 template< typename Vector2D >
 void matrixTranspose( Vector2D &matrix, Vector2D &result )
 {
    // get the input matrix's dimensions and check if they are equal or not
    const int row = matrix.size( );
    const int column = matrix[ 0 ].size( );

    if( row != column )
    {
        std::ostringstream errorMessage;
        errorMessage << std::endl;
        errorMessage << "ERROR in matrix transpose, dimension mismatch" << std::endl;
        errorMessage << std::endl;
        throw std::runtime_error( errorMessage.str( ) );
    }

    // do the transpose
    for( int i = 0; i < row; i++ )
    {
        for( int j = 0; j < column; j++ )
        {
            result[ i ][ j ] = matrix[ j ][ i ];
        }
    }
 }

} // namespace naos

#endif // BASIC_MATH_HPP
