/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <catch.hpp>

#include <vector>

#include "NAOS/basicMath.hpp"

namespace naos
{
namespace tests
{

TEST_CASE( "Test Matrix Transpose", "[matrix-transpose]" )
{
    SECTION( "Input matrix is a 3x3 matrix" )
    {
        std::vector< std::vector< double > > inputMatrix { { 5.0, 6.0, 8.0 },
                                                           { 8.0, 4.0, 9.0 },
                                                           { 10.0, 63.0, 2.0 } };

        std::vector< std::vector< double > > outputMatrix( 3, std::vector< double >( 3 ) );

        naos::matrixTranspose( inputMatrix,
                               outputMatrix );

        // check the computed values
        REQUIRE( outputMatrix[ 0 ][ 0 ] == Approx( 5.0 ).epsilon( 1.0e-15 ) );
        REQUIRE( outputMatrix[ 0 ][ 1 ] == Approx( 8.0 ).epsilon( 1.0e-15 ) );
        REQUIRE( outputMatrix[ 0 ][ 2 ] == Approx( 10.0 ).epsilon( 1.0e-15 ) );

        REQUIRE( outputMatrix[ 1 ][ 0 ] == Approx( 6.0 ).epsilon( 1.0e-15 ) );
        REQUIRE( outputMatrix[ 1 ][ 1 ] == Approx( 4.0 ).epsilon( 1.0e-15 ) );
        REQUIRE( outputMatrix[ 1 ][ 2 ] == Approx( 63.0 ).epsilon( 1.0e-15 ) );

        REQUIRE( outputMatrix[ 2 ][ 0 ] == Approx( 8.0 ).epsilon( 1.0e-15 ) );
        REQUIRE( outputMatrix[ 2 ][ 1 ] == Approx( 9.0 ).epsilon( 1.0e-15 ) );
        REQUIRE( outputMatrix[ 2 ][ 2 ] == Approx( 2.0 ).epsilon( 1.0e-15 ) );
    }
}

} // namespace tests
} // namespace naos
