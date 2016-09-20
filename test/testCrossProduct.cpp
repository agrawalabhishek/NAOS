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

TEST_CASE( "Test cross product function", "[cross-product]" )
{
    // input values were individually tested on wolfram alpha
    const std::vector< double > alpha = { 1.0, 2.0, 3.0 };
    const std::vector< double > beta = { 3.0, 4.0, 5.0 };

    // comute the cross product
    std::vector< double > result = crossProduct( alpha, beta );

    // check the computed values
    REQUIRE( result[ 0 ] == Approx( -2.0 ).epsilon( 1.0e-9 ) );
    REQUIRE( result[ 1 ] == Approx( 4.0 ).epsilon( 1.0e-9 ) );
    REQUIRE( result[ 2 ] == Approx( -2.0 ).epsilon( 1.0e-9 ) );
}

} // namespace tests
} // namespace naos
