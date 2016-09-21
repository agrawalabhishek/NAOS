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

TEST_CASE( "Test dot product function", "[dot-product]" )
{
    // input values were individually tested on wolfram alpha
    const std::vector< double > alpha = { 1.0, 2.0, 3.0 };
    const std::vector< double > beta = { 3.0, 4.0, 5.0 };

    // compute the dot product
    const double result = dotProduct( alpha, beta );

    // check the computed values
    REQUIRE( result == Approx( 26.0 ).epsilon( 1.0e-9 ) );
}

} // namespace tests
} // namespace naos
