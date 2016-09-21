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

TEST_CASE( "Test vector norm function", "[vector-norm]" )
{
    // input values were individually tested on wolfram alpha
    const std::vector< double > alpha = { 1.0, 2.0, 3.0 };

    // compute the norm
    const double result = vectorNorm( alpha );

    // check the computed values
    REQUIRE( result == Approx( 3.741657386773941 ).epsilon( 1.0e-9 ) );
}

} // namespace tests
} // namespace naos
