/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <catch.hpp>

#include "NAOS/cubicRoot.hpp"

namespace naos
{
namespace tests
{

TEST_CASE( "Test computation of max real cubic root", "[cubic-root]" )
{
    // input values were individually tested on wolfram alpha
    const double A2 = -7.0;
    const double A1 = 4.0;
    const double A0 = 12.0;

    const double B2 = 3.0 / 2.0;
    const double B1 = -11.0 / 2.0;
    const double B0 = -3.0;

    // compute the max real roots for the two cubic equations
    const double maxRealRoot1 = computeMaxRealCubicRoot( A2, A1, A0 );
    const double maxRealRoot2 = computeMaxRealCubicRoot( B2, B1, B0 );

    // check the computed values
    REQUIRE( maxRealRoot1 == Approx( 6.0 ).epsilon( 1.0e-9 ) );
    REQUIRE( maxRealRoot2 == Approx( 2.0 ).epsilon( 1.0e-9 ) );
}

} // namespace tests
} // namespace naos
