/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <catch.hpp>

#include <vector>

#include "NAOS/constants.hpp"
#include "NAOS/basicAstro.hpp"

namespace naos
{
namespace tests
{

TEST_CASE( "Test funtion to convert kepler elements to cartesian elements", "[kepler-cartesian]" )
{
    // input test values taken from MGOD course (Noomen 2014, TUDelft), example 1
    Vector6 keplerElements( 6 );
    keplerElements[ 0 ] = 6787746.891;
    keplerElements[ 1 ] = 0.000731104;
    keplerElements[ 2 ] = 51.68714486;
    keplerElements[ 3 ] = 127.5486706;
    keplerElements[ 4 ] = 74.21987137;
    keplerElements[ 5 ] = 24.10027677;
    const double gravParameter = 3.98600441e14;

    // comute the equivalent cartesian coordinates
    Vector6 cartesianElements = convertKeplerianElementsToCartesianCoordinates< Vector6 >(
                                                keplerElements,
                                                gravParameter );
    // check the computed values
    REQUIRE( cartesianElements[ xPositionIndex ] == Approx( -2700816.14 ).epsilon( 1.0e-9 ) );
    REQUIRE( cartesianElements[ yPositionIndex ] == Approx( -3314092.80 ).epsilon( 1.0e-9 ) );
    REQUIRE( cartesianElements[ zPositionIndex ] == Approx( 5266346.42 ).epsilon( 1.0e-9 ) );
    REQUIRE( cartesianElements[ xVelocityIndex ] == Approx( 5168.606550 ).epsilon( 1.0e-9 ) );
    REQUIRE( cartesianElements[ yVelocityIndex ] == Approx( -5597.546618 ).epsilon( 1.0e-9 ) );
    REQUIRE( cartesianElements[ zVelocityIndex ] == Approx( -868.878445 ).epsilon( 1.0e-9 ) );
}

} // namespace tests
} // namespace naos
