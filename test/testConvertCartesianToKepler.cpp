/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <catch.hpp>

#include <vector>

#include "NAOS/constants.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/basicMath.hpp"

namespace naos
{
namespace tests
{

TEST_CASE( "Test funtion to convert cartesian elements to kepler elements", "[cartesian-kepler]" )
{
    // input test values taken from:
    // github.com/openAstro/Astro/test/testOrbitalElementConversions.cpp
    Vector6 cartesianElements( 6 );
    cartesianElements[ 0 ] = 3.75e6;
    cartesianElements[ 1 ] = 4.24e6;
    cartesianElements[ 2 ] = -1.39e6;
    cartesianElements[ 3 ] = -4.65e3;
    cartesianElements[ 4 ] = -2.21e3;
    cartesianElements[ 5 ] = 1.66e3;

    Vector6 keplerianElements( 6 );
    keplerianElements[ 0 ] = 3.707478199246163e6;
    keplerianElements[ 1 ] = 0.949175203660321;
    keplerianElements[ 2 ] = convertRadiansToDegree( 0.334622356632438 );
    keplerianElements[ 3 ] = convertRadiansToDegree( 1.630852596545341 );
    keplerianElements[ 4 ] = convertRadiansToDegree( 2.168430616511167 );
    keplerianElements[ 5 ] = convertRadiansToDegree( 3.302032232567084 );

    const double gravParameter = 3.986004415e14;

    // comute the equivalent cartesian coordinates
    Vector6 testKeplerElements( 6 );
    convertCartesianCoordinatesToKeplerianElements< Vector6 >( cartesianElements,
                                                               gravParameter,
                                                               testKeplerElements );
    // check the computed values
    REQUIRE( testKeplerElements[ 0 ] == Approx( keplerianElements[ 0 ] ).epsilon( 1.0e-10 ) );
    REQUIRE( testKeplerElements[ 1 ] == Approx( keplerianElements[ 1 ] ).epsilon( 1.0e-10 ) );
    REQUIRE( testKeplerElements[ 2 ] == Approx( keplerianElements[ 2 ] ).epsilon( 1.0e-10 ) );
    REQUIRE( testKeplerElements[ 3 ] == Approx( keplerianElements[ 3 ] ).epsilon( 1.0e-10 ) );
    REQUIRE( testKeplerElements[ 4 ] == Approx( keplerianElements[ 4 ] ).epsilon( 1.0e-10 ) );
    REQUIRE( testKeplerElements[ 5 ] == Approx( keplerianElements[ 5 ] ).epsilon( 1.0e-10 ) );
}

} // namespace tests
} // namespace naos
