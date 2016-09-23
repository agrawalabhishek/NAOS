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

TEST_CASE( "Test funtion to convert cartesian elements to kepler elements", "[cartesian-kepler]" )
{
    // input test values taken from DA Vallado - fundamentals of astrodynamics book
    Vector6 cartesianElements( 6 );
    cartesianElements[ xPositionIndex ] = 6524.834;
    cartesianElements[ yPositionIndex ] = 6862.875;
    cartesianElements[ zPositionIndex ] = 6448.296;
    cartesianElements[ xVelocityIndex ] = 4.901327;
    cartesianElements[ yVelocityIndex ] = 5.533756;
    cartesianElements[ zVelocityIndex ] = -1.976341;

    Vector6 keplerElements( 6 );
    keplerElements[ 0 ] = 36127.3376;
    keplerElements[ 1 ] = 0.832853;
    keplerElements[ 2 ] = 87.870;
    keplerElements[ 3 ] = 227.89826;
    keplerElements[ 4 ] = 53.384930;
    keplerElements[ 5 ] = 92.335;
    const double gravParameter = 398600.4418;

    // comute the equivalent cartesian coordinates
    Vector6 testKeplerElements( 6 );
    convertCartesianCoordinatesToKeplerianElements< Vector6 >( cartesianElements,
                                                               gravParameter,
                                                               testKeplerElements );
    // check the computed values
    REQUIRE( testKeplerElements[ 0 ] == Approx( keplerElements[ 0 ] ).epsilon( 1.0e-5 ) );
    REQUIRE( testKeplerElements[ 1 ] == Approx( keplerElements[ 1 ] ).epsilon( 1.0e-5 ) );
    REQUIRE( testKeplerElements[ 2 ] == Approx( keplerElements[ 2 ] ).epsilon( 1.0e-5 ) );
    REQUIRE( testKeplerElements[ 3 ] == Approx( keplerElements[ 3 ] ).epsilon( 1.0e-5 ) );
    REQUIRE( testKeplerElements[ 4 ] == Approx( keplerElements[ 4 ] ).epsilon( 1.0e-5 ) );
    REQUIRE( testKeplerElements[ 5 ] == Approx( keplerElements[ 5 ] ).epsilon( 1.0e-5 ) );
}

} // namespace tests
} // namespace naos
