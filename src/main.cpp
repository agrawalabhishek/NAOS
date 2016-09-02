/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include "NAOS/ellipsoidSurfacePoints.hpp"
#include "NAOS/cubicRoot.hpp"
#include "NAOS/constants.hpp"
#include "NAOS/ellipsoidGravitationalAcceleration.hpp"

int main( const int numberOfInputs, const char* inputArguments[ ] )
{
    // Generate points for the surface of an ellipsoid.
    std::ostringstream ellipsoidSurfacePointsFile;
    ellipsoidSurfacePointsFile << "../../data/ellipsoidSurfacePoints.csv";
    naos::computeEllipsoidSurfacePoints( 12.0, 8.0, 5.0, 10.0, 10.0 , ellipsoidSurfacePointsFile );

    // Test functionality of the maximum real root for a cubic polynomial, compare results with
    // wolfram alpha.
    const double maxRealRootTest1 = naos::computeMaxRealCubicRoot( -7.0, 4.0, 12.0 );
    std::cout << "max Real Root Test 1 = " << maxRealRootTest1 << std::endl;
    const double maxRealRootTest2 = naos::computeMaxRealCubicRoot( ( 3.0 / 2.0 ),
                                                                   ( -11.0 / 2.0 ),
                                                                   ( -3.0 ) );
    std::cout << "max Real Root Test 2 = " << maxRealRootTest2 << std::endl;

    // Test evaluate gravitational acceleration using previously computed surface points.
    std::ifstream ellipsoidSurfacePoints;
    ellipsoidSurfacePoints.open( "../../data/ellipsoidSurfacePoints.csv" );
    std::string headers;
    std::getline( ellipsoidSurfacePoints, headers );
    std::cout << headers << std::endl;
    double xCoordinate = 0.0;
    std::string numberString;
    std::getline( ellipsoidSurfacePoints, numberString, ',' );
    xCoordinate = std::stod( numberString );
    std::cout << xCoordinate << std::endl;

    return EXIT_SUCCESS;
}
