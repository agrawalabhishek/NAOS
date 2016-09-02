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
#include <cmath>

#include "NAOS/ellipsoidSurfacePoints.hpp"
#include "NAOS/cubicRoot.hpp"
#include "NAOS/constants.hpp"
#include "NAOS/ellipsoidGravitationalAcceleration.hpp"

int main( const int numberOfInputs, const char* inputArguments[ ] )
{
    // Open file to store points data for the surface of an ellipsoid.
    std::ostringstream ellipsoidSurfacePointsFile;
    ellipsoidSurfacePointsFile << "../../data/ellipsoidSurfacePoints.csv";

    // Physical parameters for Asteroid Eros, all in SI units.
    const double alpha = 20.0 * 1.0e3;
    const double beta = 7.0 * 1.0e3;
    const double gamma = 7.0 * 1.0e3;
    const double density = 3.2 * ( 10.0e-3 ) / ( 10.0e-6 );
    const double mass = ( 4.0 * naos::PI / 3.0 ) * density * alpha * beta * gamma;
    const double gravitationalParameter = naos::GRAVITATIONAL_CONSTANT * mass;

    // Generate surface coordinates for the Asteroid
    const double stepSizeAzimuthDegree = 1.0;
    const double stepSizeElevationDegree = 1.0;
    naos::computeEllipsoidSurfacePoints( alpha, beta, gamma,
                                         stepSizeAzimuthDegree, stepSizeElevationDegree,
                                         ellipsoidSurfacePointsFile );

    // Test functionality of the maximum real root for a cubic polynomial, compare results with
    // wolfram alpha.
    const double maxRealRootTest1 = naos::computeMaxRealCubicRoot( -7.0, 4.0, 12.0 );
    std::cout << "max Real Root Test 1 = " << maxRealRootTest1 << std::endl;
    const double maxRealRootTest2 = naos::computeMaxRealCubicRoot( ( 3.0 / 2.0 ),
                                                                   ( -11.0 / 2.0 ),
                                                                   ( -3.0 ) );
    std::cout << "max Real Root Test 2 = " << maxRealRootTest2 << std::endl;

    // Test evaluation of gravitational acceleration using previously computed surface points.
    std::ifstream ellipsoidSurfacePoints;
    ellipsoidSurfacePoints.open( "../../data/ellipsoidSurfacePoints.csv" );

    // Extract the column headers first from the surface points file
    std::string headers;
    std::getline( ellipsoidSurfacePoints, headers );
    std::cout << headers << std::endl;

    // Extract numeric data, calculate the gravitational acceleration at all those points and save
    // it in a CSV file
    double xCoordinate = 0.0;
    double yCoordinate = 0.0;
    double zCoordinate = 0.0;
    double latitude = 0.0;
    double longitude = 0.0;
    double range = 0.0;
    std::string numberString;
    naos::Vector3 gravitationalAcceleration { 0.0, 0.0, 0.0 };

    std::ofstream ellipsoidSurfaceAccelerationFile;
    ellipsoidSurfaceAccelerationFile.open( "../../data/ellipsoidSurfaceAcceleration.csv" );
    ellipsoidSurfaceAccelerationFile << "Ux" << "," << "Uy" << "," << "Uz" << ",";
    ellipsoidSurfaceAccelerationFile << "U" << "," << "latitude" << "," << "longitude" << std::endl;

    while( !ellipsoidSurfacePoints.eof( ) )
    {
        std::getline( ellipsoidSurfacePoints, numberString, ',' );
        xCoordinate = std::stod( numberString );

        std::getline( ellipsoidSurfacePoints, numberString, ',' );
        yCoordinate = std::stod( numberString );

        std::getline( ellipsoidSurfacePoints, numberString, ',' );
        zCoordinate = std::stod( numberString );

        std::getline( ellipsoidSurfacePoints, numberString, ',' );
        latitude = std::stod( numberString );

        std::getline( ellipsoidSurfacePoints, numberString, ',' );
        longitude = std::stod( numberString );

        std::getline( ellipsoidSurfacePoints, numberString );
        range = std::stod( numberString );

        naos::computeEllipsoidGravitationalAcceleration( alpha, beta, gamma,
                                                         gravitationalParameter,
                                                         xCoordinate,
                                                         yCoordinate,
                                                         zCoordinate,
                                                         gravitationalAcceleration );

        ellipsoidSurfaceAccelerationFile << gravitationalAcceleration[ 0 ] << ",";
        ellipsoidSurfaceAccelerationFile << gravitationalAcceleration[ 1 ] << ",";
        ellipsoidSurfaceAccelerationFile << gravitationalAcceleration[ 2 ] << ",";

        double accelerationMagnitude =
                    std::sqrt( gravitationalAcceleration[ 0 ] * gravitationalAcceleration[ 0 ]
                             + gravitationalAcceleration[ 1 ] * gravitationalAcceleration[ 1 ]
                             + gravitationalAcceleration[ 2 ] * gravitationalAcceleration[ 2 ] );

        ellipsoidSurfaceAccelerationFile << accelerationMagnitude << ",";
        ellipsoidSurfaceAccelerationFile << latitude << "," << longitude << std::endl;
    }
    ellipsoidSurfacePoints.close( );
    ellipsoidSurfaceAccelerationFile.close( );

    return EXIT_SUCCESS;
}
