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
#include "NAOS/basicMath.hpp"

namespace naos
{

//! Compute gravitational Acceleration on the surface of the ellipsoid shaped asteroid
/*!
 * Compute absolute value of the gravitational acceleration on the surface of the ellipsoid and
 * store it in a CSV file.
 * @sa ellipsoidGravitationalAcceleration.cpp
 * @sa cubicRoot.cpp
 *
 * @param[in] alpha, beta, gamma        semi-axes of the ellipsoid such that alpha>beta>gamma
 * @param[in] gravitationalParameter    gravitational parameter of the ellipsoid
 */
const void computeEllipsoidSurfaceGravitationalAcceleration( const double alpha,
                                                             const double beta,
                                                             const double gamma,
                                                             const double gravitationalParameter )
{
    // Test evaluation of gravitational acceleration using previously computed surface points.
    std::ifstream ellipsoidSurfacePoints;
    ellipsoidSurfacePoints.open( "../../data/ellipsoidSurfacePoints.csv" );

    // Extract the column headers first from the surface points file
    std::string headers;
    std::getline( ellipsoidSurfacePoints, headers );

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

    while( std::getline( ellipsoidSurfacePoints, numberString, ',' ) )
    {
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
}

//! Compute normalized gravitational Acceleration on the surface of the ellipsoid shaped asteroid
/*!
 * Compute normalized value of the gravitational acceleration on the surface of the ellipsoid and
 * store it in a CSV file.
 * @sa ellipsoidGravitationalAcceleration.cpp
 * @sa cubicRoot.cpp
 *
 * @param[in] alpha, beta, gamma        semi-axes of the ellipsoid such that alpha>beta>gamma
 * @param[in] density                   density of the ellipsoid
 * @param[in] Wmagnitude                magnitude of the angular rotational velocity of the asteroid
 */
const void computeNonDimensionalEllipsoidSurfaceGravitationalAcceleration( const double alpha,
                                                                           const double beta,
                                                                           const double gamma,
                                                                           const double density,
                                                                           const double Wmagnitude )
{
    // Test evaluation of gravitational acceleration using previously computed surface points.
    std::ifstream ellipsoidSurfacePoints;
    ellipsoidSurfacePoints.open( "../../data/ellipsoidSurfacePoints.csv" );

    // Extract the column headers first from the surface points file
    std::string headers;
    std::getline( ellipsoidSurfacePoints, headers );

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
    std::stringstream filepath;
    filepath << "../../data/nondimensionalellipsoidSurfaceAcceleration.csv";
    ellipsoidSurfaceAccelerationFile.open( filepath.str( ) );
    ellipsoidSurfaceAccelerationFile << "Ux" << "," << "Uy" << "," << "Uz" << ",";
    ellipsoidSurfaceAccelerationFile << "U" << "," << "latitude" << "," << "longitude" << std::endl;

    while( std::getline( ellipsoidSurfacePoints, numberString, ',' ) )
    {
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

        naos::computeNonDimensionalEllipsoidGravitationalAcceleration( alpha, beta, gamma,
                                                                       density,
                                                                       Wmagnitude,
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
}

} //namespace naos
