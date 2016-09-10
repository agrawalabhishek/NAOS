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
#include "NAOS/orbiterEquationsOfMotion.hpp"
#include "NAOS/rk4.hpp"
#include "NAOS/basicMath.hpp"

int main( const int numberOfInputs, const char* inputArguments[ ] )
{
    // Open file to store points data for the surface of an ellipsoid.
    std::ostringstream ellipsoidSurfacePointsFile;
    ellipsoidSurfacePointsFile << "../../data/ellipsoidSurfacePoints.csv";

    // Physical parameters for Asteroid Eros, all in SI units, modelled as an ellipsoid.
    const double alpha = 20.0 * 1.0e3;
    const double beta = 7.0 * 1.0e3;
    const double gamma = 7.0 * 1.0e3;
    const double density = 3.2 * ( 10.0e-3 ) / ( 10.0e-6 );
    const double mass = ( 4.0 * naos::PI / 3.0 ) * density * alpha * beta * gamma;
    const double gravitationalParameter = naos::GRAVITATIONAL_CONSTANT * mass;
    const double Wx = 0.0; // rotational rate around principal x axis [rad/s]
    const double Wy = 0.0; // rotational rate around principal y axis [rad/s]
    const double Wz = 0.00033118202125129593; // rotational rate around principal z axis [rad/s]
    naos::Vector3 W { Wx, Wy, Wz };

    // Generate surface coordinates for the Asteroid
    const double stepSizeAzimuthDegree = 10.0;
    const double stepSizeElevationDegree = 10.0;
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

    // Use RK4 integrator to evaluate the equations of motion for an orbiter around a uniformly
    // rotating triaxial ellipsoid (URE). Store the solutions for the entire time of integration in
    // a CSV file.
    std::ofstream eomOrbiterUREFile;
    eomOrbiterUREFile.open( "../../data/eomOrbiterURESolution.csv" );
    eomOrbiterUREFile << "x" << "," << "y" << "," << "z" << ",";
    eomOrbiterUREFile << "vx" << "," << "vy" << "," << "vz" << "," << "t" << std::endl;

    // Specify initial values for the state vector of the orbiter in SI units.
    naos::Vector6 initialStateVector( 6 );

    // The positions given below correspond to a point on the surface of asteroid Eros.
    initialStateVector[ naos::xPositionIndex ] = 13268.2789633788;
    initialStateVector[ naos::yPositionIndex ] = 2681.15555091642;
    initialStateVector[ naos::zPositionIndex ] = 4499.51326780578;

    // Linear velocity of the surface point due to rotation of the asteroid
    naos::Vector3 particlePosition { initialStateVector[ naos::xPositionIndex ],
                                     initialStateVector[ naos::yPositionIndex ],
                                     initialStateVector[ naos::zPositionIndex ] };
    naos::Vector3 particleSurfaceVelocity( 3 );
    particleSurfaceVelocity = naos::crossProduct< naos::Vector3 > ( particlePosition, W );

    // Get the initial velocity of the particle.
    // initialStateVector[ naos::xVelocityIndex ] = 20 * particleSurfaceVelocity[ 0 ];
    // initialStateVector[ naos::yVelocityIndex ] = 20 * particleSurfaceVelocity[ 1 ];
    // initialStateVector[ naos::zVelocityIndex ] = 20 * particleSurfaceVelocity[ 2 ];
    initialStateVector[ naos::xVelocityIndex ] = 6.0;
    initialStateVector[ naos::yVelocityIndex ] = 10.0;
    initialStateVector[ naos::zVelocityIndex ] = 6.0;

    // Specify the step size value [s]
    double stepSize = 0.01;

    // Specify the integration time limits [s]
    const double tStart = 0.0;
    const double tEnd = 100000.0;
    double tCurrent = tStart;

    // Save initial values in the CSV file
    eomOrbiterUREFile << initialStateVector[ naos::xPositionIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ naos::yPositionIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ naos::zPositionIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ naos::xVelocityIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ naos::yVelocityIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ naos::zVelocityIndex ] << ",";
    eomOrbiterUREFile << tCurrent << std::endl;

    // Define a vector to store latest/current state values in the integration loop
    naos::Vector6 currentStateVector = initialStateVector;

    // Define a vector to store the integrated state vector values
    naos::Vector6 nextStateVector = initialStateVector;

    // Start the integration outer loop
    while( tCurrent != tEnd )
    {
        // calculate the new time value
        double tNext = tCurrent + stepSize;

        // check if the next time value is greater than the final time value
        if( tNext > tEnd )
        {
            // make the next time value equal to the final time value and recalculate the step size
            tNext = tEnd;
            stepSize = tEnd - tCurrent;
        }

        // Evaluate gravitational acceleration values at the current state values
        naos::Vector3 currentGravAcceleration( 3 );
        naos::computeEllipsoidGravitationalAcceleration( alpha, beta, gamma,
                                                         gravitationalParameter,
                                                         currentStateVector[ naos::xPositionIndex ],
                                                         currentStateVector[ naos::yPositionIndex ],
                                                         currentStateVector[ naos::zPositionIndex ],
                                                         currentGravAcceleration );

        // Create an object of the struct containing the equations of motion (in this case the eom
        // for an orbiter around a uniformly rotating tri-axial ellipsoid) and initialize it to the
        // values of the constant angular rate of the ellipsoidal asteroid and the current
        // gravitational accelerations
        naos::eomOrbiterURE derivatives( Wz, currentGravAcceleration );

        // Run an instance of RK4 integrator to evaluate the state at the next time value
        naos::rk4< naos::Vector6, naos::eomOrbiterURE > ( currentStateVector,
                                                          tCurrent,
                                                          stepSize,
                                                          nextStateVector,
                                                          derivatives );

        // Check if the new state is valid or not i.e. the particle/orbiter is not
        // inside the surface of the asteroid. If it is then terminate the outer loop. Evaluate the
        // function phi(x,y,z;0) and check if it has a positive sign (i.e. point is outside)
        double xCoordinateSquare = nextStateVector[ naos::xPositionIndex ]
                                    * nextStateVector[ naos::xPositionIndex ];
        double yCoordinateSquare = nextStateVector[ naos::yPositionIndex ]
                                    * nextStateVector[ naos::yPositionIndex ];
        double zCoordinateSquare = nextStateVector[ naos::zPositionIndex ]
                                    * nextStateVector[ naos::zPositionIndex ];
        double phiCheck = xCoordinateSquare / ( alpha * alpha )
                            + yCoordinateSquare / ( beta * beta )
                            + zCoordinateSquare / ( gamma * gamma )
                            - 1.0;
        if( phiCheck <= 0.0 )
        {
            // point is either inside or on the surface of the ellipsoid, so terminate the
            // simulator and exit the outer loop after saving the data.
            eomOrbiterUREFile << nextStateVector[ naos::xPositionIndex ] << ",";
            eomOrbiterUREFile << nextStateVector[ naos::yPositionIndex ] << ",";
            eomOrbiterUREFile << nextStateVector[ naos::zPositionIndex ] << ",";
            eomOrbiterUREFile << nextStateVector[ naos::xVelocityIndex ] << ",";
            eomOrbiterUREFile << nextStateVector[ naos::yVelocityIndex ] << ",";
            eomOrbiterUREFile << nextStateVector[ naos::zVelocityIndex ] << ",";
            eomOrbiterUREFile << tNext << std::endl;
            std::cout << std::endl;
            std::cout << "Houston, we've got a problem!" << std::endl;
            std::cout << "Particle at or inside the ellipsoidal surface" << std::endl;
            std::cout << "Event occurence at: " << tNext << " seconds" << std::endl;
            break;
        }

        // Store the integrated state values in the CSV file
        eomOrbiterUREFile << nextStateVector[ naos::xPositionIndex ] << ",";
        eomOrbiterUREFile << nextStateVector[ naos::yPositionIndex ] << ",";
        eomOrbiterUREFile << nextStateVector[ naos::zPositionIndex ] << ",";
        eomOrbiterUREFile << nextStateVector[ naos::xVelocityIndex ] << ",";
        eomOrbiterUREFile << nextStateVector[ naos::yVelocityIndex ] << ",";
        eomOrbiterUREFile << nextStateVector[ naos::zVelocityIndex ] << ",";
        eomOrbiterUREFile << tNext << std::endl;

        // save the new state values in the vector of current state values. these will be used in
        // the next loop iteration
        tCurrent = tNext;
        currentStateVector = nextStateVector;
    }

    return EXIT_SUCCESS;
}
