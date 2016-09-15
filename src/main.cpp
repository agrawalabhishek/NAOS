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
#include "NAOS/computeEllipsoidSurfaceGravitationalAcceleration.hpp"

int main( const int numberOfInputs, const char* inputArguments[ ] )
{
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
    const double Wmagnitude = std::sqrt( Wx * Wx + Wy * Wy + Wz * Wz );

    // Generate surface coordinates for the Asteroid
    std::ostringstream ellipsoidSurfacePointsFile;
    ellipsoidSurfacePointsFile << "../../data/ellipsoidSurfacePoints.csv";
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

    // Compute gravitational acceleration on the surface of the ellipsoidal asteroid
    naos::computeEllipsoidSurfaceGravitationalAcceleration( alpha, beta, gamma,
                                                            gravitationalParameter );

    // Compute normalized gravitational acceleration values for the surface of the asteroid
    naos::computeNonDimensionalEllipsoidSurfaceGravitationalAcceleration( alpha, beta, gamma,
                                                                          density, Wmagnitude );

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
    // initialStateVector[ naos::xVelocityIndex ] = 50.0 * particleSurfaceVelocity[ 0 ];
    // initialStateVector[ naos::yVelocityIndex ] = 50.0 * particleSurfaceVelocity[ 1 ];
    // initialStateVector[ naos::zVelocityIndex ] = 50.0 * particleSurfaceVelocity[ 2 ];
    initialStateVector[ naos::xVelocityIndex ] = 0.1;
    initialStateVector[ naos::yVelocityIndex ] = 0.1;
    initialStateVector[ naos::zVelocityIndex ] = 0.0;

    // Specify the step size value [s]
    double stepSize = 0.01;

    // Specify the integration time limits [s]
    const double tStart = 0.0;
    const double tEnd = 1000.0;
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
        naos::computeNonDimensionalEllipsoidGravitationalAcceleration(
            alpha,
            beta,
            gamma,
            density,
            Wmagnitude,
            currentStateVector[ naos::xPositionIndex ],
            currentStateVector[ naos::yPositionIndex ],
            currentStateVector[ naos::zPositionIndex ],
            currentGravAcceleration );

        // Create an object of the struct containing the equations of motion (in this case the eom
        // for an orbiter around a uniformly rotating tri-axial ellipsoid) and initialize it to the
        // values of the ellipsoidal asteroid's current gravitational accelerations
        naos::eomOrbiterURE derivatives( currentGravAcceleration );

        // Run an instance of RK4 integrator to evaluate the state at the next time value
        // The input state vector should be normalized
        naos::Vector6 currentStateVectorNormalized( 6 );
        for( int i = 0; i < 3; i++ )
        {
            currentStateVectorNormalized[ i ] = currentStateVector[ i ] / alpha;
            currentStateVectorNormalized[ i + 3 ] = currentStateVector[ i + 3 ]
                                                        / ( alpha * Wmagnitude );
        }
        naos::Vector6 nextStateVectorNormalized( 6 );
        naos::rk4< naos::Vector6, naos::eomOrbiterURE > ( currentStateVectorNormalized,
                                                          tCurrent,
                                                          stepSize,
                                                          nextStateVectorNormalized,
                                                          derivatives );
        // Dimensionalize the integrated state
        for( int i = 0; i < 3; i++ )
        {
            nextStateVector[ i ] = nextStateVectorNormalized[ i ] * alpha;
            nextStateVectorNormalized[ i + 3 ] = nextStateVectorNormalized[ i + 3 ]
                                                    * ( alpha * Wmagnitude );
        }

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
