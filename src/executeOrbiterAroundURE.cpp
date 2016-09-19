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
#include "NAOS/basicAstro.hpp"

namespace naos
{

//! Execute routine to compute orbiter/particle motion around a uniformly rotating ellipsoid (URE)
/*!
 * This routine will compute the motion of the orbiter or particle around an asteroid modelled as
 * a uniformly rotating triaxial ellipsoid. The routine is called in main for execution.
 * @sa ellipsoidGravitationalAcceleration.cpp, orbiterEquationsOfMotion.hpp, rk4.hpp
 *
 */
 const void executeOrbiterAroundURE( const double alpha,
                                     const double beta,
                                     const double gamma,
                                     const double gravParameter,
                                     const double density,
                                     std::vector< double > Wvector,
                                     const double Wmagnitude,
                                     const bool initialVectorIsCartesian,
                                     const std::vector< double > &initialVector,
                                     const double integrationStepSize,
                                     const double startTime,
                                     const double endTime,
                                     std::ostringstream &filePath )
 {
    // Use a numerical integrator to evaluate the equations of motion for an orbiter around a uniformly
    // rotating triaxial ellipsoid (URE). Store the solutions for the entire time of integration in
    // a CSV file.
    std::ofstream eomOrbiterUREFile;
    eomOrbiterUREFile.open( filePath.str( ) );
    eomOrbiterUREFile << "x" << "," << "y" << "," << "z" << ",";
    eomOrbiterUREFile << "vx" << "," << "vy" << "," << "vz" << "," << "t" << std::endl;

    // get the initial state from the input arguments
    Vector6 initialStateVector( 6 );
    if( initialVectorIsCartesian == false )
    {
        Vector6 keplerElements( 6 );
        for( int i = 0; i < 6; i++ )
        {
            keplerElements[ i ] = initialVector[ i ];
        }

        // Specify initial values for the state vector of the orbiter in SI units.
        initialStateVector = convertKeplerianElementsToCartesianCoordinates< Vector6 >(
                                keplerElements,
                                gravParameter );

        // get velocity in the body frame
        Vector3 positionVector { initialStateVector[ xPositionIndex ],
                                 initialStateVector[ yPositionIndex ],
                                 initialStateVector[ zPositionIndex ] };

        Vector3 coriolisTerm = crossProduct< Vector3 >( Wvector, positionVector );

        initialStateVector[ xVelocityIndex ] = initialStateVector[ xVelocityIndex ]
                                                - coriolisTerm[ 0 ];
        initialStateVector[ yVelocityIndex ] = initialStateVector[ yVelocityIndex ]
                                                - coriolisTerm[ 1 ];
        initialStateVector[ zVelocityIndex ] = initialStateVector[ zVelocityIndex ]
                                                - coriolisTerm[ 2 ];
    }
    else
    {
        initialStateVector = initialVector;
    }

    // Specify the step size value [s]
    double stepSize = integrationStepSize;

    // Specify the integration time limits [s]
    const double tStart = startTime;
    const double tEnd = endTime;
    double tCurrent = tStart;

    // Save initial values in the CSV file
    eomOrbiterUREFile << initialStateVector[ xPositionIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ yPositionIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ zPositionIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ xVelocityIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ yVelocityIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ zVelocityIndex ] << ",";
    eomOrbiterUREFile << tCurrent << std::endl;

    // Define a vector to store latest/current state values in the integration loop
    Vector6 currentStateVector = initialStateVector;

    // Define a vector to store the integrated state vector values
    Vector6 nextStateVector = initialStateVector;

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
        Vector3 currentGravAcceleration( 3 );
        computeNonDimensionalEllipsoidGravitationalAcceleration(
            alpha,
            beta,
            gamma,
            density,
            Wmagnitude,
            currentStateVector[ xPositionIndex ],
            currentStateVector[ yPositionIndex ],
            currentStateVector[ zPositionIndex ],
            currentGravAcceleration );

        // Create an object of the struct containing the equations of motion (in this case the eom
        // for an orbiter around a uniformly rotating tri-axial ellipsoid) and initialize it to the
        // values of the ellipsoidal asteroid's current gravitational accelerations
        eomOrbiterURE derivatives( currentGravAcceleration );

        // Run an instance of RK4 integrator to evaluate the state at the next time value
        // The input state vector should be normalized
        Vector6 currentStateVectorNormalized( 6 );
        for( int i = 0; i < 3; i++ )
        {
            currentStateVectorNormalized[ i ] = currentStateVector[ i ] / alpha;
            currentStateVectorNormalized[ i + 3 ] = currentStateVector[ i + 3 ]
                                                        / ( alpha * Wmagnitude );
        }
        Vector6 nextStateVectorNormalized( 6 );
        rk4< Vector6, eomOrbiterURE > ( currentStateVectorNormalized,
                                                          tCurrent,
                                                          stepSize,
                                                          nextStateVectorNormalized,
                                                          derivatives );
        // Dimensionalize the integrated state
        for( int i = 0; i < 3; i++ )
        {
            nextStateVector[ i ] = nextStateVectorNormalized[ i ] * alpha;
            nextStateVector[ i + 3 ] = nextStateVectorNormalized[ i + 3 ]
                                                    * ( alpha * Wmagnitude );
        }

        // Check if the new state is valid or not i.e. the particle/orbiter is not
        // inside the surface of the asteroid. If it is then terminate the outer loop. Evaluate the
        // function phi(x,y,z;0) and check if it has a positive sign (i.e. point is outside)
        double xCoordinateSquare = nextStateVector[ xPositionIndex ]
                                    * nextStateVector[ xPositionIndex ];
        double yCoordinateSquare = nextStateVector[ yPositionIndex ]
                                    * nextStateVector[ yPositionIndex ];
        double zCoordinateSquare = nextStateVector[ zPositionIndex ]
                                    * nextStateVector[ zPositionIndex ];
        double phiCheck = xCoordinateSquare / ( alpha * alpha )
                            + yCoordinateSquare / ( beta * beta )
                            + zCoordinateSquare / ( gamma * gamma )
                            - 1.0;
        if( phiCheck <= 0.0 )
        {
            // point is either inside or on the surface of the ellipsoid, so terminate the
            // simulator and exit the outer loop after saving the data.
            eomOrbiterUREFile << nextStateVector[ xPositionIndex ] << ",";
            eomOrbiterUREFile << nextStateVector[ yPositionIndex ] << ",";
            eomOrbiterUREFile << nextStateVector[ zPositionIndex ] << ",";
            eomOrbiterUREFile << nextStateVector[ xVelocityIndex ] << ",";
            eomOrbiterUREFile << nextStateVector[ yVelocityIndex ] << ",";
            eomOrbiterUREFile << nextStateVector[ zVelocityIndex ] << ",";
            eomOrbiterUREFile << tNext << std::endl;
            std::cout << std::endl;
            std::cout << "Houston, we've got a problem!" << std::endl;
            std::cout << "Particle at or inside the ellipsoidal surface" << std::endl;
            std::cout << "Event occurence at: " << tNext << " seconds" << std::endl;
            break;
        }

        // Store the integrated state values in the CSV file
        eomOrbiterUREFile << nextStateVector[ xPositionIndex ] << ",";
        eomOrbiterUREFile << nextStateVector[ yPositionIndex ] << ",";
        eomOrbiterUREFile << nextStateVector[ zPositionIndex ] << ",";
        eomOrbiterUREFile << nextStateVector[ xVelocityIndex ] << ",";
        eomOrbiterUREFile << nextStateVector[ yVelocityIndex ] << ",";
        eomOrbiterUREFile << nextStateVector[ zVelocityIndex ] << ",";
        eomOrbiterUREFile << tNext << std::endl;

        // save the new state values in the vector of current state values. these will be used in
        // the next loop iteration
        tCurrent = tNext;
        currentStateVector = nextStateVector;
    }
 }

} // namespace naos
