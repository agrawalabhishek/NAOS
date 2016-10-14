/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef EXECUTE_ORBITER_AROUND_URE_POINT_MASS_GRAVITY_HPP
#define EXECUTE_ORBITER_AROUND_URE_POINT_MASS_GRAVITY_HPP

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

#include "NAOS/constants.hpp"
#include "NAOS/orbiterEquationsOfMotion.hpp"
#include "NAOS/rk4.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/misc.hpp"

namespace naos
{

//! Execute routine to compute orbiter/particle motion around a uniformly rotating ellipsoid (URE)
/*!
 * This routine will compute the motion of the orbiter or particle around an asteroid modelled as
 * a uniformly rotating triaxial ellipsoid. The routine is called in main for execution.
 * @sa ellipsoidGravitationalAcceleration.cpp, orbiterEquationsOfMotion.hpp, rk4.hpp
 *
 * @param[in] alpha                         largest semi axis of the ellipsoid
 * @param[in] beta                          intermediate semi axis of the ellipsoid
 * @param[in] gamma                         smallest largest semi axis of the ellipsoid
 * @param[in] gravParameter                 gravitational parameter of the ellipsoid
 * @param[in] density                       density of the ellipsoid
 * @param[in] Wvector                       angular velocity vector for the ellipsoid
 * @param[in] Wmagnitude                    magnitude of the angular velocity
 * @param[in] initialVectorIsCartesian      boolean flag to indicate type of initial vector
 *                                          value is false for kepler elements, and true for
 *                                          cartesian elements.
 * @param[in] initialVector                 vector containing the initial state of the orbiter
 * @param[in] integrationStepSize           step size to be used for integration
 * @param[in] startTime                     integration start time
 * @param[in] endTime                       integration end time
 * @param[in] filePath                      path to csv file where all the results are stored
 */
 template< typename Vector >
 void executeOrbiterAroundUREPointMassGravity( const double alpha,
                                               const double beta,
                                               const double gamma,
                                               const double gravParameter,
                                               const double density,
                                               Vector Wvector,
                                               const double Wmagnitude,
                                               const bool initialVectorIsCartesian,
                                               const Vector &initialVector,
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
    eomOrbiterUREFile << "vx" << "," << "vy" << "," << "vz" << "," << "t" << ",";
    eomOrbiterUREFile << "jacobian" << "," << "semiMajor" << "," << "eccentricity" << ",";
    eomOrbiterUREFile << "inclination" << "," << "RAAN" << "," << "AOP" << "," << "TA";
    eomOrbiterUREFile << std::endl;

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
        // std::cout << "Inertial frame cartesian coordinates (at t=0):";
        // printVector( initialStateVector, 6 );

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
    // std::cout << "Body frame cartesian coordinates (at t=0):";
    // printVector( initialStateVector, 6 );

    // convert the initial state vector in body frame to inertial frame and obtain orbital elements
    double phiAngle = Wmagnitude * startTime;
    phiAngle = mod360( phiAngle );
    phiAngle = convertDegreeToRadians( phiAngle );

    double xInertial = initialStateVector[ xPositionIndex ] * std::cos( phiAngle )
                        - initialStateVector[ yPositionIndex ] * std::sin( phiAngle );
    double yInertial = initialStateVector[ xPositionIndex ] * std::sin( phiAngle )
                        + initialStateVector[ yPositionIndex ] * std::cos( phiAngle );
    double zInertial = initialStateVector[ zPositionIndex ];

    Vector3 positionVector = { initialStateVector[ xPositionIndex ],
                               initialStateVector[ yPositionIndex ],
                               initialStateVector[ zPositionIndex ] };
    Vector3 wCrossR( 3 );
    wCrossR = crossProduct( Wvector, positionVector );
    double vxInertial = initialStateVector[ xVelocityIndex ] + wCrossR[ 0 ];
    double vyInertial = initialStateVector[ yVelocityIndex ] + wCrossR[ 1 ];
    double vzInertial = initialStateVector[ zVelocityIndex ] + wCrossR[ 2 ];

    Vector6 orbitalElements( 6 );
    Vector6 inertialCartesianElements = { xInertial,
                                          yInertial,
                                          zInertial,
                                          vxInertial,
                                          vyInertial,
                                          vzInertial };
    // std::cout << "Inertial cartesian coordinates, obtained from rotating body frame coordinates";
    // printVector( inertialCartesianElements, 6 );
    convertCartesianCoordinatesToKeplerianElements( inertialCartesianElements,
                                                    gravParameter,
                                                    orbitalElements );
    // std::cout << "orbital elements from inertial coordinates at t=0:";
    // printVector( orbitalElements, 6 );

    // Specify the step size value [s]
    static double stepSize = integrationStepSize;

    // Specify the integration time limits [s]
    const double tStart = startTime;
    const double tEnd = endTime;
    double tCurrent = tStart;

    // Calculate the initial jacobian
    double xVelocitySquare = initialStateVector[ xVelocityIndex ]
                                * initialStateVector[ xVelocityIndex ];
    double yVelocitySquare = initialStateVector[ yVelocityIndex ]
                                * initialStateVector[ yVelocityIndex ];
    double zVelocitySquare = initialStateVector[ zVelocityIndex ]
                                * initialStateVector[ zVelocityIndex ];

    double xPositionSquare = initialStateVector[ xPositionIndex ]
                                * initialStateVector[ xPositionIndex ];
    double yPositionSquare = initialStateVector[ yPositionIndex ]
                                * initialStateVector[ yPositionIndex ];
    double zPositionSquare = initialStateVector[ zPositionIndex ]
                                * initialStateVector[ zPositionIndex ];

    double radialDistance = std::sqrt( xPositionSquare + yPositionSquare + zPositionSquare );

    double wSquare = Wmagnitude * Wmagnitude;

    double gravPotential = 0.0;
    gravPotential = gravParameter / radialDistance;

    double jacobian = 0.5 * ( xVelocitySquare + yVelocitySquare + zVelocitySquare )
                        - 0.5 * wSquare * ( xPositionSquare + yPositionSquare )
                        - gravPotential;

    // Save initial values in the CSV file
    eomOrbiterUREFile << initialStateVector[ xPositionIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ yPositionIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ zPositionIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ xVelocityIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ yVelocityIndex ] << ",";
    eomOrbiterUREFile << initialStateVector[ zVelocityIndex ] << ",";
    eomOrbiterUREFile << tCurrent << ",";
    eomOrbiterUREFile << jacobian << ",";
    eomOrbiterUREFile << orbitalElements[ 0 ] << ",";
    eomOrbiterUREFile << orbitalElements[ 1 ] << ",";
    eomOrbiterUREFile << orbitalElements[ 2 ] << ",";
    eomOrbiterUREFile << orbitalElements[ 3 ] << ",";
    eomOrbiterUREFile << orbitalElements[ 4 ] << ",";
    eomOrbiterUREFile << orbitalElements[ 5 ] << std::endl;

    // Define a vector to store latest/current state values in the integration loop
    Vector6 currentStateVector = initialStateVector;

    // Define a vector to store the integrated state vector values
    Vector6 nextStateVector = initialStateVector;

    // set a boolean flag to indicate start of integration
    bool startOfIntegration = true;

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

        // Integrate the EOMs using the rk4 integrator.
        Vector currentGravAcceleration( 3 );
        double radialDistanceCube = radialDistance * radialDistance * radialDistance;

        currentGravAcceleration[ xPositionIndex ]
                = -1.0 * gravParameter * currentStateVector[ xPositionIndex ] / radialDistanceCube;
        currentGravAcceleration[ yPositionIndex ]
                = -1.0 * gravParameter * currentStateVector[ yPositionIndex ] / radialDistanceCube;
        currentGravAcceleration[ zPositionIndex ]
                = -1.0 * gravParameter * currentStateVector[ zPositionIndex ] / radialDistanceCube;

        eomOrbiterURE derivatives( currentGravAcceleration, Wmagnitude );

        rk4< Vector6, eomOrbiterURE >( currentStateVector,
                                       tCurrent,
                                       stepSize,
                                       nextStateVector,
                                       derivatives );

        // convert the next state vector in body frame to inertial frame
        phiAngle = Wmagnitude * tCurrent;
        phiAngle = mod360( phiAngle );
        phiAngle = convertDegreeToRadians( phiAngle );

        xInertial = nextStateVector[ xPositionIndex ] * std::cos( phiAngle )
                            - nextStateVector[ yPositionIndex ] * std::sin( phiAngle );
        yInertial = nextStateVector[ xPositionIndex ] * std::sin( phiAngle )
                            + nextStateVector[ yPositionIndex ] * std::cos( phiAngle );
        zInertial = nextStateVector[ zPositionIndex ];

        positionVector = { nextStateVector[ xPositionIndex ],
                           nextStateVector[ yPositionIndex ],
                           nextStateVector[ zPositionIndex ] };
        wCrossR = crossProduct( Wvector, positionVector );
        vxInertial = nextStateVector[ xVelocityIndex ] + wCrossR[ 0 ];
        vyInertial = nextStateVector[ yVelocityIndex ] + wCrossR[ 1 ];
        vzInertial = nextStateVector[ zVelocityIndex ] + wCrossR[ 2 ];

        inertialCartesianElements = { xInertial,
                                      yInertial,
                                      zInertial,
                                      vxInertial,
                                      vyInertial,
                                      vzInertial };
        convertCartesianCoordinatesToKeplerianElements( inertialCartesianElements,
                                                        gravParameter,
                                                        orbitalElements );

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

        // calculate the jacobian after each integration step
        xVelocitySquare = nextStateVector[ xVelocityIndex ]
                            * nextStateVector[ xVelocityIndex ];
        yVelocitySquare = nextStateVector[ yVelocityIndex ]
                            * nextStateVector[ yVelocityIndex ];
        zVelocitySquare = nextStateVector[ zVelocityIndex ]
                            * nextStateVector[ zVelocityIndex ];

        xPositionSquare = nextStateVector[ xPositionIndex ]
                            * nextStateVector[ xPositionIndex ];
        yPositionSquare = nextStateVector[ yPositionIndex ]
                            * nextStateVector[ yPositionIndex ];
        zPositionSquare = nextStateVector[ zPositionIndex ]
                            * nextStateVector[ zPositionIndex ];
        radialDistance = std::sqrt( xPositionSquare + yPositionSquare + zPositionSquare );

        wSquare = Wmagnitude * Wmagnitude;

        gravPotential = gravParameter / radialDistance;

        jacobian = 0.5 * ( xVelocitySquare + yVelocitySquare + zVelocitySquare )
                            - 0.5 * wSquare * ( xPositionSquare + yPositionSquare )
                            - gravPotential;

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
            eomOrbiterUREFile << tNext << ",";
            eomOrbiterUREFile << jacobian << ",";
            eomOrbiterUREFile << orbitalElements[ 0 ] << ",";
            eomOrbiterUREFile << orbitalElements[ 1 ] << ",";
            eomOrbiterUREFile << orbitalElements[ 2 ] << ",";
            eomOrbiterUREFile << orbitalElements[ 3 ] << ",";
            eomOrbiterUREFile << orbitalElements[ 4 ] << ",";
            eomOrbiterUREFile << orbitalElements[ 5 ] << std::endl;
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
        eomOrbiterUREFile << tNext << ",";
        eomOrbiterUREFile << jacobian << ",";
        eomOrbiterUREFile << orbitalElements[ 0 ] << ",";
        eomOrbiterUREFile << orbitalElements[ 1 ] << ",";
        eomOrbiterUREFile << orbitalElements[ 2 ] << ",";
        eomOrbiterUREFile << orbitalElements[ 3 ] << ",";
        eomOrbiterUREFile << orbitalElements[ 4 ] << ",";
        eomOrbiterUREFile << orbitalElements[ 5 ] << std::endl;

        // save the new state values in the vector of current state values. these will be used in
        // the next loop iteration
        tCurrent = tNext;
        currentStateVector = nextStateVector;
    }
 }

} // namespace naos

#endif // EXECUTE_ORBITER_AROUND_URE_POINT_MASS_GRAVITY_HPP
