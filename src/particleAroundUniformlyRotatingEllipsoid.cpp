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
#include <limits>
#include <vector>

#include <boost/numeric/odeint.hpp>
#include <SQLiteCpp/SQLiteCpp.h>
#include <sqlite3.h>

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/misc.hpp"
#include "NAOS/ellipsoidGravitationalAcceleration.hpp"
#include "NAOS/ellipsoidPotential.hpp"

namespace naos
{

//! equations of motion (for a particle around the asteroid modelled as an Ellipsoid)
/*!
 * first order differential equations describing the motion of a particle or spacecraft around a
 * central body modeled as an ellipsoid using ellipsoid gravitational model.
 */
class equationsOfMotionParticleAroundEllipsoid
{
    // declare parameters, gravitational parameter and the semi major axes of the ellipsoid
    const double gravParameter;
    const double alpha;
    const double beta;
    const double gamma;
    const double zRotation;

public:
    // Default constructor with member initializer list
    equationsOfMotionParticleAroundEllipsoid(
               const double aGravParameter,
               const double aAlpha,
               const double aBeta,
               const double aGamma,
               const double aZRotation )
            : gravParameter( aGravParameter ),
              alpha( aAlpha ),
              beta( aBeta ),
              gamma( aGamma ),
              zRotation( aZRotation )
    { }
    void operator() ( const std::vector< double > &stateVector,
                      std::vector< double > &dXdt,
                      const double currentTime )
    {
        // calculate the gravitational accelerations first
        std::vector< double > gravAcceleration( 3, 0.0 );

        computeEllipsoidGravitationalAcceleration( alpha,
                                                   beta,
                                                   gamma,
                                                   gravParameter,
                                                   stateVector[ xPositionIndex ],
                                                   stateVector[ yPositionIndex ],
                                                   stateVector[ zPositionIndex ],
                                                   gravAcceleration );

        // now calculate the derivatives
        dXdt[ xPositionIndex ] = stateVector[ xVelocityIndex ];
        dXdt[ yPositionIndex ] = stateVector[ yVelocityIndex ];
        dXdt[ zPositionIndex ] = stateVector[ zVelocityIndex ];

        dXdt[ xVelocityIndex ] = gravAcceleration[ xPositionIndex ]
                                + 2.0 * zRotation * stateVector[ yVelocityIndex ]
                                + zRotation * zRotation * stateVector[ xPositionIndex ];

        dXdt[ yVelocityIndex ] = gravAcceleration[ yPositionIndex ]
                                - 2.0 * zRotation * stateVector[ xVelocityIndex ]
                                + zRotation * zRotation * stateVector[ yPositionIndex ];

        dXdt[ zVelocityIndex ] = gravAcceleration[ zPositionIndex ];
    }
};

//! Store intermediate state values and time( if needed )
/*!
 * This structure contains members that will save all intermediate state values and times when
 * an object of this structure is passed as an argument to the integrator function.
 */
struct pushBackStateAndTime
{
    // declare containers to store state and time
    std::vector< std::vector< double > > &stateContainer;
    std::vector< double > &timeContainer;

    //member initializer list
    pushBackStateAndTime( std::vector< std::vector< double > > &aState,
                          std::vector< double > &aTime )
                : stateContainer( aState ),
                  timeContainer( aTime )
    { }

    void operator() ( const std::vector< double > &singleStateVector, const double singleTime )
    {
        // store the intermediate state and time values in the containers
        stateContainer.push_back( singleStateVector );
        timeContainer.push_back( singleTime );
    }
};

//! particle around ellipsoid integration
/*!
 * integrate the equations of motion for a particle around an ellipsoid. The gravitational accelerations
 * calculated using the ellipsoid gravitational potential model.
 */
void executeParticleAroundEllipsoid( const double alpha,
                                     const double beta,
                                     const double gamma,
                                     const double gravParameter,
                                     std::vector< double > asteroidRotationVector,
                                     std::vector< double > &initialOrbitalElements,
                                     const double initialStepSize,
                                     const double startTime,
                                     const double endTime,
                                     std::ostringstream &outputFilePath,
                                     const int dataSaveIntervals )
{
    //! open the output csv file to save data. Declare file headers.
    std::ofstream outputFile;
    outputFile.open( outputFilePath.str( ) );
    outputFile << "x" << ",";
    outputFile << "y" << ",";
    outputFile << "z" << ",";
    outputFile << "vx" << ",";
    outputFile << "vy" << ",";
    outputFile << "vz" << ",";
    outputFile << "t" << std::endl;
    outputFile.precision( 16 );

    //! convert the initial orbital elements to cartesian state
    std::vector< double > initialStateInertial( 6, 0.0 );
    initialStateInertial = convertKeplerianElementsToCartesianCoordinates( initialOrbitalElements,
                                                                           gravParameter );

    //! account for non-zero start time value and calculate the initial state in body frame
    double phi = asteroidRotationVector[ zPositionIndex ] * startTime;

    std::vector< double > initialState( 6, 0.0 );

    // get the initial body frame position
    initialState[ xPositionIndex ]
            = initialStateInertial[ xPositionIndex ] * std::cos( phi )
            + initialStateInertial[ yPositionIndex ] * std::sin( phi );

    initialState[ yPositionIndex ]
            = -1.0 * initialStateInertial[ xPositionIndex ] * std::sin( phi )
            + initialStateInertial[ yPositionIndex ] * std::cos( phi );

    initialState[ zPositionIndex ] = initialStateInertial[ zPositionIndex ];

    // get the initial body frame velocity
    std::vector< double > inertialPositionVector = { initialStateInertial[ xPositionIndex ],
                                                     initialStateInertial[ yPositionIndex ],
                                                     initialStateInertial[ zPositionIndex ] };

    std::vector< double > omegaCrossPosition( 3, 0.0 );
    omegaCrossPosition = crossProduct( asteroidRotationVector, inertialPositionVector );

    double xbodyFrameVelocityInertialCoordinates
                = initialStateInertial[ xVelocityIndex ] - omegaCrossPosition[ 0 ];

    double ybodyFrameVelocityInertialCoordinates
                = initialStateInertial[ yVelocityIndex ] - omegaCrossPosition[ 1 ];

    double zbodyFrameVelocityInertialCoordinates
                = initialStateInertial[ zVelocityIndex ] - omegaCrossPosition[ 2 ];

    initialState[ xVelocityIndex ]
            = xbodyFrameVelocityInertialCoordinates * std::cos( phi )
            + ybodyFrameVelocityInertialCoordinates * std::sin( phi );

    initialState[ yVelocityIndex ]
            = -1.0 * xbodyFrameVelocityInertialCoordinates * std::sin( phi )
            + ybodyFrameVelocityInertialCoordinates * std::cos( phi );

    initialState[ zVelocityIndex] = zbodyFrameVelocityInertialCoordinates;

    // set up boost odeint
    const double absoluteTolerance = 1.0e-15;
    const double relativeTolerance = 1.0e-15;
    typedef boost::numeric::odeint::runge_kutta_fehlberg78< std::vector< double > > stepperType;

    // state step size guess (at each step this initial guess will be used)
    double stepSizeGuess = initialStepSize;

    // initialize the ode system
    const double zRotation = asteroidRotationVector[ zPositionIndex ];
    equationsOfMotionParticleAroundEllipsoid particleAroundEllipsoidProblem( gravParameter,
                                                                             alpha,
                                                                             beta,
                                                                             gamma,
                                                                             zRotation );

    // initialize current state vector and time
    std::vector< double > currentStateVector = initialState;
    double currentTime = startTime;
    double intermediateEndTime = currentTime + dataSaveIntervals;

    // save the initial state vector
    outputFile << currentStateVector[ xPositionIndex ] << ",";
    outputFile << currentStateVector[ yPositionIndex ] << ",";
    outputFile << currentStateVector[ zPositionIndex ] << ",";
    outputFile << currentStateVector[ xVelocityIndex ] << ",";
    outputFile << currentStateVector[ yVelocityIndex ] << ",";
    outputFile << currentStateVector[ zVelocityIndex ] << ",";
    outputFile << currentTime << std::endl;

    // start the integration outer loop
    while( intermediateEndTime <= endTime )
    {
        // perform integration, integrated result stored in currentStateVector
        size_t steps = boost::numeric::odeint::integrate_adaptive(
                            make_controlled( absoluteTolerance, relativeTolerance, stepperType( ) ),
                            particleAroundEllipsoidProblem,
                            currentStateVector,
                            currentTime,
                            intermediateEndTime,
                            stepSizeGuess );

        // update the time variables
        currentTime = intermediateEndTime;
        intermediateEndTime = currentTime + dataSaveIntervals;

        // save data
        outputFile << currentStateVector[ xPositionIndex ] << ",";
        outputFile << currentStateVector[ yPositionIndex ] << ",";
        outputFile << currentStateVector[ zPositionIndex ] << ",";
        outputFile << currentStateVector[ xVelocityIndex ] << ",";
        outputFile << currentStateVector[ yVelocityIndex ] << ",";
        outputFile << currentStateVector[ zVelocityIndex ] << ",";
        outputFile << currentTime << std::endl;

    } // end of outer while loop for integration

    outputFile.close( );
}

//! Trajectory calculation for regolith around an asteroid (modelled as ellipsoid here)
/*!
 * Same as the previous function, except that the initial conditions are now given as a cartesian
 * state. The initial cartesian state should be given in body fixed frame of the asteroid.
 */
void singleRegolithTrajectoryCalculator( const double alpha,
                                         const double beta,
                                         const double gamma,
                                         const double gravParameter,
                                         std::vector< double > asteroidRotationVector,
                                         std::vector< double > &initialCartesianStateVector,
                                         const double initialStepSize,
                                         const double startTime,
                                         const double endTime,
                                         std::ostringstream &outputFilePath,
                                         const int dataSaveIntervals )
{
    //! open the output csv file to save data. Declare file headers.
    std::ofstream outputFile;
    outputFile.open( outputFilePath.str( ) );
    outputFile << "x" << ",";
    outputFile << "y" << ",";
    outputFile << "z" << ",";
    outputFile << "vx" << ",";
    outputFile << "vy" << ",";
    outputFile << "vz" << ",";
    outputFile << "t" << std::endl;
    outputFile.precision( 16 );

    //! get the initial cartesian state vector in a seperate container
    std::vector< double > initialState = initialCartesianStateVector;

    // set up boost odeint
    const double absoluteTolerance = 1.0e-15;
    const double relativeTolerance = 1.0e-15;
    typedef boost::numeric::odeint::runge_kutta_fehlberg78< std::vector< double > > stepperType;

    // state step size guess (at each step this initial guess will be used)
    double stepSizeGuess = initialStepSize;

    // initialize the ode system
    const double zRotation = asteroidRotationVector[ zPositionIndex ];
    equationsOfMotionParticleAroundEllipsoid particleAroundEllipsoidProblem( gravParameter,
                                                                             alpha,
                                                                             beta,
                                                                             gamma,
                                                                             zRotation );

    // initialize current state vector and time
    std::vector< double > currentStateVector = initialState;
    double currentTime = startTime;
    double intermediateEndTime = currentTime + dataSaveIntervals;

    // save the initial state vector
    outputFile << currentStateVector[ xPositionIndex ] << ",";
    outputFile << currentStateVector[ yPositionIndex ] << ",";
    outputFile << currentStateVector[ zPositionIndex ] << ",";
    outputFile << currentStateVector[ xVelocityIndex ] << ",";
    outputFile << currentStateVector[ yVelocityIndex ] << ",";
    outputFile << currentStateVector[ zVelocityIndex ] << ",";
    outputFile << currentTime << std::endl;

    // start the integration outer loop
    while( intermediateEndTime <= endTime )
    {
        // save the last know state vector for when the particle is outside the asteroid
        std::vector< double > lastStateVector = currentStateVector;

        // perform integration, integrated result stored in currentStateVector
        size_t steps = boost::numeric::odeint::integrate_adaptive(
                            make_controlled( absoluteTolerance, relativeTolerance, stepperType( ) ),
                            particleAroundEllipsoidProblem,
                            currentStateVector,
                            currentTime,
                            intermediateEndTime,
                            stepSizeGuess );

        //! check if the particle is on an escape trajectory or not
        // check for energy if the particle is far away from the asteroid
        double radialDistance
            = std::sqrt( currentStateVector[ xPositionIndex ] * currentStateVector[ xPositionIndex ]
                    + currentStateVector[ yPositionIndex ] * currentStateVector[ yPositionIndex ]
                    + currentStateVector[ zPositionIndex ] * currentStateVector[ zPositionIndex ] );

        if( radialDistance >= 10.0 * alpha )
        {
            // before calculating energy, convert the body frame state vector into inertial frame state vector
            double rotationAngle = asteroidRotationVector[ zPositionIndex ] * intermediateEndTime;

            std::vector< double > inertialState( 6, 0.0 );

            // get the inertial position
            inertialState[ xPositionIndex ]
                        = currentStateVector[ xPositionIndex ] * std::cos( rotationAngle )
                        - currentStateVector[ yPositionIndex ] * std::sin( rotationAngle );

            inertialState[ yPositionIndex ]
                        = currentStateVector[ xPositionIndex ] * std::sin( rotationAngle )
                        + currentStateVector[ yPositionIndex ] * std::cos( rotationAngle );

            inertialState[ zPositionIndex ] = currentStateVector[ zPositionIndex ];

            // get the omega cross position vector for use in the transport theorem
            std::vector< double > bodyFramePositionVector = { currentStateVector[ xPositionIndex ],
                                                              currentStateVector[ yPositionIndex ],
                                                              currentStateVector[ zPositionIndex ] };

            std::vector< double > omegaCrossPosition( 3, 0.0 );
            omegaCrossPosition = crossProduct( asteroidRotationVector, bodyFramePositionVector );

            // use the transport theorem to get the body frame velocity in inertial coordinates
            double xBodyFrameVelocityInertialCoordinates
                        = currentStateVector[ xVelocityIndex ] + omegaCrossPosition[ 0 ];

            double yBodyFrameVelocityInertialCoordinates
                        = currentStateVector[ yVelocityIndex ] + omegaCrossPosition[ 1 ];

            double zBodyFrameVelocityInertialCoordinates
                        = currentStateVector[ zVelocityIndex ] + omegaCrossPosition[ 2 ];

            // use the rotation matrix to get the velocity in the inertial frame
            inertialState[ xVelocityIndex ]
                        = xBodyFrameVelocityInertialCoordinates * std::cos( rotationAngle )
                        - yBodyFrameVelocityInertialCoordinates * std::sin( rotationAngle );

            inertialState[ yVelocityIndex ]
                        = xBodyFrameVelocityInertialCoordinates * std::sin( rotationAngle )
                        + yBodyFrameVelocityInertialCoordinates * std::cos( rotationAngle );

            inertialState[ zVelocityIndex ] = zBodyFrameVelocityInertialCoordinates;

            // calculate the orbital elements and check the sign of the eccentricity
            std::vector< double > orbitalElements( 6, 0.0 );
            convertCartesianCoordinatesToKeplerianElements( inertialState,
                                                            gravParameter,
                                                            orbitalElements );
            bool eccentricityFlag = false;
            if( orbitalElements[ 1 ] < 0.0 || orbitalElements[ 1 ] >= 1.0 )
            {
                eccentricityFlag = true;
            }

            // calculate energy and check the sign
            std::vector< double > inertialPositionVector = { inertialState[ xPositionIndex ],
                                                             inertialState[ yPositionIndex ],
                                                             inertialState[ zPositionIndex ] };

            std::vector< double > inertialVelocityVector = { inertialState[ xVelocityIndex ],
                                                             inertialState[ yVelocityIndex ],
                                                             inertialState[ zVelocityIndex ] };

            double inertialVelocityMagnitude = vectorNorm( inertialVelocityVector );

            double inertialPositionMagnitude = vectorNorm( inertialPositionVector );

            double gravPotential;
            computeEllipsoidGravitationalPotential( alpha,
                                                    beta,
                                                    gamma,
                                                    gravParameter,
                                                    inertialPositionVector[ xPositionIndex ],
                                                    inertialPositionVector[ yPositionIndex ],
                                                    inertialPositionVector[ zPositionIndex ],
                                                    gravPotential );
            double particleEnergy
                    = inertialVelocityMagnitude * inertialVelocityMagnitude / 2.0
                    - gravPotential;

            bool energyFlag = false;
            if( particleEnergy > 0.0 )
            {
                energyFlag = true;
            }

            // check the eccentricity and energy flags
            if( eccentricityFlag && energyFlag )
            {
                // particle is on an escape trajectory, save data and stop integration
                // update the time variables
                currentTime = intermediateEndTime;

                // save data
                outputFile << currentStateVector[ xPositionIndex ] << ",";
                outputFile << currentStateVector[ yPositionIndex ] << ",";
                outputFile << currentStateVector[ zPositionIndex ] << ",";
                outputFile << currentStateVector[ xVelocityIndex ] << ",";
                outputFile << currentStateVector[ yVelocityIndex ] << ",";
                outputFile << currentStateVector[ zVelocityIndex ] << ",";
                outputFile << currentTime << std::endl;

                break;
            }
        }

        //! check if the particle is inside the surface of the asteroid
        double xSquare = currentStateVector[ xPositionIndex ] * currentStateVector[ xPositionIndex ];
        double ySquare = currentStateVector[ yPositionIndex ] * currentStateVector[ yPositionIndex ];
        double zSquare = currentStateVector[ zPositionIndex ] * currentStateVector[ zPositionIndex ];

        double crashCheck = xSquare / ( alpha * alpha )
                            + ySquare / ( beta * beta )
                            + zSquare / ( gamma * gamma )
                            - 1.0;

        if( crashCheck == 0.0 )
        {
            // particle is on the surface of the asteroid, save data and stop integration
            // update the time variables
            currentTime = intermediateEndTime;

            // save data
            outputFile << currentStateVector[ xPositionIndex ] << ",";
            outputFile << currentStateVector[ yPositionIndex ] << ",";
            outputFile << currentStateVector[ zPositionIndex ] << ",";
            outputFile << currentStateVector[ xVelocityIndex ] << ",";
            outputFile << currentStateVector[ yVelocityIndex ] << ",";
            outputFile << currentStateVector[ zVelocityIndex ] << ",";
            outputFile << currentTime << std::endl;

            break;
        }

        if( crashCheck < 0.0 )
        {
            double stepSize = 1.0;

            // const double machinePrecision = std::numeric_limits< double >::epsilon( );
            const double machinePrecision = 1.0e-15;

            // if the particle is on the surface of the asteroid, then the while condition
            // will be false, for all other cases it will be true. Within the while loop, the
            // inside or outside differnetiation takes place.
            while( std::fabs( crashCheck ) > machinePrecision )
            {
                // std::cout << "crash check value = " << crashCheck << std::endl;
                if( crashCheck < 0.0 ) // particle is still inside the surface
                {
                    // std::cout << "particle inside the surface" << std::endl << std::endl;
                    // particle is inside the surface of the asteroid. restart the integration from last
                    // known state external to the asteroid
                    currentStateVector = lastStateVector;
                    stepSize = 0.5 * stepSize;
                }
                else // particle is outside the asteroid at the end of last integration step
                {
                    // std::cout << "particle outside the surface" << std::endl << std::endl;
                    lastStateVector = currentStateVector;
                    currentTime = intermediateEndTime;
                }

                typedef boost::numeric::odeint::runge_kutta_fehlberg78< std::vector< double > > errorStepperType;
                intermediateEndTime =  boost::numeric::odeint::integrate_n_steps(
                                                errorStepperType( ),
                                                particleAroundEllipsoidProblem,
                                                currentStateVector,
                                                currentTime,
                                                stepSize,
                                                1 );

                xSquare = currentStateVector[ xPositionIndex ] * currentStateVector[ xPositionIndex ];
                ySquare = currentStateVector[ yPositionIndex ] * currentStateVector[ yPositionIndex ];
                zSquare = currentStateVector[ zPositionIndex ] * currentStateVector[ zPositionIndex ];

                crashCheck = xSquare / ( alpha * alpha )
                            + ySquare / ( beta * beta )
                            + zSquare / ( gamma * gamma )
                            - 1.0;
            }

            // particle is on the surface of the asteroid, save data and stop integration
            // update the time variables
            currentTime = intermediateEndTime;

            // save data
            outputFile << currentStateVector[ xPositionIndex ] << ",";
            outputFile << currentStateVector[ yPositionIndex ] << ",";
            outputFile << currentStateVector[ zPositionIndex ] << ",";
            outputFile << currentStateVector[ xVelocityIndex ] << ",";
            outputFile << currentStateVector[ yVelocityIndex ] << ",";
            outputFile << currentStateVector[ zVelocityIndex ] << ",";
            outputFile << currentTime << std::endl;

            break;
        }

        // update the time variables
        currentTime = intermediateEndTime;
        intermediateEndTime = currentTime + dataSaveIntervals;

        // save data
        outputFile << currentStateVector[ xPositionIndex ] << ",";
        outputFile << currentStateVector[ yPositionIndex ] << ",";
        outputFile << currentStateVector[ zPositionIndex ] << ",";
        outputFile << currentStateVector[ xVelocityIndex ] << ",";
        outputFile << currentStateVector[ yVelocityIndex ] << ",";
        outputFile << currentStateVector[ zVelocityIndex ] << ",";
        outputFile << currentTime << std::endl;

    } // end of outer while loop for integration

    outputFile.close( );
}

//! convert body frame state vector to inertial frame for a uniformly rotating case
/*!
 * this function cnverts a given state vector in asteroid body frame coordinates to inertial frame.
 * The body frame is rotating uniformly with respect to the inertial frame about the z axis.
 */
void convertBodyFrameVectorToInertialFrame( const std::vector< double > &asteroidRotationVector,
                                            const std::vector< double > &currentStateVector,
                                            const double currentTime,
                                            std::vector< double > &inertialState )
{
    // convert the body frame state vector into inertial frame state vector
    double rotationAngle = asteroidRotationVector[ zPositionIndex ] * currentTime;

    // get the inertial position
    inertialState[ xPositionIndex ]
                = currentStateVector[ xPositionIndex ] * std::cos( rotationAngle )
                - currentStateVector[ yPositionIndex ] * std::sin( rotationAngle );

    inertialState[ yPositionIndex ]
                = currentStateVector[ xPositionIndex ] * std::sin( rotationAngle )
                + currentStateVector[ yPositionIndex ] * std::cos( rotationAngle );

    inertialState[ zPositionIndex ] = currentStateVector[ zPositionIndex ];

    // get the omega cross position vector for use in the transport theorem
    std::vector< double > bodyFramePositionVector = { currentStateVector[ xPositionIndex ],
                                                      currentStateVector[ yPositionIndex ],
                                                      currentStateVector[ zPositionIndex ] };

    std::vector< double > omegaCrossPosition( 3, 0.0 );
    omegaCrossPosition = crossProduct( asteroidRotationVector, bodyFramePositionVector );

    // use the transport theorem to get the body frame velocity in inertial coordinates
    double xBodyFrameVelocityInertialCoordinates
                = currentStateVector[ xVelocityIndex ] + omegaCrossPosition[ 0 ];

    double yBodyFrameVelocityInertialCoordinates
                = currentStateVector[ yVelocityIndex ] + omegaCrossPosition[ 1 ];

    double zBodyFrameVelocityInertialCoordinates
                = currentStateVector[ zVelocityIndex ] + omegaCrossPosition[ 2 ];

    // use the rotation matrix to get the velocity in the inertial frame
    inertialState[ xVelocityIndex ]
                = xBodyFrameVelocityInertialCoordinates * std::cos( rotationAngle )
                - yBodyFrameVelocityInertialCoordinates * std::sin( rotationAngle );

    inertialState[ yVelocityIndex ]
                = xBodyFrameVelocityInertialCoordinates * std::sin( rotationAngle )
                + yBodyFrameVelocityInertialCoordinates * std::cos( rotationAngle );

    inertialState[ zVelocityIndex ] = zBodyFrameVelocityInertialCoordinates;
}

//! Trajectory calculation for regolith around an asteroid (modelled as ellipsoid here)
/*!
 * Same as the previous function, except that the initial conditions are now given as a cartesian
 * state and the data is not saved in a csv file.
 * The initial cartesian state should be given in body fixed frame of the asteroid.
 */
void executeSingleRegolithTrajectoryCalculation( const double alpha,
                                                 const double beta,
                                                 const double gamma,
                                                 const double gravParameter,
                                                 std::vector< double > asteroidRotationVector,
                                                 std::vector< double > &initialCartesianStateVector,
                                                 const double launchAzimuth,
                                                 const double launchDeclination,
                                                 const double initialStepSize,
                                                 const double startTime,
                                                 const double endTime,
                                                 SQLite::Statement &databaseQuery,
                                                 const int dataSaveIntervals )
{
    //! get the initial cartesian state vector in a seperate container
    std::vector< double > initialState = initialCartesianStateVector;

    //! get the initial orbital elements
    std::vector< double > initialOrbitalElements( 6, 0.0 );
    // convert the initial body frame state to inertial frame first
    std::vector< double > initialInertialState( 6, 0.0 );
    convertBodyFrameVectorToInertialFrame( asteroidRotationVector,
                                           initialState,
                                           startTime,
                                           initialInertialState );
    // convert to orbital elements
    convertCartesianCoordinatesToKeplerianElements( initialInertialState,
                                                    gravParameter,
                                                    initialOrbitalElements );

    //! get the initial energy of the particle
    double initialGravPotential;
    computeEllipsoidGravitationalPotential( alpha,
                                            beta,
                                            gamma,
                                            gravParameter,
                                            initialInertialState[ xPositionIndex ],
                                            initialInertialState[ yPositionIndex ],
                                            initialInertialState[ zPositionIndex ],
                                            initialGravPotential );

    std::vector< double > initialInertialVelocityVector = { initialInertialState[ xVelocityIndex ],
                                                            initialInertialState[ yVelocityIndex ],
                                                            initialInertialState[ zVelocityIndex ] };

    double initialInertialVelocityMagnitude = vectorNorm( initialInertialVelocityVector );

    double initialParticleEnergy
            = initialInertialVelocityMagnitude * initialInertialVelocityMagnitude / 2.0
            - initialGravPotential;

    // set up boost odeint
    const double absoluteTolerance = 1.0e-15;
    const double relativeTolerance = 1.0e-15;
    typedef boost::numeric::odeint::runge_kutta_fehlberg78< std::vector< double > > stepperType;

    // state step size guess (at each step this initial guess will be used)
    double stepSizeGuess = initialStepSize;

    // initialize the ode system
    const double zRotation = asteroidRotationVector[ zPositionIndex ];
    equationsOfMotionParticleAroundEllipsoid particleAroundEllipsoidProblem( gravParameter,
                                                                             alpha,
                                                                             beta,
                                                                             gamma,
                                                                             zRotation );

    // initialize current state vector and time
    std::vector< double > currentStateVector = initialState;
    double currentTime = startTime;
    double intermediateEndTime = currentTime + dataSaveIntervals;

    // initialize and set flags for crash, escape
    int escapeFlag = 0;
    int crashFlag = 0;

    // save the initial state vector(body frame) and orbital elements
    databaseQuery.bind( ":position_x", initialState[ xPositionIndex ] );
    databaseQuery.bind( ":position_y", initialState[ yPositionIndex ] );
    databaseQuery.bind( ":position_z", initialState[ zPositionIndex ] );
    databaseQuery.bind( ":velocity_x", initialState[ xVelocityIndex ] );
    databaseQuery.bind( ":velocity_y", initialState[ yVelocityIndex ] );
    databaseQuery.bind( ":velocity_z", initialState[ zVelocityIndex ] );
    databaseQuery.bind( ":time", currentTime );

    databaseQuery.bind( ":launch_azimuth", convertRadiansToDegree( launchAzimuth ) );
    databaseQuery.bind( ":launch_declination", convertRadiansToDegree( launchDeclination ) );

    databaseQuery.bind( ":sma", initialOrbitalElements[ 0 ] );
    databaseQuery.bind( ":eccentricity", initialOrbitalElements[ 1 ] );
    databaseQuery.bind( ":inclination", initialOrbitalElements[ 2 ] );
    databaseQuery.bind( ":raan", initialOrbitalElements[ 3 ] );
    databaseQuery.bind( ":aop", initialOrbitalElements[ 4 ] );
    databaseQuery.bind( ":ta", initialOrbitalElements[ 5 ] );

    databaseQuery.bind( ":energy", initialParticleEnergy );

    databaseQuery.bind( ":escape_flag", escapeFlag );
    databaseQuery.bind( ":crash_flag", crashFlag );

    // Execute insert query.
    databaseQuery.executeStep( );

    // Reset SQL insert query.
    databaseQuery.reset( );

    // start the integration outer loop
    while( intermediateEndTime <= endTime )
    {
        // save the last know state vector for when the particle is outside the asteroid
        std::vector< double > lastStateVector = currentStateVector;

        // perform integration, integrated result stored in currentStateVector
        size_t steps = boost::numeric::odeint::integrate_adaptive(
                            make_controlled( absoluteTolerance, relativeTolerance, stepperType( ) ),
                            particleAroundEllipsoidProblem,
                            currentStateVector,
                            currentTime,
                            intermediateEndTime,
                            stepSizeGuess );

        //! check if the particle is on an escape trajectory or not
        // check for energy if the particle is far away from the asteroid
        double radialDistance
            = std::sqrt( currentStateVector[ xPositionIndex ] * currentStateVector[ xPositionIndex ]
                    + currentStateVector[ yPositionIndex ] * currentStateVector[ yPositionIndex ]
                    + currentStateVector[ zPositionIndex ] * currentStateVector[ zPositionIndex ] );

        if( radialDistance >= 10.0 * alpha )
        {
            // convert body frame state vector to inertial frame
            std::vector< double > inertialState( 6, 0.0 );
            convertBodyFrameVectorToInertialFrame( asteroidRotationVector,
                                                   currentStateVector,
                                                   intermediateEndTime,
                                                   inertialState );

            // calculate the orbital elements and check the sign of the eccentricity
            std::vector< double > orbitalElements( 6, 0.0 );
            convertCartesianCoordinatesToKeplerianElements( inertialState,
                                                            gravParameter,
                                                            orbitalElements );
            bool eccentricityFlag = false;
            if( orbitalElements[ 1 ] < 0.0 || orbitalElements[ 1 ] >= 1.0 )
            {
                eccentricityFlag = true;
            }

            // calculate energy and check the sign
            std::vector< double > inertialPositionVector = { inertialState[ xPositionIndex ],
                                                             inertialState[ yPositionIndex ],
                                                             inertialState[ zPositionIndex ] };

            std::vector< double > inertialVelocityVector = { inertialState[ xVelocityIndex ],
                                                             inertialState[ yVelocityIndex ],
                                                             inertialState[ zVelocityIndex ] };

            double inertialVelocityMagnitude = vectorNorm( inertialVelocityVector );

            double inertialPositionMagnitude = vectorNorm( inertialPositionVector );

            double gravPotential;
            computeEllipsoidGravitationalPotential( alpha,
                                                    beta,
                                                    gamma,
                                                    gravParameter,
                                                    inertialPositionVector[ xPositionIndex ],
                                                    inertialPositionVector[ yPositionIndex ],
                                                    inertialPositionVector[ zPositionIndex ],
                                                    gravPotential );
            double particleEnergy
                    = inertialVelocityMagnitude * inertialVelocityMagnitude / 2.0
                    - gravPotential;

            bool energyFlag = false;
            if( particleEnergy > 0.0 )
            {
                energyFlag = true;
            }

            // check the eccentricity and energy flags
            if( eccentricityFlag && energyFlag )
            {
                // particle is on an escape trajectory, save data and stop integration
                // update the time variables
                currentTime = intermediateEndTime;

                // set the flags
                escapeFlag = 1;
                crashFlag = 0;

                // save data
                databaseQuery.bind( ":position_x", currentStateVector[ xPositionIndex ] );
                databaseQuery.bind( ":position_y", currentStateVector[ yPositionIndex ] );
                databaseQuery.bind( ":position_z", currentStateVector[ zPositionIndex ] );
                databaseQuery.bind( ":velocity_x", currentStateVector[ xVelocityIndex ] );
                databaseQuery.bind( ":velocity_y", currentStateVector[ yVelocityIndex ] );
                databaseQuery.bind( ":velocity_z", currentStateVector[ zVelocityIndex ] );
                databaseQuery.bind( ":time", currentTime );

                databaseQuery.bind( ":launch_azimuth", convertRadiansToDegree( launchAzimuth ) );
                databaseQuery.bind( ":launch_declination", convertRadiansToDegree( launchDeclination ) );

                databaseQuery.bind( ":sma", orbitalElements[ 0 ] );
                databaseQuery.bind( ":eccentricity", orbitalElements[ 1 ] );
                databaseQuery.bind( ":inclination", orbitalElements[ 2 ] );
                databaseQuery.bind( ":raan", orbitalElements[ 3 ] );
                databaseQuery.bind( ":aop", orbitalElements[ 4 ] );
                databaseQuery.bind( ":ta", orbitalElements[ 5 ] );

                databaseQuery.bind( ":energy", particleEnergy );

                databaseQuery.bind( ":escape_flag", escapeFlag );
                databaseQuery.bind( ":crash_flag", crashFlag );

                // Execute insert query.
                databaseQuery.executeStep( );

                // Reset SQL insert query.
                databaseQuery.reset( );

                break;
            }
        }

        //! check if the particle is inside the surface of the asteroid
        double xSquare = currentStateVector[ xPositionIndex ] * currentStateVector[ xPositionIndex ];
        double ySquare = currentStateVector[ yPositionIndex ] * currentStateVector[ yPositionIndex ];
        double zSquare = currentStateVector[ zPositionIndex ] * currentStateVector[ zPositionIndex ];

        double crashCheck = xSquare / ( alpha * alpha )
                            + ySquare / ( beta * beta )
                            + zSquare / ( gamma * gamma )
                            - 1.0;

        if( crashCheck == 0.0 )
        {
            // particle is on the surface of the asteroid, save data and stop integration
            // update the time variables
            currentTime = intermediateEndTime;

            // set the flags
            escapeFlag = 0;
            crashFlag = 1;

            // get the orbital elements and the energy
            std::vector< double > inertialState( 6, 0.0 );
            convertBodyFrameVectorToInertialFrame( asteroidRotationVector,
                                                   currentStateVector,
                                                   currentTime,
                                                   inertialState );

            std::vector< double > orbitalElements( 6, 0.0 );
            convertCartesianCoordinatesToKeplerianElements( inertialState,
                                                            gravParameter,
                                                            orbitalElements );

            std::vector< double > inertialVelocityVector = { inertialState[ xVelocityIndex ],
                                                             inertialState[ yVelocityIndex ],
                                                             inertialState[ zVelocityIndex ] };

            double inertialVelocityMagnitude = vectorNorm( inertialVelocityVector );

            double gravPotential;
            computeEllipsoidGravitationalPotential( alpha,
                                                    beta,
                                                    gamma,
                                                    gravParameter,
                                                    currentStateVector[ xPositionIndex ],
                                                    currentStateVector[ yPositionIndex ],
                                                    currentStateVector[ zPositionIndex ],
                                                    gravPotential );
            double particleEnergy
                    = inertialVelocityMagnitude * inertialVelocityMagnitude / 2.0
                    - gravPotential;

            // save data
            databaseQuery.bind( ":position_x", currentStateVector[ xPositionIndex ] );
            databaseQuery.bind( ":position_y", currentStateVector[ yPositionIndex ] );
            databaseQuery.bind( ":position_z", currentStateVector[ zPositionIndex ] );
            databaseQuery.bind( ":velocity_x", currentStateVector[ xVelocityIndex ] );
            databaseQuery.bind( ":velocity_y", currentStateVector[ yVelocityIndex ] );
            databaseQuery.bind( ":velocity_z", currentStateVector[ zVelocityIndex ] );
            databaseQuery.bind( ":time", currentTime );

            databaseQuery.bind( ":launch_azimuth", convertRadiansToDegree( launchAzimuth ) );
            databaseQuery.bind( ":launch_declination", convertRadiansToDegree( launchDeclination ) );

            databaseQuery.bind( ":sma", orbitalElements[ 0 ] );
            databaseQuery.bind( ":eccentricity", orbitalElements[ 1 ] );
            databaseQuery.bind( ":inclination", orbitalElements[ 2 ] );
            databaseQuery.bind( ":raan", orbitalElements[ 3 ] );
            databaseQuery.bind( ":aop", orbitalElements[ 4 ] );
            databaseQuery.bind( ":ta", orbitalElements[ 5 ] );

            databaseQuery.bind( ":energy", particleEnergy );

            databaseQuery.bind( ":escape_flag", escapeFlag );
            databaseQuery.bind( ":crash_flag", crashFlag );

            // Execute insert query.
            databaseQuery.executeStep( );

            // Reset SQL insert query.
            databaseQuery.reset( );

            break;
        }

        if( crashCheck < 0.0 )
        {
            double stepSize = 1.0;

            // const double machinePrecision = std::numeric_limits< double >::epsilon( );
            const double machinePrecision = 1.0e-15;

            // if the particle is on the surface of the asteroid, then the while condition
            // will be false, for all other cases it will be true. Within the while loop, the
            // inside or outside differnetiation takes place.
            while( std::fabs( crashCheck ) > machinePrecision )
            {
                // std::cout << "crash check value = " << crashCheck << std::endl;
                if( crashCheck < 0.0 ) // particle is still inside the surface
                {
                    // std::cout << "particle inside the surface" << std::endl << std::endl;
                    // particle is inside the surface of the asteroid. restart the integration from last
                    // known state external to the asteroid
                    currentStateVector = lastStateVector;
                    stepSize = 0.5 * stepSize;
                }
                else // particle is outside the asteroid at the end of last integration step
                {
                    // std::cout << "particle outside the surface" << std::endl << std::endl;
                    lastStateVector = currentStateVector;
                    currentTime = intermediateEndTime;
                }

                typedef boost::numeric::odeint::runge_kutta_fehlberg78< std::vector< double > > errorStepperType;
                intermediateEndTime =  boost::numeric::odeint::integrate_n_steps(
                                                errorStepperType( ),
                                                particleAroundEllipsoidProblem,
                                                currentStateVector,
                                                currentTime,
                                                stepSize,
                                                1 );

                xSquare = currentStateVector[ xPositionIndex ] * currentStateVector[ xPositionIndex ];
                ySquare = currentStateVector[ yPositionIndex ] * currentStateVector[ yPositionIndex ];
                zSquare = currentStateVector[ zPositionIndex ] * currentStateVector[ zPositionIndex ];

                crashCheck = xSquare / ( alpha * alpha )
                            + ySquare / ( beta * beta )
                            + zSquare / ( gamma * gamma )
                            - 1.0;
            }

            // particle is on the surface of the asteroid, save data and stop integration
            // update the time variables
            currentTime = intermediateEndTime;

            // set the flags
            escapeFlag = 0;
            crashFlag = 1;

            // get the orbital elements and the energy
            std::vector< double > inertialState( 6, 0.0 );
            convertBodyFrameVectorToInertialFrame( asteroidRotationVector,
                                                   currentStateVector,
                                                   currentTime,
                                                   inertialState );

            std::vector< double > orbitalElements( 6, 0.0 );
            convertCartesianCoordinatesToKeplerianElements( inertialState,
                                                            gravParameter,
                                                            orbitalElements );

            std::vector< double > inertialVelocityVector = { inertialState[ xVelocityIndex ],
                                                             inertialState[ yVelocityIndex ],
                                                             inertialState[ zVelocityIndex ] };

            double inertialVelocityMagnitude = vectorNorm( inertialVelocityVector );

            double gravPotential;
            computeEllipsoidGravitationalPotential( alpha,
                                                    beta,
                                                    gamma,
                                                    gravParameter,
                                                    currentStateVector[ xPositionIndex ],
                                                    currentStateVector[ yPositionIndex ],
                                                    currentStateVector[ zPositionIndex ],
                                                    gravPotential );
            double particleEnergy
                    = inertialVelocityMagnitude * inertialVelocityMagnitude / 2.0
                    - gravPotential;

            // save data
            databaseQuery.bind( ":position_x", currentStateVector[ xPositionIndex ] );
            databaseQuery.bind( ":position_y", currentStateVector[ yPositionIndex ] );
            databaseQuery.bind( ":position_z", currentStateVector[ zPositionIndex ] );
            databaseQuery.bind( ":velocity_x", currentStateVector[ xVelocityIndex ] );
            databaseQuery.bind( ":velocity_y", currentStateVector[ yVelocityIndex ] );
            databaseQuery.bind( ":velocity_z", currentStateVector[ zVelocityIndex ] );
            databaseQuery.bind( ":time", currentTime );

            databaseQuery.bind( ":launch_azimuth", convertRadiansToDegree( launchAzimuth ) );
            databaseQuery.bind( ":launch_declination", convertRadiansToDegree( launchDeclination ) );

            databaseQuery.bind( ":sma", orbitalElements[ 0 ] );
            databaseQuery.bind( ":eccentricity", orbitalElements[ 1 ] );
            databaseQuery.bind( ":inclination", orbitalElements[ 2 ] );
            databaseQuery.bind( ":raan", orbitalElements[ 3 ] );
            databaseQuery.bind( ":aop", orbitalElements[ 4 ] );
            databaseQuery.bind( ":ta", orbitalElements[ 5 ] );

            databaseQuery.bind( ":energy", particleEnergy );

            databaseQuery.bind( ":escape_flag", escapeFlag );
            databaseQuery.bind( ":crash_flag", crashFlag );

            // Execute insert query.
            databaseQuery.executeStep( );

            // Reset SQL insert query.
            databaseQuery.reset( );

            break;
        }

        // update the time variables
        currentTime = intermediateEndTime;
        intermediateEndTime = currentTime + dataSaveIntervals;

        // set the flags
        escapeFlag = 0;
        crashFlag = 0;

        // get the orbital elements and the energy
        std::vector< double > inertialState( 6, 0.0 );
        convertBodyFrameVectorToInertialFrame( asteroidRotationVector,
                                               currentStateVector,
                                               currentTime,
                                               inertialState );

        std::vector< double > orbitalElements( 6, 0.0 );
        convertCartesianCoordinatesToKeplerianElements( inertialState,
                                                        gravParameter,
                                                        orbitalElements );

        std::vector< double > inertialVelocityVector = { inertialState[ xVelocityIndex ],
                                                         inertialState[ yVelocityIndex ],
                                                         inertialState[ zVelocityIndex ] };

        double inertialVelocityMagnitude = vectorNorm( inertialVelocityVector );

        double gravPotential;
        computeEllipsoidGravitationalPotential( alpha,
                                                beta,
                                                gamma,
                                                gravParameter,
                                                currentStateVector[ xPositionIndex ],
                                                currentStateVector[ yPositionIndex ],
                                                currentStateVector[ zPositionIndex ],
                                                gravPotential );
        double particleEnergy
                = inertialVelocityMagnitude * inertialVelocityMagnitude / 2.0
                - gravPotential;

        // save data
        databaseQuery.bind( ":position_x", currentStateVector[ xPositionIndex ] );
        databaseQuery.bind( ":position_y", currentStateVector[ yPositionIndex ] );
        databaseQuery.bind( ":position_z", currentStateVector[ zPositionIndex ] );
        databaseQuery.bind( ":velocity_x", currentStateVector[ xVelocityIndex ] );
        databaseQuery.bind( ":velocity_y", currentStateVector[ yVelocityIndex ] );
        databaseQuery.bind( ":velocity_z", currentStateVector[ zVelocityIndex ] );
        databaseQuery.bind( ":time", currentTime );

        databaseQuery.bind( ":launch_azimuth", convertRadiansToDegree( launchAzimuth ) );
        databaseQuery.bind( ":launch_declination", convertRadiansToDegree( launchDeclination ) );

        databaseQuery.bind( ":sma", orbitalElements[ 0 ] );
        databaseQuery.bind( ":eccentricity", orbitalElements[ 1 ] );
        databaseQuery.bind( ":inclination", orbitalElements[ 2 ] );
        databaseQuery.bind( ":raan", orbitalElements[ 3 ] );
        databaseQuery.bind( ":aop", orbitalElements[ 4 ] );
        databaseQuery.bind( ":ta", orbitalElements[ 5 ] );

        databaseQuery.bind( ":energy", particleEnergy );

        databaseQuery.bind( ":escape_flag", escapeFlag );
        databaseQuery.bind( ":crash_flag", crashFlag );

        // Execute insert query.
        databaseQuery.executeStep( );

        // Reset SQL insert query.
        databaseQuery.reset( );
    } // end of outer while loop for integration
}

} // namespace naos
