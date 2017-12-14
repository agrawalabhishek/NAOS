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
#include <stdexcept>

#include <boost/numeric/odeint.hpp>
#include <SQLiteCpp/SQLiteCpp.h>
#include <sqlite3.h>

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/misc.hpp"
#include "NAOS/ellipsoidGravitationalAcceleration.hpp"
#include "NAOS/ellipsoidPotential.hpp"
#include "NAOS/perturbingAccelerations.hpp"
#include "NAOS/sunAsteroidKeplerProblemSolver.hpp"

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

class perturbedEquationsOfMotionParticleAroundEllipsoid
{
private:
    // declare parameters, gravitational parameter and the semi major axes of the ellipsoid
    const double gravParameter;
    const double alpha;
    const double beta;
    const double gamma;
    std::vector< double > asteroidRotationVector;
    const double initialTime;
    const double initialMeanAnomalyRadian;
    const std::vector< double > initialSunOrbitalElements;
    const double sunMeanMotion;

public:
    // Default constructor with member initializer list
    perturbedEquationsOfMotionParticleAroundEllipsoid(
               const double aGravParameter,
               const double aAlpha,
               const double aBeta,
               const double aGamma,
               std::vector< double > &anAsteroidRotationVector,
               const double anInitialTime,
               const double anInitialMeanAnomalyRadian,
               const std::vector< double > &anInitialSunOrbitalElements,
               const double aSunMeanMotion )
            : gravParameter( aGravParameter ),
              alpha( aAlpha ),
              beta( aBeta ),
              gamma( aGamma ),
              asteroidRotationVector( anAsteroidRotationVector ),
              initialTime( anInitialTime ),
              initialMeanAnomalyRadian( anInitialMeanAnomalyRadian ),
              initialSunOrbitalElements( anInitialSunOrbitalElements ),
              sunMeanMotion( aSunMeanMotion )
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

        // compute perturbing acceleration from the third body effect of Sun
        std::vector< double > regolithPositionVector = { stateVector[ xPositionIndex ],
                                                         stateVector[ yPositionIndex ],
                                                         stateVector[ zPositionIndex ] };

        std::vector< double > sunThirdBodyEffectAcceleration( 3, 0.0 );
        sunThirdBodyEffectAcceleration
            = computeSunThirdBodyEffectAcceleration( regolithPositionVector,
                                                     asteroidRotationVector,
                                                     initialTime,
                                                     initialMeanAnomalyRadian,
                                                     initialSunOrbitalElements,
                                                     currentTime,
                                                     sunMeanMotion );

        // compute perturbing acceleration from the solar radiation pressure
        std::vector< double > solarRadiationPressureAcceleration( 3, 0.0 );
        solarRadiationPressureAcceleration
            = computeSolarRadiationPressureAcceleration( regolithPositionVector,
                                                         asteroidRotationVector,
                                                         initialTime,
                                                         initialMeanAnomalyRadian,
                                                         initialSunOrbitalElements,
                                                         currentTime,
                                                         sunMeanMotion );

        double zRotation = asteroidRotationVector[ 2 ];

        // now calculate the derivatives
        dXdt[ xPositionIndex ] = stateVector[ xVelocityIndex ];
        dXdt[ yPositionIndex ] = stateVector[ yVelocityIndex ];
        dXdt[ zPositionIndex ] = stateVector[ zVelocityIndex ];

        dXdt[ xVelocityIndex ] = gravAcceleration[ xPositionIndex ]
                                + 2.0 * zRotation * stateVector[ yVelocityIndex ]
                                + zRotation * zRotation * stateVector[ xPositionIndex ]
                                + sunThirdBodyEffectAcceleration[ xPositionIndex ]
                                + solarRadiationPressureAcceleration[ xPositionIndex ];

        dXdt[ yVelocityIndex ] = gravAcceleration[ yPositionIndex ]
                                - 2.0 * zRotation * stateVector[ xVelocityIndex ]
                                + zRotation * zRotation * stateVector[ yPositionIndex ]
                                + sunThirdBodyEffectAcceleration[ yPositionIndex ]
                                + solarRadiationPressureAcceleration[ yPositionIndex ];

        dXdt[ zVelocityIndex ] = gravAcceleration[ zPositionIndex ]
                                + sunThirdBodyEffectAcceleration[ zPositionIndex ]
                                + solarRadiationPressureAcceleration[ zPositionIndex ];
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
    outputFile << "jacobi" << ",";
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

    double xBodyFrameVelocityInInertialComponents
                = initialStateInertial[ xVelocityIndex ] - omegaCrossPosition[ 0 ];

    double yBodyFrameVelocityInInertialComponents
                = initialStateInertial[ yVelocityIndex ] - omegaCrossPosition[ 1 ];

    double zBodyFrameVelocityInInertialComponents
                = initialStateInertial[ zVelocityIndex ] - omegaCrossPosition[ 2 ];

    initialState[ xVelocityIndex ]
            = xBodyFrameVelocityInInertialComponents * std::cos( phi )
            + yBodyFrameVelocityInInertialComponents * std::sin( phi );

    initialState[ yVelocityIndex ]
            = -1.0 * xBodyFrameVelocityInInertialComponents * std::sin( phi )
            + yBodyFrameVelocityInInertialComponents * std::cos( phi );

    initialState[ zVelocityIndex] = zBodyFrameVelocityInInertialComponents;

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

    //! get the initial jacobi integral
    double xVelocitySquare = initialState[ xVelocityIndex ] * initialState[ xVelocityIndex ];
    double yVelocitySquare = initialState[ yVelocityIndex ] * initialState[ yVelocityIndex ];
    double zVelocitySquare = initialState[ zVelocityIndex ] * initialState[ zVelocityIndex ];

    double xPositionSquare = initialState[ xPositionIndex ] * initialState[ xPositionIndex ];
    double yPositionSquare = initialState[ yPositionIndex ] * initialState[ yPositionIndex ];

    double omegaSquare = asteroidRotationVector[ zPositionIndex ] * asteroidRotationVector[ zPositionIndex ];

    double bodyFrameGravPotential;
    computeEllipsoidGravitationalPotential( alpha,
                                            beta,
                                            gamma,
                                            gravParameter,
                                            initialState[ xPositionIndex ],
                                            initialState[ yPositionIndex ],
                                            initialState[ zPositionIndex ],
                                            bodyFrameGravPotential );

    double jacobiIntegral
        = 0.5 * ( xVelocitySquare + yVelocitySquare + zVelocitySquare )
        - 0.5 * omegaSquare * ( xPositionSquare + yPositionSquare )
        - bodyFrameGravPotential;

    // save the initial state vector
    outputFile << currentStateVector[ xPositionIndex ] << ",";
    outputFile << currentStateVector[ yPositionIndex ] << ",";
    outputFile << currentStateVector[ zPositionIndex ] << ",";
    outputFile << currentStateVector[ xVelocityIndex ] << ",";
    outputFile << currentStateVector[ yVelocityIndex ] << ",";
    outputFile << currentStateVector[ zVelocityIndex ] << ",";
    outputFile << jacobiIntegral << ",";
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

        //! get the jacobi integral
        xVelocitySquare = currentStateVector[ xVelocityIndex ] * currentStateVector[ xVelocityIndex ];
        yVelocitySquare = currentStateVector[ yVelocityIndex ] * currentStateVector[ yVelocityIndex ];
        zVelocitySquare = currentStateVector[ zVelocityIndex ] * currentStateVector[ zVelocityIndex ];

        xPositionSquare = currentStateVector[ xPositionIndex ] * currentStateVector[ xPositionIndex ];
        yPositionSquare = currentStateVector[ yPositionIndex ] * currentStateVector[ yPositionIndex ];

        omegaSquare = asteroidRotationVector[ zPositionIndex ] * asteroidRotationVector[ zPositionIndex ];

        bodyFrameGravPotential;
        computeEllipsoidGravitationalPotential( alpha,
                                                beta,
                                                gamma,
                                                gravParameter,
                                                currentStateVector[ xPositionIndex ],
                                                currentStateVector[ yPositionIndex ],
                                                currentStateVector[ zPositionIndex ],
                                                bodyFrameGravPotential );

        jacobiIntegral
            = 0.5 * ( xVelocitySquare + yVelocitySquare + zVelocitySquare )
            - 0.5 * omegaSquare * ( xPositionSquare + yPositionSquare )
            - bodyFrameGravPotential;

        // save data
        outputFile << currentStateVector[ xPositionIndex ] << ",";
        outputFile << currentStateVector[ yPositionIndex ] << ",";
        outputFile << currentStateVector[ zPositionIndex ] << ",";
        outputFile << currentStateVector[ xVelocityIndex ] << ",";
        outputFile << currentStateVector[ yVelocityIndex ] << ",";
        outputFile << currentStateVector[ zVelocityIndex ] << ",";
        outputFile << jacobiIntegral << ",";
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

    // use the transport theorem to get the inertial velocity in body frame components
    double xinertialVelocityInBodyFrameComponents
                = currentStateVector[ xVelocityIndex ] + omegaCrossPosition[ 0 ];

    double yinertialVelocityInBodyFrameComponents
                = currentStateVector[ yVelocityIndex ] + omegaCrossPosition[ 1 ];

    double zinertialVelocityInBodyFrameComponents
                = currentStateVector[ zVelocityIndex ] + omegaCrossPosition[ 2 ];

    // use the rotation matrix to get the velocity in the inertial frame
    inertialState[ xVelocityIndex ]
                = xinertialVelocityInBodyFrameComponents * std::cos( rotationAngle )
                - yinertialVelocityInBodyFrameComponents * std::sin( rotationAngle );

    inertialState[ yVelocityIndex ]
                = xinertialVelocityInBodyFrameComponents * std::sin( rotationAngle )
                + yinertialVelocityInBodyFrameComponents * std::cos( rotationAngle );

    inertialState[ zVelocityIndex ] = zinertialVelocityInBodyFrameComponents;
}

//! Calculate particle energy
/*!
 * The routine calculates the energy of the orbiting particle, kinetic, potential and total kepler
 * energy.
 */
void computeParticleEnergy( const double alpha,
                            const double beta,
                            const double gamma,
                            const double gravParameter,
                            const std::vector< double > &inertialState,
                            const std::vector< double > &bodyFrameStateVector,
                            double &kineticEnergy,
                            double &potentialEnergy,
                            double &totalEnergy,
                            double &truePotentialEnergy,
                            double &trueTotalEnergy )
{
    double trueGravPotential;
    computeEllipsoidGravitationalPotential( alpha,
                                            beta,
                                            gamma,
                                            gravParameter,
                                            bodyFrameStateVector[ xPositionIndex ],
                                            bodyFrameStateVector[ yPositionIndex ],
                                            bodyFrameStateVector[ zPositionIndex ],
                                            trueGravPotential );

    double gravPotential;
    std::vector< double > positionVector = { inertialState[ xPositionIndex ],
                                             inertialState[ yPositionIndex ],
                                             inertialState[ zPositionIndex ] };
    double radialDistance = vectorNorm( positionVector );
    gravPotential = gravParameter / radialDistance;

    std::vector< double > velocityVector = { inertialState[ xVelocityIndex ],
                                             inertialState[ yVelocityIndex ],
                                             inertialState[ zVelocityIndex ] };

    double velocityMagnitude = vectorNorm( velocityVector );

    kineticEnergy = velocityMagnitude * velocityMagnitude / 2.0;
    potentialEnergy = -1.0 * gravPotential;
    totalEnergy = kineticEnergy + potentialEnergy;

    truePotentialEnergy = -1.0 * trueGravPotential;
    trueTotalEnergy = kineticEnergy + truePotentialEnergy;
}

//! check jacobi conservation
/*!
 * this routine compares current jacobi with the previous jacobi value at every time step of the
 * integration to see if the value is conserved.
 */
void jacobiChecker( const double currentJacobi,
                    const double previousJacobi,
                    const double eventTime,
                    const double xPosition,
                    const double yPosition,
                    const double zPosition,
                    const double xVelocity,
                    const double yVelocity,
                    const double zVelocity,
                    const double gravAcceleration,
                    const double gravPotential,
                    std::ostringstream &location )
{
    // check if jacobi value is same as the last computed value, it should remain constant
    if( std::fabs( currentJacobi - previousJacobi ) >= 1.0e-10 )
    {
        // something's went wrong and that is why we are here
        std::cout << "Jacobi integral not conserved!" << std::endl;
        std::cout << "Simulator at " << location.str( ) << std::endl;
        std::cout << "time = " << eventTime << std::endl << std::endl;

        double range = std::sqrt( xPosition * xPosition
                                + yPosition * yPosition
                                + zPosition * zPosition );

        std::cout << "x = " << xPosition << std::endl;
        std::cout << "y = " << yPosition << std::endl;
        std::cout << "z = " << zPosition << std::endl;
        std::cout << "range = " << range << std::endl;

        double velocity = std::sqrt( xVelocity * xVelocity
                                   + yVelocity * yVelocity
                                   + zVelocity * zVelocity );

        std::cout << "vx = " << xVelocity << std::endl;
        std::cout << "vy = " << yVelocity << std::endl;
        std::cout << "vz = " << zVelocity << std::endl;
        std::cout << "velocity = " << velocity << std::endl;

        std::cout << "grav. potential = " << gravPotential << std::endl;
        std::cout << "grav. acceleration = " << gravAcceleration << std::endl << std::endl;
    }
}

//! Computes the solar phase angle at a given time instance
double computeSolarPhaseAngle( const double regolithArgumentOfPeriapsis,
                               const double currentTime,
                               const double initialSunMeanAnomalyRadian,
                               const double initialTimeForSun,
                               const std::vector< double > initialSunOrbitalElements,
                               const double sunMeanMotion )
{
    // radian conversion
    // double regolithArgumentOfPeriapsisRadian = convertDegreeToRadians( regolithArgumentOfPeriapsis );

    // calculate the sun's longitude for the current time instance
    // double timeDifference = ( currentTime - initialTimeForSun );
    double timeDifference = ( currentTime );
    double meanAnomalyRadian = sunMeanMotion * timeDifference + initialSunMeanAnomalyRadian;

    double eccentricity = initialSunOrbitalElements[ 1 ];
    double eccentricitySquare = eccentricity * eccentricity;

    // get the corresponding eccentric anomaly
    double eccentricAnomalyRadian
        = convertMeanAnomalyToEccentricAnomaly( eccentricity,
                                                meanAnomalyRadian );

    // get the corresponding true anomaly
    double sineOfTrueAnomaly
        = std::sqrt( 1.0 - eccentricitySquare ) * std::sin( eccentricAnomalyRadian )
        / ( 1.0 - eccentricity * std::cos( eccentricAnomalyRadian ) );

    double cosineOfTrueAnomaly
        = ( std::cos( eccentricAnomalyRadian ) - eccentricity )
        / ( 1.0 - eccentricity * std::cos( eccentricAnomalyRadian ) );

    double trueAnomalyRadian
        = std::atan2( sineOfTrueAnomaly, cosineOfTrueAnomaly );

    double trueAnomaly = naos::convertRadiansToDegree( trueAnomalyRadian );

    if( trueAnomaly < 0.0 )
    {
        trueAnomaly = trueAnomaly + 360.0;
    }

    // double longitudeOfSun = trueAnomalyRadian;

    // double longitudeOfSun = naos::convertDegreeToRadians( initialSunOrbitalElements[ 5 ] );

    // double solarPhaseAngle = naos::PI - ( longitudeOfSun - regolithArgumentOfPeriapsisRadian );

    // double solarPhaseAngle = longitudeOfSun;

    double solarPhaseAngle = trueAnomaly;
    // solarPhaseAngle = std::fmod( solarPhaseAngle, 360.0 );

    return solarPhaseAngle;
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
                                                 const double localRegolithDirectionEscapeSpeed,
                                                 const double inertialDirectionalEscapeSpeed,
                                                 const double initialStepSize,
                                                 const double startTime,
                                                 const double endTime,
                                                 const double initialSunMeanAnomalyRadian,
                                                 const double initialTimeForSun,
                                                 const std::vector< double > initialSunOrbitalElements,
                                                 const double sunMeanMotion,
                                                 SQLite::Statement &databaseQuery,
                                                 const double dataSaveIntervals )
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
    double kineticEnergy = 0.0;
    double potentialEnergy = 0.0;
    double initialParticleEnergy = 0.0;
    double truePotentialEnergy = 0.0;
    double trueInitialParticleEnergy = 0.0;

    computeParticleEnergy( alpha,
                           beta,
                           gamma,
                           gravParameter,
                           initialInertialState,
                           initialState,
                           kineticEnergy,
                           potentialEnergy,
                           initialParticleEnergy,
                           truePotentialEnergy,
                           trueInitialParticleEnergy );

    //! get the initial jacobi integral
    double xVelocitySquare = initialState[ xVelocityIndex ] * initialState[ xVelocityIndex ];
    double yVelocitySquare = initialState[ yVelocityIndex ] * initialState[ yVelocityIndex ];
    double zVelocitySquare = initialState[ zVelocityIndex ] * initialState[ zVelocityIndex ];

    double xPositionSquare = initialState[ xPositionIndex ] * initialState[ xPositionIndex ];
    double yPositionSquare = initialState[ yPositionIndex ] * initialState[ yPositionIndex ];
    double zPositionSquare = initialState[ zPositionIndex ] * initialState[ zPositionIndex ];

    double omegaSquare = asteroidRotationVector[ zPositionIndex ] * asteroidRotationVector[ zPositionIndex ];

    double bodyFrameGravPotential;
    computeEllipsoidGravitationalPotential( alpha,
                                            beta,
                                            gamma,
                                            gravParameter,
                                            initialState[ xPositionIndex ],
                                            initialState[ yPositionIndex ],
                                            initialState[ zPositionIndex ],
                                            bodyFrameGravPotential );

    double jacobiIntegral
        = 0.5 * ( xVelocitySquare + yVelocitySquare + zVelocitySquare )
        - 0.5 * omegaSquare * ( xPositionSquare + yPositionSquare )
        - bodyFrameGravPotential;

    //! get the initial position and velocity magnitudes
    double initialPositionMagnitude = std::sqrt( xPositionSquare + yPositionSquare + zPositionSquare );
    double initialVelocityMagnitude = std::sqrt( xVelocitySquare + yVelocitySquare + zVelocitySquare );

    std::cout << std::endl << std::endl;
    std::cout << "Current Solar Phase angle = " << initialSunOrbitalElements[ 5 ];
    std::cout << std::endl;
    std::cout << "Current launch azimuth = " << naos::convertRadiansToDegree( launchAzimuth );
    std::cout << std::endl;
    std::cout << "Current launch declination = " << naos::convertRadiansToDegree( launchDeclination );
    std::cout << std::endl;
    std::cout << "Current launch velocity = " << initialVelocityMagnitude;

    // set up boost odeint
    const double absoluteTolerance = 1.0e-15;
    const double relativeTolerance = 1.0e-15;
    typedef boost::numeric::odeint::runge_kutta_fehlberg78< std::vector< double > > stepperType;

    // state step size guess (at each step this initial guess will be used)
    double stepSizeGuess = initialStepSize;

    // initialize the ode system
    // const double zRotation = asteroidRotationVector[ zPositionIndex ];
    // equationsOfMotionParticleAroundEllipsoid particleAroundEllipsoidProblem( gravParameter,
    //                                                                          alpha,
    //                                                                          beta,
    //                                                                          gamma,
    //                                                                          zRotation );

    // initialize the (perturbed) ode system
    // Note - initial time for sun's position (corresponding true anomaly for sun) is independant of
    // start time for regolith trajectory simulation
    perturbedEquationsOfMotionParticleAroundEllipsoid particleAroundEllipsoidProblem( gravParameter,
                                                                                      alpha,
                                                                                      beta,
                                                                                      gamma,
                                                                                      asteroidRotationVector,
                                                                                      initialTimeForSun,
                                                                                      initialSunMeanAnomalyRadian,
                                                                                      initialSunOrbitalElements,
                                                                                      sunMeanMotion );

    // compute the (initial) solar phase angle
    double initialArgumentOfPeriapsisForRegolith
        = naos::convertDegreeToRadians( initialOrbitalElements[ 4 ] );

    double initialLongitudeOfSun
        = naos::convertDegreeToRadians( initialSunOrbitalElements[ 5 ] );

    // double solarPhaseAngle
    //     = naos::PI - ( initialLongitudeOfSun - initialArgumentOfPeriapsisForRegolith );
    double solarPhaseAngle = initialLongitudeOfSun;

    solarPhaseAngle = naos::convertRadiansToDegree( solarPhaseAngle );
    solarPhaseAngle = std::fmod( solarPhaseAngle, 360.0 );
    double initialSolarPhaseAngle = solarPhaseAngle;

    // std::cout << std::endl;
    // std::cout << "solar longitude = " << initialSunOrbitalElements[ 5 ];
    // std::cout << std::endl;
    // std::cout << "regolith initial AOP = " << initialOrbitalElements[ 4 ];
    // std::cout << std::endl;
    // std::cout << "solar phase angle = " << solarPhaseAngle;
    // printVector( initialOrbitalElements, 6 );

    // initialize current state vector and time
    double variableDataSaveInterval = dataSaveIntervals;
    std::vector< double > currentStateVector = initialState;
    double currentTime = startTime;
    double intermediateEndTime = currentTime + variableDataSaveInterval;

    // initialize and set flags for crash, escape
    int escapeFlag = 0;
    int crashFlag = 0;
    int startFlag = 1;
    int endFlag = 0;

    // initialize the trajectory id for every unique regolith (in terms of its initial launch conditions)
    static int trajectoryID = 1;
    std::cout << std::endl;
    std::cout << "Trajectory ID = " << trajectoryID << std::endl;

    // get SRP components for the current time
    std::vector< double > currentSRPComponents( 3, 0.0 );
    std::vector< double > positionVectorForRegolith = { currentStateVector[ xPositionIndex ],
                                                        currentStateVector[ yPositionIndex ],
                                                        currentStateVector[ zPositionIndex ] };
    currentSRPComponents = computeSolarRadiationPressureAcceleration( positionVectorForRegolith,
                                                                      asteroidRotationVector,
                                                                      initialTimeForSun,
                                                                      initialSunMeanAnomalyRadian,
                                                                      initialSunOrbitalElements,
                                                                      currentTime,
                                                                      sunMeanMotion );

    // get solar tidal effect perturbations for the current time
    std::vector< double > currentSolarTidalEffectComponents( 3, 0.0 );
    currentSolarTidalEffectComponents
        = computeSunThirdBodyEffectAcceleration( positionVectorForRegolith,
                                                 asteroidRotationVector,
                                                 initialTimeForSun,
                                                 initialSunMeanAnomalyRadian,
                                                 initialSunOrbitalElements,
                                                 currentTime,
                                                 sunMeanMotion );

    // get gravitational acceleration for the current time
    std::vector< double > gravAcceleration( 3, 0.0 );
    computeEllipsoidGravitationalAcceleration( alpha,
                                               beta,
                                               gamma,
                                               gravParameter,
                                               currentStateVector[ xPositionIndex ],
                                               currentStateVector[ yPositionIndex ],
                                               currentStateVector[ zPositionIndex ],
                                               gravAcceleration );

    // save the initial state vector(body frame) and orbital elements
    databaseQuery.bind( ":trajectory_id", trajectoryID );

    databaseQuery.bind( ":initial_position_x", initialState[ xPositionIndex ] );
    databaseQuery.bind( ":initial_position_y", initialState[ yPositionIndex ] );
    databaseQuery.bind( ":initial_position_z", initialState[ zPositionIndex ] );
    databaseQuery.bind( ":initial_position_magnitude", initialPositionMagnitude );
    databaseQuery.bind( ":initial_velocity_x", initialState[ xVelocityIndex ] );
    databaseQuery.bind( ":initial_velocity_y", initialState[ yVelocityIndex ] );
    databaseQuery.bind( ":initial_velocity_z", initialState[ zVelocityIndex ] );
    databaseQuery.bind( ":initial_velocity_magnitude", initialVelocityMagnitude );

    databaseQuery.bind( ":initial_inertial_position_x", initialInertialState[ xPositionIndex ] );
    databaseQuery.bind( ":initial_inertial_position_y", initialInertialState[ yPositionIndex ] );
    databaseQuery.bind( ":initial_inertial_position_z", initialInertialState[ zPositionIndex ] );
    databaseQuery.bind( ":initial_inertial_velocity_x", initialInertialState[ xVelocityIndex ] );
    databaseQuery.bind( ":initial_inertial_velocity_y", initialInertialState[ yVelocityIndex ] );
    databaseQuery.bind( ":initial_inertial_velocity_z", initialInertialState[ zVelocityIndex ] );

    databaseQuery.bind( ":position_x", initialState[ xPositionIndex ] );
    databaseQuery.bind( ":position_y", initialState[ yPositionIndex ] );
    databaseQuery.bind( ":position_z", initialState[ zPositionIndex ] );
    databaseQuery.bind( ":velocity_x", initialState[ xVelocityIndex ] );
    databaseQuery.bind( ":velocity_y", initialState[ yVelocityIndex ] );
    databaseQuery.bind( ":velocity_z", initialState[ zVelocityIndex ] );
    databaseQuery.bind( ":time", currentTime );

    databaseQuery.bind( ":inertial_position_x", initialInertialState[ xPositionIndex ] );
    databaseQuery.bind( ":inertial_position_y", initialInertialState[ yPositionIndex ] );
    databaseQuery.bind( ":inertial_position_z", initialInertialState[ zPositionIndex ] );
    databaseQuery.bind( ":inertial_velocity_x", initialInertialState[ xVelocityIndex ] );
    databaseQuery.bind( ":inertial_velocity_y", initialInertialState[ yVelocityIndex ] );
    databaseQuery.bind( ":inertial_velocity_z", initialInertialState[ zVelocityIndex ] );

    databaseQuery.bind( ":launch_azimuth", convertRadiansToDegree( launchAzimuth ) );
    databaseQuery.bind( ":launch_declination", convertRadiansToDegree( launchDeclination ) );

    databaseQuery.bind( ":directional_escape_speed", localRegolithDirectionEscapeSpeed );
    databaseQuery.bind( ":directional_inertial_escape_speed", inertialDirectionalEscapeSpeed );

    databaseQuery.bind( ":sma", initialOrbitalElements[ 0 ] );
    databaseQuery.bind( ":eccentricity", initialOrbitalElements[ 1 ] );
    databaseQuery.bind( ":inclination", initialOrbitalElements[ 2 ] );
    databaseQuery.bind( ":raan", initialOrbitalElements[ 3 ] );
    databaseQuery.bind( ":aop", initialOrbitalElements[ 4 ] );
    databaseQuery.bind( ":ta", initialOrbitalElements[ 5 ] );

    databaseQuery.bind( ":kinetic_energy", kineticEnergy );
    databaseQuery.bind( ":potential_energy", potentialEnergy );
    databaseQuery.bind( ":total_energy", initialParticleEnergy );

    databaseQuery.bind( ":true_potential_energy", truePotentialEnergy );
    databaseQuery.bind( ":true_total_energy", trueInitialParticleEnergy );

    databaseQuery.bind( ":jacobi_integral", jacobiIntegral );

    databaseQuery.bind( ":initial_solar_phase_angle", initialSolarPhaseAngle );
    databaseQuery.bind( ":solar_phase_angle", solarPhaseAngle );
    databaseQuery.bind( ":srp_x", currentSRPComponents[ 0 ] );
    databaseQuery.bind( ":srp_y", currentSRPComponents[ 1 ] );
    databaseQuery.bind( ":srp_z", currentSRPComponents[ 2 ] );
    databaseQuery.bind( ":solarTide_x", currentSolarTidalEffectComponents[ 0 ] );
    databaseQuery.bind( ":solarTide_y", currentSolarTidalEffectComponents[ 1 ] );
    databaseQuery.bind( ":solarTide_z", currentSolarTidalEffectComponents[ 2 ] );

    databaseQuery.bind( ":gravAcc_x", gravAcceleration[ 0 ] );
    databaseQuery.bind( ":gravAcc_y", gravAcceleration[ 1 ] );
    databaseQuery.bind( ":gravAcc_z", gravAcceleration[ 2 ] );

    databaseQuery.bind( ":start_flag", startFlag );
    databaseQuery.bind( ":end_flag", endFlag );
    databaseQuery.bind( ":escape_flag", escapeFlag );
    databaseQuery.bind( ":crash_flag", crashFlag );

    startFlag = 0;

    // Execute insert query.
    databaseQuery.executeStep( );

    // Reset SQL insert query.
    databaseQuery.reset( );

    // start the integration outer loop
    while( intermediateEndTime <= endTime )
    {
        // save the last know state vector for when the particle is outside the asteroid
        std::vector< double > lastStateVector = currentStateVector;

        // save the last calculated jacobi integral value
        double previousJacobi = jacobiIntegral;

        // perform integration, integrated result stored in currentStateVector
        size_t steps = boost::numeric::odeint::integrate_adaptive(
                            make_controlled( absoluteTolerance, relativeTolerance, stepperType( ) ),
                            particleAroundEllipsoidProblem,
                            currentStateVector,
                            currentTime,
                            intermediateEndTime,
                            stepSizeGuess );

        double radialDistance
            = std::sqrt( currentStateVector[ xPositionIndex ] * currentStateVector[ xPositionIndex ]
                    + currentStateVector[ yPositionIndex ] * currentStateVector[ yPositionIndex ]
                    + currentStateVector[ zPositionIndex ] * currentStateVector[ zPositionIndex ] );

        // check if the particle is in an outer ellipse so as to reduce the data save interval
        double xCoordinate  = currentStateVector[ xPositionIndex ];
        double yCoordinate  = currentStateVector[ yPositionIndex ];
        double zCoordinate  = currentStateVector[ zPositionIndex ];
        double outerAlpha   = alpha + 1000.0;
        double outerBeta    = beta + 1000.0;
        double outerGamma   = gamma + 1000.0;

        double boundaryCheck
            = ( xCoordinate * xCoordinate ) / ( outerAlpha * outerAlpha )
            + ( yCoordinate * yCoordinate ) / ( outerBeta * outerBeta )
            + ( zCoordinate * zCoordinate ) / ( outerGamma * outerGamma );

        //! reduce the data save interval variable if the particle is near the asteroid
        if( boundaryCheck <= 1.0 )
        {
            // reduce the data save step size
            variableDataSaveInterval = 0.1;
        }
        else if( radialDistance >= 5.0 * alpha )
        {
            // increase the data save step size
            variableDataSaveInterval = 500.0;
        }
        else
        {
            // use the original high value for data save intervals
            variableDataSaveInterval = dataSaveIntervals;
        }

        //! check if the particle is on an escape trajectory or not
        // check for energy if the particle is far away from the asteroid
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
            double particleEnergy = 0.0;
            double trueParticleEnergy = 0.0;
            computeParticleEnergy( alpha,
                                   beta,
                                   gamma,
                                   gravParameter,
                                   inertialState,
                                   currentStateVector,
                                   kineticEnergy,
                                   potentialEnergy,
                                   particleEnergy,
                                   truePotentialEnergy,
                                   trueParticleEnergy );

            bool energyFlag = false;
            if( particleEnergy > 0.0 )
            {
                //! sanity check 1 for energy calculation for escape condition
                // double sanityCheckEnergy = -1.0 * gravParameter / ( 2.0 * orbitalElements[ 0 ] );
                // std::cout << std::endl << "Sanity check for energy:" << std::endl;
                // std::cout << "sanity Check Energy value = " << sanityCheckEnergy << std::endl;
                // std::cout << "computed energy = " << particleEnergy << std::endl;
                // std::cout << "eccentricity = " << orbitalElements[ 1 ] << std::endl;
                // std::cout << "semi-major axis = " << orbitalElements[ 0 ] << std::endl << std::endl;
                // if( orbitalElements[ 0 ] > 0 )
                // {
                //     std::ostringstream errorMessage;
                //     errorMessage << std::endl;
                //     errorMessage << "positive energy but semi-major axis is not negative";
                //     errorMessage << std::endl;
                //     throw std::runtime_error( errorMessage.str( ) );
                // }

                energyFlag = true;
            }

            // check the eccentricity and energy flags
            if( eccentricityFlag && energyFlag )
            {
                //! get the jacobi integral
                xVelocitySquare = currentStateVector[ xVelocityIndex ] * currentStateVector[ xVelocityIndex ];
                yVelocitySquare = currentStateVector[ yVelocityIndex ] * currentStateVector[ yVelocityIndex ];
                zVelocitySquare = currentStateVector[ zVelocityIndex ] * currentStateVector[ zVelocityIndex ];

                xPositionSquare = currentStateVector[ xPositionIndex ] * currentStateVector[ xPositionIndex ];
                yPositionSquare = currentStateVector[ yPositionIndex ] * currentStateVector[ yPositionIndex ];

                omegaSquare = asteroidRotationVector[ zPositionIndex ] * asteroidRotationVector[ zPositionIndex ];

                computeEllipsoidGravitationalPotential( alpha,
                                                        beta,
                                                        gamma,
                                                        gravParameter,
                                                        currentStateVector[ xPositionIndex ],
                                                        currentStateVector[ yPositionIndex ],
                                                        currentStateVector[ zPositionIndex ],
                                                        bodyFrameGravPotential );

                jacobiIntegral
                    = 0.5 * ( xVelocitySquare + yVelocitySquare + zVelocitySquare )
                    - 0.5 * omegaSquare * ( xPositionSquare + yPositionSquare )
                    - bodyFrameGravPotential;

                // particle is on an escape trajectory, save data and stop integration
                // update the time variables
                currentTime = intermediateEndTime;

                // set the flags
                escapeFlag = 1;
                crashFlag = 0;

                // compute the current solar phase angle
                solarPhaseAngle = computeSolarPhaseAngle( orbitalElements[ 4 ],
                                                          currentTime,
                                                          initialSunMeanAnomalyRadian,
                                                          initialTimeForSun,
                                                          initialSunOrbitalElements,
                                                          sunMeanMotion );

                // get SRP components for the current time
                positionVectorForRegolith = { currentStateVector[ xPositionIndex ],
                                              currentStateVector[ yPositionIndex ],
                                              currentStateVector[ zPositionIndex ] };
                currentSRPComponents
                    = computeSolarRadiationPressureAcceleration( positionVectorForRegolith,
                                                                 asteroidRotationVector,
                                                                 initialTimeForSun,
                                                                 initialSunMeanAnomalyRadian,
                                                                 initialSunOrbitalElements,
                                                                 currentTime,
                                                                 sunMeanMotion );

                // get solar tidal effect perturbations for the current time
                currentSolarTidalEffectComponents
                    = computeSunThirdBodyEffectAcceleration( positionVectorForRegolith,
                                                             asteroidRotationVector,
                                                             initialTimeForSun,
                                                             initialSunMeanAnomalyRadian,
                                                             initialSunOrbitalElements,
                                                             currentTime,
                                                             sunMeanMotion );

                // get gravitational acceleration for the current time
                computeEllipsoidGravitationalAcceleration( alpha,
                                                           beta,
                                                           gamma,
                                                           gravParameter,
                                                           currentStateVector[ xPositionIndex ],
                                                           currentStateVector[ yPositionIndex ],
                                                           currentStateVector[ zPositionIndex ],
                                                           gravAcceleration );

                // save data
                databaseQuery.bind( ":trajectory_id", trajectoryID );

                databaseQuery.bind( ":initial_position_x", initialState[ xPositionIndex ] );
                databaseQuery.bind( ":initial_position_y", initialState[ yPositionIndex ] );
                databaseQuery.bind( ":initial_position_z", initialState[ zPositionIndex ] );
                databaseQuery.bind( ":initial_position_magnitude", initialPositionMagnitude );
                databaseQuery.bind( ":initial_velocity_x", initialState[ xVelocityIndex ] );
                databaseQuery.bind( ":initial_velocity_y", initialState[ yVelocityIndex ] );
                databaseQuery.bind( ":initial_velocity_z", initialState[ zVelocityIndex ] );
                databaseQuery.bind( ":initial_velocity_magnitude", initialVelocityMagnitude );
                databaseQuery.bind( ":initial_inertial_position_x", initialInertialState[ xPositionIndex ] );
                databaseQuery.bind( ":initial_inertial_position_y", initialInertialState[ yPositionIndex ] );
                databaseQuery.bind( ":initial_inertial_position_z", initialInertialState[ zPositionIndex ] );
                databaseQuery.bind( ":initial_inertial_velocity_x", initialInertialState[ xVelocityIndex ] );
                databaseQuery.bind( ":initial_inertial_velocity_y", initialInertialState[ yVelocityIndex ] );
                databaseQuery.bind( ":initial_inertial_velocity_z", initialInertialState[ zVelocityIndex ] );

                databaseQuery.bind( ":position_x", currentStateVector[ xPositionIndex ] );
                databaseQuery.bind( ":position_y", currentStateVector[ yPositionIndex ] );
                databaseQuery.bind( ":position_z", currentStateVector[ zPositionIndex ] );
                databaseQuery.bind( ":velocity_x", currentStateVector[ xVelocityIndex ] );
                databaseQuery.bind( ":velocity_y", currentStateVector[ yVelocityIndex ] );
                databaseQuery.bind( ":velocity_z", currentStateVector[ zVelocityIndex ] );
                databaseQuery.bind( ":time", currentTime );

                databaseQuery.bind( ":inertial_position_x", inertialState[ xPositionIndex ] );
                databaseQuery.bind( ":inertial_position_y", inertialState[ yPositionIndex ] );
                databaseQuery.bind( ":inertial_position_z", inertialState[ zPositionIndex ] );
                databaseQuery.bind( ":inertial_velocity_x", inertialState[ xVelocityIndex ] );
                databaseQuery.bind( ":inertial_velocity_y", inertialState[ yVelocityIndex ] );
                databaseQuery.bind( ":inertial_velocity_z", inertialState[ zVelocityIndex ] );

                databaseQuery.bind( ":launch_azimuth", convertRadiansToDegree( launchAzimuth ) );
                databaseQuery.bind( ":launch_declination", convertRadiansToDegree( launchDeclination ) );

                databaseQuery.bind( ":directional_escape_speed", localRegolithDirectionEscapeSpeed );
                databaseQuery.bind( ":directional_inertial_escape_speed", inertialDirectionalEscapeSpeed );

                databaseQuery.bind( ":sma", orbitalElements[ 0 ] );
                databaseQuery.bind( ":eccentricity", orbitalElements[ 1 ] );
                databaseQuery.bind( ":inclination", orbitalElements[ 2 ] );
                databaseQuery.bind( ":raan", orbitalElements[ 3 ] );
                databaseQuery.bind( ":aop", orbitalElements[ 4 ] );
                databaseQuery.bind( ":ta", orbitalElements[ 5 ] );

                databaseQuery.bind( ":kinetic_energy", kineticEnergy );
                databaseQuery.bind( ":potential_energy", potentialEnergy );
                databaseQuery.bind( ":total_energy", particleEnergy );

                databaseQuery.bind( ":true_potential_energy", truePotentialEnergy );
                databaseQuery.bind( ":true_total_energy", trueParticleEnergy );

                databaseQuery.bind( ":jacobi_integral", jacobiIntegral );

                databaseQuery.bind( ":initial_solar_phase_angle", initialSolarPhaseAngle );
                databaseQuery.bind( ":solar_phase_angle", solarPhaseAngle );
                databaseQuery.bind( ":srp_x", currentSRPComponents[ 0 ] );
                databaseQuery.bind( ":srp_y", currentSRPComponents[ 1 ] );
                databaseQuery.bind( ":srp_z", currentSRPComponents[ 2 ] );
                databaseQuery.bind( ":solarTide_x", currentSolarTidalEffectComponents[ 0 ] );
                databaseQuery.bind( ":solarTide_y", currentSolarTidalEffectComponents[ 1 ] );
                databaseQuery.bind( ":solarTide_z", currentSolarTidalEffectComponents[ 2 ] );

                databaseQuery.bind( ":gravAcc_x", gravAcceleration[ 0 ] );
                databaseQuery.bind( ":gravAcc_y", gravAcceleration[ 1 ] );
                databaseQuery.bind( ":gravAcc_z", gravAcceleration[ 2 ] );

                databaseQuery.bind( ":start_flag", startFlag );
                databaseQuery.bind( ":end_flag", endFlag );
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

            double particleEnergy = 0.0;
            double trueParticleEnergy = 0.0;
            computeParticleEnergy( alpha,
                                   beta,
                                   gamma,
                                   gravParameter,
                                   inertialState,
                                   currentStateVector,
                                   kineticEnergy,
                                   potentialEnergy,
                                   particleEnergy,
                                   truePotentialEnergy,
                                   trueParticleEnergy );

            //! get the jacobi integral
            xVelocitySquare = currentStateVector[ xVelocityIndex ] * currentStateVector[ xVelocityIndex ];
            yVelocitySquare = currentStateVector[ yVelocityIndex ] * currentStateVector[ yVelocityIndex ];
            zVelocitySquare = currentStateVector[ zVelocityIndex ] * currentStateVector[ zVelocityIndex ];

            xPositionSquare = currentStateVector[ xPositionIndex ] * currentStateVector[ xPositionIndex ];
            yPositionSquare = currentStateVector[ yPositionIndex ] * currentStateVector[ yPositionIndex ];

            omegaSquare = asteroidRotationVector[ zPositionIndex ] * asteroidRotationVector[ zPositionIndex ];

            computeEllipsoidGravitationalPotential( alpha,
                                                    beta,
                                                    gamma,
                                                    gravParameter,
                                                    currentStateVector[ xPositionIndex ],
                                                    currentStateVector[ yPositionIndex ],
                                                    currentStateVector[ zPositionIndex ],
                                                    bodyFrameGravPotential );

            jacobiIntegral
                = 0.5 * ( xVelocitySquare + yVelocitySquare + zVelocitySquare )
                - 0.5 * omegaSquare * ( xPositionSquare + yPositionSquare )
                - bodyFrameGravPotential;

            // compute solar phase angle
            solarPhaseAngle = computeSolarPhaseAngle( orbitalElements[ 4 ],
                                                      currentTime,
                                                      initialSunMeanAnomalyRadian,
                                                      initialTimeForSun,
                                                      initialSunOrbitalElements,
                                                      sunMeanMotion );

            // get SRP components for the current time
            positionVectorForRegolith = { currentStateVector[ xPositionIndex ],
                                          currentStateVector[ yPositionIndex ],
                                          currentStateVector[ zPositionIndex ] };
            currentSRPComponents
                = computeSolarRadiationPressureAcceleration( positionVectorForRegolith,
                                                             asteroidRotationVector,
                                                             initialTimeForSun,
                                                             initialSunMeanAnomalyRadian,
                                                             initialSunOrbitalElements,
                                                             currentTime,
                                                             sunMeanMotion );

            // get solar tidal effect perturbations for the current time
            currentSolarTidalEffectComponents
                = computeSunThirdBodyEffectAcceleration( positionVectorForRegolith,
                                                         asteroidRotationVector,
                                                         initialTimeForSun,
                                                         initialSunMeanAnomalyRadian,
                                                         initialSunOrbitalElements,
                                                         currentTime,
                                                         sunMeanMotion );

            // get gravitational acceleration for the current time
            computeEllipsoidGravitationalAcceleration( alpha,
                                                       beta,
                                                       gamma,
                                                       gravParameter,
                                                       currentStateVector[ xPositionIndex ],
                                                       currentStateVector[ yPositionIndex ],
                                                       currentStateVector[ zPositionIndex ],
                                                       gravAcceleration );

            // save data
            databaseQuery.bind( ":trajectory_id", trajectoryID );

            databaseQuery.bind( ":initial_position_x", initialState[ xPositionIndex ] );
            databaseQuery.bind( ":initial_position_y", initialState[ yPositionIndex ] );
            databaseQuery.bind( ":initial_position_z", initialState[ zPositionIndex ] );
            databaseQuery.bind( ":initial_position_magnitude", initialPositionMagnitude );
            databaseQuery.bind( ":initial_velocity_x", initialState[ xVelocityIndex ] );
            databaseQuery.bind( ":initial_velocity_y", initialState[ yVelocityIndex ] );
            databaseQuery.bind( ":initial_velocity_z", initialState[ zVelocityIndex ] );
            databaseQuery.bind( ":initial_velocity_magnitude", initialVelocityMagnitude );
            databaseQuery.bind( ":initial_inertial_position_x", initialInertialState[ xPositionIndex ] );
            databaseQuery.bind( ":initial_inertial_position_y", initialInertialState[ yPositionIndex ] );
            databaseQuery.bind( ":initial_inertial_position_z", initialInertialState[ zPositionIndex ] );
            databaseQuery.bind( ":initial_inertial_velocity_x", initialInertialState[ xVelocityIndex ] );
            databaseQuery.bind( ":initial_inertial_velocity_y", initialInertialState[ yVelocityIndex ] );
            databaseQuery.bind( ":initial_inertial_velocity_z", initialInertialState[ zVelocityIndex ] );

            databaseQuery.bind( ":position_x", currentStateVector[ xPositionIndex ] );
            databaseQuery.bind( ":position_y", currentStateVector[ yPositionIndex ] );
            databaseQuery.bind( ":position_z", currentStateVector[ zPositionIndex ] );
            databaseQuery.bind( ":velocity_x", currentStateVector[ xVelocityIndex ] );
            databaseQuery.bind( ":velocity_y", currentStateVector[ yVelocityIndex ] );
            databaseQuery.bind( ":velocity_z", currentStateVector[ zVelocityIndex ] );
            databaseQuery.bind( ":time", currentTime );

            databaseQuery.bind( ":inertial_position_x", inertialState[ xPositionIndex ] );
            databaseQuery.bind( ":inertial_position_y", inertialState[ yPositionIndex ] );
            databaseQuery.bind( ":inertial_position_z", inertialState[ zPositionIndex ] );
            databaseQuery.bind( ":inertial_velocity_x", inertialState[ xVelocityIndex ] );
            databaseQuery.bind( ":inertial_velocity_y", inertialState[ yVelocityIndex ] );
            databaseQuery.bind( ":inertial_velocity_z", inertialState[ zVelocityIndex ] );

            databaseQuery.bind( ":launch_azimuth", convertRadiansToDegree( launchAzimuth ) );
            databaseQuery.bind( ":launch_declination", convertRadiansToDegree( launchDeclination ) );

            databaseQuery.bind( ":directional_escape_speed", localRegolithDirectionEscapeSpeed );
            databaseQuery.bind( ":directional_inertial_escape_speed", inertialDirectionalEscapeSpeed );

            databaseQuery.bind( ":sma", orbitalElements[ 0 ] );
            databaseQuery.bind( ":eccentricity", orbitalElements[ 1 ] );
            databaseQuery.bind( ":inclination", orbitalElements[ 2 ] );
            databaseQuery.bind( ":raan", orbitalElements[ 3 ] );
            databaseQuery.bind( ":aop", orbitalElements[ 4 ] );
            databaseQuery.bind( ":ta", orbitalElements[ 5 ] );

            databaseQuery.bind( ":kinetic_energy", kineticEnergy );
            databaseQuery.bind( ":potential_energy", potentialEnergy );
            databaseQuery.bind( ":total_energy", particleEnergy );

            databaseQuery.bind( ":true_potential_energy", truePotentialEnergy );
            databaseQuery.bind( ":true_total_energy", trueParticleEnergy );

            databaseQuery.bind( ":jacobi_integral", jacobiIntegral );

            databaseQuery.bind( ":initial_solar_phase_angle", initialSolarPhaseAngle );
            databaseQuery.bind( ":solar_phase_angle", solarPhaseAngle );
            databaseQuery.bind( ":srp_x", currentSRPComponents[ 0 ] );
            databaseQuery.bind( ":srp_y", currentSRPComponents[ 1 ] );
            databaseQuery.bind( ":srp_z", currentSRPComponents[ 2 ] );
            databaseQuery.bind( ":solarTide_x", currentSolarTidalEffectComponents[ 0 ] );
            databaseQuery.bind( ":solarTide_y", currentSolarTidalEffectComponents[ 1 ] );
            databaseQuery.bind( ":solarTide_z", currentSolarTidalEffectComponents[ 2 ] );

            databaseQuery.bind( ":gravAcc_x", gravAcceleration[ 0 ] );
            databaseQuery.bind( ":gravAcc_y", gravAcceleration[ 1 ] );
            databaseQuery.bind( ":gravAcc_z", gravAcceleration[ 2 ] );

            databaseQuery.bind( ":start_flag", startFlag );
            databaseQuery.bind( ":end_flag", endFlag );
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

            double particleEnergy = 0.0;
            double trueParticleEnergy = 0.0;
            computeParticleEnergy( alpha,
                                   beta,
                                   gamma,
                                   gravParameter,
                                   inertialState,
                                   currentStateVector,
                                   kineticEnergy,
                                   potentialEnergy,
                                   particleEnergy,
                                   truePotentialEnergy,
                                   trueParticleEnergy );

            //! get the jacobi integral
            xVelocitySquare = currentStateVector[ xVelocityIndex ] * currentStateVector[ xVelocityIndex ];
            yVelocitySquare = currentStateVector[ yVelocityIndex ] * currentStateVector[ yVelocityIndex ];
            zVelocitySquare = currentStateVector[ zVelocityIndex ] * currentStateVector[ zVelocityIndex ];

            xPositionSquare = currentStateVector[ xPositionIndex ] * currentStateVector[ xPositionIndex ];
            yPositionSquare = currentStateVector[ yPositionIndex ] * currentStateVector[ yPositionIndex ];

            omegaSquare = asteroidRotationVector[ zPositionIndex ] * asteroidRotationVector[ zPositionIndex ];

            computeEllipsoidGravitationalPotential( alpha,
                                                    beta,
                                                    gamma,
                                                    gravParameter,
                                                    currentStateVector[ xPositionIndex ],
                                                    currentStateVector[ yPositionIndex ],
                                                    currentStateVector[ zPositionIndex ],
                                                    bodyFrameGravPotential );

            jacobiIntegral
                = 0.5 * ( xVelocitySquare + yVelocitySquare + zVelocitySquare )
                - 0.5 * omegaSquare * ( xPositionSquare + yPositionSquare )
                - bodyFrameGravPotential;

            // compute solar phase angle
            solarPhaseAngle = computeSolarPhaseAngle( orbitalElements[ 4 ],
                                                      currentTime,
                                                      initialSunMeanAnomalyRadian,
                                                      initialTimeForSun,
                                                      initialSunOrbitalElements,
                                                      sunMeanMotion );

            // get SRP components for the current time
            positionVectorForRegolith = { currentStateVector[ xPositionIndex ],
                                          currentStateVector[ yPositionIndex ],
                                          currentStateVector[ zPositionIndex ] };
            currentSRPComponents
                = computeSolarRadiationPressureAcceleration( positionVectorForRegolith,
                                                             asteroidRotationVector,
                                                             initialTimeForSun,
                                                             initialSunMeanAnomalyRadian,
                                                             initialSunOrbitalElements,
                                                             currentTime,
                                                             sunMeanMotion );

            // get solar tidal effect perturbations for the current time
            currentSolarTidalEffectComponents
                = computeSunThirdBodyEffectAcceleration( positionVectorForRegolith,
                                                         asteroidRotationVector,
                                                         initialTimeForSun,
                                                         initialSunMeanAnomalyRadian,
                                                         initialSunOrbitalElements,
                                                         currentTime,
                                                         sunMeanMotion );

            // get gravitational acceleration for the current time
            computeEllipsoidGravitationalAcceleration( alpha,
                                                       beta,
                                                       gamma,
                                                       gravParameter,
                                                       currentStateVector[ xPositionIndex ],
                                                       currentStateVector[ yPositionIndex ],
                                                       currentStateVector[ zPositionIndex ],
                                                       gravAcceleration );

            // save data
            databaseQuery.bind( ":trajectory_id", trajectoryID );

            databaseQuery.bind( ":initial_position_x", initialState[ xPositionIndex ] );
            databaseQuery.bind( ":initial_position_y", initialState[ yPositionIndex ] );
            databaseQuery.bind( ":initial_position_z", initialState[ zPositionIndex ] );
            databaseQuery.bind( ":initial_position_magnitude", initialPositionMagnitude );
            databaseQuery.bind( ":initial_velocity_x", initialState[ xVelocityIndex ] );
            databaseQuery.bind( ":initial_velocity_y", initialState[ yVelocityIndex ] );
            databaseQuery.bind( ":initial_velocity_z", initialState[ zVelocityIndex ] );
            databaseQuery.bind( ":initial_velocity_magnitude", initialVelocityMagnitude );
            databaseQuery.bind( ":initial_inertial_position_x", initialInertialState[ xPositionIndex ] );
            databaseQuery.bind( ":initial_inertial_position_y", initialInertialState[ yPositionIndex ] );
            databaseQuery.bind( ":initial_inertial_position_z", initialInertialState[ zPositionIndex ] );
            databaseQuery.bind( ":initial_inertial_velocity_x", initialInertialState[ xVelocityIndex ] );
            databaseQuery.bind( ":initial_inertial_velocity_y", initialInertialState[ yVelocityIndex ] );
            databaseQuery.bind( ":initial_inertial_velocity_z", initialInertialState[ zVelocityIndex ] );

            databaseQuery.bind( ":position_x", currentStateVector[ xPositionIndex ] );
            databaseQuery.bind( ":position_y", currentStateVector[ yPositionIndex ] );
            databaseQuery.bind( ":position_z", currentStateVector[ zPositionIndex ] );
            databaseQuery.bind( ":velocity_x", currentStateVector[ xVelocityIndex ] );
            databaseQuery.bind( ":velocity_y", currentStateVector[ yVelocityIndex ] );
            databaseQuery.bind( ":velocity_z", currentStateVector[ zVelocityIndex ] );
            databaseQuery.bind( ":time", currentTime );

            databaseQuery.bind( ":inertial_position_x", inertialState[ xPositionIndex ] );
            databaseQuery.bind( ":inertial_position_y", inertialState[ yPositionIndex ] );
            databaseQuery.bind( ":inertial_position_z", inertialState[ zPositionIndex ] );
            databaseQuery.bind( ":inertial_velocity_x", inertialState[ xVelocityIndex ] );
            databaseQuery.bind( ":inertial_velocity_y", inertialState[ yVelocityIndex ] );
            databaseQuery.bind( ":inertial_velocity_z", inertialState[ zVelocityIndex ] );

            databaseQuery.bind( ":launch_azimuth", convertRadiansToDegree( launchAzimuth ) );
            databaseQuery.bind( ":launch_declination", convertRadiansToDegree( launchDeclination ) );

            databaseQuery.bind( ":directional_escape_speed", localRegolithDirectionEscapeSpeed );
            databaseQuery.bind( ":directional_inertial_escape_speed", inertialDirectionalEscapeSpeed );

            databaseQuery.bind( ":sma", orbitalElements[ 0 ] );
            databaseQuery.bind( ":eccentricity", orbitalElements[ 1 ] );
            databaseQuery.bind( ":inclination", orbitalElements[ 2 ] );
            databaseQuery.bind( ":raan", orbitalElements[ 3 ] );
            databaseQuery.bind( ":aop", orbitalElements[ 4 ] );
            databaseQuery.bind( ":ta", orbitalElements[ 5 ] );

            databaseQuery.bind( ":kinetic_energy", kineticEnergy );
            databaseQuery.bind( ":potential_energy", potentialEnergy );
            databaseQuery.bind( ":total_energy", particleEnergy );

            databaseQuery.bind( ":true_potential_energy", truePotentialEnergy );
            databaseQuery.bind( ":true_total_energy", trueParticleEnergy );

            databaseQuery.bind( ":jacobi_integral", jacobiIntegral );

            databaseQuery.bind( ":initial_solar_phase_angle", initialSolarPhaseAngle );
            databaseQuery.bind( ":solar_phase_angle", solarPhaseAngle );
            databaseQuery.bind( ":srp_x", currentSRPComponents[ 0 ] );
            databaseQuery.bind( ":srp_y", currentSRPComponents[ 1 ] );
            databaseQuery.bind( ":srp_z", currentSRPComponents[ 2 ] );
            databaseQuery.bind( ":solarTide_x", currentSolarTidalEffectComponents[ 0 ] );
            databaseQuery.bind( ":solarTide_y", currentSolarTidalEffectComponents[ 1 ] );
            databaseQuery.bind( ":solarTide_z", currentSolarTidalEffectComponents[ 2 ] );

            databaseQuery.bind( ":gravAcc_x", gravAcceleration[ 0 ] );
            databaseQuery.bind( ":gravAcc_y", gravAcceleration[ 1 ] );
            databaseQuery.bind( ":gravAcc_z", gravAcceleration[ 2 ] );

            databaseQuery.bind( ":start_flag", startFlag );
            databaseQuery.bind( ":end_flag", endFlag );
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
        intermediateEndTime = currentTime + variableDataSaveInterval;

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

        double particleEnergy = 0.0;
        double trueParticleEnergy = 0.0;
        computeParticleEnergy( alpha,
                               beta,
                               gamma,
                               gravParameter,
                               inertialState,
                               currentStateVector,
                               kineticEnergy,
                               potentialEnergy,
                               particleEnergy,
                               truePotentialEnergy,
                               trueParticleEnergy );

        //! get the jacobi integral
        xVelocitySquare = currentStateVector[ xVelocityIndex ] * currentStateVector[ xVelocityIndex ];
        yVelocitySquare = currentStateVector[ yVelocityIndex ] * currentStateVector[ yVelocityIndex ];
        zVelocitySquare = currentStateVector[ zVelocityIndex ] * currentStateVector[ zVelocityIndex ];

        xPositionSquare = currentStateVector[ xPositionIndex ] * currentStateVector[ xPositionIndex ];
        yPositionSquare = currentStateVector[ yPositionIndex ] * currentStateVector[ yPositionIndex ];

        omegaSquare = asteroidRotationVector[ zPositionIndex ] * asteroidRotationVector[ zPositionIndex ];

        computeEllipsoidGravitationalPotential( alpha,
                                                beta,
                                                gamma,
                                                gravParameter,
                                                currentStateVector[ xPositionIndex ],
                                                currentStateVector[ yPositionIndex ],
                                                currentStateVector[ zPositionIndex ],
                                                bodyFrameGravPotential );

        jacobiIntegral
            = 0.5 * ( xVelocitySquare + yVelocitySquare + zVelocitySquare )
            - 0.5 * omegaSquare * ( xPositionSquare + yPositionSquare )
            - bodyFrameGravPotential;

        // compute solar phase angle
        solarPhaseAngle = computeSolarPhaseAngle( orbitalElements[ 4 ],
                                                  currentTime,
                                                  initialSunMeanAnomalyRadian,
                                                  initialTimeForSun,
                                                  initialSunOrbitalElements,
                                                  sunMeanMotion );

        // get SRP components for the current time
        positionVectorForRegolith = { currentStateVector[ xPositionIndex ],
                                      currentStateVector[ yPositionIndex ],
                                      currentStateVector[ zPositionIndex ] };
        currentSRPComponents
            = computeSolarRadiationPressureAcceleration( positionVectorForRegolith,
                                                         asteroidRotationVector,
                                                         initialTimeForSun,
                                                         initialSunMeanAnomalyRadian,
                                                         initialSunOrbitalElements,
                                                         currentTime,
                                                         sunMeanMotion );

        // get solar tidal effect perturbations for the current time
        currentSolarTidalEffectComponents
            = computeSunThirdBodyEffectAcceleration( positionVectorForRegolith,
                                                     asteroidRotationVector,
                                                     initialTimeForSun,
                                                     initialSunMeanAnomalyRadian,
                                                     initialSunOrbitalElements,
                                                     currentTime,
                                                     sunMeanMotion );

        // get gravitational acceleration for the current time
        computeEllipsoidGravitationalAcceleration( alpha,
                                                   beta,
                                                   gamma,
                                                   gravParameter,
                                                   currentStateVector[ xPositionIndex ],
                                                   currentStateVector[ yPositionIndex ],
                                                   currentStateVector[ zPositionIndex ],
                                                   gravAcceleration );

        // save data
        databaseQuery.bind( ":trajectory_id", trajectoryID );

        databaseQuery.bind( ":initial_position_x", initialState[ xPositionIndex ] );
        databaseQuery.bind( ":initial_position_y", initialState[ yPositionIndex ] );
        databaseQuery.bind( ":initial_position_z", initialState[ zPositionIndex ] );
        databaseQuery.bind( ":initial_position_magnitude", initialPositionMagnitude );
        databaseQuery.bind( ":initial_velocity_x", initialState[ xVelocityIndex ] );
        databaseQuery.bind( ":initial_velocity_y", initialState[ yVelocityIndex ] );
        databaseQuery.bind( ":initial_velocity_z", initialState[ zVelocityIndex ] );
        databaseQuery.bind( ":initial_velocity_magnitude", initialVelocityMagnitude );
        databaseQuery.bind( ":initial_inertial_position_x", initialInertialState[ xPositionIndex ] );
        databaseQuery.bind( ":initial_inertial_position_y", initialInertialState[ yPositionIndex ] );
        databaseQuery.bind( ":initial_inertial_position_z", initialInertialState[ zPositionIndex ] );
        databaseQuery.bind( ":initial_inertial_velocity_x", initialInertialState[ xVelocityIndex ] );
        databaseQuery.bind( ":initial_inertial_velocity_y", initialInertialState[ yVelocityIndex ] );
        databaseQuery.bind( ":initial_inertial_velocity_z", initialInertialState[ zVelocityIndex ] );

        databaseQuery.bind( ":position_x", currentStateVector[ xPositionIndex ] );
        databaseQuery.bind( ":position_y", currentStateVector[ yPositionIndex ] );
        databaseQuery.bind( ":position_z", currentStateVector[ zPositionIndex ] );
        databaseQuery.bind( ":velocity_x", currentStateVector[ xVelocityIndex ] );
        databaseQuery.bind( ":velocity_y", currentStateVector[ yVelocityIndex ] );
        databaseQuery.bind( ":velocity_z", currentStateVector[ zVelocityIndex ] );
        databaseQuery.bind( ":time", currentTime );

        databaseQuery.bind( ":inertial_position_x", inertialState[ xPositionIndex ] );
        databaseQuery.bind( ":inertial_position_y", inertialState[ yPositionIndex ] );
        databaseQuery.bind( ":inertial_position_z", inertialState[ zPositionIndex ] );
        databaseQuery.bind( ":inertial_velocity_x", inertialState[ xVelocityIndex ] );
        databaseQuery.bind( ":inertial_velocity_y", inertialState[ yVelocityIndex ] );
        databaseQuery.bind( ":inertial_velocity_z", inertialState[ zVelocityIndex ] );

        databaseQuery.bind( ":launch_azimuth", convertRadiansToDegree( launchAzimuth ) );
        databaseQuery.bind( ":launch_declination", convertRadiansToDegree( launchDeclination ) );

        databaseQuery.bind( ":directional_escape_speed", localRegolithDirectionEscapeSpeed );
        databaseQuery.bind( ":directional_inertial_escape_speed", inertialDirectionalEscapeSpeed );

        databaseQuery.bind( ":sma", orbitalElements[ 0 ] );
        databaseQuery.bind( ":eccentricity", orbitalElements[ 1 ] );
        databaseQuery.bind( ":inclination", orbitalElements[ 2 ] );
        databaseQuery.bind( ":raan", orbitalElements[ 3 ] );
        databaseQuery.bind( ":aop", orbitalElements[ 4 ] );
        databaseQuery.bind( ":ta", orbitalElements[ 5 ] );

        databaseQuery.bind( ":kinetic_energy", kineticEnergy );
        databaseQuery.bind( ":potential_energy", potentialEnergy );
        databaseQuery.bind( ":total_energy", particleEnergy );

        databaseQuery.bind( ":true_potential_energy", truePotentialEnergy );
        databaseQuery.bind( ":true_total_energy", trueParticleEnergy );

        databaseQuery.bind( ":jacobi_integral", jacobiIntegral );

        databaseQuery.bind( ":initial_solar_phase_angle", initialSolarPhaseAngle );
        databaseQuery.bind( ":solar_phase_angle", solarPhaseAngle );
        databaseQuery.bind( ":srp_x", currentSRPComponents[ 0 ] );
        databaseQuery.bind( ":srp_y", currentSRPComponents[ 1 ] );
        databaseQuery.bind( ":srp_z", currentSRPComponents[ 2 ] );
        databaseQuery.bind( ":solarTide_x", currentSolarTidalEffectComponents[ 0 ] );
        databaseQuery.bind( ":solarTide_y", currentSolarTidalEffectComponents[ 1 ] );
        databaseQuery.bind( ":solarTide_z", currentSolarTidalEffectComponents[ 2 ] );

        databaseQuery.bind( ":gravAcc_x", gravAcceleration[ 0 ] );
        databaseQuery.bind( ":gravAcc_y", gravAcceleration[ 1 ] );
        databaseQuery.bind( ":gravAcc_z", gravAcceleration[ 2 ] );

        databaseQuery.bind( ":start_flag", startFlag );
        databaseQuery.bind( ":escape_flag", escapeFlag );
        databaseQuery.bind( ":crash_flag", crashFlag );

        if( currentTime == endTime || intermediateEndTime >= endTime )
        {
            endFlag = 1;
            databaseQuery.bind( ":end_flag", endFlag );
        }
        else
        {
            endFlag = 0;
            databaseQuery.bind( ":end_flag", endFlag );
        }

        // Execute insert query.
        databaseQuery.executeStep( );

        // Reset SQL insert query.
        databaseQuery.reset( );
    } // end of outer while loop for integration

    // sanity check for trajectory id increment
    // std::cout << std::endl;
    // std::cout << "Trajectory ID = " << trajectoryID << std::endl;
    trajectoryID++;
}

} // namespace naos
