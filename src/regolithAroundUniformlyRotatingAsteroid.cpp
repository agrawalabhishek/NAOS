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
class regolithAroundAsteroidNoPerturbations
{
    // declare parameters, gravitational parameter, asteroid rotation and the semi major axes
    // of the ellipsoid
    const double gravParameter;
    const double alpha;
    const double beta;
    const double gamma;
    const double zRotation;

public:
    // Default constructor with member initializer list
    regolithAroundAsteroidNoPerturbations(
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

        // point mass gravitational acceleration
        // const std::vector< double > positionVector = { stateVector[ xPositionIndex ],
        //                                                stateVector[ yPositionIndex ],
        //                                                stateVector[ zPositionIndex ] };
        // const double range = vectorNorm( positionVector );
        // double rangeCube = range * range * range;
        // gravAcceleration[ 0 ] = -1.0 * gravParameter * stateVector[ xPositionIndex ] / rangeCube;
        // gravAcceleration[ 1 ] = -1.0 * gravParameter * stateVector[ yPositionIndex ] / rangeCube;
        // gravAcceleration[ 2 ] = -1.0 * gravParameter * stateVector[ zPositionIndex ] / rangeCube;

        // gravitaitonal acceleration from ellipsoid potential model
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

//! Calculate Jacobi integral
/*!
 * This function calculates the jacobi integral of the asteroid-regolith problem at hand.
 */
double calculateJacobiIntegral( const double alpha,
                                const double beta,
                                const double gamma,
                                const double gravParameter,
                                const std::vector< double > &stateVector,
                                const std::vector< double > &asteroidRotationVector )
{
    const double xVelocitySquare = stateVector[ xVelocityIndex ] * stateVector[ xVelocityIndex ];
    const double yVelocitySquare = stateVector[ yVelocityIndex ] * stateVector[ yVelocityIndex ];
    const double zVelocitySquare = stateVector[ zVelocityIndex ] * stateVector[ zVelocityIndex ];

    const double xPositionSquare = stateVector[ xPositionIndex ] * stateVector[ xPositionIndex ];
    const double yPositionSquare = stateVector[ yPositionIndex ] * stateVector[ yPositionIndex ];

    const double zRotationSquare = asteroidRotationVector[ 2 ] * asteroidRotationVector[ 2 ];

    double gravPotential = 0.0;

    // compute point mass potential
    // const std::vector< double > positionVector = { stateVector[ xPositionIndex ],
    //                                                stateVector[ yPositionIndex ],
    //                                                stateVector[ zPositionIndex ] };
    // const double range = vectorNorm( positionVector );
    // gravPotential = gravParameter / range;

    // compute ellipsoid potential
    computeEllipsoidGravitationalPotential( alpha,
                                            beta,
                                            gamma,
                                            gravParameter,
                                            stateVector[ xPositionIndex ],
                                            stateVector[ yPositionIndex ],
                                            stateVector[ zPositionIndex ],
                                            gravPotential );

    double jacobi = ( 1.0 / 2.0 ) * ( xVelocitySquare + yVelocitySquare + zVelocitySquare )
                    - ( 1.0 / 2.0 ) * zRotationSquare * ( xPositionSquare + yPositionSquare )
                    - gravPotential;

    return jacobi;
}

//! Store data in SQLite database
/*!
 * This function stores the data in a SQLite database
 */
void storeData( const std::vector< double > &bodyFrameStateVector,
                const std::vector< double > &inertialFrameStateVector,
                const std::vector< double > &orbitalElements,
                const double currentTime,
                const double launchAzimuth,
                const double launchDeclination,
                const double kineticEnergy,
                const double potentialEnergy,
                const double totalEnergy,
                const double jacobi,
                const int startFlag,
                const int escapeFlag,
                const int crashFlag,
                SQLite::Statement &query )
{
    query.bind( ":position_x",          bodyFrameStateVector[ xPositionIndex ] );
    query.bind( ":position_y",          bodyFrameStateVector[ yPositionIndex ] );
    query.bind( ":position_z",          bodyFrameStateVector[ zPositionIndex ] );
    query.bind( ":velocity_x",          bodyFrameStateVector[ xVelocityIndex ] );
    query.bind( ":velocity_y",          bodyFrameStateVector[ yVelocityIndex ] );
    query.bind( ":velocity_z",          bodyFrameStateVector[ zVelocityIndex ] );

    query.bind( ":time",                currentTime );

    query.bind( ":inertial_position_x", inertialFrameStateVector[ xPositionIndex ] );
    query.bind( ":inertial_position_y", inertialFrameStateVector[ yPositionIndex ] );
    query.bind( ":inertial_position_z", inertialFrameStateVector[ zPositionIndex ] );
    query.bind( ":inertial_velocity_x", inertialFrameStateVector[ xVelocityIndex ] );
    query.bind( ":inertial_velocity_y", inertialFrameStateVector[ yVelocityIndex ] );
    query.bind( ":inertial_velocity_z", inertialFrameStateVector[ zVelocityIndex ] );

    query.bind( ":launch_azimuth",      convertRadiansToDegree( launchAzimuth ) );
    query.bind( ":launch_declination",  convertRadiansToDegree( launchDeclination ) );

    query.bind( ":sma",                 orbitalElements[ 0 ] );
    query.bind( ":eccentricity",        orbitalElements[ 1 ] );
    query.bind( ":inclination",         orbitalElements[ 2 ] );
    query.bind( ":raan",                orbitalElements[ 3 ] );
    query.bind( ":aop",                 orbitalElements[ 4 ] );
    query.bind( ":ta",                  orbitalElements[ 5 ] );

    query.bind( ":kinetic_energy",      kineticEnergy );
    query.bind( ":potential_energy",    potentialEnergy );
    query.bind( ":energy",              totalEnergy );

    query.bind( ":jacobi_integral",     jacobi );

    query.bind( ":start_flag",          startFlag );
    query.bind( ":escape_flag",         escapeFlag );
    query.bind( ":crash_flag",          crashFlag );

    query.executeStep( );
    query.reset( );
}

//! Calculate trajectory of a single regolith around the asteroid
/*!
 * Given initial conditions in the asteroid's principle axis frame, the function performs numerical
 * integration of the EOMs and stores the ephemeris in a SQLite database.
 */
void calculateRegolithTrajectoryAroundAsteroid( const double alpha,
                                                const double beta,
                                                const double gamma,
                                                const double gravParameter,
                                                const std::vector< double > &asteroidRotationVector,
                                                const std::vector< double > &initialState,
                                                const double launchAzimuth,
                                                const double launchDeclination,
                                                const double initialStepSize,
                                                const double startTime,
                                                const double endTime,
                                                const double dataSaveIntervals,
                                                SQLite::Statement &query )
{
    // Start by saving the initial data points
    int startFlag = 1;
    int escapeFlag = 0;
    int crashFlag = 0;

    std::vector< double > bodyFrameStateVector = initialState;
    double currentTime = startTime;

    double jacobi = calculateJacobiIntegral( alpha,
                                             beta,
                                             gamma,
                                             gravParameter,
                                             bodyFrameStateVector,
                                             asteroidRotationVector );

    std::vector< double > inertialFrameStateVector( 6, 0.0 );
    std::vector< double > orbitalElements( 6, 0.0 );
    double kineticEnergy = 0.0;
    double potentialEnergy = 0.0;
    double totalEnergy = 0.0;

    storeData( bodyFrameStateVector,
               inertialFrameStateVector,
               orbitalElements,
               currentTime,
               launchAzimuth,
               launchDeclination,
               kineticEnergy,
               potentialEnergy,
               totalEnergy,
               jacobi,
               startFlag,
               escapeFlag,
               crashFlag,
               query );

    // reset the start flag to 0
    startFlag = 0;

    // set up boost odeint
    const double absoluteTolerance = 1.0e-15;
    const double relativeTolerance = 1.0e-15;
    typedef boost::numeric::odeint::runge_kutta_fehlberg78< std::vector< double > > stepperType;

    // initialize the ODE system with the constant parameters
    regolithAroundAsteroidNoPerturbations classObjectForEOM( gravParameter,
                                                             alpha,
                                                             beta,
                                                             gamma,
                                                             asteroidRotationVector[ 2 ] );

    // initialize time values
    currentTime = startTime;
    double intermediateEndTime = currentTime + dataSaveIntervals;

    // integration while loop
    while( intermediateEndTime <= endTime )
    {
        // perform integration with automatic step size and error control
        // The integrated state is stored in the input state vector
        // The same initial step size guess is used for every run of the loop, this step
        // size is internally altered to meet the accuracy requirements.
        size_t steps = boost::numeric::odeint::integrate_adaptive(
                            make_controlled( absoluteTolerance, relativeTolerance, stepperType( ) ),
                            classObjectForEOM,
                            bodyFrameStateVector,
                            currentTime,
                            intermediateEndTime,
                            initialStepSize );

        // update time variables
        currentTime = intermediateEndTime;
        intermediateEndTime = intermediateEndTime + dataSaveIntervals;

        // calculate the jacobi
        jacobi = calculateJacobiIntegral( alpha,
                                          beta,
                                          gamma,
                                          gravParameter,
                                          bodyFrameStateVector,
                                          asteroidRotationVector );

        // data storage
        storeData( bodyFrameStateVector,
                   inertialFrameStateVector,
                   orbitalElements,
                   currentTime,
                   launchAzimuth,
                   launchDeclination,
                   kineticEnergy,
                   potentialEnergy,
                   totalEnergy,
                   jacobi,
                   startFlag,
                   escapeFlag,
                   crashFlag,
                   query );

    } // end of outer integration while loop
}

} // namespace naos
