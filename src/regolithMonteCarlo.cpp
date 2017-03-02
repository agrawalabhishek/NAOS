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
#include <cstring>
#include <stdexcept>

#include <SQLiteCpp/SQLiteCpp.h>
#include <sqlite3.h>

#include "NAOS/ellipsoidSurfacePoints.hpp"
#include "NAOS/cubicRoot.hpp"
#include "NAOS/constants.hpp"
#include "NAOS/ellipsoidGravitationalAcceleration.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/computeEllipsoidSurfaceGravitationalAcceleration.hpp"
#include "NAOS/ellipsoidPotential.hpp"
#include "NAOS/misc.hpp"
#include "NAOS/particleAroundUniformlyRotatingEllipsoid.hpp"
#include "NAOS/regolithTrajectoryCalculator.hpp"
#include "NAOS/regolithAroundUniformlyRotatingAsteroid.hpp"

namespace naos
{

//! Function to perform sanity check on the regolith trajectory propagator
/*!
 * This function bypasses the function that calculates the initial position and velocity of regolith
 * on the surface of the asteroid, and directly calls the trajectory calculator function. This
 * function is meant to be used only as a test bed.
 */
void testBedForRegolithTrajectoryCalculator( const double alpha,
                                             const double beta,
                                             const double gamma,
                                             const double gravitationalParameter,
                                             const std::vector< double > &angularVelocityVector,
                                             const double integrationStepSize,
                                             const double startTime,
                                             const double dataSaveIntervals,
                                             SQLite::Statement &databaseQuery )
{
    //! Sanity check for the use of "executeSingleRegolithTrajectoryCalculation" function
    std::vector< double > initialKeplerianState = { 28000.0,
                                                    0.01,
                                                    20.0,
                                                    10.0,
                                                    10.0,
                                                    0.0 };

    std::vector< double > initialStateInertial
        = convertKeplerianElementsToCartesianCoordinates( initialKeplerianState,
                                                          gravitationalParameter );

    // account for non-zero start time value and calculate the initial state in body frame
    double phi = angularVelocityVector[ zPositionIndex ] * startTime;

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
    omegaCrossPosition = crossProduct( angularVelocityVector, inertialPositionVector );

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

    // executeSingleRegolithTrajectoryCalculation( alpha,
    //                                             beta,
    //                                             gamma,
    //                                             gravitationalParameter,
    //                                             angularVelocityVector,
    //                                             initialState,
    //                                             0.0,
    //                                             0.0,
    //                                             0.0,
    //                                             0.0,
    //                                             integrationStepSize,
    //                                             startTime,
    //                                             50000.0,
    //                                             databaseQuery,
    //                                             dataSaveIntervals );
}

//! Function that performs sanity check for ellipsoid acceleration and potential
/*!
 * This function performs a simple sanity check for the ellipsoid gravitation and potential functions
 * by comparing the values of acceleration and potential for a particlur point outside the ellipsoid
 * with the ones computed by 'Nicola' from his own routine in matlab.
 */
void sanityCheckForEllipsoidAccelerationAndPotential( const double alpha,
                                                      const double beta,
                                                      const double gamma )
{
    //! sanity check for ellipsoid grav acceleration and potential
    // (match values with ones given by nicola)
    const double testGravitationalParameter = 446382.0;
    std::vector< double > testLocation = { 10000.0,
                                           13000.0,
                                           8000.0 };

    std::vector< double > testAcceleration( 3, 0.0 );
    computeEllipsoidGravitationalAcceleration( alpha,
                                               beta,
                                               gamma,
                                               testGravitationalParameter,
                                               testLocation[ xPositionIndex ],
                                               testLocation[ yPositionIndex ],
                                               testLocation[ zPositionIndex ],
                                               testAcceleration );

    double testGravPotential = 0.0;
    computeEllipsoidGravitationalPotential( alpha,
                                            beta,
                                            gamma,
                                            testGravitationalParameter,
                                            testLocation[ xPositionIndex ],
                                            testLocation[ yPositionIndex ],
                                            testLocation[ zPositionIndex ],
                                            testGravPotential );

    const double potentialValueCheck = 2.37100525543964e-05 * 1.0e6;
    const std::vector< double > gravAccelerationValuesCheck = { -4.47629167383408e-07 * 1.0e3,
                                                                -9.6233888139995e-07 * 1.0e3,
                                                                -5.92208542399969e-07 * 1.0e3 };

    std::cout.precision( 15 );

    if( std::fabs( testGravPotential - potentialValueCheck ) >= 1.0e-12 )
    {
        std::cout << std::endl;
        std::cout << "ERROR!: Ellipsoid potential sanity check failed!" << std::endl;
        std::cout << "Computed Value = " << testGravPotential << std::endl;
        std::cout << "Test Value = " << potentialValueCheck << std::endl << std::endl;
    }

    if( std::fabs( testAcceleration[ 0 ] - gravAccelerationValuesCheck[ 0 ] ) >= 1.0e-12 )
    {
        std::cout << std::endl;
        std::cout << "ERROR!: Ellipsoid grav. acceleration 'X' component sanity check failed!" << std::endl;
        std::cout << "Computed Value = " << testAcceleration[ 0 ] << std::endl;
        std::cout << "Test Value = " << gravAccelerationValuesCheck[ 0 ] << std::endl << std::endl;
    }

    if( std::fabs( testAcceleration[ 1 ] - gravAccelerationValuesCheck[ 1 ] ) >= 1.0e-12 )
    {
        std::cout << std::endl;
        std::cout << "ERROR!: Ellipsoid grav. acceleration 'Y' component sanity check failed!" << std::endl;
        std::cout << "Computed Value = " << testAcceleration[ 1 ] << std::endl;
        std::cout << "Test Value = " << gravAccelerationValuesCheck[ 1 ] << std::endl << std::endl;
    }

    if( std::fabs( testAcceleration[ 2 ] - gravAccelerationValuesCheck[ 2 ] ) >= 1.0e-12 )
    {
        std::cout << std::endl;
        std::cout << "ERROR!: Ellipsoid grav. acceleration 'Z' component sanity check failed!" << std::endl;
        std::cout << "Computed Value = " << testAcceleration[ 2 ] << std::endl;
        std::cout << "Test Value = " << gravAccelerationValuesCheck[ 2 ] << std::endl << std::endl;
    }
}

//! create table to store the regolith trajectory results
/*!
 * this function creates a table in the database everytime it is called
 */
void createDatabaseTable( SQLite::Database &database )
{
    // Drop table from database if it exists.
    database.exec( "DROP TABLE IF EXISTS regolith_trajectory_results;" );

    // Set up SQL command to create table to store regolith_trajectory_results.
    std::ostringstream regolithTableCreate;
    regolithTableCreate
        << "CREATE TABLE regolith_trajectory_results ("
        << "\"sql_id\"                                          INTEGER PRIMARY KEY AUTOINCREMENT,"

        << "\"trajectory_id\"                                   INTEGER,"

        << "\"initial_position_x\"                              REAL,"
        << "\"initial_position_y\"                              REAL,"
        << "\"initial_position_z\"                              REAL,"
        << "\"initial_position_magnitude\"                      REAL,"
        << "\"initial_velocity_x\"                              REAL,"
        << "\"initial_velocity_y\"                              REAL,"
        << "\"initial_velocity_z\"                              REAL,"
        << "\"initial_velocity_magnitude\"                      REAL,"
        << "\"initial_inertial_position_x\"                     REAL,"
        << "\"initial_inertial_position_y\"                     REAL,"
        << "\"initial_inertial_position_z\"                     REAL,"
        << "\"initial_inertial_velocity_x\"                     REAL,"
        << "\"initial_inertial_velocity_y\"                     REAL,"
        << "\"initial_inertial_velocity_z\"                     REAL,"

        << "\"position_x\"                                      REAL,"
        << "\"position_y\"                                      REAL,"
        << "\"position_z\"                                      REAL,"
        << "\"velocity_x\"                                      REAL,"
        << "\"velocity_y\"                                      REAL,"
        << "\"velocity_z\"                                      REAL,"

        << "\"time\"                                            REAL,"

        << "\"inertial_position_x\"                             REAL,"
        << "\"inertial_position_y\"                             REAL,"
        << "\"inertial_position_z\"                             REAL,"
        << "\"inertial_velocity_x\"                             REAL,"
        << "\"inertial_velocity_y\"                             REAL,"
        << "\"inertial_velocity_z\"                             REAL,"

        << "\"launch_azimuth\"                                  REAL,"
        << "\"launch_declination\"                              REAL,"

        << "\"directional_escape_speed\"                        REAL,"
        << "\"directional_inertial_escape_speed\"               REAL,"

        << "\"sma\"                                             REAL,"
        << "\"eccentricity\"                                    REAL,"
        << "\"inclination\"                                     REAL,"
        << "\"raan\"                                            REAL,"
        << "\"aop\"                                             REAL,"
        << "\"ta\"                                              REAL,"

        << "\"kinetic_energy\"                                  REAL,"
        << "\"potential_energy\"                                REAL,"
        << "\"total_energy\"                                    REAL,"

        << "\"jacobi_integral\"                                 REAL,"

        << "\"start_flag\"                                      INTEGER,"
        << "\"end_flag\"                                        INTEGER,"
        << "\"escape_flag\"                                     INTEGER,"
        << "\"crash_flag\"                                      INTEGER"
        <<                                                      ");";

    // Execute command to create table.
    database.exec( regolithTableCreate.str( ).c_str( ) );

    if ( !database.tableExists( "regolith_trajectory_results" ) )
    {
        throw std::runtime_error( "ERROR: Creating table 'regolith_trajectory_results' failed!" );
    }
}

//! monte carlo simulation for regolith
/*!
 * This function performs monte carlo simulation to simulate trajectories for multiple regolith
 * lofted from the surface of an asteroid.
 *
 */
void executeRegolithMonteCarlo( const double alpha,
                                const double beta,
                                const double gamma,
                                const double gravitationalParameter,
                                const std::vector< double > &angularVelocityVector,
                                const double angularVelocityMagnitude,
                                const double integrationStepSize,
                                const double startTime,
                                const double endTime,
                                const double dataSaveIntervals,
                                std::ostringstream &databaseFilePath )
{
    double aXValue = 1.0;
    double aYValue = 0.0;
    double aZValue = 0.0;
    const double velocityMagnitudeFactor = 0.9;
    const double coneAngleDeclination = naos::convertDegreeToRadians( 45.0 );
    const double coneAngleAzimuthFactor = 1.0;

    // Open database in read/write mode.
    SQLite::Database database( databaseFilePath.str( ),
                               SQLITE_OPEN_READWRITE|SQLITE_OPEN_CREATE );

    // Create table for trajectory results in SQLite database.
    createDatabaseTable( database );
    std::cout << "SQLite database set up successfully!" << std::endl;

    // Start SQL transaction
    SQLite::Transaction transaction( database );

    // set up insert query
    std::ostringstream regolithTrajectoryTableInsert;
    regolithTrajectoryTableInsert
        << "INSERT INTO regolith_trajectory_results VALUES ("
        << "NULL,"
        << ":trajectory_id,"
        << ":initial_position_x,"
        << ":initial_position_y,"
        << ":initial_position_z,"
        << ":initial_position_magnitude,"
        << ":initial_velocity_x,"
        << ":initial_velocity_y,"
        << ":initial_velocity_z,"
        << ":initial_velocity_magnitude,"
        << ":initial_inertial_position_x,"
        << ":initial_inertial_position_y,"
        << ":initial_inertial_position_z,"
        << ":initial_inertial_velocity_x,"
        << ":initial_inertial_velocity_y,"
        << ":initial_inertial_velocity_z,"
        << ":position_x,"
        << ":position_y,"
        << ":position_z,"
        << ":velocity_x,"
        << ":velocity_y,"
        << ":velocity_z,"
        << ":time,"
        << ":inertial_position_x,"
        << ":inertial_position_y,"
        << ":inertial_position_z,"
        << ":inertial_velocity_x,"
        << ":inertial_velocity_y,"
        << ":inertial_velocity_z,"
        << ":launch_azimuth,"
        << ":launch_declination,"
        << ":directional_escape_speed,"
        << ":directional_inertial_escape_speed,"
        << ":sma,"
        << ":eccentricity,"
        << ":inclination,"
        << ":raan,"
        << ":aop,"
        << ":ta,"
        << ":kinetic_energy,"
        << ":potential_energy,"
        << ":total_energy,"
        << ":jacobi_integral,"
        << ":start_flag,"
        << ":end_flag,"
        << ":escape_flag,"
        << ":crash_flag"
        << ");";

    SQLite::Statement databaseQuery( database, regolithTrajectoryTableInsert.str( ) );

    sanityCheckForEllipsoidAccelerationAndPotential( alpha, beta, gamma );

    // testBedForRegolithTrajectoryCalculator( alpha,
    //                                         beta,
    //                                         gamma,
    //                                         gravitationalParameter,
    //                                         angularVelocityVector,
    //                                         integrationStepSize,
    //                                         startTime,
    //                                         dataSaveIntervals,
    //                                         databaseQuery );

    // azimuth angle iterator begins here
    for( int azimuthIterator = 0; azimuthIterator < 360; azimuthIterator = azimuthIterator + 1 )
    {
        // calculate the azimuth angle
        const double coneAngleAzimuth
                = naos::convertDegreeToRadians( coneAngleAzimuthFactor * azimuthIterator );
        // const double coneAngleAzimuth = naos::convertDegreeToRadians( 16.0 );

        executeRegolithTrajectoryCalculation( alpha,
                                              beta,
                                              gamma,
                                              gravitationalParameter,
                                              angularVelocityVector,
                                              angularVelocityMagnitude,
                                              aXValue,
                                              aYValue,
                                              aZValue,
                                              coneAngleAzimuth,
                                              coneAngleDeclination,
                                              velocityMagnitudeFactor,
                                              integrationStepSize,
                                              startTime,
                                              endTime,
                                              dataSaveIntervals,
                                              databaseQuery );
    } // end of azimuth angle iterator loop

    transaction.commit( );
}

} // namespace naos
