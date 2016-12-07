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
#include "NAOS/springMassIntegratorTest.hpp"
#include "NAOS/postAnalysis.hpp"
#include "NAOS/particleAroundUniformlyRotatingEllipsoid.hpp"
#include "NAOS/regolithTrajectoryCalculator.hpp"

namespace naos
{

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
        << "\"trajectory_id\"                           INTEGER PRIMARY KEY AUTOINCREMENT,"

        << "\"position_x\"                              REAL,"
        << "\"position_y\"                              REAL,"
        << "\"position_z\"                              REAL,"
        << "\"velocity_x\"                              REAL,"
        << "\"velocity_y\"                              REAL,"
        << "\"velocity_z\"                              REAL,"
        << "\"time\"                                    REAL,"

        << "\"inertial_position_x\"                     REAL,"
        << "\"inertial_position_y\"                     REAL,"
        << "\"inertial_position_z\"                     REAL,"
        << "\"inertial_velocity_x\"                     REAL,"
        << "\"inertial_velocity_y\"                     REAL,"
        << "\"inertial_velocity_z\"                     REAL,"

        << "\"launch_azimuth\"                          REAL,"
        << "\"launch_declination\"                      REAL,"

        << "\"sma\"                                     REAL,"
        << "\"eccentricity\"                            REAL,"
        << "\"inclination\"                             REAL,"
        << "\"raan\"                                    REAL,"
        << "\"aop\"                                     REAL,"
        << "\"ta\"                                      REAL,"

        << "\"kinetic_energy\"                          REAL,"
        << "\"potential_energy\"                        REAL,"
        << "\"total_energy\"                            REAL,"

        << "\"start_flag\"                              INTEGER,"
        << "\"escape_flag\"                             INTEGER,"
        << "\"crash_flag\"                              INTEGER"
        <<                                              ");";

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
    double aYValue = 1.0;
    double aZValue = 1.0;
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
        << ":sma,"
        << ":eccentricity,"
        << ":inclination,"
        << ":raan,"
        << ":aop,"
        << ":ta,"
        << ":kinetic_energy,"
        << ":potential_energy,"
        << ":energy,"
        << ":start_flag,"
        << ":escape_flag,"
        << ":crash_flag"
        << ");";

    SQLite::Statement databaseQuery( database, regolithTrajectoryTableInsert.str( ) );

    // azimuth angle iterator begins here
    for( int azimuthIterator = 0; azimuthIterator < 360; azimuthIterator++ )
    {
        // calculate the azimuth angle
        const double coneAngleAzimuth
                = naos::convertDegreeToRadians( coneAngleAzimuthFactor * azimuthIterator );

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
