/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef REGOLITH_MONTE_CARLO
#define REGOLITH_MONTE_CARLO

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
                                             SQLite::Statement &databaseQuery );

//! function that performs sanity check for ellipsoid acceleration and potential
/*!
 * This function performs a simple sanity check for the ellipsoid gravitation and potential functions
 * by comparing the values of acceleration and potential for a particlur point outside the ellipsoid
 * with the ones computed by 'Nicola' from his own routine in matlab.
 */
void sanityCheckForEllipsoidAccelerationAndPotential( const double alpha,
                                                      const double beta,
                                                      const double gamma );

//! create table to store the regolith trajectory results
/*!
 * this function creates a table in the database everytime it is called
 */
void createDatabaseTable( SQLite::Database &database );

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
                                std::ostringstream &databaseFilePath );

} // namespace naos

#endif // REGOLITH_MONTE_CARLO
