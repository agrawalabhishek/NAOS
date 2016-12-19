/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef REGOLITH_AROUND_UNIFORMLY_ROTATING_ASTEROID_HPP
#define REGOLITH_AROUND_UNIFORMLY_ROTATING_ASTEROID_HPP

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

//! Calculate Jacobi integral
/*!
 * This function calculates the jacobi integral of the asteroid-regolith problem at hand.
 */
double calculateJacobiIntegral( const double alpha,
                                const double beta,
                                const double gamma,
                                const double gravParameter,
                                const std::vector< double > &stateVector,
                                const std::vector< double > &asteroidRotationVector );

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
                SQLite::Statement &query );

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
                                                SQLite::Statement &query );

} // namespace naos

#endif // REGOLITH_AROUND_UNIFORMLY_ROTATING_ASTEROID_HPP
