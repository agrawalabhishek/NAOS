/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef PARTICLE_AROUND_UNIFORMLY_ROTATING_ELLIPSOID
#define PARTICLE_AROUND_UNIFORMLY_ROTATING_ELLIPSOID

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

namespace naos
{

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
                                     const int dataSaveIntervals );

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
                                         const int dataSaveIntervals );

//! Trajectory calculation for regolith around an asteroid (data saved in SQL db)
/*!
 * Same as the previous function, except that the initial conditions are now given as a cartesian
 * state. The initial cartesian state should be given in body fixed frame of the asteroid.
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
                                                 const int dataSaveIntervals );

//! convert body frame state vector to inertial frame for a uniformly rotating case
/*!
 * this function cnverts a given state vector in asteroid body frame coordinates to inertial frame.
 * The body frame is rotating uniformly with respect to the inertial frame about the z axis.
 */
void convertBodyFrameVectorToInertialFrame( const std::vector< double > &asteroidRotationVector,
                                            const std::vector< double > &currentStateVector,
                                            const double currentTime,
                                            std::vector< double > &inertialState );

//! Calculate particle energy
/*!
 * The routine calculates the energy of the orbiting particle, kinetic, potential and total energy.
 */
void computeParticleEnergy( const double alpha,
                            const double beta,
                            const double gamma,
                            const double gravParameter,
                            std::vector< double > inertialState,
                            double &kineticEnergy,
                            double &potentialEnergy,
                            double &totalEnergy );

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
                    std::ostringstream &location );

} // namespace naos

#endif // PARTICLE_AROUND_UNIFORMLY_ROTATING_ELLIPSOID
