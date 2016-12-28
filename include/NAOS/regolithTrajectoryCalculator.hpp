/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef REGOLITH_TRAJECTORY_CALCULATOR
#define REGOLITH_TRAJECTORY_CALCULATOR

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

#include <SQLiteCpp/SQLiteCpp.h>
#include <sqlite3.h>

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/particleAroundUniformlyRotatingEllipsoid.hpp"

namespace naos
{

//! Compute the velocity vector of lofted regolith
/*!
 * This sub-routine computes the velocity vector, given a magnitude, a unit normal vector from the
 * point on the surface where the regolith is lofted from, and the two conic angles which describe
 * the velocity vector's direction relative to the normal vector. A backwards approach is used
 * to go from the velocity vector in the final rotated intermediate frame back to the body fixed
 * frame. More details are given in the thesis report and author's personal notes.
 *
 */
void computeRegolithVelocityVector( std::vector< double > regolithPositionVector,
                                    const double velocityMagnitude,
                                    const double coneAngleAzimuth,
                                    const double coneAngleDeclination,
                                    std::vector< double > &unitNormalVector,
                                    std::vector< double > &regolithVelocityVector );

void computeRegolithVelocityVector2( std::vector< double > regolithPositionVector,
                                     const double velocityMagnitude,
                                     const double coneAngleAzimuth,
                                     const double coneAngleDeclination,
                                     std::vector< double > &unitNormalVector,
                                     std::vector< double > &regolithVelocityVector,
                                     std::vector< double > &regolithDirectionUnitVector );

//! Compute trajectory for regolith ejected from the surface of an asteroid
/*!
 * This routine computes the trajectory for a single regolith ejected from the surface of an
 * asteroid by first computing the appropriate initial conditions (IC) and then integrating the
 * equations of motion using those ICs.
 *
 */
void calculateRegolithTrajectory( const double alpha,
                                  const double beta,
                                  const double gamma,
                                  const double gravitationalParameter,
                                  const double density,
                                  const std::vector< double > W,
                                  const double Wmagnitude,
                                  const double aXValue,
                                  const double aYValue,
                                  const double aZValue,
                                  const double coneAngleAzimuth,
                                  const double coneAngleDeclination,
                                  const double velocityMagnitudeFactor,
                                  const double integrationStepSize,
                                  const double startTime,
                                  const double endTime,
                                  const double dataSaveIntervals,
                                  std::ostringstream &filePath );

//! Compute trajectory for regolith ejected from the surface of an asteroid (data saved in SQL db)
/*!
 * This routine computes the trajectory for a single regolith ejected from the surface of an
 * asteroid by first computing the appropriate initial conditions (IC) and then integrating the
 * equations of motion using those ICs.
 *
 */
void executeRegolithTrajectoryCalculation( const double alpha,
                                           const double beta,
                                           const double gamma,
                                           const double gravitationalParameter,
                                           const std::vector< double > W,
                                           const double Wmagnitude,
                                           const double aXValue,
                                           const double aYValue,
                                           const double aZValue,
                                           const double coneAngleAzimuth,
                                           const double coneAngleDeclination,
                                           const double velocityMagnitudeFactor,
                                           const double integrationStepSize,
                                           const double startTime,
                                           const double endTime,
                                           const double dataSaveIntervals,
                                           SQLite::Statement &databaseQuery );
} // namespace naos

#endif // REGOLITH_TRAJECTORY_CALCULATOR
