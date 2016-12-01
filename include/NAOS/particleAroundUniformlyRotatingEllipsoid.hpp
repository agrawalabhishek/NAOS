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

} // namespace naos

#endif // PARTICLE_AROUND_UNIFORMLY_ROTATING_ELLIPSOID
