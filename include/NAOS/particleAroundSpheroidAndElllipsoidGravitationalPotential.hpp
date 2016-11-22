/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef PARTICLE_AROUND_SPHEROID_AND_ELLIPSOID_GRAVITATIONAL_POTENTIAL_HPP
#define PARTICLE_AROUND_SPHEROID_AND_ELLIPSOID_GRAVITATIONAL_POTENTIAL_HPP

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
//! particle around spheroid integration
/*!
 * integrate the equations of motion for a particle around a spheroid. The gravitational accelerations
 * calculated using the ellipsoid gravitational potential model.
 */
void executeParticleAroundSpheroid( const double alpha,
                                    const double gravParameter,
                                    std::vector< double > asteroidRotationVector,
                                    std::vector< double > &initialOrbitalElements,
                                    const double initialStepSize,
                                    const double startTime,
                                    const double endTime,
                                    std::ostringstream &outputFilePath,
                                    const int dataSaveIntervals );
} // namespace naos

#endif // PARTICLE_AROUND_SPHEROID_AND_ELLIPSOID_GRAVITATIONAL_POTENTIAL_HPP
