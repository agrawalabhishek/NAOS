/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef ASTEROID_SURFACE_GRAVITATIONAL_AND_PERTURBING_ACCELERATIONS_HPP
#define ASTEROID_SURFACE_GRAVITATIONAL_AND_PERTURBING_ACCELERATIONS_HPP

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <vector>
#include <stdexcept>

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/misc.hpp"
#include "NAOS/ellipsoidGravitationalAcceleration.hpp"
#include "NAOS/perturbingAccelerations.hpp"

namespace naos
{

//! compute the gravitational and perturbing accelerations
/*!
 * this function computes and stores the gravitational and perturbing accelerations on the surface
 * an asteroid. The data can be plotted to analyze the relative strengths of the perturbing
 * accelerations.
 */
void computeGravityAndPerturbingAccelerations( const double alpha,
                                               const double beta,
                                               const double gamma,
                                               const double gravParameter,
                                               std::vector< double > &asteroidRotationVector,
                                               std::ostringstream &filePath );

} // namespace naos

#endif // ASTEROID_SURFACE_GRAVITATIONAL_AND_PERTURBING_ACCELERATIONS_HPP
