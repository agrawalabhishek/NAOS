/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef COMPUTE_PERTURBATIONS_FOR_SPECIFIC_COORDINATES_HPP
#define COMPUTE_PERTURBATIONS_FOR_SPECIFIC_COORDINATES_HPP

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
void computeGravitationalAndPerturbingAccelerationsForSpecificCoordinates( const double alpha,
                                                                           const double beta,
                                                                           const double gamma,
                                                                           const double gravParameter,
                                                                           std::vector< double > &asteroidRotationVector,
                                                                           std::ostringstream &filePath );

} // namespace naos

#endif // COMPUTE_PERTURBATIONS_FOR_SPECIFIC_COORDINATES_HPP
