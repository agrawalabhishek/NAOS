/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef PERTURBING_ACCELERATIONS_HPP
#define PERTURBING_ACCELERATIONS_HPP

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <vector>
#include <stdexcept>

#include "NAOS/basicAstro.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/constants.hpp"
#include "NAOS/misc.hpp"

namespace naos
{

//! open sun ephemeris and extract data for a specified time value
/*!
 * This file opens the ephemeris data for the sun around the asteroid at specified time values in
 * multiples of ten. If a finer ephemeris is required then the sunAsteroidTwoBodyProblem must be
 * re-run prior to using this mode, with a finer data save interval value.
 */
std::vector< double > extractSunEphemeris( const double timeValue,
                                           std::ostringstream &ephemerisFilePath );

//! Third body effect function based only on position vectors as input - general function
/*!
 * Compute the perturbing acceleration from the thrid body effect of any perturbing body
 */
std::vector< double > computeThirdBodyEffect( std::vector< double > perturberPositionVector,
                                              std::vector< double > targetPositionVector,
                                              const double perturberGravParameter );

//! Third-body effect
/*!
 * Compute the perturbing acceleration from the thrid body effect of the sun.
 */
std::vector< double > computeSunThirdBodyEffectAcceleration( const std::vector< double > &regolithPositionVector,
                                                             std::vector< double > &asteroidRotationVector,
                                                             const std::vector< double > &initialSunOrbitalElements,
                                                             const double timeValue );

//! Solar radiation pressure general function
/*!
 * General computation of the perturbing acceleration from solar radition pressure.
 */
std::vector< double > computeSRPAcceleration( std::vector< double > &positionVectorTargetToSource,
                                              const double targetAlbedo,
                                              const double incidentArea,
                                              const double targetMass,
                                              const double solarConstant );

//! Solar radiation pressure
/*!
 * comupte the perturbing acceleration from solar radition pressure.
 */
std::vector< double > computeSolarRadiationPressureAcceleration( const std::vector< double > &regolithPositionVector,
                                                                 std::vector< double > &asteroidRotationVector,
                                                                 const std::vector< double > &initialSunOrbitalElements,
                                                                 const double timeValue );

} // namespace naos

#endif // PERTURBING_ACCELERATIONS_HPP
