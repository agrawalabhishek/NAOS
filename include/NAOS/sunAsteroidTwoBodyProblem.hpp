/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef SUN_ASTEROID_TWO_BODY_PROBLEM
#define SUN_ASTEROID_TWO_BODY_PROBLEM

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <vector>
#include <stdexcept>

#include <boost/numeric/odeint.hpp>

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/misc.hpp"

namespace naos
{

//! convert state vector from inertial frame to body frame
/*!
 * function converts state vector from inertial frame centred at the asteroid to a rotating frame
 * which is alligned with the asteroids principle axes.
 */
void convertInertialFrameVectorToBodyFrame( const std::vector< double > &asteroidRotationVector,
                                            const std::vector< double > &inertialStateVector,
                                            const double currentTime,
                                            std::vector< double > &bodyFrameStateVector );

//! compute jacobi constant
/*!
 * compute the jacobi integral for the two body problem. It should be conserved since eom's have no
 * explicit time dependance.
 */
double computeSunAsteroidJacobiConstant( const std::vector< double > &bodyFrameStateVector,
                                         const std::vector< double > &asteroidRotationVector,
                                         const double gravParameter );

//! compute energy
/*!
 * computes the energy for sun orbiting the asteroid
 */
double computeSunAsteroidEnergy( const std::vector< double > &inertialStateVector,
                                 const double gravParameter );

//! sun-asteroid Restricted two body problem integration
/*!
 * integrate the equations of motion for the sun around the asteroid.
 */
void executeSunAsteroidTwoBodyProblem( const double gravParameter,
                                       const std::vector< double > &asteroidRotationVector,
                                       std::vector< double > &initialOrbitalElements,
                                       const double initialStepSize,
                                       const double startTime,
                                       const double endTime,
                                       std::ostringstream &outputFilePath,
                                       const int dataSaveIntervals );

} // namespace naos

#endif // SUN_ASTEROID_TWO_BODY_PROBLEM
