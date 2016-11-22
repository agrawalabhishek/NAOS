/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef POST_ANALYSIS_HPP
#define POST_ANALYSIS_HPP

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <vector>

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/misc.hpp"
#include "NAOS/ellipsoidPotential.hpp"

namespace naos
{

//! Calculate the jacobian
/*!
 * calculate the jacobian for the point mass gravity potential case.
 */
void calculateJacobianPointMassGravity( const std::vector< double > &angularVelocityVector,
                                        const double gravParameter,
                                        std::ostringstream &inputFilePath,
                                        std::ostringstream &outputFilePath );

//! Calculate the jacobian
/*!
 * calculate the jacobian for ellipsoid gravity potential case.
 */
void calculateJacobianEllipsoidPotentialModel( const double alpha,
                                               const double beta,
                                               const double gamma,
                                               const std::vector< double > &angularVelocityVector,
                                               const double gravParameter,
                                               std::ostringstream &inputFilePath,
                                               std::ostringstream &outputFilePath );

//! Convert cartesian states to orbital elements post simulation.
/*!
 * Convert the cartesian states stored in a csv file to orbital elements and store it in a seperate
 * csv file. The cartesian states are to be in body fixed frame. And the central body is assumed to
 * be in a unoform rotation mode. This information is important since the cartesian states will be
 * transformed from the body fixed frame to the inertial frame before converting it to the orbital
 * elements.
 */
void postSimulationOrbitalElementsConversion( const std::vector< double > &angularVelocityVector,
                                              const double gravParameter,
                                              std::ostringstream &inputFilePath,
                                              std::ostringstream &outputFilePath );

} // namespace naos

#endif // POST_ANALYSIS_HPP
