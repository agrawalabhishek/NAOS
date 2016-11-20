/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef BOOST_INTEGRATOR_RESTRICTED_TWO_BODY_PROBLEM_HPP
#define BOOST_INTEGRATOR_RESTRICTED_TWO_BODY_PROBLEM_HPP

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

namespace naos
{

//! Restricted two body problem integration
/*!
 * integrate the equations of motion for a particle around a point mass gravitational centre.
 */
void boostIntegratorRestrictedTwoBodyProblem( const double gravParameter,
                                              std::vector< double > &initialOrbitalElements,
                                              const double initialStepSize,
                                              const double startTime,
                                              const double endTime,
                                              std::ostringstream &outputFilePath,
                                              const int dataSaveIntervals );

} // namespace naos

#endif // BOOST_INTEGRATOR_RESTRICTED_TWO_BODY_PROBLEM_HPP
