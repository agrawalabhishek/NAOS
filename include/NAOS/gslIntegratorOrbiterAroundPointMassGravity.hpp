/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef GSL_INTEGRATOR_ORBITER_AROUND_POINT_MASS_GRAVITY_HPP
#define GSL_INTEGRATOR_ORBITER_AROUND_POINT_MASS_GRAVITY_HPP

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/misc.hpp"

namespace naos
{

//! Equations of motion for spacecraft or particle around a uniformly rotating ellipsoid (URE)
/*!
 * function defining the equation of motion which will be used in the integrator routine of GSL
 *
 */
 struct asteroidParameters
 {
     double xGravAcceleration;
     double yGravAcceleration;
     double zGravAcceleration;
     double zRotation;
 };

 int equationsOfMotion( double t, const double X[], double dXdt[], void *parameters );

//! Execute routine to compute orbiter/particle motion around a uniformly rotating ellipsoid (URE)
/*!
 * This routine will compute the motion of the orbiter or particle around an asteroid modelled as
 * a uniformly rotating triaxial ellipsoid. The routine is called in main for execution.
 * @sa ellipsoidGravitationalAcceleration.cpp, orbiterEquationsOfMotion.hpp, rk4.hpp
 *
 * @param[in] alpha                         largest semi axis of the ellipsoid
 * @param[in] beta                          intermediate semi axis of the ellipsoid
 * @param[in] gamma                         smallest largest semi axis of the ellipsoid
 * @param[in] gravParameter                 gravitational parameter of the ellipsoid
 * @param[in] density                       density of the ellipsoid
 * @param[in] Wvector                       angular velocity vector for the ellipsoid
 * @param[in] Wmagnitude                    magnitude of the angular velocity
 * @param[in] initialVectorIsCartesian      boolean flag to indicate type of initial vector
 *                                          value is false for kepler elements, and true for
 *                                          cartesian elements.
 * @param[in] initialVector                 vector containing the initial state of the orbiter
 * @param[in] integrationStepSize           step size to be used for integration
 * @param[in] startTime                     integration start time
 * @param[in] endTime                       integration end time
 * @param[in] filePath                      path to csv file where all the results are stored
 */
 void executeOrbiterAroundPointMassGravity( const double alpha,
                                            const double beta,
                                            const double gamma,
                                            const double gravParameter,
                                            const double density,
                                            std::vector< double > Wvector,
                                            const double Wmagnitude,
                                            const bool initialVectorIsCartesian,
                                            const std::vector< double > &initialVector,
                                            const double integrationStepSize,
                                            const double startTime,
                                            const double endTime,
                                            std::ostringstream &filePath );

} // namespace naos

#endif // GSL_INTEGRATOR_ORBITER_AROUND_POINT_MASS_GRAVITY_HPP
