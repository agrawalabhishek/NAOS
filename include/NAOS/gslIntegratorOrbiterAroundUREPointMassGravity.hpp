/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef GSL_INTEGRATOR_ORBITER_AROUND_URE_POINT_MASS_GRAVITY_HPP
#define GSL_INTEGRATOR_ORBITER_AROUND_URE_POINT_MASS_GRAVITY_HPP

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/misc.hpp"

namespace naos
{

//! equations of motion
/*!
 * first order differential equations describing the motion of a particle or spacecraft around a
 * central body modeled as a point mass gravity. Look at GSL odeiv2 examples to see the format
 * of how the equations are written for a c++ application.
 */
 struct asteroidParameters
 {
     double xGravAcceleration;
     double yGravAcceleration;
     double zGravAcceleration;
     double zRotation;
 };

 int equationsOfMotion( double currentTime,
                        const double stateVector[ ],
                        double stateTimeDerivative[ ],
                        void *parameters );

//! Get the jacobian for the integrator
int jacobianForIntegrator( double currentTime,
                           const double stateVector[ ],
                           double *dfdy[ ],
                           double dfdt[ ],
                           void *parameters );

//! Calculate the jacobian
/*!
 * calculate the jacobian for the point mass gravity potential case.
 */
double calculateJacobianPointMassGravity( const double stateVector[ ],
                                          const double angularVelocity,
                                          const double gravParameter );

//! Execute orbiter around an asteroid (URE) modeled as a point mass
/*!
 * integrate the equations of motion for a particle around an asteroid using the GSL odeiv2 library
 * for a given set of initial conditions given in the form of classical orbital elements. The
 * cartesian state vector is store in a csv file.
 */
void gslIntegratorOrbiterAroundUREPointMassGravity( const double gravParameter,
                                                    std::vector< double > &asteroidRotationVector,
                                                    std::vector< double > &initialOrbitalElements,
                                                    const double initialStepSize,
                                                    const double startTime,
                                                    const double endTime,
                                                    std::ostringstream &outputFilePath );

} // namespace naos

#endif // GSL_INTEGRATOR_ORBITER_AROUND_URE_POINT_MASS_GRAVITY_HPP
