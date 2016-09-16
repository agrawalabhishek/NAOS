/*
 * Copyright (c) 2016 Abhishek Agrawal abhishek.agrawal@protonmail.com
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef EXECUTE_ORBITER_AROUND_URE_HPP
#define EXECUTE_ORBITER_AROUND_URE_HPP

namespace naos
{

//! Execute routine to compute orbiter/particle motion around a uniformly rotating ellipsoid (URE)
/*!
 * This routine will compute the motion of the orbiter or particle around an asteroid modelled as
 * a uniformly rotating triaxial ellipsoid. The routine is called in main for execution.
 * @sa ellipsoidGravitationalAcceleration.cpp, orbiterEquationsOfMotion.hpp, rk4.hpp
 *
 */
 const void executeOrbiterAroundURE( const double alpha,
                                     const double beta,
                                     const double gamma,
                                     const double gravParameter,
                                     const double density,
                                     std::vector< double > Wvector,
                                     const double Wmagnitude,
                                     const double semiMajor,
                                     const double eccentricity,
                                     const double inclination,
                                     const double RAAN,
                                     const double AOP,
                                     const double TA,
                                     const double integrationStepSize,
                                     const double startTime,
                                     const double endTime,
                                     std::ostringstream &filePath );

} // namespace naos

#endif // EXECUTE_ORBITER_AROUND_URE_HPP
