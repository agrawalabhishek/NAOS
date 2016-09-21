/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef COMPUTE_ELLIPSOID_SURFACE_GRAVITATIONAL_ACCELERATION_HPP
#define COMPUTE_ELLIPSOID_SURFACE_GRAVITATIONAL_ACCELERATION_HPP

namespace naos
{

//! Compute gravitational Acceleration on the surface of the ellipsoid shaped asteroid
/*!
 * Compute absolute value of the gravitational acceleration on the surface of the ellipsoid and
 * store it in a CSV file.
 * @sa ellipsoidGravitationalAcceleration.cpp
 * @sa cubicRoot.cpp
 *
 * @param[in] alpha, beta, gamma        semi-axes of the ellipsoid such that alpha>beta>gamma
 * @param[in] gravitationalParameter    gravitational parameter of the ellipsoid
 */
const void computeEllipsoidSurfaceGravitationalAcceleration( const double alpha,
                                                             const double beta,
                                                             const double gamma,
                                                             const double gravitationalParameter );

} // namespace naos

#endif // COMPUTE_ELLIPSOID_SURFACE_GRAVITATIONAL_ACCELERATION_HPP
