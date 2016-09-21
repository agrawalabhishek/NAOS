/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef NAOS_ELLIPSOID_GRAVITATIONAL_ACCELERATION_HPP
#define NAOS_ELLIPSOID_GRAVITATIONAL_ACCELERATION_HPP

#include "NAOS/constants.hpp"

namespace naos
{
//! Compute ellipsoid gravitational Acceleration
/*!
 * Compute the gravitational acceleration i.e. partial derivative of the external gravitational
 * potential of a tri-axial ellipsoid. The (x,y,z) coordinates for this should be provided in the
 * body fixed principal axis frame [1]. The accelerations are in the body frame.
 *
 * @param[in] alpha                 Largest semi-axis of ellipsoid
 * @param[in] beta                  Intermediate semi-axis of ellipsoid
 * @param[in] gamma                 Smallest semi-axis of ellipsoid
 * @param[in] gravParameter         gravitational parameter of the ellipsoid
 * @param[in] xCoordinate           X coordinate where the gravitational acceleration has to be
 *                                  calculated
 * @param[in] yCoordinate           Y coordinate where the gravitational acceleration has to be
 *                                  calculated
 * @param[in] zCoordinate           Z coordinate where the gravitational acceleration has to be
 *                                  calculated
 * @param[in/out] gravAcceleration  Output vector containing the gravitational acceleration
 *                                  components
 */
 const void computeEllipsoidGravitationalAcceleration( const double alpha,
                                                       const double beta,
                                                       const double gamma,
                                                       const double gravParameter,
                                                       const double xCoordinate,
                                                       const double yCoordinate,
                                                       const double zCoordinate,
                                                       Vector3 &gravAcceleration );

} // namespace naos

#endif // NAOS_ELLIPSOID_GRAVITATIONAL_ACCELERATION_HPP
