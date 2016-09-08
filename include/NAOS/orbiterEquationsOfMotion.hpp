/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef ORBITER_EQUATIONS_OF_MOTION
#define ORBITER_EQUATIONS_OF_MOTION

#include <vector>

#include "NAOS/constants.hpp"

namespace naos
{
//! Equations of motion for spacecraft or particle around a uniformly rotating ellipsoid (URE)
/*!
 * A struct defining the discretized equations of motion (EOM) for an orbiter around a single
 * uniformly rotating triaxial ellipsoid. The equations can be integrated numerically in a seperate
 * routine. The system is modular in the sense that the user can provide the value of the constant
 * i.e. the rotation rate (about the z axis in this case), in the main without editing the
 * equations here in the struct. The EOMs are defined in the body fixed principal axis frame.
 * 't' denotes time, 'X' is the state vector, 'dXdt' is the time derivative of the state vector.
 *
 * Numerical integration provides the solution to the differential equation for a user given value
 * for the rotation rate of the ellipsoid and all three components of gravitational acceleration.

 * An object of this struct, returns the value of dXdt at time t.
 * @param[in] t         time
 * @param[in] X         State vector containing the values defined at time t
 * @param[out] dXdt     Vector of derivatives of the state vector at time t
 */
struct eomOrbiterURE
{
    const double Wz; // rotation rate about the principal axis Z
    const double Ux;
    const double Uy;
    const double Uz;
    // Default constructor with member initializer list
    eomOrbiterURE( const double wValue,
                   Vector3 gravAcceleration )
                    : Wz( wValue ),
                      Ux( gravAcceleration[ naos::xPositionIndex ] ),
                      Uy( gravAcceleration[ naos::yPositionIndex ] ),
                      Uz( gravAcceleration[ naos::zPositionIndex ] )
    { }
    void operator() ( const double t, Vector6 &X, Vector6 &dXdt )
    {
        dXdt[ naos::xPositionIndex ] = X[ naos::xVelocityIndex ];
        dXdt[ naos::yPositionIndex ] = X[ naos::yVelocityIndex ];
        dXdt[ naos::zPositionIndex ] = X[ naos::zVelocityIndex ];
        dXdt[ naos::xVelocityIndex ] = Ux + Wz * X[ naos::yVelocityIndex ];
        dXdt[ naos::yVelocityIndex ] = Uy - Wz * X[ naos::xVelocityIndex ];
        dXdt[ naos::zVelocityIndex ] = Uz;
    }
};

} // namespace naos

#endif // ORBITER_EQUATIONS_OF_MOTION