/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef SPRING_MASS_EQUATIONS_OF_MOTION
#define SPRING_MASS_EQUATIONS_OF_MOTION

#include <vector>

namespace naos
{
//! Equations of motion for a basic spring mass system
/*!
 * A struct defining the equations of motion for a basic autonomous spring mass system. This system
 * is used to verify the routine for the numerical integrator as the spring mass system already has
 * an analytical solution.

 * An object of this struct, returns the value of dXdt at time t.
 * @param[in] t         time
 * @param[in] X         State vector containing the values defined at time t
 * @param[out] dXdt     Vector of derivatives of the state vector at time t
 */
struct eomSpringMass
{
    const double mass;
    const double springConstant;

    const int positionIndex = 0;
    const int velocityIndex = 1;

    // Default constructor with member initializer list
    eomSpringMass( double aMassValue, double aSpringConstantValue )
                    : mass( aMassValue ),
                      springConstant( aSpringConstantValue )
    { }
    void operator() ( const double t,
                      std::vector< double > &X,
                      std::vector< double > &dXdt )
    {
        dXdt[ positionIndex ] = X[ velocityIndex ];
        dXdt[ velocityIndex ] = -1.0 * ( springConstant / mass ) * X[ positionIndex ];
    }
};

} // namespace naos

#endif // SPRING_MASS_EQUATIONS_OF_MOTION
