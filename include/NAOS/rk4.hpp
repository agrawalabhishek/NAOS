/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef RK4_HPP
#define RK4_HPP

#include <vector>

namespace naos
{

//! RK4 integrator
/*!
 * The classic RK4 integrator routine to integrate a set of
 * first order differential equations. The equations have to be provided in a seperate header file
 * in a struct. See orbiterEquationsOfMotion.hpp for an example. Addaptive step size control will
 * be applied from the outside, before calling the rk4 function.
 *
 * @param[in] Xcurrent      currently known state vector value i.e. before the integration
 * @param[in] t             current time value for which the state vector is known
 * @param[in] stepSize      the step size of integration
 * @param[out] Xnext        the state vector value after integration i.e. at next time step
 * @param[in] derivatives   An object of the struct containing the differential equations. The
 *                          object should have already gotten the initializing values before
 *                          passing it to the argument of the function rk4.
 */
template< typename Vector, class orbiterEquationsOfMotion >
void rk4(
    Vector &Xcurrent,
    const double t,
    const double stepSize,
    Vector &Xnext,
    orbiterEquationsOfMotion &derivatives )
{
    const double h = stepSize;
    const int numberOfElements = Xcurrent.size( );

    // Evaluate K1 step in RK4 algorithm.
    Vector dXdtCurrent( numberOfElements );
    Vector K1( numberOfElements );
    derivatives( t, Xcurrent, dXdtCurrent );
    for( int i = 0; i < numberOfElements; i++ )
    {
        K1[ i ] = h * dXdtCurrent[ i ];
    }

    // Evaluate K2 step in RK4 algorithm.
    Vector dXdt_K2( numberOfElements );
    Vector K2( numberOfElements );
    const double t_K2 = t + 0.5 * h;
    Vector X_K2( numberOfElements );
    for( int i = 0; i < numberOfElements; i++ )
    {
        X_K2[ i ] = Xcurrent[ i ] + 0.5 * K1[ i ];
    }
    derivatives( t_K2, X_K2, dXdt_K2 );
    for( int i = 0; i < numberOfElements; i++ )
    {
        K2[ i ] = h * dXdt_K2[ i ];
    }

    // Evaluate K3 step in RK4 algorithm.
    Vector dXdt_K3( numberOfElements );
    Vector K3( numberOfElements );
    const double t_K3 = t + 0.5 * h;
    Vector X_K3( numberOfElements );
    for( int i = 0; i < numberOfElements; i++ )
    {
        X_K3[ i ] = Xcurrent[ i ] + 0.5 * K2[ i ];
    }
    derivatives( t_K3, X_K3, dXdt_K3 );
    for( int i = 0; i < numberOfElements; i++ )
    {
        K3[ i ] = h * dXdt_K3[ i ];
    }

    // Evaluate K4 step in RK4 algorithm.
    Vector dXdt_K4( numberOfElements );
    Vector K4( numberOfElements );
    const double t_K4 = t + h;
    Vector X_K4( numberOfElements );
    for( int i = 0; i < numberOfElements; i++ )
    {
        X_K4[ i ] = Xcurrent[ i ] + K3[ i ];
    }
    derivatives( t_K4, X_K4, dXdt_K4 );
    for( int i = 0; i < numberOfElements; i++ )
    {
        K4[ i ] = h * dXdt_K4[ i ];
    }

    // Final step, evaluate the weighted summation.
    for( int i = 0; i < numberOfElements; i++ )
    {
        Xnext[ i ] = Xcurrent[ i ]
                        + ( 1.0 / 6.0 ) * K1[ i ]
                        + ( 1.0 / 3.0 ) * K2[ i ]
                        + ( 1.0 / 3.0 ) * K3[ i ]
                        + ( 1.0 / 6.0 ) * K4[ i ];
    }
}

} // namespace naos

#endif // RK4_HPP
