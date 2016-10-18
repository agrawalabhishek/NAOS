/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef RK54_HPP
#define RK54_HPP

#include <vector>
#include <limits>
#include <stdexcept>

#include "NAOS/constants.hpp"
#include "NAOS/ellipsoidGravitationalAcceleration.hpp"

namespace naos
{

//! RK5(4) integrator
/*!
 * The RK5(4) integrator routine to integrate a set of
 * first order differential equations. The equations have to be provided in a seperate header file
 * in a struct. See orbiterEquationsOfMotion.hpp for an example.
 *
 * @param[in]  Xcurrent      currently known state vector value i.e. before the integration
 * @param[in]  t             current time value for which the state vector is known
 * @param[in]  stepSize      the step size of integration
 * @param[out] Xnext         the state vector value after integration i.e. at next time step
 * @param[in]  derivatives   An object of the struct containing the differential equations. The
 *                           object should have already gotten the initializing values before
 *                           passing it to the argument of the function rk4.
 */
template< typename Vector, class orbiterEquationsOfMotion >
void rk54(
    Vector &Xcurrent,
    const double t,
    const double stepSize,
    Vector &Xnext,
    Vector &delta,
    orbiterEquationsOfMotion &derivatives )
{
    const double h = stepSize;
    const int numberOfElements = Xcurrent.size( );

    // declare the coefficients for the RK5(4) algorithm
    const std::vector< double > c { 0.0,
                                    1.0 / 5.0,
                                    3.0 / 10.0,
                                    4.0 / 5.0,
                                    8.0 / 9.0,
                                    1.0,
                                    1.0 };

    const std::vector< std::vector< double > > a {
        { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { 1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { 3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0 },
        { 44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0 },
        { 19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0 },
        { 9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0 },
        { 35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0 } };

    const std::vector< double > b { 35.0 / 384.0,
                                    0.0,
                                    500.0 / 1113.0,
                                    125.0 / 192.0,
                                    -2187.0 / 6784.0,
                                    11.0 / 84.0,
                                    0.0 };

    const std::vector< double > bStar { 5179.0 / 57600.0,
                                        0.0,
                                        7571.0 / 16695.0,
                                        393.0 / 640.0,
                                        -92097.0 / 339200.0,
                                        187.0 / 2100.0,
                                        1.0 / 40.0 };

    // declare the vector that stores the values for the right hand side of the EOMs
    Vector dXdt( numberOfElements );

    // Evaluate K1 step in RK5(4) algorithm.
    Vector K1( numberOfElements );
    derivatives( t, Xcurrent, dXdt );
    for( int i = 0; i < numberOfElements; i++ )
    {
        K1[ i ] = h * dXdt[ i ];
    }

    // Evaluate K2 step in RK5(4) algorithm.
    Vector K2( numberOfElements );
    const double t_K2 = t + c[ 1 ] * h;
    Vector X_K2( numberOfElements );
    for( int i = 0; i < numberOfElements; i++ )
    {
        X_K2[ i ] = Xcurrent[ i ] + a[ 1 ][ 0 ] * K1[ i ];
    }
    derivatives( t_K2, X_K2, dXdt );
    for( int i = 0; i < numberOfElements; i++ )
    {
        K2[ i ] = h * dXdt[ i ];
    }

    // Evaluate K3 step in RK5(4) algorithm.
    Vector K3( numberOfElements );
    const double t_K3 = t + c[ 2 ] * h;
    Vector X_K3( numberOfElements );
    for( int i = 0; i < numberOfElements; i++ )
    {
        X_K3[ i ] = Xcurrent[ i ] + a[ 2 ][ 0 ] * K1[ i ]
                                  + a[ 2 ][ 1 ] * K2[ i ];
    }
    derivatives( t_K3, X_K3, dXdt );
    for( int i = 0; i < numberOfElements; i++ )
    {
        K3[ i ] = h * dXdt[ i ];
    }

    // Evaluate K4 step in RK5(4) algorithm.
    Vector K4( numberOfElements );
    const double t_K4 = t + c[ 3 ] * h;
    Vector X_K4( numberOfElements );
    for( int i = 0; i < numberOfElements; i++ )
    {
        X_K4[ i ] = Xcurrent[ i ] + a[ 3 ][ 0 ] * K1[ i ]
                                  + a[ 3 ][ 1 ] * K2[ i ]
                                  + a[ 3 ][ 2 ] * K3[ i ];
    }
    derivatives( t_K4, X_K4, dXdt );
    for( int i = 0; i < numberOfElements; i++ )
    {
        K4[ i ] = h * dXdt[ i ];
    }

    // Evaluate K5 step in RK5(4) algorithm.
    Vector K5( numberOfElements );
    const double t_K5 = t + c[ 4 ] * h;
    Vector X_K5( numberOfElements );
    for( int i = 0; i < numberOfElements; i++ )
    {
        X_K5[ i ] = Xcurrent[ i ] + a[ 4 ][ 0 ] * K1[ i ]
                                  + a[ 4 ][ 1 ] * K2[ i ]
                                  + a[ 4 ][ 2 ] * K3[ i ]
                                  + a[ 4 ][ 3 ] * K4[ i ];
    }
    derivatives( t_K5, X_K5, dXdt );
    for( int i = 0; i < numberOfElements; i++ )
    {
        K5[ i ] = h * dXdt[ i ];
    }

    // Evaluate K6 step in RK5(4) algorithm.
    Vector K6( numberOfElements );
    const double t_K6 = t + c[ 5 ] * h;
    Vector X_K6( numberOfElements );
    for( int i = 0; i < numberOfElements; i++ )
    {
        X_K6[ i ] = Xcurrent[ i ] + a[ 5 ][ 0 ] * K1[ i ]
                                  + a[ 5 ][ 1 ] * K2[ i ]
                                  + a[ 5 ][ 2 ] * K3[ i ]
                                  + a[ 5 ][ 3 ] * K4[ i ]
                                  + a[ 5 ][ 4 ] * K5[ i ];
    }
    derivatives( t_K6, X_K6, dXdt );
    for( int i = 0; i < numberOfElements; i++ )
    {
        K6[ i ] = h * dXdt[ i ];
    }

    // Final step, evaluate the weighted summation.
    for( int i = 0; i < numberOfElements; i++ )
    {
        Xnext[ i ] = Xcurrent[ i ]
                        + b[ 0 ] * K1[ i ]
                        + b[ 1 ] * K2[ i ]
                        + b[ 2 ] * K3[ i ]
                        + b[ 3 ] * K4[ i ]
                        + b[ 4 ] * K5[ i ]
                        + b[ 5 ] * K6[ i ];
    }

    // calculate the embedded fourth order solution
    Vector XnextStar = Xnext;
    for( int i = 0; i < numberOfElements; i++ )
    {
        XnextStar[ i ] = Xcurrent[ i ]
                        + bStar[ 0 ] * K1[ i ]
                        + bStar[ 1 ] * K2[ i ]
                        + bStar[ 2 ] * K3[ i ]
                        + bStar[ 3 ] * K4[ i ]
                        + bStar[ 4 ] * K5[ i ]
                        + bStar[ 5 ] * K6[ i ];
    }

    // calculate delta, the difference vector
    for( int i = 0; i < numberOfElements; i++ )
    {
        delta[ i ] = Xnext[ i ] - XnextStar[ i ];
    }
}

//! RK54 integrator wrapper
/*!
 * to compute gravitational acceleration values for the current state
 * vector and then perform the step integration using RK54 routine.
 *
 * @param[in] alpha             largest semi axis of the ellipsoid
 * @param[in] beta              intermediate semi axis of the ellipsoid
 * @param[in] gamma             smallest semi axis of the ellipsoid
 * @param[in] gravParameter     gravitational parameter of the ellipsoid
 * @param[in] Wmagnitude        angular velcotiy magnitude of the ellipsoid
 * @param[in] Xcurrent          current known state vector value
 * @param[in] t                 time
 * @param[in] stepSize          step size of integration
 * @param[in] Xnext             next or integrated state vector
 */
template< typename Vector, class orbiterEquationsOfMotion >
void rk54Integrator(
    const double alpha,
    const double beta,
    const double gamma,
    const double gravParameter,
    const double Wmagnitude,
    Vector &Xcurrent,
    const double t,
    double &stepSize,
    Vector &Xnext )
{
    // specify tolerance values
    const double absoluteTolerance = 10.0e-5;
    const double relativeTolerance = 10.0e-5;

    // Evaluate gravitational acceleration values at the current state values
    Vector currentGravAcceleration( 3 );
    computeEllipsoidGravitationalAcceleration(
        alpha,
        beta,
        gamma,
        gravParameter,
        Xcurrent[ xPositionIndex ],
        Xcurrent[ yPositionIndex ],
        Xcurrent[ zPositionIndex ],
        currentGravAcceleration );

    // Create an object of the struct containing the equations of motion (in this case the eom
    // for an orbiter around a uniformly rotating tri-axial ellipsoid) and initialize it to the
    // values of the ellipsoidal asteroid's current gravitational accelerations
    orbiterEquationsOfMotion derivatives( currentGravAcceleration, Wmagnitude );

    // run the step integrator for using rk54 routine
    bool coldStart = true;
    bool stepReject = false;

    // Declare the output variables
    const double numberOfElements = Xcurrent.size( );
    Vector integratedState( numberOfElements );
    Vector delta( numberOfElements );

    // set safety factor, and max and min scaling factors for scaling the step size
    const double safetyFactor = 0.9;
    // new step size should not increase or decrease by a factor of 10 or 0.5 respectively.
    const double maxScale = 2.0;
    const double minScale = 0.5;
    double scalingFactor = 0.0;

    while( coldStart || stepReject )
    {
        coldStart = false;

        rk54< Vector, orbiterEquationsOfMotion >( Xcurrent,
                                                  t,
                                                  stepSize,
                                                  integratedState,
                                                  delta,
                                                  derivatives );

        // compute the error
        double errorValue = 0.0;
        Vector scale( numberOfElements );
        double sum = 0.0;
        for( int i = 0; i < numberOfElements; i++ )
        {
            scale[ i ] = absoluteTolerance
                        + relativeTolerance * std::max( Xcurrent[ i ], Xnext[ i ] );
            sum += ( delta[ i ] / scale[ i ] ) * ( delta[ i ] / scale[ i ] );
        }
        errorValue = std::sqrt( sum / numberOfElements );

        // check if error is <= 1.0
        if( errorValue <= 1.0 )
        {
            // accept the step
            stepReject = false;

            // calculate the stepsize value for the next integration step
            if( errorValue == 0.0 )
            {
                // choose the maximum scaling factor for the next step size
                scalingFactor = maxScale;
            }
            else
            {
                // choose the scaling factor as per equation 17.2.12 in numerical recipes in C book
                scalingFactor = safetyFactor * std::pow( ( 1.0 / errorValue ), ( 1.0 / 5.0 ) );
                if( scalingFactor < minScale )
                {
                    scalingFactor = minScale;
                }
                if( scalingFactor > maxScale )
                {
                    scalingFactor = maxScale;
                }
            }
            // if the previous step was rejected, the new step size for the next integration step
            // should not be increased
            if( stepReject )
            {
                stepSize = stepSize * std::min( scalingFactor, 1.0 );
                stepReject = false;
                break; // exit the while loop and return the integrated state
            }
            else
            {
                // choose the next step size
                stepSize = stepSize * scalingFactor;
                stepReject = false;
                break; // exit the while loop and return the integrated state
            }
        }
        else
        {
            // truncation error is large, redo the integration step with a reduced step size value
            stepReject = true;
            scalingFactor = safetyFactor * std::pow( ( 1.0 / errorValue ), ( 1.0 / 5.0 ) );
            // scaling factor of the step size shouldn't go below the min value
            scalingFactor = std::min( scalingFactor, minScale );
            stepSize = stepSize * scalingFactor;
        }
    } // break in while loop leads here

    // give back the final result
    Xnext = integratedState;
}

//! RK54 integrator general wrapper
/*!
 * perform the step integration for any set of differential equations using RK54 routine.
 *
 * @param[in] Xcurrent          current known state vector value
 * @param[in] t                 time
 * @param[in] stepSize          step size of integration
 * @param[in] Xnext             next or integrated state vector
 */
template< typename Vector, class equationsOfMotion >
void rk54GeneralIntegrator(
    Vector &Xcurrent,
    const double t,
    double &stepSize,
    Vector &Xnext,
    equationsOfMotion derivatives,
    bool &stepReject,
    double &previousErrorValue )
{
    // specify tolerance values
    const double absoluteTolerance = 10.0e-3;
    const double relativeTolerance = 10.0e-9;

    // run the step integrator for using rk54 routine
    bool coldStart = true;

    // Declare the output variables
    const double numberOfElements = Xcurrent.size( );
    Vector integratedState( numberOfElements );
    Vector delta( numberOfElements );

    // declare the error storing variable
    double errorValue = 0.0;

    // set safety factor, and max and min scaling factors for scaling the step size
    const double safetyFactor = 0.8;
    // new step size should not increase or decrease by a factor of 10 or 0.5 respectively.
    const double maxScale = 4.0;
    const double minScale = 0.1;
    double scalingFactor = 0.0;

    // terms to be used in the PI adaptive step size control algorithm
    const double alphaExponent = 1.0 / 5.0;
    const double betaExponent = 0.04;
    double oldErrorValue = previousErrorValue; // from previous integration step

    const double eps = 10.0 * std::numeric_limits< double >::epsilon( );

    while( coldStart || stepReject )
    {
        coldStart = false;

        if( stepSize < eps )
        {
            std::ostringstream errorMessage;
            errorMessage << std::endl;
            errorMessage << "ERROR: step size underflow!" << std::endl;
            errorMessage << std::endl;
            throw std::runtime_error( errorMessage.str( ) );
        }

        rk54< Vector, equationsOfMotion >( Xcurrent,
                                           t,
                                           stepSize,
                                           integratedState,
                                           delta,
                                           derivatives );

        // compute the error
        Vector scale( numberOfElements );
        double sum = 0.0;
        for( int i = 0; i < numberOfElements; i++ )
        {
            scale[ i ] = absoluteTolerance
                        + relativeTolerance * std::max( Xcurrent[ i ], integratedState[ i ] );
            sum += ( delta[ i ] / scale[ i ] ) * ( delta[ i ] / scale[ i ] );
        }
        errorValue = std::sqrt( sum / numberOfElements );

        // check if error is <= 1.0
        if( errorValue <= 1.0 )
        {
            // accept the step
            stepReject = false;

            // calculate the stepsize value for the next integration step
            if( errorValue == 0.0 )
            {
                // choose the maximum scaling factor for the next step size
                scalingFactor = maxScale;
            }
            else
            {
                // choose the scaling factor as per equation 17.2.12 in numerical recipes in C book
                scalingFactor = safetyFactor * std::pow( ( 1.0 / errorValue ), alphaExponent )
                                * std::pow( oldErrorValue, betaExponent );
                if( scalingFactor < minScale )
                {
                    scalingFactor = minScale;
                }
                if( scalingFactor > maxScale )
                {
                    scalingFactor = maxScale;
                }
            }
            // if the previous step was rejected, the new step size for the next integration step
            // should not be increased
            if( stepReject )
            {
                scalingFactor = std::min( scalingFactor, 1.0 );
                stepSize = stepSize * scalingFactor;
                stepReject = false;
                break; // exit the while loop and return the integrated state
            }
            else
            {
                // choose the next step size
                stepSize = stepSize * scalingFactor;
                stepReject = false;
                break; // exit the while loop and return the integrated state
            }
        }
        else
        {
            // truncation error is large, redo the integration step with a reduced step size value
            stepReject = true;
            scalingFactor = safetyFactor * std::pow( ( 1.0 / errorValue ), alphaExponent );
            // scaling factor of the step size shouldn't go below the min value
            scalingFactor = std::max( scalingFactor, minScale );
            stepSize = stepSize * scalingFactor;
        }
    } // break in while loop leads here

    previousErrorValue = errorValue;

    // give back the final result
    Xnext = integratedState;
}

} // namespace naos

#endif // RK54_HPP
