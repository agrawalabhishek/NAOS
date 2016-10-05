/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef RK4_HPP
#define RK4_HPP

#include <vector>

#include "NAOS/computeEllipsoidSurfaceGravitationalAcceleration.hpp"
#include "NAOS/constants.hpp"

namespace naos
{

//! Bulirsch-Stoer integrator
/*!
 * The bulirsch-stoer integrator routine to integrate a set of
 * first order differential equations. The equations have to be provided in a seperate header file
 * in a struct. See orbiterEquationsOfMotion.hpp for an example. This function is the ultimate wrap.
 * It makes use of the mid point rule method for actual integration, a routine for adaptive step size
 * and order control, and a routine for richardson extrapolation which is used for step size control.
 *
 */
template< typename Vector, class orbiterEquationsOfMotion >
void bulirschStoer(
    const double alpha,
    const double beta,
    const double gamma,
    const double gravParameter,
    const double Wmagnitude,
    Vector &Xcurrent,
    const double t,
    double &stepSize,
    const Vector &subStepSequence,
    double &k,
    const double maxK,
    Vector &Xnext,
    orbiterEquationsOfMotion &derivatives )
{
    double H = stepSize;
    const int numberOfElements = Xcurrent.size( );
}

//! modified mid point rule (MMR) evaluation
/*!
 * This function will advance the state vector from a point 't' in time to 't + H' in a sequence
 * of 'n' substeps each of same size. For more details, refer to "numerical recipes in C" book.
 * This routine does not have the authority to modify the stepSize or the value 'k'.
 */
template< typename Vector, class orbiterEquationsOfMotion >
void modifiedMidPointRule(
    const double alpha,
    const double beta,
    const double gamma,
    const double gravParameter,
    const double Wmagnitude,
    Vector &Xcurrent,
    const double t,
    const double stepSize,
    const Vector &subStepSequence,
    const double k,
    Vector &Xnext )
{
    const int numberOfSubSteps = subStepSequence[ k ];
    const int h = stepSize / numberOfSubSteps;

    // intermediate variable declaration
    const double sizeOfVector = Xcurrent.size( );
    Vector dXdt( sizeOfVector );
    // Z is a 2D matrix, rows = number of substeps+1, columns = state vector size
    std::vector< std::vector< double > > Z( numberOfSubSteps + 1,
                                            std::vector< double >( sizeOfVector ) );

    // Zeroth step of MMR
    Vector state = Xcurrent;
    for( int i = 0; i < sizeOfVector; i++ )
    {
        Z[ 0 ][ i ] = state[ i ];
    }

    // Performing computations for the first step of MMR

    // Evaluate gravitational acceleration values at the current state values
    Vector3 currentGravAcceleration( 3 );
    computeEllipsoidGravitationalAcceleration(
        alpha,
        beta,
        gamma,
        gravParameter,
        Z[ 0 ][ xPositionIndex ],
        Z[ 0 ][ yPositionIndex ],
        Z[ 0 ][ zPositionIndex ],
        currentGravAcceleration );

    // Create an object of the struct containing the equations of motion (in this case the eom
    // for an orbiter around a uniformly rotating tri-axial ellipsoid) and initialize it to the
    // values of the ellipsoidal asteroid's current gravitational accelerations
    orbiterEquationsOfMotion derivatives( currentGravAcceleration, Wmagnitude );

    // evaluate the right hand side of the differential equations and get the values for the
    // derivatives
    derivatives( t, Z[ 0 ], dXdt );

    // finally, compute the first step of MMR
    for( int i = 0; i < sizeOfVector; i++ )
    {
        Z[ 1 ][ i ] = Z[ 0 ][ i ] + h * dXdt[ i ];
    }

    // Now compute all the general steps until, except the last step

    // Again following the notation in the book, "numerical recipes in C", the general MMR steps
    // are calculated
    for( int m = 1; m <= numberOfSubSteps - 1; m++ )
    {
        // get the gravitational acceleration for the intermediate state vector computations
        computeEllipsoidGravitationalAcceleration(
            alpha,
            beta,
            gamma,
            gravParameter,
            Z[ m ][ xPositionIndex ],
            Z[ m ][ yPositionIndex ],
            Z[ m ][ zPositionIndex ],
            currentGravAcceleration );

        // get the value of the derivatives or EOMs
        orbiterEquationsOfMotion IntermediateDerivatives( currentGravAcceleration, Wmagnitude );
        double timeNow = t + m * h;
        IntermediateDerivatives( timeNow, Z[ m ], dXdt );

        // compute the intermediate step now
        for( int i = 0; i < sizeOfVector; i++ )
        {
            Z[ m + 1 ][ i ] = Z[ m - 1 ][ i ] + 2.0 * h * dXdt[ i ];
        }
    }

    // now calculate the final step of the MMR, which is also the final state value
    Vector finalStateVector( sizeOfVector );

    // get the gravitational acceleration values
    computeEllipsoidGravitationalAcceleration(
        alpha,
        beta,
        gamma,
        gravParameter,
        Z[ numberOfSubSteps ][ xPositionIndex ],
        Z[ numberOfSubSteps ][ yPositionIndex ],
        Z[ numberOfSubSteps ][ zPositionIndex ],
        currentGravAcceleration );

    // compute the final derivatives
    orbiterEquationsOfMotion finalDerivatives( currentGravAcceleration, Wmagnitude );
    double finalTime = t + stepSize;
    finalDerivatives( finalTime, Z[ numberOfSubSteps ], dXdt );

    // compute the final step now
    for( int i = 0; i < sizeOfVector; i++ )
    {
        finalStateVector[ i ] = 0.5 * ( Z[ numberOfSubSteps ][ i ]
                                        + Z[ numberOfSubSteps - 1 ][ i ]
                                        + h * dXdt[ i ] );
    }

    // return the final state vector that was computed
    Xnext = finalStateVector;
}

} // namespace naos

#endif // RK4_HPP
