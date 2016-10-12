/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef BULIRSCH_STOER_HPP
#define BULIRSCH_STOER_HPP

#include <algorithm>
#include <vector>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <sstream>

#include "NAOS/ellipsoidGravitationalAcceleration.hpp"
#include "NAOS/constants.hpp"
#include "NAOS/orbiterEquationsOfMotion.hpp"
#include "NAOS/misc.hpp"

namespace naos
{

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
    std::vector< int > subStepSequence,
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

//! Polynomial extrapolation routine
/*!
 * This routine performs the richardson polynomial extrapolation to extrapolate the values computed
 * by the modified mid point rule, to a point where the number of substeps are infinity and h=0.
 * refer to numerical recipes in C book for more details. Also refer to the aitken-neville algorithm
 * equation 9.10 in the book "solving ordinary differential equations 1" by springer publications.
 *
 */
struct stateVector
{
    double xPosition;
    double yPosition;
    double zPosition;
    double xVelocity;
    double yVelocity;
    double zVelocity;
};

template< typename Vector >
void polynomialExtrapolation(
        const double k,
        std::vector< std::vector< stateVector > > &table,
        const Vector subStepSequence )
{
    double nK = subStepSequence[ k ];
    for( int j = 0; j <= k - 1; j++ )
    {
        // slight change in the computation of the nKminusJ term since the index j starts with 0,
        // and not 1.
        double nKminusJ = subStepSequence[ k - ( j + 1 ) ];
        double tempSquareTerm = ( nK / nKminusJ ) * ( nK / nKminusJ );

        // Perform extrapolation, for each element of the state vector
        table[ k ][ j + 1 ].xPosition =
                table[ k ][ j ].xPosition
                    + ( ( table[ k ][ j ].xPosition - table[ k - 1 ][ j ].xPosition )
                        / ( tempSquareTerm - 1.0 ) );

        table[ k ][ j + 1 ].yPosition =
                table[ k ][ j ].yPosition
                    + ( ( table[ k ][ j ].yPosition - table[ k - 1 ][ j ].yPosition )
                        / ( tempSquareTerm - 1.0 ) );

        table[ k ][ j + 1 ].zPosition =
                table[ k ][ j ].zPosition
                    + ( ( table[ k ][ j ].zPosition - table[ k - 1 ][ j ].zPosition )
                        / ( tempSquareTerm - 1.0 ) );

        table[ k ][ j + 1 ].xVelocity =
                table[ k ][ j ].xVelocity
                    + ( ( table[ k ][ j ].xVelocity - table[ k - 1 ][ j ].xVelocity )
                        / ( tempSquareTerm - 1.0 ) );

        table[ k ][ j + 1 ].yVelocity =
                table[ k ][ j ].yVelocity
                    + ( ( table[ k ][ j ].yVelocity - table[ k - 1 ][ j ].yVelocity )
                        / ( tempSquareTerm - 1.0 ) );

        table[ k ][ j + 1 ].zVelocity =
                table[ k ][ j ].zVelocity
                    + ( ( table[ k ][ j ].zVelocity - table[ k - 1 ][ j ].zVelocity )
                        / ( tempSquareTerm - 1.0 ) );
    }
}

//! Bulirsch-Stoer integrator
/*!
 * The bulirsch-stoer integrator routine to integrate a set of
 * first order differential equations. The equations have to be provided in a seperate header file
 * in a struct. See orbiterEquationsOfMotion.hpp for an example. This function is the ultimate wrap.
 * It makes use of the mid point rule method for actual integration, a routine for adaptive step size
 * and order control, and a routine for richardson extrapolation which is used for step size control.
 *
 */
template< typename Vector >
Vector bulirschStoerIntegrator(
    const double alpha,
    const double beta,
    const double gamma,
    const double gravParameter,
    const double Wmagnitude,
    Vector &Xcurrent,
    const double t,
    double &stepSize,
    double &targetK )
{
    const int numberOfElements = Xcurrent.size( );

    // declare the output state vector
    Vector Xnext( numberOfElements );

    // some intermediate boolean flags used in step size and order control
    static bool stepReject = false;

    // machine precision value, difference between 1.0 and the next big value, upto 10 significant
    // digits after the decimal point.
    const double eps = 10.0 * std::numeric_limits< double >::epsilon( );

    // declare absolute and relative tolerances for the adaptive step size and order control method
    const double absoluteTolerance = 10.0e-5;
    const double relativeTolerance = 10.0e-5;

    // substep parameters. The sub step series is the deuflhard sequence. Refer to numerical
    // recipes in C for more details.
    // note k value starts from 0 and ends at 8, so totally 9 elements.
    const int maxK = 8;
    std::vector< int > subStepSequence( maxK + 1 );
    for( int k = 0; k <= maxK; k++ )
    {
        subStepSequence[ k ] = 2 * ( k + 1 );
    }

    // set a cold start flag to start the while loop below for the first time
    bool coldStart = true;

    // initialize the table to be used by the polynomialExtrapolation routine later on.
    // the table is a general 2D matrix of size maxK+1 by maxK+1 but not all columns will be used.
    // Only the lower triangle will be utilized.
    std::vector< std::vector< stateVector > > table(
                                                maxK + 1, std::vector< stateVector >( maxK + 1 ) );

    // declare a counter to count the number of times the order value for a integration step was
    // overridden
    int countTargetKOverride = 0;

    // calculate the cost and the work for obtaining column k of the table
    std::vector< double > cost( maxK + 1 );
    cost[ 0 ] = subStepSequence[ 0 ] + 1;
    for( int i = 0; i <= maxK; i++ )
    {
        // compute the cost using the recursive formulas given in the numerical recipes in C book
        cost[ i + 1 ] = cost[ i ] + subStepSequence[ i + 1 ];
    }

    // start the while loop that integrates, checks for local error and modifies the stepsize
    while( coldStart || stepReject )
    {
        // make sure step size (H) is positive
        stepSize = std::abs( stepSize );

        // reset the flags to appropriate values
        stepReject = false;
        coldStart = false;

        // check whether step size underflows, i.e. its value becomes less than machine precision
        if( stepSize < eps )
        {
            std::ostringstream errorMessage;
            errorMessage << std::endl;
            errorMessage << "bulirsch stoer integrator step size underflow" << std::endl;
            errorMessage << std::endl;
            throw std::runtime_error( errorMessage.str( ) );
        }

        // ensure that the target k is in the range: 2 <= targetK <= maxK-1, this was taken from
        // the book solving ordinary differential equations 1 by springer publications.
        if( targetK < 2 )
        {
            // computation of error for the k-1 order will not be possible if target order k = 1
            targetK = 2;
            countTargetKOverride++;
        }
        else if( targetK > maxK - 1 )
        {
            // computation of error for the k+1 order will not be possible if target k > maxK-1
            targetK = maxK - 1;
            countTargetKOverride++;
        }

        // check the count value for targetK override, to avoid the integrator from getting stuck
        if( countTargetKOverride > 10 )
        {
            std::ostringstream errorMessage;
            errorMessage << std::endl;
            errorMessage << "ERROR: Order target value override exceeded set threshold limit ";
            errorMessage << "in Bulirsch stoer interator routine. Manual intervention needed. ";
            errorMessage << "Check routine bulirschStoer.hpp" << std::endl;
            errorMessage << std::endl;
            throw std::runtime_error( errorMessage.str( ) );
        }

        // declare the integration error variable
        std::vector< double > errorValue( maxK + 1 );
        errorValue[ 0 ] = 0.0;

        // declare constant factors to be used in the calculation of the work done to obtain a
        // column k of the table
        // notation of numerical recipes in C book
        const double S1 = 0.65;
        const double S2 = 0.94;
        const double S3 = 0.02;
        const double S4 = 4.0;
        std::vector< double > work( maxK + 1 );
        work[ 0 ] = 0;

        // declare the vector to store estimated optimal values for the step sizes for various
        // orders k
        std::vector< double > optimalStepSize( maxK + 1 );
        optimalStepSize[ 0 ] = stepSize; // this is equivalent to storing a dummy value

        // declare the variable to store the optimal order values, determined later on based on the
        // value of error in extension with the value of the work required for that particular order
        double optimalOrderK = targetK;
        double newStepSize = stepSize;

        // start the order or 'k' loop, k value can have a maximum cycle range from 0 through 8.
        for( int k = 0; k <= targetK + 1; k++ )
        {
            // perform the modified mid point rule evaluation for the values of k in the loop
            modifiedMidPointRule< Vector, eomOrbiterURE >( alpha,
                                                           beta,
                                                           gamma,
                                                           gravParameter,
                                                           Wmagnitude,
                                                           Xcurrent,
                                                           t,
                                                           stepSize,
                                                           subStepSequence,
                                                           k,
                                                           Xnext );

            // store the result in the first column of the table and the corresponding k'th row.
            // index value followed from the numerical recipes in C book.
            table[ k ][ 0 ].xPosition = Xnext[ xPositionIndex ];
            table[ k ][ 0 ].yPosition = Xnext[ yPositionIndex ];
            table[ k ][ 0 ].zPosition = Xnext[ zPositionIndex ];
            table[ k ][ 0 ].xVelocity = Xnext[ xVelocityIndex ];
            table[ k ][ 0 ].yVelocity = Xnext[ yVelocityIndex ];
            table[ k ][ 0 ].zVelocity = Xnext[ zVelocityIndex ];

            // check if the value of k is greater than 0, then perform polynomial extrapolation.
            // if k=0, then the next step in the for loop is performed following which the if loop
            // is activated.
            if( k != 0 )
            {
                // use the extrapolation function to calculate and store the terms in the k'th row
                // of the table
                polynomialExtrapolation( k, table, subStepSequence );

                // declare the scaling variable used in normalized error computation
                std::vector< double > scale( numberOfElements );

                // calculate the normalized error value for each order k
                for( int index = 0; index < numberOfElements; index++ )
                {
                    // compute the scaling factors to perform the normalized norm later
                    scale[ index ] = absoluteTolerance
                        + relativeTolerance * std::max( Xcurrent[ index ], Xnext[ index ] );
                }
                // compute the summation part of the normalized error formula, eqn 17.2.9
                // of the numerical recipes in C book
                double xPositionScaledError = ( Xnext[ xPositionIndex ]
                                        - table[ 0 ][ 0 ].xPosition ) / scale[ xPositionIndex ];

                double yPositionScaledError = ( Xnext[ yPositionIndex ]
                                        - table[ 0 ][ 0 ].yPosition ) / scale[ yPositionIndex ];

                double zPositionScaledError = ( Xnext[ zPositionIndex ]
                                        - table[ 0 ][ 0 ].zPosition ) / scale[ zPositionIndex ];

                double xVelocityScaledError = ( Xnext[ xVelocityIndex ]
                                        - table[ 0 ][ 0 ].xVelocity ) / scale[ xVelocityIndex ];

                double yVelocityScaledError = ( Xnext[ yVelocityIndex ]
                                        - table[ 0 ][ 0 ].yVelocity ) / scale[ yVelocityIndex ];

                double zVelocityScaledError = ( Xnext[ zVelocityIndex ]
                                        - table[ 0 ][ 0 ].zVelocity ) / scale[ zVelocityIndex ];

                errorValue[ k ] = xPositionScaledError * xPositionScaledError
                                + yPositionScaledError * yPositionScaledError
                                + zPositionScaledError * zPositionScaledError
                                + xVelocityScaledError * xVelocityScaledError
                                + yVelocityScaledError * yVelocityScaledError
                                + zVelocityScaledError * zVelocityScaledError;

                errorValue[ k ] = std::sqrt( errorValue[ k ] / numberOfElements );

                // compute optimal step size for the current non zero order value, following which
                // the work for the current order will be calculated based on the optimal order
                // value. The formulas and procedure are taken from numerical recipes in C book
                double exponentValue = 1.0 / ( 2.0 * k + 1.0 );
                // equation 17.3.22 in numerical recipes in C book
                double minimumFactor = std::pow( S3, exponentValue );
                double factor = 0.0;
                if( errorValue[ k ] == 0.0 )
                {
                    // the upper limit in 17.3.22
                    factor = 1.0 / minimumFactor;
                }
                else
                {
                    // equation 17.3.11 in numerical recipes in C
                    factor = S1 * std::pow( ( S2 / errorValue[ k ] ), exponentValue );
                    // evaluating the value for the lower limit in equation 17.3.22
                    factor = std::max( ( minimumFactor / S4 ),
                                         std::min( ( 1.0 / minimumFactor ), factor ) );
                }
                // now obtain the optimal step size for current order k, numerical recipes in C
                optimalStepSize[ k ] = std::abs( stepSize * factor );
                // and finally obtain the work done for order k
                work[ k ] = cost[ k ] / optimalStepSize[ k ];

                // check for convergence in order targetk - 1
                if( k == targetK - 1 )
                {
                    if( errorValue[ targetK - 1 ] <= 1.0 )
                    {
                        // accept table[k-1][k-1] as the numerical solution for the integration step
                        Xnext[ xPositionIndex ] = table[ targetK - 1 ][ targetK - 1 ].xPosition;
                        Xnext[ yPositionIndex ] = table[ targetK - 1 ][ targetK - 1 ].yPosition;
                        Xnext[ zPositionIndex ] = table[ targetK - 1 ][ targetK - 1 ].zPosition;
                        Xnext[ xVelocityIndex ] = table[ targetK - 1 ][ targetK - 1 ].xVelocity;
                        Xnext[ yVelocityIndex ] = table[ targetK - 1 ][ targetK - 1 ].yVelocity;
                        Xnext[ zVelocityIndex ] = table[ targetK - 1 ][ targetK - 1 ].zVelocity;

                        // set boolean flag for step reject to false
                        stepReject = false;

                        // choose a new order and step value for the next integration step
                        optimalOrderK = targetK - 1;
                        if( work[ targetK - 2 ] < 0.8 * work[ targetK - 1 ] )
                        {
                            optimalOrderK = targetK - 2;
                        }
                        if( work[ targetK - 1 ] < 0.9 * work[ targetK - 2 ] )
                        {
                            optimalOrderK = targetK;
                        }

                        if( optimalOrderK == targetK - 1 || optimalOrderK == targetK - 2 )
                        {
                            newStepSize = optimalStepSize[ optimalOrderK ];
                        }
                        else
                        {
                            newStepSize = optimalStepSize[ targetK - 1 ]
                                            * ( cost[ targetK ] / cost[ targetK - 1 ] );
                        }
                        break; // end the for loop and enter the final steps of while loop
                    }
                    else
                    {
                        // this is the part if convergence wasn't achieved in order targetK - 1
                        double tempTerm1 = ( subStepSequence[ targetK ] / subStepSequence[ 0 ] );
                        double tempSquareTerm1 = tempTerm1 * tempTerm1;
                        double tempTerm2 = ( subStepSequence[ targetK + 1 ] / subStepSequence[ 0 ] );
                        double tempSquareTerm2 = tempTerm2 * tempTerm2;
                        if( errorValue[ targetK - 1 ] > ( tempSquareTerm1 * tempSquareTerm2 ) )
                        {
                            // set the step rejection flag to true
                            stepReject = true;

                            // choose an appropriate step and order value to restart this
                            // integration step. note after step rejection, order and step size
                            // can't be increased. follow equation 17.3.14 and 17.3.15
                            optimalOrderK = targetK - 1;
                            if( work[ targetK - 2 ] < 0.8 * work[ targetK - 1 ] )
                            {
                                optimalOrderK = targetK - 2;
                            }
                            newStepSize = optimalStepSize[ optimalOrderK ];

                            // break the for loop and enter into the final steps of the while loop
                            break;
                        }
                        else
                        {
                            // continue with the computation of the table and go to the check
                            // for convergence in order targetK.
                            // The for loop continues with k=targetK
                            continue;
                        }
                    }
                } // end of if loop that checks for convergence in order targetk - 1

                // check for convergence in order targetk
                if( k == targetK )
                {
                    if( errorValue[ targetK ] <= 1.0 )
                    {
                        // step integration is accepted
                        stepReject = false;

                        // accept table[targetK][targetK] as the numerical solution
                        // for the integration step
                        Xnext[ xPositionIndex ] = table[ targetK ][ targetK ].xPosition;
                        Xnext[ yPositionIndex ] = table[ targetK ][ targetK ].yPosition;
                        Xnext[ zPositionIndex ] = table[ targetK ][ targetK ].zPosition;
                        Xnext[ xVelocityIndex ] = table[ targetK ][ targetK ].xVelocity;
                        Xnext[ yVelocityIndex ] = table[ targetK ][ targetK ].yVelocity;
                        Xnext[ zVelocityIndex ] = table[ targetK ][ targetK ].zVelocity;

                        // compute the optimal order and step size for the next integration step
                        optimalOrderK = targetK;
                        if( work[ targetK - 1 ] < 0.8 * work[ targetK ] )
                        {
                            optimalOrderK = targetK - 1;
                        }
                        if( work[ targetK ] < 0.9 * work[ targetK - 1 ] )
                        {
                            optimalOrderK = targetK + 1;
                        }

                        if( optimalOrderK == targetK || optimalOrderK == targetK - 1 )
                        {
                            newStepSize = optimalStepSize[ optimalOrderK ];
                        }
                        else
                        {
                            newStepSize = optimalStepSize[ targetK ]
                                            * ( cost[ targetK + 1 ] / cost[ targetK ] );
                        }

                        // break the for loop and go to the final steps of the while loop
                        break;
                    }
                    else
                    {
                        // convergence in target k was not achieved
                        double tempTerm = subStepSequence[ targetK + 1 ] / subStepSequence[ 0 ];
                        double tempSquareTerm = tempTerm * tempTerm;
                        if( errorValue[ targetK ] > tempSquareTerm )
                        {
                            stepReject = true;

                            // choose new order for repeat of while loop
                            optimalOrderK = targetK - 1;

                            // choose new step size
                            newStepSize = optimalStepSize[ optimalOrderK ];

                            // break out of the for loop
                            break;
                        }
                        else
                        {
                            // continue with the next iteration of the for loop
                            continue;
                        }
                    }
                } // end of if loop that checks for convergence in the order targetK

                // check for convergence in order targetk + 1
                if( k == targetK + 1 )
                {
                    if( errorValue[ targetK + 1 ] <= 1.0 )
                    {
                        // step is accepted
                        stepReject = false;

                        // accept table[targetK + 1][targetK + 1] as the numerical solution
                        // for the integration step
                        Xnext[ xPositionIndex ] = table[ targetK + 1 ][ targetK + 1 ].xPosition;
                        Xnext[ yPositionIndex ] = table[ targetK + 1 ][ targetK + 1 ].yPosition;
                        Xnext[ zPositionIndex ] = table[ targetK + 1 ][ targetK + 1 ].zPosition;
                        Xnext[ xVelocityIndex ] = table[ targetK + 1 ][ targetK + 1 ].xVelocity;
                        Xnext[ yVelocityIndex ] = table[ targetK + 1 ][ targetK + 1 ].yVelocity;
                        Xnext[ zVelocityIndex ] = table[ targetK + 1 ][ targetK + 1 ].zVelocity;

                        // select optimal step size and order for the next integration step
                        optimalOrderK = targetK;
                        if( work[ targetK - 1 ] < 0.8 * work[ targetK ] )
                        {
                            optimalOrderK = targetK - 1;
                        }
                        if( work[ targetK + 1 ] < 0.9 * work[ optimalOrderK ] )
                        {
                            optimalOrderK = targetK + 1;
                        }

                        if( optimalOrderK == targetK || optimalOrderK == targetK - 1 )
                        {
                            newStepSize = optimalStepSize[ optimalOrderK ];
                        }
                        else
                        {
                            newStepSize = optimalStepSize[ targetK ]
                                            * ( cost[ targetK + 1 ] / cost[ targetK ] );
                        }

                        // break the for loop and go to final steps of while loop
                        break;
                    }
                    else
                    {
                        // step rejected
                        stepReject = true;

                        // choose new order for repeat of while loop
                        optimalOrderK = targetK - 1;

                        // choose new step size
                        newStepSize = optimalStepSize[ optimalOrderK ];

                        // break out of the for loop
                        break;
                    }
                } // end of if loop that checks for convergence in the order targetK+1
            } // end of the if loop which performs polynomial extrapolation, error computation, order and step size control
        } // end of the for loop which runs over all desired orders

        // all break statements lead to this point.
        //the following values for order and step size will be used for the next integration step
        targetK = optimalOrderK;
        stepSize = newStepSize;
    } // end of while loop for single integration step

    return Xnext;
} // scope end of the burlirschStoerIntegrator function

} // namespace naos

#endif // BULIRSCH_STOER_HPP
