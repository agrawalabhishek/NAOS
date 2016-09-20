/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef NORMALIZED_POTENTIAL_URE_HPP
#define NORMALIZED_POTENTIAL_URE_HPP

#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <vector>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_errno.h>

#include "NAOS/cubicRoot.hpp"
#include "NAOS/constants.hpp"

namespace naos
{

//! compute the normalized gravitational potential for a URE
/*!
 * Given a set of position coordinates, calculate the gravitational potential for a uniformly
 * rotating ellipsoid.
 *
 * @param[in] positionVector        dimensional position coordinates
 * @param[in] alpha, beta, gamma    largest to smallest semi axes of ellipsoid
 * @param[in] density               density of the ellipsoid
 * @param[in] Wmagnitude            magnitude of angular velocity of ellipsoid
 * @return    U                     the value of the normalized potential at a given point
 */
template < typename Vector3 >
double computeNormalizedPotentialURE(
    const Vector3 &positionVector,
    const double alpha,
    const double beta,
    const double gamma,
    const double density,
    const double Wmagnitude )
{
    const double xCoordinate = positionVector[ xPositionIndex ];
    const double yCoordinate = positionVector[ yPositionIndex ];
    const double zCoordinate = positionVector[ zPositionIndex ];

    // Calculate the lower limit for the elliptical integrals by solving the cubic polynomial
    // in Lambda. See [3] for more details. All terms with "star" in the name are normalized terms.

    // Compute the normalized terms.
    const double betaStar   = beta / alpha;
    const double gammaStar  = gamma / alpha;
    const double xStar      = xCoordinate / alpha;
    const double yStar      = yCoordinate / alpha;
    const double zStar      = zCoordinate / alpha;

    // Compute delta, an intermediate term used for potential calculation
    const double alphaCube = alpha * alpha * alpha;
    const double mass = ( 4.0 * PI / 3.0 ) * density * alpha * beta * gamma;
    const double gravParameter = GRAVITATIONAL_CONSTANT * mass;
    const double delta = gravParameter / ( Wmagnitude * Wmagnitude * alphaCube );

    const double betaStarSquare     = betaStar * betaStar;
    const double gammaStarSquare    = gammaStar * gammaStar;
    const double xStarSquare        = xStar * xStar;
    const double yStarSquare        = yStar * yStar;
    const double zStarSquare        = zStar * zStar;

    // Get the coefficients for the polynomila to solve for Lambda (lower integration limit)
    const double coefficientA2 = 1.0 + betaStarSquare + gammaStarSquare
                                    - xStarSquare - yStarSquare - zStarSquare;

    const double coefficientA1 = betaStarSquare * gammaStarSquare
                                    + ( betaStarSquare + gammaStarSquare )
                                    - ( betaStarSquare + gammaStarSquare ) * xStarSquare
                                    - ( 1.0 + gammaStarSquare ) * yStarSquare
                                    - ( 1.0 + betaStarSquare ) * zStarSquare;

    const double coefficientA0 = ( betaStarSquare * gammaStarSquare )
                                    - xStarSquare * betaStarSquare * gammaStarSquare
                                    - gammaStarSquare * yStarSquare
                                    - betaStarSquare * zStarSquare;

    // Perform a check on the sign of phi(r,0) to determine the value of lambda [3].
    int phi = ( xStarSquare / 1.0 )
                 + ( yStarSquare / betaStarSquare )
                 + ( zStarSquare / gammaStarSquare ) - 1.0;
    double lambda;
    if( phi <= 0 )
    {
        // The case where the point is on the surface or inside the ellipsoid.
        lambda = 0;
    }
    else
    {
        // The case where the point is outside the ellipsoid.
        lambda = computeMaxRealCubicRoot( coefficientA2,
                                          coefficientA1,
                                          coefficientA0 );
    }

    // Note: For the sake of clarity between the code and the equations in the reference documents,
    // all vectors will have their 0th element as zero, so that the usable indices start from 1
    // until 4.

    // Calculate the first three parts of the elliptic integral in the normalized potential
    // equation.

    // Define the coefficients "a" and "b" for the elliptical integral
    std::vector< double > a = { 0.0, 1.0, betaStarSquare, gammaStarSquare, 1.0 };
    std::vector< double > b = { 0.0, 0.0, 1.0, 1.0, 1.0 };

    // evaluate the intermediate variables "Y"
    std::vector< double > Y( 5 );
    for( int i = 1; i <= 4; i++ )
    {
        Y[ i ] = std::sqrt( a[ i ] + ( b[ i ] * lambda ) );
    }

    // evaluate the intermediate variables U12, U13, U14
    double U12 = std::sqrt( b[ 1 ] * b[ 2 ] ) * Y[ 3 ] * Y[ 4 ]
                            + Y[ 1 ] * Y[ 2 ] * std::sqrt( b[ 3 ] * b[ 4 ] );

    double U13 = std::sqrt( b[ 1 ] * b[ 3 ] ) * Y[ 2 ] * Y[ 4 ]
                            + Y[ 1 ] * Y[ 3 ] * std::sqrt( b[ 2 ] * b[ 4 ] );

    double U14 = std::sqrt( b[ 1 ] * b[ 4 ] ) * Y[ 3 ] * Y[ 2 ]
                            + Y[ 1 ] * Y[ 4 ] * std::sqrt( b[ 3 ] * b[ 2 ] );

    // evaluate intermediate variables d12, d13
    double d12 = a[ 1 ] * b[ 2 ] - a[ 2 ] * b[ 1 ];
    double d13 = a[ 1 ] * b[ 3 ] - a[ 3 ] * b[ 1 ];

    // Set default error handler of GSL off to manually handle GSL errors.
    gsl_set_error_handler_off( );

    // Evaluate the elliptical integral.
    gsl_sf_result result;
    int status = gsl_sf_ellint_RD_e( U12 * U12,
                                     U13 * U13,
                                     U14 * U14,
                                     GSL_PREC_DOUBLE,
                                     &result );
    if( status != GSL_SUCCESS )
    {
        std::string errorMessage;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "!!!!!!!ERROR!!!!!!!:" << std::endl;
        std::cout << "elliptical integral (RD) evaluation failed in normalizedPotentialURE.hpp";
        std::cout << std::endl;
        errorMessage = gsl_strerror( status );
        throw std::runtime_error( errorMessage );
    }

    // evaluate the first three parts of the elliptic integral in the potential expression
    const double V1 = ( 3.0 / 4.0 ) * xStarSquare * ( 2.0 / 3.0 ) * d12 * d13 * result.val;
    const double V2 = ( 3.0 / 4.0 ) * yStarSquare * ( 2.0 / 3.0 ) * d12 * d13 * result.val;
    const double V3 = ( 3.0 / 4.0 ) * zStarSquare * ( 2.0 / 3.0 ) * d12 * d13 * result.val;

    // calculate the last part of the elliptic integral

    // Define the coefficients "a" and "b" for the elliptical integral
    a  = { 0.0, 1.0, 1.0, betaStarSquare, gammaStarSquare };
    b  = { 0.0, 0.0, 1.0, 1.0, 1.0 };

    // evaluate the intermediate variables "Y"
    for( int i = 1; i <= 4; i++ )
    {
        Y[ i ] = std::sqrt( a[ i ] + ( b[ i ] * lambda ) );
    }

    // evaluate the intermediate variables U12, U13, U14
    U12 = std::sqrt( b[ 1 ] * b[ 2 ] ) * Y[ 3 ] * Y[ 4 ]
                            + Y[ 1 ] * Y[ 2 ] * std::sqrt( b[ 3 ] * b[ 4 ] );

    U13 = std::sqrt( b[ 1 ] * b[ 3 ] ) * Y[ 2 ] * Y[ 4 ]
                            + Y[ 1 ] * Y[ 3 ] * std::sqrt( b[ 2 ] * b[ 4 ] );

    U14 = std::sqrt( b[ 1 ] * b[ 4 ] ) * Y[ 3 ] * Y[ 2 ]
                            + Y[ 1 ] * Y[ 4 ] * std::sqrt( b[ 3 ] * b[ 2 ] );

    // Set default error handler of GSL off to manually handle GSL errors.
    gsl_set_error_handler_off( );

    // Evaluate the elliptical integral.
    status = gsl_sf_ellint_RF_e( U12 * U12,
                                 U13 * U13,
                                 U14 * U14,
                                 GSL_PREC_DOUBLE,
                                 &result );
    if( status != GSL_SUCCESS )
    {
        std::string errorMessage;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "!!!!!!!ERROR!!!!!!!" << std::endl;
        std::cout << "elliptical integral (RF) evaluation failed in normalizedPotentialURE.hpp";
        std::cout << std::endl;
        errorMessage = gsl_strerror( status );
        throw std::runtime_error( errorMessage );
    }

    // evaluate the first three parts of the elliptic integral in the potential expression
    const double V4 = ( -3.0 / 4.0 ) * ( 2.0 ) * result.val;

    // get the entire elliptical integral expression now
    const double V = V1 + V2 + V3 + V4;

    // get the value for the normalized potential
    const double U = 0.5 * ( xStarSquare + yStarSquare ) - delta * V;

    // return the result
    return U;
}

} // namespace naos

#endif // NORMALIZED_POTENTIAL_URE_HPP

//! References
/*!
 * [1] Orbital motion around strongly perturbed environments. Author: Daniel J. Scheeres.
 * [2] A Table of elliptic integrals of the second kind. Author: B.C. Carlson.
 * [3] Dynamics about uniformly rotating tri-axial ellipsoids. Author: Daniel J. Scheeres.
 */
