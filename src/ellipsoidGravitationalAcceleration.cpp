/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_errno.h>

#include "NAOS/cubicRoot.hpp"
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
                                                       Vector3 &gravAcceleration )
 {
    // Calculate the lower limit for the elliptical integrals by solving the cubic polynomial
    // in Lambda. See [1] for more details.
    const double alphaSquare = alpha * alpha;
    const double betaSquare = beta * beta;
    const double gammaSquare = gamma * gamma;
    const double xCoordinateSquare = xCoordinate * xCoordinate;
    const double yCoordinateSquare = yCoordinate * yCoordinate;
    const double zCoordinateSquare = zCoordinate * zCoordinate;

    const double coefficientA2 = alphaSquare + betaSquare + gammaSquare
                                    - xCoordinateSquare - yCoordinateSquare - zCoordinateSquare;

    const double coefficientA1 = ( alphaSquare * betaSquare ) + ( alphaSquare * gammaSquare )
                                    + ( betaSquare * gammaSquare )
                                    - ( xCoordinateSquare * ( betaSquare + gammaSquare ) )
                                    - ( yCoordinateSquare * ( alphaSquare + gammaSquare ) )
                                    - ( zCoordinateSquare * ( alphaSquare + betaSquare ) );

    const double coefficientA0 = ( alphaSquare * betaSquare * gammaSquare )
                                    - ( xCoordinateSquare * gammaSquare * betaSquare )
                                    - ( yCoordinateSquare * alphaSquare * gammaSquare )
                                    - ( zCoordinateSquare * alphaSquare * betaSquare );

    // Perform a check on the sign of phi(r,0) to determine the value of lambda [3].
    int phi = ( xCoordinateSquare / alphaSquare )
                 + ( yCoordinateSquare / betaSquare )
                 + ( zCoordinateSquare / gammaSquare ) - 1.0;
    double lambda;
    if( phi <= 0 )
    {
        // The case where the point is on the surface or inside the ellipsoid.
        lambda = 0;
    }
    else
    {
        // The case where the point is outside the ellipsoid.
        lambda = naos::computeMaxRealCubicRoot( coefficientA2,
                                                coefficientA1,
                                                coefficientA0 );
    }

    // Note: For the sake of clarity between the code and the equations in the reference documents,
    // all vectors will have their 0th element as zero, so that the usable indices start from 1
    // until 4.

    // Define intermediate coefficients "a_i" and "b_i" for Ux, Uy, Uz (the partial derivatives of
    // gravitational potential). See section 2 of [2] for usage of these coefficients.
    const std::vector< double > a_Ux { 0.0, 1.0, betaSquare, gammaSquare, alphaSquare };
    const std::vector< double > a_Uy { 0.0, 1.0, alphaSquare, gammaSquare, betaSquare };
    const std::vector< double > a_Uz { 0.0, 1.0, alphaSquare, betaSquare, gammaSquare };

    const std::vector< double > b_Ux { 0.0, 0.0, 1.0, 1.0, 1.0 };
    const std::vector< double > b_Uy { 0.0, 0.0, 1.0, 1.0, 1.0 };
    const std::vector< double > b_Uz { 0.0, 0.0, 1.0, 1.0, 1.0 };

    // Evaluate the intermediate variable "Y_i" for Ux, Uy, Uz.
    std::vector< double > Y_Ux { 0.0, 0.0, 0.0, 0.0, 0.0 };
    std::vector< double > Y_Uy { 0.0, 0.0, 0.0, 0.0, 0.0 };
    std::vector< double > Y_Uz { 0.0, 0.0, 0.0, 0.0, 0.0 };

    for( int i = 1; i <= 4; i++ )
    {
        Y_Ux[ i ] = std::sqrt( a_Ux[ i ] + ( b_Ux[ i ] * lambda ) );
        Y_Uy[ i ] = std::sqrt( a_Uy[ i ] + ( b_Uy[ i ] * lambda ) );
        Y_Uz[ i ] = std::sqrt( a_Uz[ i ] + ( b_Uz[ i ] * lambda ) );
    }

    // Evaluate the intermediate variable "U_ij" for Ux, Uy, Uz.
    const double U12_Ux = std::sqrt( b_Ux[ 1 ] * b_Ux[ 2 ] ) * Y_Ux[ 3 ] * Y_Ux[ 4 ]
                            + Y_Ux[ 1 ] * Y_Ux[ 2 ] * std::sqrt( b_Ux[ 3 ] * b_Ux[ 4 ] );

    const double U13_Ux = std::sqrt( b_Ux[ 1 ] * b_Ux[ 3 ] ) * Y_Ux[ 2 ] * Y_Ux[ 4 ]
                            + Y_Ux[ 1 ] * Y_Ux[ 3 ] * std::sqrt( b_Ux[ 2 ] * b_Ux[ 4 ] );


    const double U14_Ux = std::sqrt( b_Ux[ 1 ] * b_Ux[ 4 ] ) * Y_Ux[ 3 ] * Y_Ux[ 2 ]
                            + Y_Ux[ 1 ] * Y_Ux[ 4 ] * std::sqrt( b_Ux[ 3 ] * b_Ux[ 2 ] );

    const double U12_Uy = std::sqrt( b_Uy[ 1 ] * b_Uy[ 2 ] ) * Y_Uy[ 3 ] * Y_Uy[ 4 ]
                            + Y_Uy[ 1 ] * Y_Uy[ 2 ] * std::sqrt( b_Uy[ 3 ] * b_Uy[ 4 ] );

    const double U13_Uy = std::sqrt( b_Uy[ 1 ] * b_Uy[ 3 ] ) * Y_Uy[ 2 ] * Y_Uy[ 4 ]
                            + Y_Uy[ 1 ] * Y_Uy[ 3 ] * std::sqrt( b_Uy[ 2 ] * b_Uy[ 4 ] );

    const double U14_Uy = std::sqrt( b_Uy[ 1 ] * b_Uy[ 4 ] ) * Y_Uy[ 3 ] * Y_Uy[ 2 ]
                            + Y_Uy[ 1 ] * Y_Uy[ 4 ] * std::sqrt( b_Uy[ 3 ] * b_Uy[ 2 ] );

    const double U12_Uz = std::sqrt( b_Uz[ 1 ] * b_Uz[ 2 ] ) * Y_Uz[ 3 ] * Y_Uz[ 4 ]
                            + Y_Uz[ 1 ] * Y_Uz[ 2 ] * std::sqrt( b_Uz[ 3 ] * b_Uz[ 4 ] );

    const double U13_Uz = std::sqrt( b_Uz[ 1 ] * b_Uz[ 3 ] ) * Y_Uz[ 2 ] * Y_Uz[ 4 ]
                            + Y_Uz[ 1 ] * Y_Uz[ 3 ] * std::sqrt( b_Uz[ 2 ] * b_Uz[ 4 ] );

    const double U14_Uz = std::sqrt( b_Uz[ 1 ] * b_Uz[ 4 ] ) * Y_Uz[ 3 ] * Y_Uz[ 2 ]
                            + Y_Uz[ 1 ] * Y_Uz[ 4 ] * std::sqrt( b_Uz[ 3 ] * b_Uz[ 2 ] );

    // Evaluate intermediate variables "d_ij" for Ux, Uy, Uz.
    const double d12_Ux = a_Ux[ 1 ] * b_Ux[ 2 ] - a_Ux[ 2 ] * b_Ux[ 1 ];
    const double d13_Ux = a_Ux[ 1 ] * b_Ux[ 3 ] - a_Ux[ 3 ] * b_Ux[ 1 ];

    const double d12_Uy = a_Uy[ 1 ] * b_Uy[ 2 ] - a_Uy[ 2 ] * b_Uy[ 1 ];
    const double d13_Uy = a_Uy[ 1 ] * b_Uy[ 3 ] - a_Uy[ 3 ] * b_Uy[ 1 ];

    const double d12_Uz = a_Uz[ 1 ] * b_Uz[ 2 ] - a_Uz[ 2 ] * b_Uz[ 1 ];
    const double d13_Uz = a_Uz[ 1 ] * b_Uz[ 3 ] - a_Uz[ 3 ] * b_Uz[ 1 ];

    // Set default error handler of GSL off to manually handle GSL errors.
    gsl_set_error_handler_off( );

    // Evaluate Ux.
    gsl_sf_result result_Ux;
    int status_Ux = gsl_sf_ellint_RD_e( U12_Ux * U12_Ux,
                                        U13_Ux * U13_Ux,
                                        U14_Ux * U14_Ux,
                                        GSL_PREC_DOUBLE,
                                        &result_Ux );
    // Checking the status of evaluation, a non-zero status is an error
    if( status_Ux != GSL_SUCCESS )
    {
        std::cout << std::endl;
        std::cout << "Error in ellipsoid gravitational acceleration computation!" << std::endl;
        std::cout << "Parameteric values that resulted in the error:" << std::endl;
        std::cout << "x Coordinate = " << xCoordinate << std::endl;
        std::cout << "y Coordinate = " << yCoordinate << std::endl;
        std::cout << "z Coordinate = " << zCoordinate << std::endl;
        std::cout << "a_Ux = ";
        for( int i = 1; i <=4; i++ )
        {
            std::cout << a_Ux[ i ] << ",";
        }
        std::cout << std::endl;
        std::cout << "b_Ux = ";
        for( int i = 1; i <=4; i++ )
        {
            std::cout << b_Ux[ i ] << ",";
        }
        std::cout << std::endl;
        std::cout << "Y_Ux = ";
        for( int i = 1; i <=4; i++ )
        {
            std::cout << Y_Ux[ i ] << ",";
        }
        std::cout << std::endl;
        std::cout << "U12_Ux = " << U12_Ux << std::endl;
        std::cout << "U13_Ux = " << U13_Ux << std::endl;
        std::cout << "U14_Ux = " << U14_Ux << std::endl;
        std::cout << "lambda = " << lambda << std::endl;
        std::cout << "Phi =  " << phi << std::endl;

        std::string errorMessage;
        errorMessage = gsl_strerror( status_Ux );
        throw std::runtime_error( errorMessage );
    }
    const double Ux = ( -3.0 * gravParameter * xCoordinate / 2.0 )
                        * ( 2.0 / 3.0 ) * d12_Ux * d13_Ux * result_Ux.val;

    // Evaluate Uy.
    const double Uy = ( -3.0 * gravParameter * yCoordinate / 2.0 )
                        * ( 2.0 / 3.0 ) * d12_Uy * d13_Uy
                        * gsl_sf_ellint_RD( U12_Uy * U12_Uy,
                                            U13_Uy * U13_Uy,
                                            U14_Uy * U14_Uy,
                                            GSL_PREC_DOUBLE );

    // Evaluate Uz.
    const double Uz = ( -3.0 * gravParameter * zCoordinate / 2.0 )
                        * ( 2.0 / 3.0 ) * d12_Uz * d13_Uz
                        * gsl_sf_ellint_RD( U12_Uz * U12_Uz,
                                            U13_Uz * U13_Uz,
                                            U14_Uz * U14_Uz,
                                            GSL_PREC_DOUBLE );

    // Store the final gravitational acceleration values.
    gravAcceleration[ 0 ] = Ux;
    gravAcceleration[ 1 ] = Uy;
    gravAcceleration[ 2 ] = Uz;

 }

//! Compute nondimensionalized ellipsoid gravitational Acceleration
/*!
 * Compute the gravitational acceleration i.e. partial derivative of the external gravitational
 * potential of a tri-axial ellipsoid. The (x,y,z) coordinates for this should be provided in the
 * body fixed principal axis frame [1]. The accelerations are in the body frame. All distance units
 * will be normalized by the largest semi-axes of the ellipsoid.
 *
 * @param[in] alpha                 Largest semi-axis of ellipsoid
 * @param[in] beta                  Intermediate semi-axis of ellipsoid
 * @param[in] gamma                 Smallest semi-axis of ellipsoid
 * @param[in] density               density of the ellipsoid in SI units
 * @param[in] w                     angular velocity magnitude of the ellipsoid
 * @param[in] xCoordinate           X coordinate where the gravitational acceleration has to be
 *                                  calculated
 * @param[in] yCoordinate           Y coordinate where the gravitational acceleration has to be
 *                                  calculated
 * @param[in] zCoordinate           Z coordinate where the gravitational acceleration has to be
 *                                  calculated
 * @param[in/out] gravAcceleration  Output vector containing the gravitational acceleration
 *                                  components
 */
 const void computeNonDimensionalEllipsoidGravitationalAcceleration( const double alpha,
                                                                     const double beta,
                                                                     const double gamma,
                                                                     const double density,
                                                                     const double w,
                                                                     const double xCoordinate,
                                                                     const double yCoordinate,
                                                                     const double zCoordinate,
                                                                     Vector3 &gravAcceleration )
 {
    // Calculate the lower limit for the elliptical integrals by solving the cubic polynomial
    // in Lambda. See [3] for more details. All terms with "star" in the name are normalized terms.
    // Compute the normalized terms.
    const double betaStar   = beta / alpha;
    const double gammaStar  = gamma / alpha;
    const double xStar      = xCoordinate / alpha;
    const double yStar      = yCoordinate / alpha;
    const double zStar      = zCoordinate / alpha;

    // cube of alpha, intermediate term used in some places
    const double alphaCube = alpha * alpha * alpha;

    const double mStar             = ( 4.0 * naos::PI / 3.0 ) * density * 1.0
                                        * betaStar * gammaStar * alphaCube;
    const double gravParameterStar = naos::GRAVITATIONAL_CONSTANT * mStar;

    const double betaStarSquare     = betaStar * betaStar;
    const double gammaStarSquare    = gammaStar * gammaStar;
    const double xStarSquare        = xStar * xStar;
    const double yStarSquare        = yStar * yStar;
    const double zStarSquare        = zStar * zStar;

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
        lambda = naos::computeMaxRealCubicRoot( coefficientA2,
                                                coefficientA1,
                                                coefficientA0 );
    }

    // Note: For the sake of clarity between the code and the equations in the reference documents,
    // all vectors will have their 0th element as zero, so that the usable indices start from 1
    // until 4.

    // Define intermediate coefficients "a_i" and "b_i" for Ux, Uy, Uz (the partial derivatives of
    // gravitational potential). See section 2 of [2] for usage of these coefficients.
    const std::vector< double > a_Ux { 0.0, 1.0, betaStarSquare, gammaStarSquare, 1.0 };
    const std::vector< double > a_Uy { 0.0, 1.0, 1.0, gammaStarSquare, betaStarSquare };
    const std::vector< double > a_Uz { 0.0, 1.0, 1.0, betaStarSquare, gammaStarSquare };

    const std::vector< double > b_Ux { 0.0, 0.0, 1.0, 1.0, 1.0 };
    const std::vector< double > b_Uy { 0.0, 0.0, 1.0, 1.0, 1.0 };
    const std::vector< double > b_Uz { 0.0, 0.0, 1.0, 1.0, 1.0 };

    // Evaluate the intermediate variable "Y_i" for Ux, Uy, Uz.
    // Y_i declaration only:
    std::vector< double > Y_Ux { 0.0, 0.0, 0.0, 0.0, 0.0 };
    std::vector< double > Y_Uy { 0.0, 0.0, 0.0, 0.0, 0.0 };
    std::vector< double > Y_Uz { 0.0, 0.0, 0.0, 0.0, 0.0 };

    for( int i = 1; i <= 4; i++ )
    {
        Y_Ux[ i ] = std::sqrt( a_Ux[ i ] + ( b_Ux[ i ] * lambda ) );
        Y_Uy[ i ] = std::sqrt( a_Uy[ i ] + ( b_Uy[ i ] * lambda ) );
        Y_Uz[ i ] = std::sqrt( a_Uz[ i ] + ( b_Uz[ i ] * lambda ) );
    }

    // Evaluate the intermediate variable "U_ij" for Ux, Uy, Uz.
    const double U12_Ux = std::sqrt( b_Ux[ 1 ] * b_Ux[ 2 ] ) * Y_Ux[ 3 ] * Y_Ux[ 4 ]
                            + Y_Ux[ 1 ] * Y_Ux[ 2 ] * std::sqrt( b_Ux[ 3 ] * b_Ux[ 4 ] );

    const double U13_Ux = std::sqrt( b_Ux[ 1 ] * b_Ux[ 3 ] ) * Y_Ux[ 2 ] * Y_Ux[ 4 ]
                            + Y_Ux[ 1 ] * Y_Ux[ 3 ] * std::sqrt( b_Ux[ 2 ] * b_Ux[ 4 ] );


    const double U14_Ux = std::sqrt( b_Ux[ 1 ] * b_Ux[ 4 ] ) * Y_Ux[ 3 ] * Y_Ux[ 2 ]
                            + Y_Ux[ 1 ] * Y_Ux[ 4 ] * std::sqrt( b_Ux[ 3 ] * b_Ux[ 2 ] );

    const double U12_Uy = std::sqrt( b_Uy[ 1 ] * b_Uy[ 2 ] ) * Y_Uy[ 3 ] * Y_Uy[ 4 ]
                            + Y_Uy[ 1 ] * Y_Uy[ 2 ] * std::sqrt( b_Uy[ 3 ] * b_Uy[ 4 ] );

    const double U13_Uy = std::sqrt( b_Uy[ 1 ] * b_Uy[ 3 ] ) * Y_Uy[ 2 ] * Y_Uy[ 4 ]
                            + Y_Uy[ 1 ] * Y_Uy[ 3 ] * std::sqrt( b_Uy[ 2 ] * b_Uy[ 4 ] );

    const double U14_Uy = std::sqrt( b_Uy[ 1 ] * b_Uy[ 4 ] ) * Y_Uy[ 3 ] * Y_Uy[ 2 ]
                            + Y_Uy[ 1 ] * Y_Uy[ 4 ] * std::sqrt( b_Uy[ 3 ] * b_Uy[ 2 ] );

    const double U12_Uz = std::sqrt( b_Uz[ 1 ] * b_Uz[ 2 ] ) * Y_Uz[ 3 ] * Y_Uz[ 4 ]
                            + Y_Uz[ 1 ] * Y_Uz[ 2 ] * std::sqrt( b_Uz[ 3 ] * b_Uz[ 4 ] );

    const double U13_Uz = std::sqrt( b_Uz[ 1 ] * b_Uz[ 3 ] ) * Y_Uz[ 2 ] * Y_Uz[ 4 ]
                            + Y_Uz[ 1 ] * Y_Uz[ 3 ] * std::sqrt( b_Uz[ 2 ] * b_Uz[ 4 ] );

    const double U14_Uz = std::sqrt( b_Uz[ 1 ] * b_Uz[ 4 ] ) * Y_Uz[ 3 ] * Y_Uz[ 2 ]
                            + Y_Uz[ 1 ] * Y_Uz[ 4 ] * std::sqrt( b_Uz[ 3 ] * b_Uz[ 2 ] );

    // Evaluate intermediate variables "d_ij" for Ux, Uy, Uz.
    const double d12_Ux = a_Ux[ 1 ] * b_Ux[ 2 ] - a_Ux[ 2 ] * b_Ux[ 1 ];
    const double d13_Ux = a_Ux[ 1 ] * b_Ux[ 3 ] - a_Ux[ 3 ] * b_Ux[ 1 ];

    const double d12_Uy = a_Uy[ 1 ] * b_Uy[ 2 ] - a_Uy[ 2 ] * b_Uy[ 1 ];
    const double d13_Uy = a_Uy[ 1 ] * b_Uy[ 3 ] - a_Uy[ 3 ] * b_Uy[ 1 ];

    const double d12_Uz = a_Uz[ 1 ] * b_Uz[ 2 ] - a_Uz[ 2 ] * b_Uz[ 1 ];
    const double d13_Uz = a_Uz[ 1 ] * b_Uz[ 3 ] - a_Uz[ 3 ] * b_Uz[ 1 ];

    // Set default error handler of GSL off to manually handle GSL errors.
    gsl_set_error_handler_off( );

    // Evaluate Ux.
    gsl_sf_result result_Ux;
    int status_Ux = gsl_sf_ellint_RD_e( U12_Ux * U12_Ux,
                                        U13_Ux * U13_Ux,
                                        U14_Ux * U14_Ux,
                                        GSL_PREC_DOUBLE,
                                        &result_Ux );
    // Checking the status of evaluation, a non-zero status is an error
    if( status_Ux != GSL_SUCCESS )
    {
        std::cout << std::endl;
        std::cout << "Error in ellipsoid gravitational acceleration computation!" << std::endl;
        std::cout << "Parameteric values that resulted in the error:" << std::endl;
        std::cout << "x Coordinate normalized = " << xStar << std::endl;
        std::cout << "y Coordinate normalized = " << yStar << std::endl;
        std::cout << "z Coordinate normalized = " << zStar << std::endl;
        std::cout << "a_Ux = ";
        for( int i = 1; i <=4; i++ )
        {
            std::cout << a_Ux[ i ] << ",";
        }
        std::cout << std::endl;
        std::cout << "b_Ux = ";
        for( int i = 1; i <=4; i++ )
        {
            std::cout << b_Ux[ i ] << ",";
        }
        std::cout << std::endl;
        std::cout << "Y_Ux = ";
        for( int i = 1; i <=4; i++ )
        {
            std::cout << Y_Ux[ i ] << ",";
        }
        std::cout << std::endl;
        std::cout << "U12_Ux = " << U12_Ux << std::endl;
        std::cout << "U13_Ux = " << U13_Ux << std::endl;
        std::cout << "U14_Ux = " << U14_Ux << std::endl;
        std::cout << "lambda = " << lambda << std::endl;
        std::cout << "Phi =  " << phi << std::endl;

        std::string errorMessage;
        errorMessage = gsl_strerror( status_Ux );
        throw std::runtime_error( errorMessage );
    }

    // Evaluate Ux normalized
    const double UxStar = xStar - ( 3.0 * gravParameterStar * xStar / ( 2.0 * alphaCube * w * w ) )
                            * ( 2.0 / 3.0 ) * d12_Ux * d13_Ux * result_Ux.val;

    // Evaluate Uy normalized
    const double UyStar = yStar - ( 3.0 * gravParameterStar * yStar / ( 2.0 * alphaCube * w * w ) )
                            * ( 2.0 / 3.0 ) * d12_Uy * d13_Uy
                            * gsl_sf_ellint_RD( U12_Uy * U12_Uy,
                                                U13_Uy * U13_Uy,
                                                U14_Uy * U14_Uy,
                                                GSL_PREC_DOUBLE );

    // Evaluate Uz normalized
    const double UzStar = ( -3.0 * gravParameterStar * zStar / ( 2.0 * alphaCube * w * w ) )
                            * ( 2.0 / 3.0 ) * d12_Uz * d13_Uz
                            * gsl_sf_ellint_RD( U12_Uz * U12_Uz,
                                                U13_Uz * U13_Uz,
                                                U14_Uz * U14_Uz,
                                                GSL_PREC_DOUBLE );

    // Store the final gravitational acceleration values.
    gravAcceleration[ 0 ] = UxStar;
    gravAcceleration[ 1 ] = UyStar;
    gravAcceleration[ 2 ] = UzStar;
 }

} // namespace naos

//! References
/*!
 * [1] Orbital motion around strongly perturbed environments. Author: Daniel J. Scheeres.
 * [2] A Table of elliptic integrals of the second kind. Author: B.C. Carlson.
 * [3] Dynamics about uniformly rotating tri-axial ellipsoids. Author: Daniel J. Scheeres.
 */
