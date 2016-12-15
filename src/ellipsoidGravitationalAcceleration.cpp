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

#include <boost/math/special_functions/ellint_rd.hpp>

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

    // Evaluate Ux.
    double xResult = boost::math::ellint_rd( ( betaSquare + lambda ),
                                             ( gammaSquare + lambda ),
                                             ( alphaSquare + lambda ) );

    const double Ux = -1.0 * gravParameter * xCoordinate * xResult;

    // Evaluate Uy.
    double yResult = boost::math::ellint_rd( ( gammaSquare + lambda ),
                                             ( alphaSquare + lambda ),
                                             ( betaSquare + lambda ) );

    const double Uy = -1.0 * gravParameter * yCoordinate * yResult;

    // Evaluate Uz.
    double zResult = boost::math::ellint_rd( ( alphaSquare + lambda ),
                                             ( betaSquare + lambda ),
                                             ( gammaSquare + lambda ) );

    const double Uz = -1.0 * gravParameter * zCoordinate * zResult;

    // Store the final gravitational acceleration values.
    gravAcceleration[ 0 ] = Ux;
    gravAcceleration[ 1 ] = Uy;
    gravAcceleration[ 2 ] = Uz;

 }

} // namespace naos

//! References
/*!
 * [1] Orbital motion around strongly perturbed environments. Author: Daniel J. Scheeres.
 * [2] A Table of elliptic integrals of the second kind. Author: B.C. Carlson.
 * [3] Dynamics about uniformly rotating tri-axial ellipsoids. Author: Daniel J. Scheeres.
 */
