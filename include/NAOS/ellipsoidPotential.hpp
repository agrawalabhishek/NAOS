/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef ELLIPSOID_POTENTIAL_HPP
#define ELLIPSOID_POTENTIAL_HPP

#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iostream>

#include <boost/math/special_functions/ellint_rd.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>

#include "NAOS/cubicRoot.hpp"
#include "NAOS/constants.hpp"

namespace naos
{

//! Compute gravitational potential of an ellipsoid
/*!
 * Compute the gravitational potential of an ellipsoid. The x,y,z coordinate are provided in the
 * body fixed fame. The potential values are defined in the body frame.
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
 * @param[in/out] gravPotential     Output which is the gravitational potential value
 */
template< typename Real >
void computeEllipsoidGravitationalPotential(
    const Real alpha,
    const Real beta,
    const Real gamma,
    const Real gravParameter,
    const Real xCoordinate,
    const Real yCoordinate,
    const Real zCoordinate,
    Real &gravPotential )
{
    // Calculate the lower limit for the elliptical integrals by solving the cubic polynomial
    // in Lambda. See [1] for more details.
    const Real alphaSquare = alpha * alpha;
    const Real betaSquare = beta * beta;
    const Real gammaSquare = gamma * gamma;
    const Real xCoordinateSquare = xCoordinate * xCoordinate;
    const Real yCoordinateSquare = yCoordinate * yCoordinate;
    const Real zCoordinateSquare = zCoordinate * zCoordinate;

    const Real coefficientA2 = alphaSquare + betaSquare + gammaSquare
                                    - xCoordinateSquare - yCoordinateSquare - zCoordinateSquare;

    const Real coefficientA1 = ( alphaSquare * betaSquare ) + ( alphaSquare * gammaSquare )
                                    + ( betaSquare * gammaSquare )
                                    - ( xCoordinateSquare * ( betaSquare + gammaSquare ) )
                                    - ( yCoordinateSquare * ( alphaSquare + gammaSquare ) )
                                    - ( zCoordinateSquare * ( alphaSquare + betaSquare ) );

    const Real coefficientA0 = ( alphaSquare * betaSquare * gammaSquare )
                                    - ( xCoordinateSquare * gammaSquare * betaSquare )
                                    - ( yCoordinateSquare * alphaSquare * gammaSquare )
                                    - ( zCoordinateSquare * alphaSquare * betaSquare );

    // Perform a check on the sign of phi(r,0) to determine the value of lambda [3].
    double phi = ( xCoordinateSquare / alphaSquare )
                 + ( yCoordinateSquare / betaSquare )
                 + ( zCoordinateSquare / gammaSquare );
    Real lambda;
    if( ( phi - 1.0 ) <= 1.0e-10 || ( phi - 1.0 ) == 0.0 )
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

    // Evaluate part 1 of the potential equation, V1
    Real result_V1 = boost::math::ellint_rd( ( betaSquare + lambda ),
                                             ( gammaSquare + lambda ),
                                             ( alphaSquare + lambda ) );

    const Real V1 = -1.0 * gravParameter * xCoordinateSquare * result_V1 / 2.0;

    // Evaluate part 2 of the potential equation, V2
    Real result_V2 = boost::math::ellint_rd( ( gammaSquare + lambda ),
                                             ( alphaSquare + lambda ),
                                             ( betaSquare + lambda ) );

    const Real V2 = -1.0 * gravParameter * yCoordinateSquare * result_V2 / 2.0;

    // Evaluate part 3 of the potential equation, V3
    Real result_V3 = boost::math::ellint_rd( ( alphaSquare + lambda ),
                                             ( betaSquare + lambda ),
                                             ( gammaSquare + lambda ) );

    const Real V3 = -1.0 * gravParameter * zCoordinateSquare * result_V3 / 2.0;

    // Evaluate part 4 of the potential equation, V4
    Real result_V4 = boost::math::ellint_rf( ( alphaSquare + lambda ),
                                             ( betaSquare + lambda ),
                                             ( gammaSquare + lambda ) );

    const Real V4 = 3.0 * gravParameter * result_V4 / 2.0;

    // final potential value
    gravPotential = V1 + V2 + V3 + V4;
}

} // namespace naos

#endif // ELLIPSOID_POTENTIAL_HPP
