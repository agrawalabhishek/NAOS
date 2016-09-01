/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef NAOS_CUBIC_ROOT_HPP
#define NAOS_CUBIC_ROOT_HPP

namespace naos
{
//! Compute maximum real root of a cubic polynomial
/*!
 * Compute the real roots of any general cubic polynomial. The equations are taken from [1].
 * To use this function, the equation should be reduced in the standard form such that the
 * coefficient of the cubed variable in the polynomial is unity.
 *
 * @param[in] a2    Coefficient of the squared variable in the cubic polynomial
 * @param[in] a1    Coefficient of the variable powered 1, in the cubic polynomial
 * @param[in] a0    The constant coefficient in the cubic polynomial
 * @return          The maximum real root of the cubic polynomial
 */
 const double computeMaxRealCubicRoot( const double a2,
                                       const double a1,
                                       const double a0 );
} // namespace naos

#endif // NAOS_CUBIC_ROOT_HPP

//! References
/*!
 * [1] Wolfram Alpha online math resource. URL: http://mathworld.wolfram.com/Ellipsoid.html
       Accessed: 29th August 2016
 *
 */
