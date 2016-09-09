/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream>

#include "NAOS/cubicRoot.hpp"
#include "NAOS/constants.hpp"

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
                                       const double a0 )
 {
    double maxRealRoot = 0.0;

    // Compute intermediate values. Same variables used as given in [1].
    const double Q = ( 3.0 * a1 - a2 * a2 ) / 9.0;
    const double R = ( 9.0 * a2 * a1 - 27.0 * a0 - ( 2.0 * a2 * a2 * a2 ) ) / 54.0;
    const double D = ( Q * Q * Q ) + ( R * R );

    if( D >= 0.0 )
    {
        // Some other intermediate values as defined in [1].
        const double S = std::cbrt( R + std::sqrt( D ) );
        const double T = std::cbrt( R - std::sqrt( D ) );

        if( D > 0.0 )
        {
            // Only one unique real root exists [1].
            maxRealRoot = ( -1.0 / 3.0 ) * a2 + ( S + T );
        }
        else if( D == 0.0 )
        {
            // All roots are real but Atmost only two unique real roots exist [1].
            double realRoot1 = ( -1.0 / 3.0 ) * a2 + ( S + T );
            double realRoot2 = ( -1.0 / 3.0 ) * a2 - ( ( 1.0 / 2.0 ) * ( S + T ) );
            // Get the maximum real root.
            maxRealRoot = std::max( realRoot1, realRoot2 );
        }
    }
    else
    {
        // Some other intermediate values [1].
        double nQ = -1.0 * Q;
        double theta = std::acos( R / std::sqrt( nQ * nQ * nQ ) );

        // All three roots are real for D < 0 as stated in [1].
        double realRoot1 = 2.0 * std::sqrt( nQ ) * std::cos( theta / 3.0 ) - ( a2 / 3.0 );
        double realRoot2 = 2.0 * std::sqrt( nQ ) * std::cos( ( theta + 2.0 * naos::PI ) / 3.0 )
                            - ( a2 / 3.0 );
        double realRoot3 = 2.0 * std::sqrt( nQ ) * std::cos( ( theta + 4.0 * naos::PI ) / 3.0 )
                            - ( a2 / 3.0 );

        // Find the maximum real root from all three.
        maxRealRoot = std::max( realRoot1, realRoot2 );
        maxRealRoot = std::max( maxRealRoot, realRoot3 );
    }

    // Check sign of maxRealRoot
    if( maxRealRoot < 0.0 )
    {
        // throw an exception
        std::ostringstream errorMessage;
        std::cout << std::endl;
        errorMessage << "ERROR: Negative Lambda value obtained!" << std::endl;
        std::cout << errorMessage.str( ) << std::endl;
        std::cout << "Coefficient A2 = " << a2 << std::endl;
        std::cout << "Coefficient A1 = " << a1 << std::endl;
        std::cout << "Coefficient A0 = " << a0 << std::endl;
        std::cout << "Max Real Root = " << maxRealRoot << std::endl;
        // throw std::runtime_error( errorMessage.str( ) );
    }

    // Return the max real root.
    return maxRealRoot;
 }

} // namespace naos

//! References
/*!
 * [1] Wolfram Alpha online math resource. URL: http://mathworld.wolfram.com/CubicFormula.html
       Accessed: 30th August 2016
 *
 */
