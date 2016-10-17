/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef MISC_HPP
#define MISC_HPP

#include <iostream>

namespace naos
{

//! Print the elements of a vector on console
/*!
 * Print the elements of a vector on the console.
 *
 * @tparam Vector       Template paramter for a 3 element vector
 * @param[in] vector    A vector of any size
 * @param[in] size      size of the vector
 */
template< typename Vector >
void printVector( const Vector &vector, const double size )
{
    std::cout << std::endl;
    for( int i = 0; i < size; i++ )
    {
        std::cout << vector[ i ] << std::endl;
    }
    std::cout << std::endl;
}

//! Progress bar
/*!
 * function that will display the progress bar when executed within a certain c++ function.
 *
 */
template< typename Real >
void displayProgressBar( Real counterValue, Real endValue )
{
    int barWidth = 70;

    std::cout << "[";
    int position = barWidth * ( counterValue / endValue );
    for ( int barIndex = 0; barIndex < barWidth; ++barIndex )
    {
        if (barIndex < position)
            std::cout << "=";
        else if (barIndex == position)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int( ( counterValue / endValue ) * 100.0 ) << "%\r";
    std::cout.flush( );
}

} // namespace naos

#endif // MISC_HPP
