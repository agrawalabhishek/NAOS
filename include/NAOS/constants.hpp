/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef NAOS_CONSTANTS_HPP
#define NAOS_CONSTANTS_HPP

#include <vector>

namespace naos
{
    typedef std::vector< double > Vector6;
    typedef std::vector< double > Vector3;

    //! Gravitational constant [m^3 s^-2]
    const static double GRAVITATIONAL_CONSTANT = 6.67259e-11;

    //! Pi
    const static double PI = 3.14159265358979323846;
} // namespace naos

#endif // NAOS_CONSTANTS_HPP
