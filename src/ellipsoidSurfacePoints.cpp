/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <cstdlib>
#include <fstream>
#include <cmath>

#include "NAOS/ellipsoidSurfacePoints.hpp"

namespace naos
{

//! Compute tri-axial ellipsoid surface points.
/*!
 * Computes the cartesian coordinates for surface points of a tri-axial ellipsoid by performing
 * a grid search over all possible coordinates and checking if they satisfy the ellipsoid equation.
 * If they satisfy the equation, then the points are stored in a CSV file.
 *
 * @param[in]
 */
const void computeEllipsoidSurfacePoints( const double alpha,
                                          const double beta,
                                          const double gamma,
                                          const double stepSize,
                                          std::ofstream &fileHandle )
{
    double latitude = 0.0;
    double longitude = 0.0;
    for( latitude = 0.0; latitude <=  )
}


} // namespace naos
