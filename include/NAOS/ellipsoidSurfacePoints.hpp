/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef NAOS_ELLIPSOID_SURFACE_POINTS_HPP
#define NAOS_ELLIPSOID_SURFACE_POINTS_HPP

namespace naos
{
//! Compute tri-axial ellipsoid surface points.
/*!
* Computes the cartesian coordinates for surface points of a tri-axial ellipsoid. The parametric
* equations are taken from [1]. The computed coordinates are stored in a CSV file.
*
* @param[in] alpha                 Largest semi-axis of the ellipsoid [km]
* @param[in] beta                  Intermediate semi-axis of the ellipsoid [km]
* @param[in] gamma                 Smallest semi-axis of the ellipsoid [km]
* @param[in] stepSizeAzimuth       Step size for succession in azimuth [deg]
* @param[in] stepSizeElevation     Step size for succession in elevation [deg]
* @param[in] outputFile            CSV file name that will store the calculated surface
*                                  coordinates. Specify file name with its location relative to
*                                  the 'bin' folder. Also specify the '.csv' extension.
*/
const void computeEllipsoidSurfacePoints( const double alpha,
                                          const double beta,
                                          const double gamma,
                                          const double stepSizeAzimuth,
                                          const double stepSizeElevation,
                                          const std::ostringstream &outputFile );
} // namespace naos

#endif // NAOS_ELLIPSOID_SURFACE_POINTS_HPP

//! References
/*!
 * [1] Wolfram Alpha online math resource. URL: http://mathworld.wolfram.com/Ellipsoid.html
       Accessed: 29th August 2016
 *
 */
