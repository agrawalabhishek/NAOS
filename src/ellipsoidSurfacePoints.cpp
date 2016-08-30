/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <stdexcept>

#include "NAOS/ellipsoidSurfacePoints.hpp"
#include "NAOS/constants.hpp"

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
 * @param[in] stepSizeAzimuth       Step size for succession in Azimuth [deg]
 * @param[in] stepSizeElevation     Step size for succession in Elevation [deg]
 * @param[in] outputFile            CSV file name that will store the calculated surface
 *                                  coordinates. Specify file name with its location relative to
 *                                  the 'bin' folder. Also specify the '.csv' extension.
 */
const void computeEllipsoidSurfacePoints( const double alpha,
                                          const double beta,
                                          const double gamma,
                                          const double stepSizeAzimuth,
                                          const double stepSizeElevation,
                                          const std::ostringstream &outputFile )
{
    std::ofstream fileHandle;
    fileHandle.precision( 15 );
    fileHandle.open( outputFile.str( ) );
    if( !fileHandle.is_open( ) )
    {
        std::ostringstream errorMessage;
        errorMessage << "ERROR: file unable to open in function computeEllipsoidSurfacePoints()";
        throw std::runtime_error( errorMessage.str( ) );
    }
    else
    {
        fileHandle << "X" << ",";
        fileHandle << "Y" << ",";
        fileHandle << "Z" << std::endl;
    }


    double azimuth = 0.0;
    double elevation = 0.0;
    double xCoordinate = 0.0;
    double yCoordinate = 0.0;
    double zCoordinate = 0.0;
    const double tolerance = 1.0e-10;

    for( azimuth = 0.0; azimuth < 360.0; azimuth = azimuth + stepSizeAzimuth )
    {
        for( elevation = 0.0; elevation <= 180.0; elevation = elevation + stepSizeElevation )
        {
            xCoordinate = alpha * std::cos( azimuth * naos::PI / 180.0 )
                                    * std::sin( elevation * naos::PI / 180.0 );
            yCoordinate = beta * std::sin( azimuth * naos::PI / 180.0 )
                                    * std::sin( elevation * naos::PI / 180.0 );
            zCoordinate = gamma * std::cos( elevation * naos::PI / 180.0 );

            // Perform a sanity-check here using the ellipsoidal equation from [1]
            double sum = ( xCoordinate * xCoordinate ) / ( alpha * alpha )
                            + ( yCoordinate * yCoordinate ) / ( beta * beta )
                                + ( zCoordinate * zCoordinate ) / ( gamma * gamma );
            if( std::abs( sum - 1.0 ) > tolerance )
            {
                std::ostringstream errorMessage;
                errorMessage << "ERROR: ellipsoidal equation sanity check failed!" << std::endl;
                throw std::runtime_error( errorMessage.str( ) );
            }

            fileHandle << xCoordinate << ",";
            fileHandle << yCoordinate << ",";
            fileHandle << zCoordinate << std::endl;
        }
    }
    fileHandle.close( );
}


} // namespace naos

//! References
/*!
 * [1] Wolfram Alpha online math resource. URL: http://mathworld.wolfram.com/Ellipsoid.html
       Accessed: 29th August 2016
 *
 */
