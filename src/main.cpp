/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fsream>

#include "NAOS/ellipsoidSurfacePoints.hpp"

int main( const int numberOfInputs, const char* inputArguments[ ] )
{
    std::ostringstream ellipsoidSurfacePointsFile;
    ellipsoidSurfacePointsFile << "../../data/ellipsoidSurfacePoints.csv";
    naos::computeEllipsoidSurfacePoints( 12.0, 8.0, 5.0, 1.0, 1.0 , ellipsoidSurfacePointsFile );

    return EXIT_SUCCESS;
}
