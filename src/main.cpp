/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>

#include "NAOS/ellipsoidSurfacePoints.hpp"

int main( const int numberOfInputs, const char* inputArguments[ ] )
{
    std::ostringstream ellipsoidSurfacePointsFile;
    ellipsoidSurfacePointsFile << "../../data/ellipsoidSurfacePoints.csv";
    naos::computeEllipsoidSurfacePoints( 10.0, 10.0, 10.0, 1.0, 1.0 , ellipsoidSurfacePointsFile );

    return EXIT_SUCCESS;
}
