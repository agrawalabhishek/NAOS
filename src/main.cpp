/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

#include "NAOS/ellipsoidSurfacePoints.hpp"
#include "NAOS/cubicRoot.hpp"
#include "NAOS/constants.hpp"
#include "NAOS/ellipsoidGravitationalAcceleration.hpp"
#include "NAOS/orbiterEquationsOfMotion.hpp"
#include "NAOS/rk4.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/computeEllipsoidSurfaceGravitationalAcceleration.hpp"
#include "NAOS/executeOrbiterAroundURE.hpp"
#include "NAOS/ellipsoidPotential.hpp"

int main( const int numberOfInputs, const char* inputArguments[ ] )
{
    // Physical parameters for Asteroid Eros, all in SI units, modelled as an ellipsoid.
    const double alpha = 20.0 * 1.0e3;
    const double beta = 7.0 * 1.0e3;
    const double gamma = 7.0 * 1.0e3;
    const double density = 3.2 * ( 10.0e-3 ) / ( 10.0e-6 );
    const double mass = ( 4.0 * naos::PI / 3.0 ) * density * alpha * beta * gamma;
    const double gravitationalParameter = naos::GRAVITATIONAL_CONSTANT * mass;
    const double Wx = 0.0; // rotational rate around principal x axis [rad/s]
    const double Wy = 0.0; // rotational rate around principal y axis [rad/s]
    const double Wz = 0.00033118202125129593; // rotational rate around principal z axis [rad/s]
    naos::Vector3 W { Wx, Wy, Wz };
    const double Wmagnitude = std::sqrt( Wx * Wx + Wy * Wy + Wz * Wz );

    // Generate surface coordinates for the Asteroid
    std::ostringstream ellipsoidSurfacePointsFile;
    ellipsoidSurfacePointsFile << "../../data/ellipsoidSurfacePoints.csv";
    const double stepSizeAzimuthDegree = 10.0;
    const double stepSizeElevationDegree = 10.0;
    naos::computeEllipsoidSurfacePoints( alpha, beta, gamma,
                                         stepSizeAzimuthDegree, stepSizeElevationDegree,
                                         ellipsoidSurfacePointsFile );

    // Compute gravitational acceleration on the surface of the ellipsoidal asteroid
    naos::computeEllipsoidSurfaceGravitationalAcceleration( alpha, beta, gamma,
                                                            gravitationalParameter );

    // compute orbit trajectory around a URE for given initial conditions
    std::ostringstream orbiterAroundUREFilePath;
    orbiterAroundUREFilePath << "../../data/eomOrbiterURESolution.csv";
    const double semiMajor = 24000.0;
    const double eccentricity = 0.1;
    const double inclination = 10.0;
    const double RAAN = 50.0;
    const double AOP = 100.0;
    const double TA = 5.0;
    naos::Vector6 initialVector { semiMajor, eccentricity, inclination, RAAN, AOP, TA };
    const bool initialVectorIsCartesian = false;
    const double integrationStepSize = 0.01;
    const double startTime = 0.0;
    const double endTime = 10000.0;

    naos::executeOrbiterAroundURE( alpha,
                                   beta,
                                   gamma,
                                   gravitationalParameter,
                                   density,
                                   W,
                                   Wmagnitude,
                                   initialVectorIsCartesian,
                                   initialVector,
                                   integrationStepSize,
                                   startTime,
                                   endTime,
                                   orbiterAroundUREFilePath );

    return EXIT_SUCCESS;
}
