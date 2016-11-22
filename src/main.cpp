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
#include <cstring>

#include "NAOS/ellipsoidSurfacePoints.hpp"
#include "NAOS/cubicRoot.hpp"
#include "NAOS/constants.hpp"
#include "NAOS/ellipsoidGravitationalAcceleration.hpp"
#include "NAOS/orbiterEquationsOfMotion.hpp"
#include "NAOS/rk4.hpp"
#include "NAOS/rk54.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/computeEllipsoidSurfaceGravitationalAcceleration.hpp"
#include "NAOS/executeOrbiterAroundURE.hpp"
#include "NAOS/ellipsoidPotential.hpp"
#include "NAOS/springMassIntegratorTest.hpp"
#include "NAOS/executeOrbiterAroundUREPointMassGravity.hpp"
#include "NAOS/postAnalysis.hpp"
#include "NAOS/particleAroundSpheroidAndElllipsoidGravitationalPotential.hpp"
// #include "NAOS/boostIntegratorRestrictedTwoBodyProblem.hpp"
// #include "NAOS/gslIntegratorOrbiterAroundUREPointMassGravity.hpp"
// #include "NAOS/gslIntegratorOrbiterAroundURE.hpp"
// #include "NAOS/regolithTrajectoryCalculator.hpp"
// #include "NAOS/bulirschStoer.hpp"

int main( const int numberOfInputs, const char* inputArguments[ ] )
{
    // get the mode specified by the user
    std::string userMode = inputArguments[ 1 ];
    std::cout << std::endl << "Executing the following mode: ";
    std::cout << userMode << std::endl << std::endl;

    // Physical parameters for Asteroid Eros, all in SI units, modelled as an ellipsoid.
    const double alpha = 20.0 * 1.0e3;
    const double beta = 20.0 * 1.0e3;
    const double gamma = 20.0 * 1.0e3;
    const double density = 3.2 * ( 10.0e-3 ) / ( 10.0e-6 );
    const double mass = ( 4.0 * naos::PI / 3.0 ) * density * alpha * beta * gamma;
    const double gravitationalParameter = naos::GRAVITATIONAL_CONSTANT * mass;
    const double Wx = 0.0; // rotational rate around principal x axis [rad/s]
    const double Wy = 0.0; // rotational rate around principal y axis [rad/s]
    const double Wz = 0.00033118202125129593; // rotational rate around principal z axis [rad/s]
    // const double Wz = 0.0;
    naos::Vector3 W { Wx, Wy, Wz };
    const double Wmagnitude = std::sqrt( Wx * Wx + Wy * Wy + Wz * Wz );

    if( userMode.compare( "generateAsteroidSurfaceCoordinates" ) == 0 )
    {
        // Generate surface coordinates for the Asteroid
        std::ostringstream ellipsoidSurfacePointsFile;
        ellipsoidSurfacePointsFile << "../../data/ellipsoidSurfacePoints.csv";
        const double stepSizeAzimuthDegree = 10.0;
        const double stepSizeElevationDegree = 10.0;
        naos::computeEllipsoidSurfacePoints( alpha, beta, gamma,
                                             stepSizeAzimuthDegree, stepSizeElevationDegree,
                                             ellipsoidSurfacePointsFile );
    }

    else if( userMode.compare( "computeAsteroidSurfaceAcceleration" ) == 0 )
    {
        // Compute gravitational acceleration on the surface of the ellipsoidal asteroid
        naos::computeEllipsoidSurfaceGravitationalAcceleration( alpha, beta, gamma,
                                                                gravitationalParameter );
    }

    else if( userMode.compare( "springMassIntegratorSanityCheck" ) == 0 )
    {
        // Perform sanity check for the numerical integrator routine using the spring mass system
        const double springConstant = 1.0;
        const double blockMass = 10.0;
        const double springMassStepSize = 0.01;
        const double springMassStartTime = 0.0;
        const double springMassEndTime = 1000.0;

        const std::vector< double > springMassInitialState { 0.2, 0.0 };

        std::ostringstream springMassFilePath;
        springMassFilePath << "../../data/springMassIntegrationSolution.csv";

        naos::executeSpringMassIntegration< std::vector< double > >( blockMass,
                                                                     springConstant,
                                                                     springMassInitialState,
                                                                     springMassStepSize,
                                                                     springMassStartTime,
                                                                     springMassEndTime,
                                                                     springMassFilePath );
    }

    // else if( userMode.compare( "executeRestricted2BP" ) == 0 )
    // {
    //     // Point mass gravity R2BP orbiter problem solution
    //     std::ostringstream pointMassFilePath;
    //     pointMassFilePath << "../../data/solutionrestricted2BP.csv";
    //     const double semiMajor = 35000.0;
    //     const double eccentricity = 0.1;
    //     const double inclination = 10.0;
    //     const double RAAN = 50.0;
    //     const double AOP = 100.0;
    //     const double TA = 0.0;

    //     naos::Vector6 initialVector { semiMajor,
    //                                   eccentricity,
    //                                   inclination,
    //                                   RAAN,
    //                                   AOP,
    //                                   TA };

    //     const double integrationStepSize = 0.01;
    //     const double startTime = 0.0;
    //     const double endTime = 24.0 * 30.0 * 24.0 * 60.0 * 60.0;
    //     const int dataSaveIntervals = 1000;

    //     double wallTimeStart = naos::getWallTime< double >( );
    //     double cpuTimeStart = naos::getCPUTime< double >( );

    //     naos::boostIntegratorRestrictedTwoBodyProblem( gravitationalParameter,
    //                                                    initialVector,
    //                                                    integrationStepSize,
    //                                                    startTime,
    //                                                    endTime,
    //                                                    pointMassFilePath,
    //                                                    dataSaveIntervals );

    //     double wallTimeEnd = naos::getWallTime< double >( );
    //     double cpuTimeEnd = naos::getCPUTime< double >( );

    //     std::cout << "Total wall time for execution = " << wallTimeEnd - wallTimeStart << std::endl;
    //     std::cout << "Total CPU time for execution = " << cpuTimeEnd - cpuTimeStart << std::endl;
    // }

    else if( userMode.compare( "executeParticleAroundSpheroid" ) == 0 )
    {
        // Particle around spheroid orbiter problem solution
        std::ostringstream filePath;
        filePath << "../../data/solutionParticleAroundRotatingSpheroid.csv";
        const double semiMajor = 35000.0;
        const double eccentricity = 0.1;
        const double inclination = 10.0;
        const double RAAN = 50.0;
        const double AOP = 10.0;
        const double TA = 0.0;

        naos::Vector6 initialVector { semiMajor,
                                      eccentricity,
                                      inclination,
                                      RAAN,
                                      AOP,
                                      TA };

        const double integrationStepSize = 0.01;
        const double startTime = 0.0;
        const double endTime = 2.0 * 365.0 * 24.0 * 60.0 * 60.0;
        const int dataSaveIntervals = 1000;

        double wallTimeStart = naos::getWallTime< double >( );
        double cpuTimeStart = naos::getCPUTime< double >( );

        naos::executeParticleAroundSpheroid( alpha,
                                             gravitationalParameter,
                                             W,
                                             initialVector,
                                             integrationStepSize,
                                             startTime,
                                             endTime,
                                             filePath,
                                             dataSaveIntervals );

        double wallTimeEnd = naos::getWallTime< double >( );
        double cpuTimeEnd = naos::getCPUTime< double >( );

        std::cout << "Total wall time for execution = " << wallTimeEnd - wallTimeStart << std::endl;
        std::cout << "Total CPU time for execution = " << cpuTimeEnd - cpuTimeStart << std::endl;
    }

    // compute orbit trajectory around a URE for given initial conditions
    // std::ostringstream orbiterAroundUREFilePath;
    // orbiterAroundUREFilePath << "../../data/RK5_OrbiterAroundSphericalEros_longTermSimulation.csv";
    // const double semiMajor = 35000.0;
    // const double eccentricity = 0.1;
    // const double inclination = 10.0;
    // const double RAAN = 60.0;
    // const double AOP = 100.0;
    // const double TA = 0.0;
    // naos::Vector6 initialVector { semiMajor, eccentricity, inclination, RAAN, AOP, TA };
    // const bool initialVectorIsCartesian = false;
    // const double integrationStepSize = 0.01;
    // const double startTime = 0.0;
    // const double endTime = 10000.0;

    // naos::executeOrbiterAroundURE( alpha,
    //                                beta,
    //                                gamma,
    //                                gravitationalParameter,
    //                                density,
    //                                W,
    //                                Wmagnitude,
    //                                initialVectorIsCartesian,
    //                                initialVector,
    //                                integrationStepSize,
    //                                startTime,
    //                                endTime,
    //                                orbiterAroundUREFilePath );

    // std::ostringstream gslOrbiterAroundUREFilePath;
    // gslOrbiterAroundUREFilePath << "../../data/gsl_RK8_OrbiterAroundSphericalEros.csv";
    // naos::gslIntegratorOrbiterAroundURE( alpha,
    //                                      beta,
    //                                      gamma,
    //                                      gravitationalParameter,
    //                                      density,
    //                                      W,
    //                                      Wmagnitude,
    //                                      initialVectorIsCartesian,
    //                                      initialVector,
    //                                      integrationStepSize,
    //                                      startTime,
    //                                      endTime,
    //                                      gslOrbiterAroundUREFilePath );

    // compute trajectory evolution for single regolith lofted from the surface of an asteroid
    // std::ostringstream regolithAroundUREFilePath;
    // regolithAroundUREFilePath << "../../data/singleRegolithEjectaURESolution.csv";

    // const double integrationStepSize = 0.01;
    // const double startTime = 0.0;
    // const double endTime = 100.0;

    // const double aXValue = 1.0;
    // const double aYValue = 2.0;
    // const double aZValue = 0.0;
    // const double velocityMagnitudeFactor = 0.7;
    // const double coneAngleAzimuth = naos::convertDegreeToRadians( 90.0 );
    // const double coneAngleDeclination = naos::convertDegreeToRadians( 45.0 );

    // naos::calculateRegolithTrajectory( alpha,
    //                                    beta,
    //                                    gamma,
    //                                    gravitationalParameter,
    //                                    density,
    //                                    W,
    //                                    Wmagnitude,
    //                                    aXValue,
    //                                    aYValue,
    //                                    aZValue,
    //                                    coneAngleAzimuth,
    //                                    coneAngleDeclination,
    //                                    velocityMagnitudeFactor,
    //                                    integrationStepSize,
    //                                    startTime,
    //                                    endTime,
    //                                    regolithAroundUREFilePath );

    else if( userMode.compare( "postAnalysisRestricted2BP" ) == 0 )
    {
        // for point mass analysis (for R2BP)
        std::vector< double > centralBodyRotationRate { 0.0, 0.0, 0.0 };

        // specify the input file path
        std::ostringstream inputFilePath;
        inputFilePath << "../../data/solutionrestricted2BP.csv";

        // specify the output file path
        std::ostringstream pointMassOrbitalElementsOutputFilePath;
        pointMassOrbitalElementsOutputFilePath << "../../data/solutionrestricted2BP_orbitalElements.csv";

        // post analysis, conversion of cartesian data into orbital elements
        naos::postSimulationOrbitalElementsConversion( centralBodyRotationRate,
                                                       gravitationalParameter,
                                                       inputFilePath,
                                                       pointMassOrbitalElementsOutputFilePath );

        //! Calculate the jacobian
        std::ostringstream pointMassJacobianOutputFilePath;
        pointMassJacobianOutputFilePath << "../../data/solutionrestricted2BP_jacobian.csv";
        naos::calculateJacobianPointMassGravity( centralBodyRotationRate,
                                                 gravitationalParameter,
                                                 inputFilePath,
                                                 pointMassJacobianOutputFilePath );
    }

    else if( userMode.compare( "postAnalysisParticleAroundSpheroid" ) == 0 )
    {
        // specify the input file path
        std::ostringstream inputFilePath;
        inputFilePath << "../../data/solutionParticleAroundRotatingSpheroid.csv";

        // specify the output file path
        std::ostringstream OrbitalElementsOutputFilePath;
        OrbitalElementsOutputFilePath << "../../data/solutionParticleAroundRotatingSpheroid_orbitalElements.csv";

        // post analysis, conversion of cartesian data into orbital elements
        naos::postSimulationOrbitalElementsConversion( W,
                                                       gravitationalParameter,
                                                       inputFilePath,
                                                       OrbitalElementsOutputFilePath );

        //! Calculate the jacobian
        std::ostringstream JacobianOutputFilePath;
        JacobianOutputFilePath << "../../data/solutionParticleAroundRotatingSpheroid_jacobian.csv";
        naos::calculateJacobianEllipsoidPotentialModel( alpha,
                                                        beta,
                                                        gamma,
                                                        W,
                                                        gravitationalParameter,
                                                        inputFilePath,
                                                        JacobianOutputFilePath );
    }

    return EXIT_SUCCESS;
}
