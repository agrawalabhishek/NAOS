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
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/computeEllipsoidSurfaceGravitationalAcceleration.hpp"
#include "NAOS/ellipsoidPotential.hpp"
#include "NAOS/springMassIntegratorTest.hpp"
#include "NAOS/postAnalysis.hpp"
#include "NAOS/particleAroundUniformlyRotatingEllipsoid.hpp"
#include "NAOS/regolithTrajectoryCalculator.hpp"
#include "NAOS/regolithMonteCarlo.hpp"
#include "NAOS/perturbingAccelerations.hpp"
#include "NAOS/asteroidSurfaceGravitationalAndPerturbingAccelerations.hpp"
#include "NAOS/computePerturbationsForSpecificCoordinates.hpp"
#include "NAOS/sunAsteroidKeplerProblemSolver.hpp"
// #include "NAOS/sunAsteroidTwoBodyProblem.hpp"
// #include "NAOS/particleAroundSpheroidAndElllipsoidGravitationalPotential.hpp"
// #include "NAOS/boostIntegratorRestrictedTwoBodyProblem.hpp"
// #include "NAOS/executeOrbiterAroundUREPointMassGravity.hpp"
// #include "NAOS/executeOrbiterAroundURE.hpp"
// #include "NAOS/orbiterEquationsOfMotion.hpp"

int main( const int numberOfInputs, const char* inputArguments[ ] )
{
    // get the mode specified by the user
    std::string userMode = inputArguments[ 1 ];
    std::cout << std::endl << "Executing the following mode: ";
    std::cout << userMode << std::endl << std::endl;

    // Physical parameters for Asteroid Eros, all in SI units, modelled as an ellipsoid.
    const double alpha = 20.0 * 1.0e3;
    const double beta = 7.0 * 1.0e3;
    const double gamma = 7.0 * 1.0e3;
    const double density = 3.2 * ( 1.0e-3 ) / ( 1.0e-6 );
    const double mass = ( 4.0 * naos::PI / 3.0 ) * density * alpha * beta * gamma;
    const double gravitationalParameter = naos::GRAVITATIONAL_CONSTANT * mass;
    std::cout << std::endl << std::endl << "Gravitational Parameter = " << gravitationalParameter;
    std::cout << std::endl << std::endl;
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

    else if( userMode.compare( "executeSunAsteroidTwoBodyProblem" ) == 0 )
    {
        const double integrationStepSize = 0.01;
        const double dataSaveIntervals = 10.0;

        double wallTimeStart = naos::getWallTime< double >( );
        double cpuTimeStart = naos::getCPUTime< double >( );

        std::ostringstream sunAsteroidFilePath;
        sunAsteroidFilePath << "../../data/regolith_launched_from_longest_edge";
        sunAsteroidFilePath << "/multiple_launch_velocity_with_perturbations";
        sunAsteroidFilePath << "/simulation_time_9_months";
        sunAsteroidFilePath << "/3.2Density_1cmSize/";
        sunAsteroidFilePath << "sunEphemeris_phase45.csv";

        std::ostringstream integratedSunAsteroidFilePath;
        integratedSunAsteroidFilePath << "../../data/sun_asteroid_2BP/sunEphemeris_phase45.csv";

        // data taken from url:
        // http://ssd.jpl.nasa.gov/sbdb.cgi?sstr=433
        // accessed 3 jan 2017
        // std::vector< double > initialOrbitalElements = { 1.457945652635353 * naos::oneAstronomicalUnit,
        //                                                  0.2225680937603629,
        //                                                  10.82771477612614,
        //                                                  304.3265065906873,
        //                                                  178.8050095729968,
        //                                                  0.0 };

        std::vector< double > initialOrbitalElements = { 1.0 * naos::oneAstronomicalUnit,
                                                         0.0,
                                                         0.0,
                                                         0.0,
                                                         0.0,
                                                         45.0 };

        // calculate the starting time based on the true anomaly
        // double semiMajorAxis = initialOrbitalElements[ 0 ];
        // double semiMajorAxisCube = semiMajorAxis * semiMajorAxis * semiMajorAxis;
        // double sunMeanMotion = std::sqrt( naos::sunGravParameter / semiMajorAxisCube );
        // double trueAnomalyRadian = naos::convertDegreeToRadians( initialOrbitalElements[ 5 ] );
        // const double startTime = trueAnomalyRadian / sunMeanMotion;
        // const double endTime = startTime + ( 10.0 * 24.0 * 60.0 * 60.0 );

        const double startTime = 0.0;
        const double endTime = 270.0 * 24.0 * 60.0 * 60.0;

        // naos::executeSunAsteroidTwoBodyProblem( naos::sunGravParameter,
        //                                         W,
        //                                         initialOrbitalElements,
        //                                         integrationStepSize,
        //                                         startTime,
        //                                         endTime,
        //                                         integratedSunAsteroidFilePath,
        //                                         dataSaveIntervals );

        naos::executeSunAsteroidKeplerSolver( naos::sunGravParameter,
                                              W,
                                              initialOrbitalElements,
                                              startTime,
                                              endTime,
                                              sunAsteroidFilePath,
                                              dataSaveIntervals );

        double wallTimeEnd = naos::getWallTime< double >( );
        double cpuTimeEnd = naos::getCPUTime< double >( );

        std::cout << "Total wall time for execution = " << wallTimeEnd - wallTimeStart << std::endl;
        std::cout << "Total CPU time for execution = " << cpuTimeEnd - cpuTimeStart << std::endl;
    }

    else if( userMode.compare( "testPerturbations" ) == 0 )
    {
        double wallTimeStart = naos::getWallTime< double >( );
        double cpuTimeStart = naos::getCPUTime< double >( );

        std::vector< double > testSunThirdBodyPerturbingAcceleration( 3, 0.0 );
        std::vector< double > testSolarRadiationPerturbingAcceleration( 3, 0.0 );

        // Test data taken from TUDAT - general 3D case
        std::vector< double > targetPositionVector = { -40000000.0,
                                                       9000000.0,
                                                       -9500000.0 };

        std::vector< double > perturberPositionVector = { 25000000.0,
                                                          -380000000.0,
                                                          -55000000.0 };
        const double testGravParameter = 4900.0e9;

        std::vector< double > expectedRealisticAcceleration =
                                    { 2.9394638802044552864534126789629e-6,
                                      2.2253948636344714462843917227786e-6,
                                      1.1680066249314399399805985890762e-6 };

        std::vector< double > thirdBodyEffectAcceleration( 3, 0.0 );
        thirdBodyEffectAcceleration = naos::computeThirdBodyEffect( perturberPositionVector,
                                                                    targetPositionVector,
                                                                    testGravParameter );

        std::cout << "Test for third body effect from TUDAT data:";
        naos::printVector( thirdBodyEffectAcceleration, 3 );

        std::cout << "\nExpected acceleration values for TUDAT test data:";
        naos::printVector( expectedRealisticAcceleration, 3 );

        // Regolith sun specific test for third body effect
        std::vector< double > testRegolithPositionVector = { 25000.0,
                                                             0.0,
                                                             0.0 };

        std::vector< double > sunPositionVector = { -1.0 * naos::oneAstronomicalUnit,
                                                    0.0,
                                                    0.0 };

        thirdBodyEffectAcceleration = naos::computeThirdBodyEffect( sunPositionVector,
                                                                    testRegolithPositionVector,
                                                                    naos::sunGravParameter );

        // expected value for the acceleration when grav parameter = 1 in above computation
        // expectedRealisticAcceleration = { 1.4934602208384173e-29, 0.0, 0.0 };

        // expected value of acceleration when grav parameter is the Sun's grav parameter
        expectedRealisticAcceleration = { 1.9820074997728747e-09, 0.0, 0.0 };

        std::cout << "third body perturbing acceleration test with Sun position vector:";
        naos::printVector( thirdBodyEffectAcceleration, 3 );

        std::cout << "\nHand calculated expected acceleration values:";
        naos::printVector( expectedRealisticAcceleration, 3 );

        // Regolith sun test with sun orbital elements
        // set initial conditions for the Sun (to compute perturbations)
        std::vector< double > initialSunOrbitalElements = { 1.0 * naos::oneAstronomicalUnit,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            180.0 };

        // set the time at which the perturbation has to be calculated
        double perturbationTime = 0.0;

        testSunThirdBodyPerturbingAcceleration
            = naos::computeSunThirdBodyEffectAcceleration( testRegolithPositionVector,
                                                           W,
                                                           initialSunOrbitalElements,
                                                           perturbationTime );

        std::cout << "third body perturbing acceleration test with Sun orbital elements:";
        naos::printVector( testSunThirdBodyPerturbingAcceleration, 3 );

        // SRP test based on TUDAT data 1 (spacecraft around venus)
        std::vector< double > positionVectorTargetToSun
                        = { 0.732 * naos::oneAstronomicalUnit / std::sqrt( 2.0 ),
                            0.732 * naos::oneAstronomicalUnit / std::sqrt( 2.0 ),
                            0.0 };

        double targetAlbedo = 0.5;
        double incidentArea = 0.005;
        double targetMass = 0.0022;
        const double solarConstantTUDAT = 1.0205062450596109e+17;

        testSolarRadiationPerturbingAcceleration = naos::computeSRPAcceleration(
                                                        positionVectorTargetToSun,
                                                        targetAlbedo,
                                                        incidentArea,
                                                        targetMass,
                                                        solarConstantTUDAT );

        std::vector< double > expectedSRPAcceleration = { -2.05147517201883e-05,
                                                          -2.05147517201883e-05,
                                                          0.0 };

        std::cout << "SRP test acceleration values from TUDAT data 1:";
        naos::printVector( testSolarRadiationPerturbingAcceleration, 3 );

        std::cout << "SRP expected acceleration values for TUDAT data 1:";
        naos::printVector( expectedSRPAcceleration, 3 );

        // SRP test based on TUDAT data 2 (random point away from Sun)
        positionVectorTargetToSun = { 94359740.25,
                                      90831886.1,
                                      14668782.92 };

        targetAlbedo = 0.4058;
        incidentArea = 514701.9505;
        targetMass = 1.0;

        expectedSRPAcceleration = { -3043733.21422537 / targetMass,
                                    -2929936.30441141 / targetMass,
                                    -473166.433773283 / targetMass };

        testSolarRadiationPerturbingAcceleration = naos::computeSRPAcceleration(
                                                        positionVectorTargetToSun,
                                                        targetAlbedo,
                                                        incidentArea,
                                                        targetMass,
                                                        solarConstantTUDAT );

        std::cout << "SRP test acceleration values from TUDAT data 2:";
        naos::printVector( testSolarRadiationPerturbingAcceleration, 3 );

        std::cout << "SRP expected acceleration values for TUDAT data 2:";
        naos::printVector( expectedSRPAcceleration, 3 );

        // Conducts test(s) for SRP now based on initial sun's orbital elements
        testRegolithPositionVector = { 0.0,
                                       0.9 * naos::oneAstronomicalUnit,
                                       0.0 };

       initialSunOrbitalElements = { 1.0 * naos::oneAstronomicalUnit,
                                     0.0,
                                     0.0,
                                     0.0,
                                     0.0,
                                     90.0 };

        perturbationTime = 0.0;

        testSolarRadiationPerturbingAcceleration
            = naos::computeSolarRadiationPressureAcceleration( testRegolithPositionVector,
                                                               W,
                                                               initialSunOrbitalElements,
                                                               perturbationTime );

        std::cout << "SRP acceleration test with Sun orbital elements:";
        naos::printVector( testSolarRadiationPerturbingAcceleration, 3 );

        double wallTimeEnd = naos::getWallTime< double >( );
        double cpuTimeEnd = naos::getCPUTime< double >( );

        std::cout << "Total wall time for execution = " << wallTimeEnd - wallTimeStart << std::endl;
        std::cout << "Total CPU time for execution = " << cpuTimeEnd - cpuTimeStart << std::endl;
    }

    else if( userMode.compare( "executeRegolithMonteCarlo" ) == 0 )
    {
        const double integrationStepSize = 0.01;
        double phaseAngle = naos::convertDegreeToRadians( 0.0 );
        double Wz = W[ 2 ];
        double timeCorrespondingToPhaseAngle = phaseAngle / Wz;
        const double startTime = timeCorrespondingToPhaseAngle;
        const double endTime = 9.0 * 30.0 * 24.0 * 60.0 * 60.0;
        const double dataSaveIntervals = 10.0;

        double wallTimeStart = naos::getWallTime< double >( );
        double cpuTimeStart = naos::getCPUTime< double >( );

        std::ostringstream databaseFilePath;
        databaseFilePath << "../../data/regolith_launched_from_trailing_edge";
        databaseFilePath << "/multiple_launch_velocity_with_perturbations";
        databaseFilePath << "/simulation_time_9_months";
        // databaseFilePath << "/7.5Density_5cmRadius";
        databaseFilePath << "/trailingEdge_3P2Density_1cmRadius.db";
        // databaseFilePath << "test.db";

        // databaseFilePath << "../../data/VandV";
        // databaseFilePath << "/new";
        // databaseFilePath << "/orbital_dynamics";
        // databaseFilePath << "/finalFateValidation.db";

        // databaseFilePath << "../../data/regolith_launched_from_longest_edge";
        // databaseFilePath << "/multiple_launch_velocity_with_perturbations";
        // databaseFilePath << "/simulation_time_9_months";
        // databaseFilePath << "/captureTest2.db";

        naos::executeRegolithMonteCarlo( alpha,
                                         beta,
                                         gamma,
                                         gravitationalParameter,
                                         W,
                                         Wmagnitude,
                                         integrationStepSize,
                                         startTime,
                                         endTime,
                                         dataSaveIntervals,
                                         databaseFilePath );

        double wallTimeEnd = naos::getWallTime< double >( );
        double cpuTimeEnd = naos::getCPUTime< double >( );

        std::cout << "Total wall time for execution = " << wallTimeEnd - wallTimeStart << std::endl;
        std::cout << "Total CPU time for execution = " << cpuTimeEnd - cpuTimeStart << std::endl;
    }

    return EXIT_SUCCESS;
}
