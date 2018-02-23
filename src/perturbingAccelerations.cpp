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
#include <limits>
#include <vector>
#include <stdexcept>

#include "NAOS/basicAstro.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/constants.hpp"
#include "NAOS/misc.hpp"
#include "NAOS/sunAsteroidKeplerProblemSolver.hpp"

namespace naos
{

//! open sun ephemeris and extract data for a specified time value
/*!
 * This file opens the ephemeris data for the sun around the asteroid at specified time values in
 * multiples of ten. If a finer ephemeris is required then the sunAsteroidTwoBodyProblem must be
 * re-run prior to using this mode, with a finer data save interval value.
 */
std::vector< double > extractSunEphemeris( const double timeValue,
                                           std::ostringstream &ephemerisFilePath )
{
    // open the file
    std::ifstream data;
    data.open( ephemerisFilePath.str( ) );

    // initialize the variables to read the data
    std::string line;
    std::getline( data, line ); // reads the header of the input csv file

    const int timeIndex = 6;

    std::vector< double > bodyFrameStateVectorAndTime( 7, 0.0 );

    bool valueFoundFlag = false;

    // read data line by line and convert to orbital elements
    while( std::getline( data, line ) )
    {
        std::stringstream linestream( line );
        std::string cell;
        std::vector< double > dataVector( 21, 0.0 );
        int index = 0;

        while( std::getline( linestream, cell, ',' ) )
        {
            std::stringstream stringData( cell );
            double floatingPointData = 1.0;
            stringData >> floatingPointData;
            dataVector[ index ] = floatingPointData;
            index++;
        }

        if( std::fabs( dataVector[ timeIndex ] - timeValue ) <= 1.0e-10 )
        {
            valueFoundFlag = true;
            for( int i = 0; i < 7; i++ )
            {
                bodyFrameStateVectorAndTime[ i ] = dataVector[ i ];
            }
            // break out of outer most while loop
            break;
        }
    }

    if( !valueFoundFlag )
    {
        std::ostringstream errorMessage;
        errorMessage << std::endl;
        errorMessage << "Error in calculating perturbing accelerations cpp file" << std::endl;
        errorMessage << "Sun's ephemeris could not be obtained because of illegal time value" << std::endl;
        errorMessage << "Time value = " << timeValue << std::endl;
        errorMessage << std::endl;
        throw std::runtime_error( errorMessage.str( ) );
    }
    // close the csv file
    data.close( );

    // return the data
    return bodyFrameStateVectorAndTime;
}

//! Third body effect function based only on position vectors as input - general function
/*!
 * Compute the perturbing acceleration from the thrid body effect of any perturbing body
 */
std::vector< double > computeThirdBodyEffect( std::vector< double > perturberPositionVector,
                                              std::vector< double > targetPositionVector,
                                              const double perturberGravParameter )
{
    // gravitational parameter of the perturbing body
    const double mu = perturberGravParameter;

    double perturberPositionVectorMagnitude = vectorNorm( perturberPositionVector );

    double perturberPositionVectorMagnitudeCube
        = perturberPositionVectorMagnitude * perturberPositionVectorMagnitude * perturberPositionVectorMagnitude;

    // compute the perturbing acceleration
    std::vector< double > targetPositionFromPerturber( 3, 0.0 );
    targetPositionFromPerturber[ xPositionIndex ]
        = targetPositionVector[ xPositionIndex ] - perturberPositionVector[ xPositionIndex ];

    targetPositionFromPerturber[ yPositionIndex ]
        = targetPositionVector[ yPositionIndex ] - perturberPositionVector[ yPositionIndex ];

    targetPositionFromPerturber[ zPositionIndex ]
        = targetPositionVector[ zPositionIndex ] - perturberPositionVector[ zPositionIndex ];

    double targetPositionFromPerturberMagnitude = vectorNorm( targetPositionFromPerturber );

    double targetPositionFromPerturberMagnitudeCube
        = targetPositionFromPerturberMagnitude * targetPositionFromPerturberMagnitude * targetPositionFromPerturberMagnitude;

    std::vector< double > accelerationTermOne( 3, 0.0 );
    std::vector< double > accelerationTermTwo( 3, 0.0 );

    accelerationTermOne[ xPositionIndex ]
        = ( -1.0 * mu / targetPositionFromPerturberMagnitudeCube ) * targetPositionFromPerturber[ xPositionIndex ];

    accelerationTermOne[ yPositionIndex ]
        = ( -1.0 * mu / targetPositionFromPerturberMagnitudeCube ) * targetPositionFromPerturber[ yPositionIndex ];

    accelerationTermOne[ zPositionIndex ]
        = ( -1.0 * mu / targetPositionFromPerturberMagnitudeCube ) * targetPositionFromPerturber[ zPositionIndex ];

    accelerationTermTwo[ xPositionIndex ]
        = ( -1.0 * mu / perturberPositionVectorMagnitudeCube ) * perturberPositionVector[ xPositionIndex ];

    accelerationTermTwo[ yPositionIndex ]
        = ( -1.0 * mu / perturberPositionVectorMagnitudeCube ) * perturberPositionVector[ yPositionIndex ];

    accelerationTermTwo[ zPositionIndex ]
        = ( -1.0 * mu / perturberPositionVectorMagnitudeCube ) * perturberPositionVector[ zPositionIndex ];

    std::vector< double > thirdBodyPerturbingAcceleration( 3, 0.0 );
    thirdBodyPerturbingAcceleration[ xPositionIndex ]
        = accelerationTermOne[ xPositionIndex ] + accelerationTermTwo[ xPositionIndex ];

    thirdBodyPerturbingAcceleration[ yPositionIndex ]
        = accelerationTermOne[ yPositionIndex ] + accelerationTermTwo[ yPositionIndex ];

    thirdBodyPerturbingAcceleration[ zPositionIndex ]
        = accelerationTermOne[ zPositionIndex ] + accelerationTermTwo[ zPositionIndex ];

    return thirdBodyPerturbingAcceleration;
}

//! Third-body effect
/*!
 * Compute the perturbing acceleration from the thrid body effect of the sun.
 */
std::vector< double > computeSunThirdBodyEffectAcceleration( const std::vector< double > &regolithPositionVector,
                                                             std::vector< double > &asteroidRotationVector,
                                                             const std::vector< double > &initialSunOrbitalElements,
                                                             const double timeValue )
{
    std::vector< double > sunStateVectorInertialFrame( 6, 0.0 );
    std::vector< double > sunStateVectorInRotatingFrame( 6, 0.0 );
    computeSunStateVector( asteroidRotationVector,
                           initialSunOrbitalElements,
                           timeValue,
                           sunStateVectorInertialFrame,
                           sunStateVectorInRotatingFrame );

    // form the body frame position vector of the sun
    std::vector< double > sunPositionVector = { sunStateVectorInRotatingFrame[ xPositionIndex ],
                                                sunStateVectorInRotatingFrame[ yPositionIndex ],
                                                sunStateVectorInRotatingFrame[ zPositionIndex ] };

    // grav parameter of Sun
    const double mu = naos::sunGravParameter;

    std::vector< double > thirdBodyPerturbingAcceleration( 3, 0.0 );
    thirdBodyPerturbingAcceleration = computeThirdBodyEffect( sunPositionVector,
                                                              regolithPositionVector,
                                                              mu );

    return thirdBodyPerturbingAcceleration;
}

//! Solar radiation pressure general function
/*!
 * General computation of the perturbing acceleration from solar radition pressure.
 */
std::vector< double > computeSRPAcceleration( std::vector< double > &positionVectorTargetToSource,
                                              const double targetAlbedo,
                                              const double incidentArea,
                                              const double targetMass,
                                              const double solarConstant )
{
    double targetSunDistance = vectorNorm( positionVectorTargetToSource );

    double targetSunDistanceCube
        = targetSunDistance * targetSunDistance * targetSunDistance;

    const double areaToMassRatio = incidentArea / targetMass;

    // compute the acceleration
    const double multiplicationConstant = -1.0 * ( 1.0 + targetAlbedo ) * solarConstant * areaToMassRatio;

    std::vector< double > solarRadiationPerturbingAcceleration( 3, 0.0 );

    solarRadiationPerturbingAcceleration[ xPositionIndex ]
        = multiplicationConstant * positionVectorTargetToSource[ xPositionIndex ] / targetSunDistanceCube;

    solarRadiationPerturbingAcceleration[ yPositionIndex ]
        = multiplicationConstant * positionVectorTargetToSource[ yPositionIndex ] / targetSunDistanceCube;

    solarRadiationPerturbingAcceleration[ zPositionIndex ]
        = multiplicationConstant * positionVectorTargetToSource[ zPositionIndex ] / targetSunDistanceCube;

    return solarRadiationPerturbingAcceleration;
}

//! Solar radiation pressure
/*!
 * comupte the perturbing acceleration from solar radition pressure.
 */
std::vector< double > computeSolarRadiationPressureAcceleration( const std::vector< double > &regolithPositionVector,
                                                                 std::vector< double > &asteroidRotationVector,
                                                                 const std::vector< double > &initialSunOrbitalElements,
                                                                 const double timeValue )
{
    std::vector< double > sunStateVectorInertialFrame( 6, 0.0 );
    std::vector< double > sunStateVectorInRotatingFrame( 6, 0.0 );
    computeSunStateVector( asteroidRotationVector,
                           initialSunOrbitalElements,
                           timeValue,
                           sunStateVectorInertialFrame,
                           sunStateVectorInRotatingFrame );

    // form the body frame position vector of the sun
    std::vector< double > sunPositionVector = { sunStateVectorInRotatingFrame[ xPositionIndex ],
                                                sunStateVectorInRotatingFrame[ yPositionIndex ],
                                                sunStateVectorInRotatingFrame[ zPositionIndex ] };

    // position vector regolith to Sun
    std::vector< double > positionVectorRegolithToSun( 3, 0.0 );
    positionVectorRegolithToSun[ xPositionIndex ]
        = -1.0 * regolithPositionVector[ xPositionIndex ] + sunPositionVector[ xPositionIndex ];

    positionVectorRegolithToSun[ yPositionIndex ]
        = -1.0 * regolithPositionVector[ yPositionIndex ] + sunPositionVector[ yPositionIndex ];

    positionVectorRegolithToSun[ zPositionIndex ]
        = -1.0 * regolithPositionVector[ zPositionIndex ] + sunPositionVector[ zPositionIndex ];

    // specify the solar constant value in SI units [kg m/s^2]
    const double solarConstant = 1.0e17;

    // specify the albedo value, for max acceleration from SRP, set rho=1
    const double rho = 1.0;

    // specify the area to mass ratio for the regolith
    const double regolithGrainDensity = 3.2 * 1.0e3;
    const double regolithGrainRadius = 1.0 * 1.0e-2;

    const double regolithCrossSectionalArea = naos::PI *
                                                ( regolithGrainRadius * regolithGrainRadius );

    const double regolithMass = ( 4.0 / 3.0 ) * naos::PI *
        ( regolithGrainRadius * regolithGrainRadius * regolithGrainRadius ) * regolithGrainDensity;

    std::vector< double > solarRadiationPerturbingAcceleration( 3, 0.0 );

    solarRadiationPerturbingAcceleration = computeSRPAcceleration( positionVectorRegolithToSun,
                                                                   rho,
                                                                   regolithCrossSectionalArea,
                                                                   regolithMass,
                                                                   solarConstant );

    return solarRadiationPerturbingAcceleration;
}

} // namespace naos
