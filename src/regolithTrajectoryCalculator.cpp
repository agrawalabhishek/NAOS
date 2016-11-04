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
#include <vector>
#include <limits>
#include <stdexcept>

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/gslIntegratorOrbiterAroundURE.hpp"
#include "NAOS/regolithTrajectoryCalculator.hpp"
#include "NAOS/misc.hpp"

namespace naos
{

//! Compute the velocity vector of lofted regolith
/*!
 * This sub-routine computes the velocity vector, given a magnitude, a unit normal vector from the
 * point on the surface where the regolith is lofted from, and the two conic angles which describe
 * the velocity vector's direction relative to the normal vector. A backwards approach is used
 * to go from the velocity vector in the final rotated intermediate frame back to the body fixed
 * frame. More details are given in the thesis report and author's personal notes.
 *
 */
void computeRegolithVelocityVector( std::vector< double > regolithPositionVector,
                                    const double velocityMagnitude,
                                    const double coneAngleAzimuth,
                                    const double coneAngleDeclination,
                                    std::vector< double > &unitNormalVector,
                                    std::vector< double > &regolithVelocityVector )
{
    // form the velocity vector, assuming that the intermediate frame's z-axis, on the surface of
    // the asteroid, is along the final velocity vector
    regolithVelocityVector[ 0 ] = 0.0;
    regolithVelocityVector[ 1 ] = 0.0;
    regolithVelocityVector[ 2 ] = velocityMagnitude;

    // get the rotatin matrix to go from the rotated intermediate frame back to the initial
    // frame where the z-axis was along the normal axis and the x-axis was pointing towards north
    // direction
    std::vector< std::vector< double > > zBasicRotationMatrix
            { { std::cos( coneAngleAzimuth ), std::sin( coneAngleAzimuth ), 0.0 },
              { -std::sin( coneAngleAzimuth ), std::cos( coneAngleAzimuth ), 0.0 },
              { 0.0, 0.0, 1.0 } };

    std::vector< std::vector< double > > yBasicRotationMatrix
            { { std::cos( coneAngleDeclination ), 0.0, -std::sin( coneAngleDeclination ) },
              { 0.0, 1.0, 0.0 },
              { std::sin( coneAngleDeclination ), 0.0, std::cos( coneAngleDeclination ) } };

    std::vector< std::vector< double > > nonTransposedRotationMatrix( 3, std::vector< double > ( 3 ) );

    matrixMultiplication( yBasicRotationMatrix,
                          zBasicRotationMatrix,
                          nonTransposedRotationMatrix,
                          3,
                          3,
                          3,
                          3 );

    std::vector< std::vector< double > > intermediateFrameRotationMatrix( 3, std::vector< double > ( 3 ) );

    matrixTranspose( nonTransposedRotationMatrix, intermediateFrameRotationMatrix );

    std::vector< std::vector< double > > rotatedIntermediateFrameVelocityVector
                { { regolithVelocityVector[ 0 ] },
                  { regolithVelocityVector[ 1 ] },
                  { regolithVelocityVector[ 2 ] } };

    std::vector< std::vector< double > > intermediateFrameVelocityVector( 3, std::vector< double > ( 1 ) );

    matrixMultiplication( intermediateFrameRotationMatrix,
                          rotatedIntermediateFrameVelocityVector,
                          intermediateFrameVelocityVector,
                          3,
                          3,
                          3,
                          1 );

    // now obtain the basis vectors for the intermediate frame expressed in the
    // body fixed frame coordinates
    std::vector< double > xUnitVector( 3 );
    std::vector< double > yUnitVector( 3 );
    std::vector< double > zUnitVector( 3 );

    zUnitVector = unitNormalVector;

    std::vector< double > bodyFrameZUnitVector { 0.0, 0.0, 1.0 };

    // get the intermediate RTN frame at the surface point
    std::vector< double > unitR = normalize( regolithPositionVector );

    std::vector< double > unitT = normalize( crossProduct( unitR, bodyFrameZUnitVector ) );

    // get the x basis vector, pointing to north
    xUnitVector = normalize( crossProduct( unitT, zUnitVector ) );

    // get the y basis vector
    yUnitVector = normalize( crossProduct( zUnitVector, xUnitVector ) );

    std::vector< double > zPrincipalAxisBodyFrame { 0.0, 0.0, 1.0 };
    std::vector< double > zNegativePrincipalAxisBodyFrame { 0.0, 0.0, -1.0 };

    // check if the position vector is along the poles
    const double positionDotPrincipalZ
            = dotProduct( normalize( regolithPositionVector ), zPrincipalAxisBodyFrame );
    const double positionDotNegativePrincipalZ
            = dotProduct( normalize( regolithPositionVector ), zNegativePrincipalAxisBodyFrame );
    if( positionDotPrincipalZ == 1.0 )
    {
        // the position vector is pointing to the poles, hence x basis vector pointing to the north
        // direction wouldn't work
        xUnitVector = { 1.0, 0.0, 0.0 };
        yUnitVector = { 0.0, 1.0, 0.0 };
    }
    else if( positionDotNegativePrincipalZ == 1.0 )
    {
        xUnitVector = { -1.0, 0.0, 0.0 };
        yUnitVector = { 0.0, 1.0, 0.0 };
    }

    // yUnitVector = normalize( crossProduct( unitNormalVector, bodyFrameZUnitVector ) );

    // xUnitVector = normalize( crossProduct( yUnitVector, unitNormalVector ) );

    // put the basis vectors in a 3x3 matrix
    std::vector< std::vector< double > > intermediateFrameBasisMatrix
            { { xUnitVector[ 0 ], yUnitVector[ 0 ], zUnitVector[ 0 ] },
              { xUnitVector[ 1 ], yUnitVector[ 1 ], zUnitVector[ 1 ] },
              { xUnitVector[ 2 ], yUnitVector[ 2 ], zUnitVector[ 2 ] } };

    std::vector< std::vector< double > > bodyFrameVelocityVector( 3, std::vector< double >( 1 ) );

    matrixMultiplication( intermediateFrameBasisMatrix,
                          intermediateFrameVelocityVector,
                          bodyFrameVelocityVector,
                          3,
                          3,
                          3,
                          1 );

    // return the final regolith velocity vector, expressed in body frame coordinates
    regolithVelocityVector[ 0 ] = bodyFrameVelocityVector[ 0 ][ 0 ];
    regolithVelocityVector[ 1 ] = bodyFrameVelocityVector[ 1 ][ 0 ];
    regolithVelocityVector[ 2 ] = bodyFrameVelocityVector[ 2 ][ 0 ];
}

void computeRegolithVelocityVector2( std::vector< double > regolithPositionVector,
                                     const double velocityMagnitude,
                                     const double coneAngleAzimuth,
                                     const double coneAngleDeclination,
                                     std::vector< double > &unitNormalVector,
                                     std::vector< double > &regolithVelocityVector )
{
    // get the z basis vector of the surface frame
    std::vector< double > zUnitVector = unitNormalVector;

    std::vector< double > xUnitVector( 3 );
    std::vector< double > yUnitVector( 3 );

    // body frame principal axis Z
    std::vector< double > zPrincipalAxisBodyFrame { 0.0, 0.0, 1.0 };
    std::vector< double > zNegativePrincipalAxisBodyFrame { 0.0, 0.0, -1.0 };

    // check if the position vector is along the poles
    const double positionDotPrincipalZ
            = dotProduct( normalize( regolithPositionVector ), zPrincipalAxisBodyFrame );
    const double positionDotNegativePrincipalZ
            = dotProduct( normalize( regolithPositionVector ), zNegativePrincipalAxisBodyFrame );

    if( positionDotPrincipalZ == 1.0 )
    {
        // the position vector is pointing to the poles, hence x basis vector pointing to the north
        // direction wouldn't work
        xUnitVector = { 1.0, 0.0, 0.0 };
        yUnitVector = { 0.0, 1.0, 0.0 };
    }
    else if( positionDotNegativePrincipalZ == 1.0 )
    {
        xUnitVector = { -1.0, 0.0, 0.0 };
        yUnitVector = { 0.0, 1.0, 0.0 };
    }
    else
    {
        // do the regular thing where the x basis is pointing to north
        // get the intermediate RTN frame at the surface point
        std::vector< double > unitR = normalize( regolithPositionVector );

        std::vector< double > unitT = normalize( crossProduct( unitR, zPrincipalAxisBodyFrame ) );

        // get the x basis vector, pointing to north
        xUnitVector = normalize( crossProduct( unitT, zUnitVector ) );

        // get the y basis vector
        yUnitVector = normalize( crossProduct( zUnitVector, xUnitVector ) );
    }

    const double cosDelta = std::cos( coneAngleDeclination );
    const double sinDelta = std::sin( coneAngleDeclination );
    const double cosGamma = std::cos( coneAngleAzimuth );
    const double sinGamma = std::sin( coneAngleAzimuth );

    regolithVelocityVector[ 0 ] = velocityMagnitude * ( cosDelta * zUnitVector[ 0 ]
                                                        + sinDelta * cosGamma * xUnitVector[ 0 ]
                                                        + sinDelta * sinGamma * yUnitVector[ 0 ] );

    regolithVelocityVector[ 1 ] = velocityMagnitude * ( cosDelta * zUnitVector[ 1 ]
                                                        + sinDelta * cosGamma * xUnitVector[ 1 ]
                                                        + sinDelta * sinGamma * yUnitVector[ 1 ] );

    regolithVelocityVector[ 2 ] = velocityMagnitude * ( cosDelta * zUnitVector[ 2 ]
                                                        + sinDelta * cosGamma * xUnitVector[ 2 ]
                                                        + sinDelta * sinGamma * yUnitVector[ 2 ] );
}

//! Compute trajectory for regolith ejected from the surface of an asteroid
/*!
 * This routine computes the trajectory for a single regolith ejected from the surface of an
 * asteroid by first computing the appropriate initial conditions (IC) and then integrating the
 * equations of motion using those ICs.
 *
 */
void calculateRegolithTrajectory( const double alpha,
                                  const double beta,
                                  const double gamma,
                                  const double gravitationalParameter,
                                  const double density,
                                  const std::vector< double > W,
                                  const double Wmagnitude,
                                  const double aXValue,
                                  const double aYValue,
                                  const double aZValue,
                                  const double coneAngleAzimuth,
                                  const double coneAngleDeclination,
                                  const double velocityMagnitudeFactor,
                                  const double integrationStepSize,
                                  const double startTime,
                                  const double endTime,
                                  std::ostringstream &filePath )
{
    // get the initial position coordinates for the particle on the surface
    std::vector< double > aVector { aXValue, aYValue, aZValue };
    std::vector< double > positionUnitVector = naos::normalize( aVector );

    double xTempTerm
        = ( positionUnitVector[ xPositionIndex ] * positionUnitVector[ xPositionIndex ] )
            / ( alpha * alpha );
    double yTempTerm
        = ( positionUnitVector[ yPositionIndex ] * positionUnitVector[ yPositionIndex ] )
            / ( beta * beta );
    double zTempTerm
        = ( positionUnitVector[ zPositionIndex ] * positionUnitVector[ zPositionIndex ] )
            / ( gamma * gamma );
    double positionMagnitude = std::sqrt( 1.0 / ( xTempTerm + yTempTerm + zTempTerm ) );

    std::vector< double > regolithPositionVector( 3 );
    regolithPositionVector[ xPositionIndex ]
                            = positionMagnitude * positionUnitVector[ xPositionIndex ];
    regolithPositionVector[ yPositionIndex ]
                            = positionMagnitude * positionUnitVector[ yPositionIndex ];
    regolithPositionVector[ zPositionIndex ]
                            = positionMagnitude * positionUnitVector[ zPositionIndex ];

    // get the normal vector for the point on the surface from where the particle is launched
    // this normal vector is not in the same direction as the position vector unless the central
    // body is a sphere.
    std::vector< double > regolithNormalVector
                               = { regolithPositionVector[ xPositionIndex ] / ( alpha * alpha ),
                                   regolithPositionVector[ yPositionIndex ] / ( beta * beta ),
                                   regolithPositionVector[ zPositionIndex ] / ( gamma * gamma ) };

    // get the unit normal vector
    std::vector< double > regolithNormalUnitVector = naos::normalize( regolithNormalVector );

    // get the velocity vector
    positionMagnitude = vectorNorm( regolithPositionVector );
    const double velocityMagnitude = velocityMagnitudeFactor
                                    * std::sqrt( 2.0 * gravitationalParameter / positionMagnitude );
    std::vector< double > regolithVelocityVector( 3 );

    computeRegolithVelocityVector2( regolithPositionVector,
                                    velocityMagnitude,
                                    coneAngleAzimuth,
                                    coneAngleDeclination,
                                    regolithNormalUnitVector,
                                    regolithVelocityVector );

    naos::Vector6 initialVector { regolithPositionVector[ xPositionIndex ],
                                  regolithPositionVector[ yPositionIndex ],
                                  regolithPositionVector[ zPositionIndex ],
                                  regolithVelocityVector[ 0 ],
                                  regolithVelocityVector[ 1 ],
                                  regolithVelocityVector[ 2 ] };
    const bool initialVectorIsCartesian = true;

    naos::gslIntegratorOrbiterAroundURE( alpha,
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
                                         filePath );
}

} // namespace naos
