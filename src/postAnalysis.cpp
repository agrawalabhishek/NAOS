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

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/misc.hpp"

namespace naos
{

//! Calculate the jacobian
/*!
 * calculate the jacobian for the point mass gravity potential case.
 */
void calculateJacobianPointMassGravity( const std::vector< double > &angularVelocityVector,
                                        const double gravParameter,
                                        std::ostringstream &inputFilePath,
                                        std::ostringstream &outputFilePath )
{
    // open the input csv file
    std::ifstream data;
    data.open( inputFilePath.str( ) );

    // open and setup the headers for the output csv file
    std::ofstream outputData;
    outputData.open( outputFilePath.str( ) );
    outputData << "time" << ",";
    outputData << "jacobian" << std::endl;

    // initialize the variables to read the data
    std::string line;
    std::getline( data, line ); // reads the header of the input csv file

    // read data line by line and convert to orbital elements
    while( std::getline( data, line ) )
    {
        std::stringstream linestream( line );
        std::string cell;
        std::vector< double > stateVector( 8 );
        int index = 0;

        while( std::getline( linestream, cell, ',' ) )
        {
            std::stringstream stringData( cell );
            double floatingPointData = 1.0;
            stringData >> floatingPointData;
            stateVector[ index ] = floatingPointData;
            index++;
        }

        double xVelocitySquare = stateVector[ xVelocityIndex ] * stateVector[ xVelocityIndex ];

        double yVelocitySquare = stateVector[ yVelocityIndex ] * stateVector[ yVelocityIndex ];

        double zVelocitySquare = stateVector[ zVelocityIndex ] * stateVector[ zVelocityIndex ];

        double xPositionSquare = stateVector[ xPositionIndex ] * stateVector[ xPositionIndex ];

        double yPositionSquare = stateVector[ yPositionIndex ] * stateVector[ yPositionIndex ];

        double zPositionSquare = stateVector[ zPositionIndex ] * stateVector[ zPositionIndex ];

        double radialDistance = std::sqrt( xPositionSquare + yPositionSquare + zPositionSquare );

        double angularVelocitySquare
            = angularVelocityVector[ zPositionIndex ] * angularVelocityVector[ zPositionIndex ];

        double gravPotential = gravParameter / radialDistance;

        double jacobian = 0.5 * ( xVelocitySquare + yVelocitySquare + zVelocitySquare )
                            - 0.5 * angularVelocitySquare * ( xPositionSquare + yPositionSquare )
                            - gravPotential;

        outputData << stateVector[ 6 ] << ",";
        outputData << jacobian << std::endl;
    } // end of outer while loop

    // close all files
    data.close( );
    outputData.close( );
}

//! Convert cartesian states to orbital elements post simulation.
/*!
 * Convert the cartesian states stored in a csv file to orbital elements and store it in a seperate
 * csv file. The cartesian states are to be in body fixed frame. And the central body is assumed to
 * be in a unoform rotation mode. This information is important since the cartesian states will be
 * transformed from the body fixed frame to the inertial frame before converting it to the orbital
 * elements.
 */
void postSimulationOrbitalElementsConversion( const std::vector< double > &angularVelocityVector,
                                              const double gravParameter,
                                              std::ostringstream &inputFilePath,
                                              std::ostringstream &outputFilePath )
{
    // open the input csv file
    std::ifstream data;
    data.open( inputFilePath.str( ) );

    // open and setup the headers for the output csv file
    std::ofstream outputData;
    outputData.open( outputFilePath.str( ) );
    outputData << "semiMajor" << ",";
    outputData << "eccentricity" << ",";
    outputData << "inclination" << ",";
    outputData << "raan" << ",";
    outputData << "aop" << ",";
    outputData << "ta" << ",";
    outputData << "time" << std::endl;

    // initialize the variables to read the data
    std::string line;
    std::getline( data, line ); // reads the header of the input csv file

    // read data line by line and convert to orbital elements
    while( std::getline( data, line ) )
    {
        std::stringstream linestream( line );
        std::string cell;
        std::vector< double > stateVector( 8 );
        int index = 0;

        while( std::getline( linestream, cell, ',' ) )
        {
            std::stringstream stringData( cell );
            double floatingPointData = 1.0;
            stringData >> floatingPointData;
            stateVector[ index ] = floatingPointData;
            index++;
        }

        // note: only elements 0 to 5 of the stateVector are the actual cartesian data points
        std::vector< double > inertialState( 6, 0.0 );

        // calculate the frame rotation angle, the 7th element of the stateVector is time, for the
        // corresponding state values
        double phi = angularVelocityVector[ zPositionIndex ] * stateVector[ 6 ];

        // get the inertial position
        inertialState[ xPositionIndex ]
                    = stateVector[ xPositionIndex ] * std::cos( phi )
                    - stateVector[ yPositionIndex ] * std::sin( phi );

        inertialState[ yPositionIndex ]
                    = stateVector[ xPositionIndex ] * std::sin( phi )
                    + stateVector[ yPositionIndex ] * std::cos( phi );

        inertialState[ zPositionIndex ] = stateVector[ zPositionIndex ];

        // get the omega cross position vector for use in the transport theorem
        std::vector< double > bodyFramePositionVector = { stateVector[ xPositionIndex ],
                                                          stateVector[ yPositionIndex ],
                                                          stateVector[ zPositionIndex ] };

        std::vector< double > omegaCrossPosition( 3, 0.0 );
        omegaCrossPosition = crossProduct( angularVelocityVector, bodyFramePositionVector );

        // use the transport theorem to get the body frame velocity in inertial coordinates
        double xBodyFrameVelocityInertialCoordinates
                    = stateVector[ xVelocityIndex ] + omegaCrossPosition[ xPositionIndex ];

        double yBodyFrameVelocityInertialCoordinates
                    = stateVector[ yVelocityIndex ] + omegaCrossPosition[ yPositionIndex ];

        double zBodyFrameVelocityInertialCoordinates
                    = stateVector[ zVelocityIndex ] + omegaCrossPosition[ zPositionIndex ];

        // use the rotation matrix to get the velocity in the inertial frame
        inertialState[ xVelocityIndex ]
                    = xBodyFrameVelocityInertialCoordinates * std::cos( phi )
                    - yBodyFrameVelocityInertialCoordinates * std::sin( phi );

        inertialState[ yVelocityIndex ]
                    = xBodyFrameVelocityInertialCoordinates * std::sin( phi )
                    + yBodyFrameVelocityInertialCoordinates * std::cos( phi );

        inertialState[ zVelocityIndex ] = zBodyFrameVelocityInertialCoordinates;

        // get the orbital elements for the inertial state vector
        std::vector< double > orbitalElements( 6, 0.0 );
        convertCartesianCoordinatesToKeplerianElements( inertialState,
                                                        gravParameter,
                                                        orbitalElements );

        // store the orbital elements in the output csv file
        outputData << orbitalElements[ 0 ] << ",";
        outputData << orbitalElements[ 1 ] << ",";
        outputData << orbitalElements[ 2 ] << ",";
        outputData << orbitalElements[ 3 ] << ",";
        outputData << orbitalElements[ 4 ] << ",";
        outputData << orbitalElements[ 5 ] << ",";
        outputData << stateVector[ 6 ] << std::endl;
    } // end of outer while loop

    // close all files
    data.close( );
    outputData.close( );
}

} // namespace naos
