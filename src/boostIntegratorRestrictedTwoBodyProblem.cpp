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

#include <boost/numeric/odeint.hpp>

#include "NAOS/constants.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/basicAstro.hpp"
#include "NAOS/misc.hpp"

namespace naos
{

//! equations of motion (restricted two body problem)
/*!
 * first order differential equations describing the motion of a particle or spacecraft around a
 * central body modeled as a point mass gravity.
 */
class equationsOfMotion
{
    // declare parameters (gravitational parameter - used for calculating accelerations)
    const double gravParameter;

public:
    // Default constructor with member initializer list, get the gravitational parameter
    equationsOfMotion( const double aGravParameter )
                    : gravParameter( aGravParameter )
    { }
    void operator() ( const std::vector< double > &stateVector,
                      std::vector< double > &dXdt,
                      const double currentTime )
    {
        // calculate the gravitational accelerations first
        std::vector< double > positionVector = { stateVector[ xPositionIndex ],
                                                 stateVector[ yPositionIndex ],
                                                 stateVector[ zPositionIndex ] };

        double radialDistance = vectorNorm( positionVector );
        double radialDistanceCube = radialDistance * radialDistance * radialDistance;

        double Ux = -1.0 * gravParameter * positionVector[ xPositionIndex ] / radialDistanceCube;

        double Uy = -1.0 * gravParameter * positionVector[ yPositionIndex ] / radialDistanceCube;

        double Uz = -1.0 * gravParameter * positionVector[ zPositionIndex ] / radialDistanceCube;

        // now calculate the derivatives
        dXdt[ xPositionIndex ] = stateVector[ xVelocityIndex ];
        dXdt[ yPositionIndex ] = stateVector[ yVelocityIndex ];
        dXdt[ zPositionIndex ] = stateVector[ zVelocityIndex ];
        dXdt[ xVelocityIndex ] = Ux;
        dXdt[ yVelocityIndex ] = Uy;
        dXdt[ zVelocityIndex ] = Uz;
    }
};

//! Store intermediate state values and time( if needed )
/*!
 * This structure contains members that will save all intermediate state values and times when
 * an object of this structure is passed as an argument to the integrator function.
 */
struct pushBackStateAndTime
{
    // declare containers to store state and time
    std::vector< std::vector< double > > &stateContainer;
    std::vector< double > &timeContainer;

    //member initializer list
    pushBackStateAndTime( std::vector< std::vector< double > > &aState,
                          std::vector< double > &aTime )
                : stateContainer( aState ),
                  timeContainer( aTime )
    { }

    void operator() ( const std::vector< double > &singleStateVector, const double singleTime )
    {
        // store the intermediate state and time values in the containers
        stateContainer.push_back( singleStateVector );
        timeContainer.push_back( singleTime );
    }
};

//! Restricted two body problem integration
/*!
 * integrate the equations of motion for a particle around a point mass gravitational centre.
 */
void boostIntegratorRestrictedTwoBodyProblem( const double gravParameter,
                                              std::vector< double > &initialOrbitalElements,
                                              const double initialStepSize,
                                              const double startTime,
                                              const double endTime,
                                              std::ostringstream &outputFilePath,
                                              const int dataSaveIntervals )
{
    //! open the output csv file to save data. Declare file headers.
    std::ofstream outputFile;
    outputFile.open( outputFilePath.str( ) );
    outputFile << "x" << ",";
    outputFile << "y" << ",";
    outputFile << "z" << ",";
    outputFile << "vx" << ",";
    outputFile << "vy" << ",";
    outputFile << "vz" << ",";
    outputFile << "t" << std::endl;
    outputFile.precision( 16 );

    //! convert the initial orbital elements to cartesian state
    std::vector< double > initialState( 6, 0.0 );
    initialState = convertKeplerianElementsToCartesianCoordinates( initialOrbitalElements,
                                                                   gravParameter );

    // set up boost odeint
    const double absoluteTolerance = 1.0e-10;
    const double relativeTolerance = 1.0e-10;
    typedef boost::numeric::odeint::runge_kutta_fehlberg78< std::vector< double > > stepperType;

    // state step size guess (at each step this initial guess will be used)
    double stepSizeGuess = initialStepSize;

    // initialize the ode system
    equationsOfMotion restrictedTwoBodyProblem( gravParameter );

    // initialize current state vector and time
    std::vector< double > currentStateVector = initialState;
    double currentTime = startTime;
    double intermediateEndTime = currentTime + dataSaveIntervals;

    // save the initial state vector
    outputFile << currentStateVector[ xPositionIndex ] << ",";
    outputFile << currentStateVector[ yPositionIndex ] << ",";
    outputFile << currentStateVector[ zPositionIndex ] << ",";
    outputFile << currentStateVector[ xVelocityIndex ] << ",";
    outputFile << currentStateVector[ yVelocityIndex ] << ",";
    outputFile << currentStateVector[ zVelocityIndex ] << ",";
    outputFile << currentTime << std::endl;

    // start the integration outer loop
    while( intermediateEndTime <= endTime )
    {
        // perform integration, integrated result stored in currentStateVector
        size_t steps = boost::numeric::odeint::integrate_adaptive(
                            make_controlled( absoluteTolerance, relativeTolerance, stepperType( ) ),
                            restrictedTwoBodyProblem,
                            currentStateVector,
                            currentTime,
                            intermediateEndTime,
                            stepSizeGuess );

        // update the time variables
        currentTime = intermediateEndTime;
        intermediateEndTime = currentTime + dataSaveIntervals;

        // save data
        outputFile << currentStateVector[ xPositionIndex ] << ",";
        outputFile << currentStateVector[ yPositionIndex ] << ",";
        outputFile << currentStateVector[ zPositionIndex ] << ",";
        outputFile << currentStateVector[ xVelocityIndex ] << ",";
        outputFile << currentStateVector[ yVelocityIndex ] << ",";
        outputFile << currentStateVector[ zVelocityIndex ] << ",";
        outputFile << currentTime << std::endl;

    } // end of outer while loop for integration

    outputFile.close( );
}

} // namespace naos
