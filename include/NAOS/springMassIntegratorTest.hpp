/*
 * Copyright (c) 2016 Abhishek Agrawal (abhishek.agrawal@protonmail.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef SPRING_MASS_INTEGRATOR_TEST_HPP
#define SPRING_MASS_INTEGRATOR_TEST_HPP

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

#include "NAOS/constants.hpp"
#include "NAOS/rk4.hpp"
#include "NAOS/rk54.hpp"
#include "NAOS/basicMath.hpp"
#include "NAOS/misc.hpp"
#include "NAOS/springMassEquationsOfMotion.hpp"

namespace naos
{

//! Execute routine to test numerical integrators performance
/*!
 * This routine will perform numerical integration of the autonomous spring mass system and compare
 * the integrated values for each time step with the analytical solution to test the integrity of the
 * numerical integrator.
 */
 template< typename Vector >
 void executeSpringMassIntegration( const double mass,
                                    const double springConstant,
                                    const Vector &initialState,
                                    const double integrationStepSize,
                                    const double startTime,
                                    const double endTime,
                                    std::ostringstream &filePath )
 {
    // open the file to store the analytical and numerically integrated results
    std::ofstream outputFile;
    outputFile.open( filePath.str( ) );
    outputFile << "integrated_distance" << "," << "integrated_velocity" << ",";
    outputFile << "analytical_distance" << "," << "analytical_velocity" << ",";
    outputFile << "stepSize" << "," << "t" << std::endl;

    // time values
    double tCurrent = startTime;
    double tEnd = endTime;
    static double stepSize = integrationStepSize;

    // declare and initialize the result storing vectors
    Vector currentStateVector = initialState;
    Vector nextStateVector = initialState;
    Vector analyticVector = initialState;

    // declare and initialize an object containing the EOMs
    eomSpringMass derivatives( mass, springConstant );

    // store the initial results in the outputFile
    outputFile << currentStateVector[ 0 ] << ",";
    outputFile << currentStateVector[ 1 ] << ",";
    outputFile << analyticVector[ 0 ] << ",";
    outputFile << analyticVector[ 1 ] << ",";
    outputFile << stepSize << ",";
    outputFile << tCurrent << std::endl;

    // calculate the natural frequency of the dynamical system
    const double omega = std::sqrt( springConstant / mass );

    // Start the integration outer loop
    while( tCurrent != tEnd )
    {
        // display the progress bar
        displayProgressBar< double >( tCurrent, tEnd );

        // calculate the new time value
        double tNext = tCurrent + stepSize;

        // check if the next time value is greater than the final time value
        if( tNext > tEnd )
        {
            // make the next time value equal to the final time value and recalculate the step size
            tNext = tEnd;
            stepSize = tEnd - tCurrent;
        }

        // perform the RK4 integration
        // rk4< Vector6, eomSpringMass >( currentStateVector,
        //                                tCurrent,
        //                                stepSize,
        //                                nextStateVector,
        //                                derivatives );

        // Perform the RK5(4) integration
        static bool stepReject = false;
        static double previousError = 10.0e-4;
        rk54GeneralIntegrator< Vector6, eomSpringMass >( currentStateVector,
                                                         tCurrent,
                                                         stepSize,
                                                         nextStateVector,
                                                         derivatives,
                                                         stepReject,
                                                         previousError );

        // get the analytical solution
        analyticVector[ 0 ] = ( initialState[ 1 ] / omega ) * std::sin( omega * tNext )
                                + initialState[ 0 ] * std::cos( omega * tNext );

        analyticVector[ 1 ] = initialState[ 1 ] * std::cos( omega * tNext )
                                - initialState[ 0 ] * omega * std::sin( omega * tNext );

        // store the results
        outputFile << nextStateVector[ 0 ] << ",";
        outputFile << nextStateVector[ 1 ] << ",";
        outputFile << analyticVector[ 0 ] << ",";
        outputFile << analyticVector[ 1 ] << ",";
        outputFile << stepSize << ",";
        outputFile << tNext << std::endl;

        // save the new state values in the vector of current state values. these will be used in
        // the next loop iteration
        tCurrent = tNext;
        currentStateVector = nextStateVector;
    }
    outputFile.close( );
 }

} // namespace naos

#endif // SPRING_MASS_INTEGRATOR_TEST_HPP
