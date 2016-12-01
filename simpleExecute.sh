#!/bin/bash/
cd build
make
cd bin
./NAOS_main executeRegolithTrajectoryCalculation
cd ../../
