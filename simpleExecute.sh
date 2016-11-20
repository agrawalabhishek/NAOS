#!/bin/bash/
cd build
make
cd bin
./NAOS_main executePointMassGravityOrbiter
./NAOS_main postAnalysis
cd ../../
