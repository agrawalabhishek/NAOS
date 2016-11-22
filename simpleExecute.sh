#!/bin/bash/
cd build
make
cd bin
./NAOS_main executeParticleAroundUniformlyRotatingEllipsoid
./NAOS_main postAnalysisParticleAroundUniformlyRotatingEllipsoid
cd ../../
