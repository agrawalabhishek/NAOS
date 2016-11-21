#!/bin/bash/
cd build
make
cd bin
./NAOS_main executeParticleAroundSpheroid
./NAOS_main postAnalysisParticleAroundSpheroid
cd ../../
