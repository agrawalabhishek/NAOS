#!/bin/bash/
cd build
make
cd bin
./NAOS_main executeRestricted2BP
./NAOS_main postAnalysisRestricted2BP
cd ../../
