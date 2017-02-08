#!/bin/bash/
cd build
make
cd bin
./NAOS_main executeRegolithMonteCarlo
# ./NAOS_main testPerturbations
# ./NAOS_main executeSunAsteroidTwoBodyProblem
cd ../../
