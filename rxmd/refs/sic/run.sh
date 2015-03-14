#!/bin/sh

source /usr/usc/intel/14.0/setup.sh
source /usr/usc/openmpi/1.8.4/setup.sh.intel

cd init/
ifort geninit.F90 -o geninit
./geninit sicnp_o_norm.xyz 
cp rxff000000 ../DAT/
cd ..
cp ../../../rxmd . 
mpirun -np 1 ./rxmd | tee std.out
