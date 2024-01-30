#!/bin/sh
module purge
module load spack-pe-gcc
module load cmake
module load PrgEnv-gnu
module load frameworks/2023.10.15.001
 
rm -rfv build
mkdir build
cd build

basedir=/lus/gecko/projects/NAQMC_RMD_aesp_CNDA/knomura/libtorch/share/cmake/
cmake \
	-DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
	-DIPEX_DIR=${basedir}/IPEX -DTorch_DIR=${basedir}/Torch \
	.. && make

#cmake \
#	-DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx \
#	-DCMAKE_CXX_FLAGS="-stdlib=libc++" \
#	-DIPEX_DIR=${basedir}/IPEX -DTorch_DIR=${basedir}/Torch \
#	.. && make
