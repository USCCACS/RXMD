#!/bin/sh
module load cmake
module load oneapi/release/2024.2.1
 
rm -rfv build
mkdir build
cd build

basedir=/lus/flare/projects/NAQMC_RMD_aesp_CNDA/knomura/libtorch/share/cmake/
cmake \
	-DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
	-DIPEX_DIR=${basedir}/IPEX -DTorch_DIR=${basedir}/Torch \
	.. && make
