#!/bin/sh

module load spack
module load cmake
module load PrgEnv-gnu
module load frameworks/2023.05.15.001
 
rm -rfv build
mkdir build
cd build

basedir=/lus/gila/projects/NAQMC_RMD_aesp_CNDA/knomura/libs/libtorch/share/cmake/
cmake -Wno-dev \
	-DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx \
	-DIPEX_DIR=${basedir}/IPEX -DTorch_DIR=${basedir}/Torch \
	.. && make
