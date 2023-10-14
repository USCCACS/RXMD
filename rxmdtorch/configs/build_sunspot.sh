#!/bin/sh

module use /soft/modulefiles/
module load spack
module load cmake
module load PrgEnv-gnu
module load frameworks/2023.05.15.001
 
rm -rfv build
mkdir build
cd build

cmake \
        -DCMAKE_CXX_COMPILER=icpx \
	-DTorch_DIR=/soft/datascience/aurora_models_frameworks-2023.2/lib/python3.9/site-packages/torch/share/cmake/Torch/ \
	.. && make
