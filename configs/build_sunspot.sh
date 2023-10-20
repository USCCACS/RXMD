#!/bin/sh

module load spack
module load frameworks/2023.05.15.001
module load cmake

rm -rfv build
mkdir build
cd build

cmake \
    -DCMAKE_CXX_COMPILER=icpx \
    -DTorch_ROOT=$CONDA_PREFIX/lib/python3.9/site-packages/torch/ \
    -DIPEX_ROOT=$CONDA_PREFIX/lib/python3.9/site-packages/intel_extension_for_pytorch/ \
    -DCMAKE_Fortran_COMPILER=mpif90 -DRXMD_ENABLE_TORCH=ON -DRXMD_ENABLE_IPEX=ON \
    .. && make
