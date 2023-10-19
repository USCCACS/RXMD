#!/bin/sh

module load spack
module load cmake
module load oneapi/eng-compiler/2023.05.15.007
 
rm -rfv build
mkdir build
cd build

cmake \
    -DCMAKE_CXX_COMPILER=icpx \
    -DTorch_ROOT=/soft/datascience/aurora_models_frameworks-2023.2/lib/python3.9/site-packages/torch/ \
    -DIPEX_ROOT=/soft/datascience/aurora_models_frameworks-2023.2/lib/python3.9/site-packages/intel_extension_for_pytorch/ \
    -DCMAKE_Fortran_COMPILER=mpif90 -DRXMD_ENABLE_TORCH=ON -DRXMD_ENABLE_IPEX=ON \
    .. && make
