#!/bin/sh

module load spack
module load cmake
module load PrgEnv-gnu
#module load frameworks/2022.12.30.001
module load frameworks/2023.05.15.001
 
rm -rfv build
mkdir build
cd build

cmake \
	-DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx \
	-DCMAKE_PREFIX_PATH="/soft/datascience/aurora_models_frameworks-2023.1/lib/python3.9/site-packages/torch/share/cmake/Torch" \
	-DINTEL_EXTENSION_FOR_PYTORCH_PATH="/soft/datascience/aurora_models_frameworks-2023.1/lib/python3.9/site-packages/intel_extension_for_pytorch" \
	.. && make
