#!/bin/sh

module load spack
module load cmake
module load PrgEnv-gnu
module load frameworks/2022.12.30.001
 
rm -rfv build
mkdir build
cd build

#sitepkg=/soft/datascience/aurora_models_frameworks-2023.0/lib/python3.9/site-packages
sitepkg=/home/knomura/ipex

cmake \
	-DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx \
	-DCMAKE_PREFIX_PATH=${sitepkg}/torch/share/cmake/Torch \
	-DINTEL_EXTENSION_FOR_PYTORCH_PATH=${sitepkg}/intel_extension_for_pytorch \
	.. && make
