#!/bin/sh

module load spack-pe-base
module load cmake
module load cudatoolkit-standalone cudnn cuda-PrgEnv-nvidia

rm -rfv build
mkdir build
cd build

cmake -DCMAKE_CXX_COMPILER=gcc-12 \
  -DCMAKE_PREFIX_PATH="/lus/eagle/projects/NAQMC_RMD_aesp/knomura/libtorch" \
  .. && make
