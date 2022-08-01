#!/bin/sh

module load PrgEnv-gnu
module unload craype-accel-nvidia80
module load cudatoolkit-standalone/11.6.2

rm -rfv build
mkdir build
cd build

cmake -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_PREFIX_PATH="/soft/datascience/conda/2022-07-19/pytorch/torch;/soft/datascience/cuda/cudnn-11.6-linux-x64-v8.4.0.27" \
  .. && make
