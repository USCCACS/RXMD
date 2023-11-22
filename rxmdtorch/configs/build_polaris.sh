#!/bin/sh

module unload craype-accel-nvidia80
module load cudatoolkit-standalone/11.6.2
module load gcc/11.2.0

rm -rfv build
mkdir build
cd build

cmake -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_PREFIX_PATH="/soft/datascience/conda/2023-10-04/pytorch/torch;/soft/libraries/cudnn/cudnn-11-linux-x64-v8.7.0.84/" \
  .. && make
