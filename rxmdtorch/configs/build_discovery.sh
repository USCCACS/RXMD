#!/bin/sh

module purge
module load usc
module load gcc/13.3.0 cuda/12.6.3 cudnn/8.9.7.29-12-cuda
module load cmake

rm -rfv build
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=/scratch1/knomura/libtorch .. && make
