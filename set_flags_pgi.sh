#!/bin/bash

module purge
module load pgi-nvhpc

export OMP_TARGET_OFFLOAD=MANDATORY
export CXX='nvc++'
export FC='nvfortran'
#export COMMON_FLAGS='-fopenmp'
export CXXFLAGS="-mp -target=gpu -gpu=cc60,cuda10.1"
export FFLAGS="-mp -target=gpu -gpu=cc60,cuda10.1"

