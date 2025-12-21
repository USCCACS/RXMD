#!/bin/sh

cmake_prefix_path=`python -c 'import torch; print(torch.utils.cmake_prefix_path)'`

#Polaris Config
#sitepackage_path=`python -c 'import site; print(*site.getsitepackages())'`
#LD_LIBRARY_PATH=${sitepackage_path}/nvidia/nccl/lib:${LD_LIBRARY_PATH}

##CARC
#module purge
#module load usc openblas gcc cuda cmake openmpi cudnn intel-oneapi-mkl
#module unload python

#Common Part
echo $cmake_prefix_path
echo $nccl_path

rm -rf build && mkdir build && cd build

#Aurora
cmake .. -DCMAKE_PREFIX_PATH=${cmake_prefix_path} -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DMPI_Fortran_COMPILER=mpifort

make -j && make install

