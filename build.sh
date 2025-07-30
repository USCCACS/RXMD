#!/bin/sh

cmake_prefix_path=`python -c 'import torch; print(torch.utils.cmake_prefix_path)'`

#Polaris Config
sitepackage_path=`python -c 'import site; print(*site.getsitepackages())'`
LD_LIBRARY_PATH=${sitepackage_path}/nvidia/nccl/lib:${LD_LIBRARY_PATH}

##Aurora Config
#cmake_prefix_path=/lus/flare/projects/NAQMC_RMD_aesp_CNDA/knomura/libtorch-2.6.10/share/cmake

##CARC
#module purge
#module load usc openblas gcc cuda cmake openmpi cudnn intel-oneapi-mkl
#module unload python


#Common Part
echo $cmake_prefix_path
echo $nccl_path

rm -rf build && mkdir build && cd build

cmake .. -DCMAKE_PREFIX_PATH=${cmake_prefix_path} 
make -j && make install

