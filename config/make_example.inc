# MPIF90 is the MPI wrapper for Fortran
# FFLAG is the flag for FC
# LINKER, optional, linker
# LIBS, optional, extra libraries needed

MPIF90 = mpif90
#FC = gfortran
FC = xlf

MPIF90_FLAGS = -O3 -qhot
#MPIF90_FLAGS = -qopenmp -O3 -xMIC-AVX512 -g  -heap-arrays 0
#MPIF90_FLAGS = -qopenmp -O3 -g -qopt-report=3  -heap-arrays 0
#MPIF90_FLAGS = -qopenmp -heap-arrays 0 -check bounds -traceback

FC_FLAGS = -O3 -g

# LINKER = 
# LIBS =
