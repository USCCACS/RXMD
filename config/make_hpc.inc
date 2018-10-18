# MPIF90 is the MPI wrapper for Fortran
# FFLAG is the flag for FC
# LINKER, optional, linker
# LIBS, optional, extra libraries needed

MPIF90 = mpif90
FC = ifort

# gfortran
MPIF90_FLAGS = -O3 -ffast-math
#MPIF90_FLAGS = -fopenmp -O3 -ffast-math
#MPIF90_FLAGS = -fcheck=all -g

FC_FLAGS = -O3 
#FC_FLAGS = -O3 -g

# LINKER = 
# LIBS =
