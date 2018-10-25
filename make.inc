# An example make.inc file
# MPIF90 is the MPI wrapper for Fortran and it is used to compile reaxFF code in src directory
# FC is fortran compiler and it is used to compite geninit.F90 code in init directory so as to create input configuration for reaxFF simulation
# MPIF90_FLAGS is flag for MPIF90
# FC_FLAGS is the flag for FC
# LINKER, optional, linker
# LIBS, optional, extra libraries needed

MPIF90 = mpif90
FC = gfortran

MPIF90_FLAGS = -O3 -qhot
#MPIF90_FLAGS = -qopenmp -O3 -xMIC-AVX512 -g  -heap-arrays 0
#MPIF90_FLAGS = -qopenmp -O3 -g -qopt-report=3  -heap-arrays 0
#MPIF90_FLAGS = -qopenmp -heap-arrays 0 -check bounds -traceback

FC_FLAGS = -O3 -g

# LINKER = 
# LIBS =
