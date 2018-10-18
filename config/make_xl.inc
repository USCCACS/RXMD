# MPIF90 is the MPI wrapper for Fortran
# FFLAG is the flag for FC
# LINKER, optional, linker
# LIBS, optional, extra libraries needed

MPIF90 = mpif90
FC = xlf

# xl fortran
MPIF90_FLAGS = -O3 -qhot
#MPIF90_FLAGS = -O3 -qhot -qsmp=omp
#MPIF90_FLAGS = -C -g

FC_FLAGS = -O3 
#FC_FLAGS = -O3 -g

# LINKER = 
# LIBS =
