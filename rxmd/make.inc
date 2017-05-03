# FC is the MPI wrapper for Fortran
# FFLAG is the flag for FC
# LINKER, optional, linker
# LIBS, optional, extra libraries needed

MPIF90 = ftn
FC = ifort

MPIF90_FLAGS = -qopenmp -O3 -g
FC_FLAGS = -O3 -g

# LINKER = 
# LIBS =
