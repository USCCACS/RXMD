# An example make.inc file
# MPIF90 is a MPI wrapper for Fortran and is used to compile reaxFF code in src directory
# MPIF90_FLAGS is the flag for MPIF90

MPIF90 = mpif90
MPIF90_FLAGS = -O3

# FC is a fortran compiler and is used to compite geninit.F90 code in init directory so as to create the input configuration for reaxFF simulation
# FC_FLAGS is the  flag for FC

FC = gfortran
FC_FLAGS = -O3

# LINKER is linker (optional)
# LIBS are extra needed libraries (optional)

# LINKER = 
# LIBS =
