# suitable for Linux Computer with Intel Compiler
# MPIF90 is a MPI wrapper for Fortran and is used to compile reaxFF code in src directory
# FC is a fortran compiler and is used to compite geninit.F90 code in init directory so as to create the input configuration for reaxFF simulation
# MPIF90_FLAGS is the flag for MPIF90
# FC_FLAGS is the  flag for FC
# LINKER is linker (optional)
# LIBS are extra needed libraries (optional) 

#MPIF90 = mpiifort -check all -traceback
MPIF90 = mpiifort -DRXMDNN 
FC = ifx

MPIF90_FLAGS = -O3

FC_FLAGS = -O3 

# LINKER = 
# LIBS =
