# suitable for Linux Computer with Intel Compiler
# MPIF90 is a MPI wrapper for Fortran and is used to compile reaxFF code in src directory
# FC is a fortran compiler and is used to compite geninit.F90 code in init directory so as to create the input configuration for reaxFF simulation
# MPIF90_FLAGS is the flag for MPIF90
# FC_FLAGS is the  flag for FC
# LINKER is linker (optional)
# LIBS are extra needed libraries (optional) 

# Sophia
MPIF90 = mpifort
FC = mpifort
MPIF90_FLAGS = -C -O2 -fallow-argument-mismatch -mcmodel=medium
FC_FLAGS = -C -O2 -fallow-argument-mismatch -mcmodel=medium -DRXMDNN 

# CARC
#MPIF90 = mpiifort
#FC = ifx
#MPIF90_FLAGS = -O2
#FC_FLAGS = -O2 -C -DRXMDNN

# Sunspot/Aurora
#MPIF90 = mpifort
#FC = mpifort
#MPIF90_FLAGS = -O2 -traceback
#FC_FLAGS = -O2 -traceback -DRXMDNN

# LINKER = 
# LIBS =
