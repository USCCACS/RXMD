#FC=gfortran
#FFLAGS=-fcheck=all -fbacktrace 
#FC=ifort
#FC=mpiifort -qopenmp -check all -traceback
FC=mpiifort  -qopenmp 
#FC=mpiifort -check all -traceback
#FFLAGS=-check all -traceback -check noarg_temp_created
FFLAGS=-O3

all:
	${FC} ${FFLAGS} stat.f90 voxels.f90 main.f90 

clean:
	rm -f *.o *.mod a.out *.dat *.xyz *.png
