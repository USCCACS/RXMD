#FC=gfortran
#FFLAGS=-fcheck=all -fbacktrace 
FC=ifort
FFLAGS=-O3

all:
	${FC} ${FFLAGS} sum.f90 

clean:
	rm -f *.o *.mod a.out *.dat *.xyz *.png
