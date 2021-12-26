CXX=icpc
FXX=ifort
MPICXX=mpiicpc
MPIF90=mpiifort
#OMPFLAGS=-fiopenmp -fopenmp-targets=spir64 -D__STRICT_ANSI__ 
#OMPFLAGS=-fiopenmp -fopenmp-targets=spir64 
OMPFLAGS=

# hybrid build
all:
	$(CXX) $(OMPFLAGS) -std=c++17 -c nn.cpp
	$(FXX) -c rxmdnn.f90
	$(FXX) $(OMPFLAGS) -lstdc++ -cxxlib rxmdnn.o nn.o 

gxx:
	g++ -c nn.cpp
	gfortran -c rxmdnn.f90
	gfortran -lstdc++ rxmdnn.o nn.o

ixx:
	icpc -c nn.cpp
	ifort -c rxmdnn.f90
	ifort -lstdc++ rxmdnn.o nn.o

cxx:
	$(CXX) $(OMPFLAGS) nn.cpp


clean:
	rm -fv *.o *.mod a.out

# C++ code test
#icpx -std=c++17 -g -fiopenmp -fopenmp-targets=spir64 -D__STRICT_ANSI__ nn.cpp && iprof ./a.out pto.xyz params.in
#icpx -std=c++17 -g nn.cpp && ./a.out pto.xyz params.in 


