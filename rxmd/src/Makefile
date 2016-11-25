# suffix rules
.SUFFIXES: .o .F90 .cpp
.F90.o: 
	$(FC) -c $(CPP) $(OPT) $<

%.o: %.mod

# Intel Compiler
#FC = mpif90 -check all -traceback
#FC = mpif90 -openmp
FC = mpif90 

# gfortran
#FC = mpif90 -fopenmp
#FC = mpif90 -fcheck=all

# xl fortran
#FC = mpif90 -O3 -qreport -qhot
#FC = mpif90 -C -g

SRCS := main.F90 pot.F90 fileio.F90 comm.F90 init.F90 bo.F90 qeq.F90 cg.F90 param.F90 stress.F90
OBJS := $(SRCS:.F90=.o) 
CPP = 

#OPT = -O3 -qhot=level=0 -qimd=auto
#LDFLAGS = -L/soft/perftools/hpctw -lmpihpm -L/soft/perftools/bgpm/lib -lbgpm 

EXE = rxmd

$(EXE): $(OBJS)
	$(FC) -o $(EXE) $(LDFLAGS) $(OBJS) module.o
	mv $(EXE) ../$(EXE)

$(OBJS) : base.mod

base.mod: module.F90
	$(FC) -c $(CPP) $(OPT) module.F90
clean:
	rm -f PI* *.o *.mod *.MOD mpif.h
