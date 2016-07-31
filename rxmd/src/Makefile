# suffix rules
.SUFFIXES: .o .F90 .cpp
.F90.o: 
	$(FC) -c $(CPP) $(OPT) $<
.cpp.o: 
	$(CC) -c $(CPP) $(OPT) $<
%.o: %.mod

#FC = mpif90 -check all -traceback
FC = mpif90 -openmp
CC = mpic++
OPT = -O3 -g

SRCS := module.F90 main.F90 pot.F90 fileio.F90 comm.F90 init.F90 bo.F90 qeq.F90 cg.F90 param.F90 stress.F90
OBJS := $(SRCS:.F90=.o) 
CPP = 

#OPT = -O3 -qhot=level=0 -qimd=auto
#LDFLAGS = -L/soft/perftools/hpctw -lmpihpm -L/soft/perftools/bgpm/lib -lbgpm 

EXE = rxmd

ifeq ($(EXE),rxmd-cpp)
  CPPSRCS := Basic.cpp cpp_interface.cpp
  CPPOBJS := $(CPPSRCS:.cpp=.o)
  LDFLAGS = -L/usr/lib64/ -lstdc++
  SRCS += $(CPPSRCS)
  OBJS += $(CPPOBJS)
  CPP += -DINTEROP
endif

$(EXE): $(OBJS) 
	$(FC) -o $(EXE) $(LDFLAGS) $(OBJS) 
	mv $(EXE) ../$(EXE)

$(OBJS) : atoms.mod

atoms.mod: module.F90
	$(FC) -c $(CPP) $(OPT) module.F90


clean:
	rm -f PI* *.o *.mod *.MOD mpif.h
