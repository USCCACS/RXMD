# suffix rules
.SUFFIXES: .o .F90 .cpp
.F90.o: 
	$(FC) -c $(CPP) $(OPT) $<

%.o: %.mod

include ../../Makefile.inc

SRCS := cg.F90 pot.F90 fileio.F90 comm.F90 init.F90 bo.F90 qeq.F90 param.F90 stress.F90 main.F90
OBJS := $(SRCS:.F90=.o) 
CPP = 

EXE = rxmd

BUILDDIR=../../src/

$(EXE): 
	@rm -f *.o *.mod
	@cp -v ${BUILDDIR}/*.mod ${BUILDDIR}/*.o . 
	$(FC) -c $(CPP) $(OPT) pot.F90
	$(FC) -o $(EXE) $(LDFLAGS) $(OBJS) module.o
#	mv $(EXE) ../$(EXE)

$(OBJS) : base.mod

base.mod: module.F90
	$(FC) -c $(CPP) $(OPT) module.F90

main.F90 : cg.mod

cg.mod : cg.o

fileio.o : param.o

clean:
	rm -f PI* *.o *.mod *.MOD mpif.h rxmd
