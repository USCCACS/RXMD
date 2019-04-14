include ../make.inc

# suffix rules
.SUFFIXES: .o .F90 .cpp
.F90.o:
	$(MPIF90) $(MPIF90_FLAGS) $(LIBS) -c $<

%.o: %.mod

SRCS := velocity.F90 pot.F90 fileio.F90 comm.F90 init.F90 bo.F90 qeq.F90 param.F90 main.F90 pqeq.F90 cmdline.F90
OBJS := $(SRCS:.F90=.o)

EXE := ../rxmd

$(EXE): $(OBJS)
	$(MPIF90) $(MPIF90_FLAGS) -o $(EXE) $(OBJS) module.o $(LIBS)

$(OBJS) : base.o

base.o: module.F90
	$(MPIF90) $(MPIF90_FLAGS) -c module.F90


init.o : fileio.o param.o
cg.o : fileio.o
main.F90 : cg.o

cmdline.o : base.o
param.o : cmdline.o

fileio.o : param.o
bo.o : param.o
pot.o : param.o
qeq.o : param.o
pqeq.o : param.o
velocity.o : param.o

main.o : init.o


nompi: nompi.f90
	$(FC) $(FC_FLAGS) -DNOMPI -c nompi.f90 module.F90
	make MPIF90=$(FC) MPIF90_FLAGS="$(FC_FLAGS) -DNOMPI" OBJS="$(OBJS) nompi.o" $(EXE) 

clean:
	rm -f PI* *.o *.mod *.MOD mpif.h ../rxmd

rebuild:
	make clean
	make -j
