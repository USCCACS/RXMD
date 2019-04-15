include ../make.inc

# suffix rules
.SUFFIXES: .o .F90 .cpp
.F90.o:
	$(MPIF90) $(MPIF90_FLAGS) $(LIBS) -c $<

%.o: %.mod

SRCS := utils.F90 memory.F90 velocity.F90 pot.F90 fileio.F90 comm.F90 init.F90 bo.F90 qeq.F90 param.F90 main.F90 pqeq.F90 cmdline.F90 fnn.F90
OBJS := $(SRCS:.F90=.o)

EXE := ../rxmd

$(EXE): $(OBJS)
	$(MPIF90) $(MPIF90_FLAGS) -o $(EXE) $(OBJS) module.o $(LIBS)

init.o : fileio.o param.o pot.o fnn.o
cg.o : fileio.o
main.F90 : cg.o

module.o : utils.o memory.o

cmdline.o : module.o memory.o fnn.o
param.o : cmdline.o utils.o

fileio.o : param.o
bo.o : param.o
pot.o : param.o bo.o comm.o pqeq.o qeq.o
qeq.o : param.o comm.o
pqeq.o : param.o comm.o
velocity.o : param.o
comm.o : param.o

main.o : init.o


nompi: nompi.f90
	$(FC) $(FC_FLAGS) -DNOMPI -c nompi.f90 utils.F90
	make MPIF90=$(FC) MPIF90_FLAGS="$(FC_FLAGS) -DNOMPI" OBJS="$(OBJS) nompi.o" $(EXE) 

clean:
	rm -f PI* *.o *.mod *.MOD mpif.h ../rxmd

rebuild:
	make clean
	make -j
