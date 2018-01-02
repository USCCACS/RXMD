include ../Makefile.inc

$(EXE): $(OBJS)
	$(FC) -o $(EXE) $(LDFLAGS) $(OBJS) module.o
	mv $(EXE) ../$(EXE)

$(OBJS) : base.mod

base.mod: module.F90
	$(FC) -c $(CPP) $(OPT) module.F90

main.F90 : cg.mod

cg.mod : cg.o

fileio.o : param.o

clean:
	rm -f PI* *.o *.mod *.MOD mpif.h ../rxmd

rebuild:
	make clean
	make -j 4
