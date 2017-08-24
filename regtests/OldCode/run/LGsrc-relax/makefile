# suffix rules
.SUFFIXES: .o .f90
.f90.o:
	$(FC) -c $(CPP) $(OPT) $<

# macros    
# available options -DTIMING, -DSTRESS, -DDEBUG
#FC = mpif90 -check all -traceback
FC = mpif90 
OPT = -ffast-math -O3
#OPT = -O3 -static
CPP = -cpp -DRXFFLG 
##CPP = -cpp -DSTRESS -DRXFFLG 
#CPP = -Mpreprocess
#CPP = -Mpreprocess -DSTRESS
#OPT = -prof_gen
#OPT = -O0 -C

EXE = rxmd
SRC = module.f90 rxmd.f90 init.f90 bo.f90 pot.f90 qeq.f90 param.f90 cg.f90 stress.f90
OBJ = module.o init.o bo.o pot.o qeq.o param.o cg.o  stress.o

# explicit rules
$(EXE): $(OBJ) $(EXE).o
	$(FC) -o $(EXE) $(OBJ) $(EXE).o
	mv $(EXE) ../rxmd_old

clean:
	rm -f PI* *.o *.mod *.MOD mpif.h


# implicit rules
rxmd.o: $(OBJ)
init.o: module.o
bo.o: module.o
pot.o: module.o
qeq.o: module.o
cg.o: module.o
	$(FC) -c $(CPP) $(OPT) cg.f90
ensemble.o: module.o
	$(FC) -c $(CPP) $(OPT) ensemble.f90
