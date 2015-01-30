# suffix rules
.SUFFIXES: .o .F90
.F90.o:
	$(FC) -c $(CPP) $(OPT) $<

# macros    
# available options -DTIMING, -DSTRESS, -DDEBUG

FC = mpif90 
#FC = mpixlf2003_r
OPT = -O2 -g
#CPP = -DSTRESS
CPP =

#OPT = -O3 -qhot=level=0 -qimd=auto
#LDFLAGS = -L/soft/perftools/hpctw -lmpihpm -L/soft/perftools/bgpm/lib -lbgpm 
#BLASLIB = -L/soft/libraries/essl/current/lib64 -lesslbg

EXE = rxmd
SRC = module.F90 rxmd.F90 init.F90 bo.F90 pot.F90 qeq.F90 param.F90 cg.F90 stress.F90
OBJ = module.o init.o bo.o pot.o qeq.o param.o cg.o stress.o

# explicit rules
$(EXE): $(OBJ) $(EXE).o
	$(FC) -o $(EXE) $(OBJ) $(EXE).o
	mv $(EXE) ../rxmd

clean:
	rm -f PI* *.o *.mod *.MOD mpif.h

# implicit rules
rxmd.o: $(OBJ)
init.o: module.o
bo.o: module.o
pot.o: module.o
qeq.o: module.o
cg.o: module.o
	$(FC) -c $(CPP) $(OPT) cg.F90
