# suffix rules
.SUFFIXES: .o .F90 .cpp
.F90.o: 
	$(FC) -c $(CPP) $(OPT) $<
.cpp.o: 
	$(CC) -c $(CPP) $(OPT) $<
%.o: %.mod

FC = mpif90 
CC = mpic++
OPT = -O3 -g
#OPT = -check all -traceback
#CPP = -DSTRESS
CPP = -DINTEROP
#CPP = 

#OPT = -O3 -qhot=level=0 -qimd=auto
#LDFLAGS = -L/soft/perftools/hpctw -lmpihpm -L/soft/perftools/bgpm/lib -lbgpm 
#BLASLIB = -L/soft/libraries/essl/current/lib64 -lesslbg

EXE = rxmd
SRC = module.F90 main.F90 init.F90 bo.F90 pot.F90 qeq.F90 param.F90 cg.F90 stress.F90
OBJ = $(SRC:.F90=.o) 

$(EXE): $(OBJ) 
	$(FC) -o $(EXE) $(OBJ) 
	mv $(EXE) ../rxmd

interop: $(OBJ) cpp_interface.o Basic.o
	$(FC) -o $(EXE) $(OBJ) cpp_interface.o Basic.o -lstdc++
	mv $(EXE) ../rxmd

clean:
	rm -f PI* *.o *.mod *.MOD mpif.h

# implicit rules
init.o: module.o
bo.o: module.o
pot.o: module.o
qeq.o: module.o
param.o: module.o
