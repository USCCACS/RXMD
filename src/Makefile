include ../make.inc

# suffix rules
.SUFFIXES: .o .F90 .cpp
.F90.o:
	$(FC) $(FC_FLAGS) -c $<

%.o: %.mod

SRCS = module.F90 utils.F90 memory.F90 velocity.F90 pot.F90 fileio.F90 comm.F90 init.F90 bo.F90 qeq.F90 reaxff.F90 main.F90 pqeq.F90 cmdline.F90 lists.F90 msd.F90 ensembles.F90 constant.F90 remd_pot_ANN_original.F90 fnn.F90 shortrep.F90 symmfunc.F90 rxmdnn.F90 nnmm.F90 stepmod.F90 harmopot.F90
OBJS = $(SRCS:.F90=.o)


#CARC
#install_root=/scratch1/knomura/torchapp/libtorch
#TORCHLIB_DIRS=-L../rxmdtorch/build -L${install_root}/lib -Wl,-R${install_root}/lib
#TORCHLIBS=-ltorch_cpu -ltorch_cuda -lc10 -lrxmdtorch

#Polaris
#install_root=/soft/datascience/conda/2023-10-04/pytorch/torch
#TORCHLIB_DIRS=-L../rxmdtorch/build -L${install_root}/lib -Wl,-R${install_root}/lib -Wl,-R/opt/cray/pe/gcc/11.2.0/snos/lib64
#TORCHLIBS=-ltorch_cpu -ltorch_cuda -lc10 -lrxmdtorch

#Sunspot
install_root=/lus/gila/projects/NAQMC_RMD_aesp_CNDA/knomura/libs/libtorch
TORCHLIB_DIRS=-L../rxmdtorch/build -L${install_root}/lib -Wl,-R${install_root}/lib -L${install_root}/lib -Wl,-R${install_root}/lib
TORCHLIBS=-ltorch_cpu -ltorch  -lc10 -lrxmdtorch -lintel-ext-pt-gpu -lintel-ext-pt-cpu

#TORCHLIBS=-ltorch_cpu -ltorch  -lc10 -lrxmdtorch -lintel-ext-pt-gpu -lintel-ext-pt-cpu 

INCLUDE_DIRS=-I${install_root}/include -I${install_root}/include/torch/csrc/api/include 

MACROS=


EXE := ../rxmd

$(EXE): $(OBJS)
	$(MPIF90) $(MPIF90_FLAGS) -o $(EXE) $(OBJS) $(LIBS) 

rxmdnn: $(OBJS)
	icpc -c ../nn/nn.cpp 
	$(MPIF90) $(MPIF90_FLAGS) -o $(EXE) $(OBJS) $(LIBS) -lstdc++ nn.o

rxmdtorch: $(OBJS) 
	make -C ../rxmdtorch/build
	$(MPIF90) $(MPIF90_FLAGS) $(TORCHLIB_DIRS) $(MACROS) -o $(EXE) $(OBJS) $(LIBS) $(TORCHLIBS) -lstdc++

#cg.o : fileio.o
#main.F90 : cg.o

nnmm.o : module.o fileio.o utils.o memory.o stepmod.o

stepmod.o : utils.o module.o

harmopot.o : utils.o module.o

shortrep.o : module.o fileio.o utils.o

remd_pot_ANN_original.o : constant.o symmfunc.o

module.o : utils.o memory.o remd_pot_ANN_original.o

fileio.o : module.o harmopot.o
comm.o : module.o
velocity.o : module.o

lists.o : module.o utils.o
fnn.o : module.o comm.o lists.o fileio.o shortrep.o nnmm.o
rxmdnn.o : module.o cmdline.o utils.o fileio.o comm.o velocity.o msd.o qeq.o

cmdline.o : utils.o module.o memory.o nnmm.o
msd.o : utils.o module.o
reaxff.o : cmdline.o utils.o fileio.o comm.o velocity.o msd.o
ensembles.o : module.o

bo.o : reaxff.o
pot.o : reaxff.o bo.o comm.o pqeq.o qeq.o harmopot.o
init.o : reaxff.o fileio.o pot.o fnn.o shortrep.o rxmdnn.o nnmm.o
qeq.o : reaxff.o comm.o lists.o
pqeq.o : reaxff.o comm.o lists.o

main.o : init.o module.o


## build No-MPI executable

nompi: MPIF90_FLAGS := $(FC_FLAGS) -DNOMPI
nompi: MPIF90 := $(FC)

nompi_o: nompi.o
	$(MPIF90) $(MPIF90_FLAGS) -c nompi.F90

nompi: nompi_o $(OBJS) 
	$(MPIF90) $(MPIF90_FLAGS) -o $(EXE) $(OBJS) nompi.o $(LIBS)


## build librxmd

librxmd: MPIF90_FLAGS := $(FC_FLAGS) -DNOMPI -DLIBRXMD -fpic
librxmd: MPIF90 := $(FC)

librxmd: nompi_o $(OBJS)
	$(FC) *.o  -shared -fpic -o librxmd.so



clean:
	rm -f PI* *.o *.mod *.MOD mpif.h *.a *.so ../rxmd

rebuild:
	make clean
	make -j
