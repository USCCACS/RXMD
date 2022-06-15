#module purge  
#module load usc intel-oneapi gcc/8.3.0 cuda cudnn

install_root=/scratch1/knomura/torchapp/libtorch

FC=ifx

LIB_DIRS=-L../build -L${install_root}/lib -Wl,-R${install_root}/lib 
INCLUDE_DIRS=-I${install_root}/include -I${install_root}/include/torch/csrc/api/include 

LIBS=-lstdc++ -ltorch -ltorch_cpu -ltorch_cuda -lc10 -lrxmdtorch
MACROS=-D_GLIBCXX_USE_CXX11_ABI=0
OBJS=rxmdnn.o

all:
	${FC} -c rxmdnn.f90
	${FC} ${MACROS} ${INCLUDE_DIRS} ${LIB_DIRS} ${OBJS} ${LIBS} -o exe 

clean:
	rm -f *.mod *.o exe
