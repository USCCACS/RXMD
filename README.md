# rxmd
ReaxFF MD code repository 

# Currently tested env
source /usr/usc/intel/14.0/setup.sh

source /usr/usc/openmpi/1.8.7/setup.sh.intel


# multithreaded non-bonding pair 
this version requires OpenMP enabled as default. If you are going to use Infiniband nodes in priya queue,

example)
- run a 16 MPIrank rxmd job on two IB nodes

qsub -I -d . -l nodes=2:ppn=16:priya_IB,walltime=2:00:00

mpirun -x OMP_NUM_THREADS=2 --bind-to none -npernode 8 -np 16 ./rxmd
