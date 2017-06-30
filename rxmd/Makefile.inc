#compiler settings

# Intel Compiler
#FC = mpif90 -check all -traceback
#FC = mpif90 -qopenmp -O3
#FC = mpif90 

# gfortran
FC = mpif90 -fopenmp -O3 -ffast-math
#FC = mpif90 -fopenmp -O3 -ffast-math -DQEQDUMP -DRFDUMP
#FC = mpif90 -O3 -ffast-math -DQEQDUMP -DRFDUMP
#FC = mpif90 -O3 -ffast-math 
#FC = mpif90 -fcheck=all -g

# xl fortran
#FC = mpif90 -O3 -qhot -qsmp=omp
#FC = mpif90 -O3 -qhot
#FC = mpif90 -C -g
