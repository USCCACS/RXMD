# rxmd : Linear Scalable Parallel ReaxFF Molecular Dynamics Simulator

**rxmd** has been developed to simulate large-scale Reactive Force Field molecular dynamics (MD) simulations on from commodity laptops to high-end supercomputing platforms. **rxmd** has been used in a various class of material studies, such as shock-induced chemical reactions, stress corrosion cracking, underwater bubble collapse, fracture of self-healing ceramics and oxidation of nanoparticles. 

## Prerequisites

**rxmd** is designed to be simple, portable and minimally dependent on 3rd party library. You will need only a Fortran compiler that supports OpenMP and MPI (Message Passing Interface) library to compile the code and run it. Modern Fortran compilers natively support OpenMP, and you can find many freely distributed MPI libraries online. Please refer to MPI library developer site about how to install their library. 

**rxmd** has been tested on following environments.

### - Fortan Compiler: 
```
GNU Fortran (GCC) 6.1.0
Intel Fortran (IFORT) 17.0.4
IBM XL Fortran V14.1
```
### - MPI library: 
```
OpenMPI 1.8.8
MPICH2
MVAPICH2 
Cray Mpich 7.6.0
```

## Getting Started

To get started,  clone this repository to your computer. 
```
~$ git clone https://github.com/USCCACS/rxmd.git
```

## How to Compile

Frist, change working directory to **rxmd/**
```
~$ cd rxmd
```
you will see following files and directories.

```
rxmd $ ls
DAT/          conf/         ffield        regtests/     src/          util/
Makefile.inc  doc/          init/         rxmd.in       unittests/
```

Among them, two directories, **src/** and **init/**, are especially important to get started here. **src/** contains all source codes and **init/** has program and input files to generate an initial configurations for MD simulation. 

To build the **rxmd** executable, first we need to make sure what compiler we want to use. There are two files **Makefile.inc** and **init/Makefile** that you would need to modify according to your computing environment, i.e. compiler, MPI library etc. 

- **Makefile.inc** is to specify what compiler you like to use to build the executable. We have several predefined compiler settings in the file. Please select and uncomment the macro **FC** you want to use. 

Example) mpif90 with gfortran optimization flags. 
```
# gfortran
FC = mpif90 -fopenmp -O3 -ffast-math
```

- **init/Makefile** defines how to build a standalone software to generate intial configuration. Any Fortran or MPI compiler that supports [the stream I/O](https://docs.oracle.com/cd/E19205-01/819-5262/aeuca/index.html) can be used here. 

Example)
```
# macros    
FC = mpif90
```

Now we have the compiler setting done! Next step is to generate initial MD configuration and the **rxmd** executable. From the working directory type the make command below, which compiles the standalone application **geninit** to read a geometry file (init.xyz by default) in the directory, replicate it to cover the entire initial MD geometry (rxff.bin), and place it in **DAT/** directory. 

```
rxmd $ make -C init/
```

Then type the command below to compile the **rxmd** executable.

```
rxmd $ make -C src/
```

Make sure you have **rxmd** and **DAT/rxff.bin** in place, then you are ready to run a simulation.

```
rxmd $ ls
DAT/          conf/         ffield        regtests/     rxmd.in       unittests/
Makefile.inc  doc/          init/         rxmd*         src/          util/
```
```
rxmd $ ls DAT/
rxff.bin
```

## How to run rxmd

Default input parameters are set to run a single process job. Type
```
rxmd $ ./rxmd
```

How to run a multi process job depends on which MPI library you use, but most likely **mpirun** just works for you. 

```
rxmd $ mpirun -np nprocessors ./rxmd
```

If you see following outputs, congratulations! You have everything working.

```
rxmd $ ./rxmd 
              rxmd has started
----------------------------------------------------------------
         req/alloc # of procs:        1  /        1
         req proc arrengement:        1        1        1
                parameter set:Reactive MD-force field: nitramines (RDX/HMX/TATB/PETN)               
                time step[fs]:    2.50E-01
 MDMODE CURRENTSTEP NTIMESTPE:  1         0       100
  isQEq,QEq_tol,NMAXQEq,qstep:     1   1.0E-07   500     1
                Lex_fqs,Lex_k:   1.000   2.000
            treq,vsfact,sstep:     300.000   1.000      100
                  fstep,pstep:   100    10
               NATOMS GNATOMS:                     168                     168
                         LBOX:       1.000       1.000       1.000
                  Hmatrix [A]:         13.180          0.000          0.000
                  Hmatrix [A]:          0.000         11.570          0.000
                  Hmatrix [A]:          0.000          0.000         10.710
               lata,latb,latc:      13.180      11.570      10.710
          lalpha,lbeta,lgamma:      90.000      90.000      90.000
               density [g/cc]:    1.8061
         # of linkedlist cell:     4     3     3
            maxrc, lcsize [A]:     3.160        3.29      3.86      3.57
    # of linkedlist cell (NB):     4     3     3
              lcsize [A] (NB):      3.29      3.86      3.57
     MAXNEIGHBS, MAXNEIGHBS10:    30   700
            NMINCELL, NBUFFER:     3    30000
    FFPath, DataDir, ParmPath:      ffield          DAT      rxmd.in
          # of atoms per type:          24 - 1          48 - 2          48 - 3          48 - 4
----------------------------------------------------------------
nstep  TE  PE  KE: 1-Ebond 2-(Elnpr,Eover,Eunder) 3-(Eval,Epen,Ecoa) 4-(Etors,Econj) 5-Ehbond 6-(Evdw,EClmb,Echarge)
        0 -9.82464E+01 -9.82464E+01  0.00000E+00 -1.369E+02  1.287E+00 -1.362E+00  5.208E-01 -1.398E-03  3.821E+01     0.00    0.00    0.00  41    0.36    0.23
       10 -9.82465E+01 -9.82467E+01  2.32025E-04 -1.369E+02  1.290E+00 -1.364E+00  5.214E-01 -1.397E-03  3.821E+01     0.08    0.00   -0.00  32    0.36    0.27
       20 -9.82466E+01 -9.82471E+01  4.80178E-04 -1.369E+02  1.287E+00 -1.366E+00  5.202E-01 -1.408E-03  3.821E+01     0.16    0.00   -0.00   4    0.36    0.25

...


       total (sec):       2.9980         2.9980
----------------------------------------------
    rxmd successfully finished

```

To learn more about **rxmd**, please refer to [RXMD Manual](https://github.com/USCCACS/rxmd/blob/master/doc/ReaxFF/RXMDManual.md).

## License

This project is licensed under the GPU 3.0 license - see the [LICENSE.md](https://github.com/USCCACS/rxmd/blob/master/LICENSE.md) file for details


## Publications
* Mechanochemistry of shock-induced nanobubble collapse near silica in water
K. Nomura, R. K. Kalia, A. Nakano, and P. Vashishta,
[Applied Physics Letters 101, 073108: 1-4  (2012)](http://aip.scitation.org/doi/10.1063/1.4746270)

* Structure and dynamics of shock-induced nanobubble collapse in water
M. Vedadi, A. Choubey, K. Nomura, R. K. Kalia, A. Nakano, P. Vashishta, and A. C. T. van Duin,
[Physical Review Letters 105, 014503: 1-4  (2010)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.105.014503)

* Embrittlement of metal by solute segregation-induced amorphization
H. Chen,R. K. Kalia, E. Kaxiras, G. Lu, A. Nakano, K. Nomura, A. C. T. van Duin, P. Vashishta, and Z. Yuan,
[Physical Review Letters 104, 155502: 1-4  (2010)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.104.155502)

* Metascalable molecular dynamics simulation of nano-mechano-chemistry
F. Shimojo, R. K. Kalia, A. Nakano, K. Nomura, and P. Vashishta,
[Journal of Physics: Condensed Matter 20, 294204: 1-9  (2008)](http://iopscience.iop.org/article/10.1088/0953-8984/20/29/294204)

* A scalable parallel algorithm for large-scale reactive force-field molecular dynamics simulations
K. Nomura, R. K. Kalia, A. Nakano, and P. Vashishta,
[Computer Physics Communications 178, 73-87  (2008)](http://www.sciencedirect.com/science/article/pii/S0010465507003748)
