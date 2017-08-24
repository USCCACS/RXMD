# rxmd : Linear Scalable Parallel ReaxFF Molecular Dynamics Simulator

**rxmd** has been developed to simulate large-scale Reactive Force Field molecular dynamcis (MD) simulations on from commodity laptops to high-end supercomputing platforms. **rxmd** has been used in various class of material researches, such as shock-induced chemical reactions, stress corrosion cracking, underwater bubble collapse, facture of self-healing ceramics and oxidation of nanoparticles. 

## Prerequisites

**rxmd** is designed to be simple, portable and minimally dependent on 3rd party library. You will only a Fortran compiler that supoprt OpenMP and MPI (Message Passing Interface) library to compile the code and run it. Modern Fortran compilers natively support OpenMP and you can find many freely distributed MPI libraries online. Please refer to MPI library developper site about how to install their library. 

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
MPICH 
MVAPCHI 
Cray Mpich 7.6.0
```

## Getting Started

To get started, simply clone this repository to your computer. 
```
~$ git clone https://github.com/USCCACS/rxmd.git
```

## How to install

Then change working directory.
```
~$ cd rxmd/rxmd
```
you will see following files and directories.

```
rxmd $ ls
DAT/          conf/         ffield        regtests/     src/          util/
Makefile.inc  doc/          init/         rxmd.in       unittests/
```

Two directories are important for you to get started here, **src/** and **init/**. **src/** to compile rxmd executable and **init/** to generate intial configuration for MD simulation. Type

```
rxmd $ make -C init/
```
to generate an initial configuration and place the input file **rxff.bin** in **DAT/** directory. Then type

```
rxmd $ make -C src/
```
to compile the rxmd exectuable. You will have **rxmd** plus **DAT/rxff.bin**, then you are ready to run a simulation. 

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
Default input parameters are set for single process job where tou can simply start simulation by typing
```
rxmd $ ./rxmd
```

How to run a multi process job depends on which MPI library you use, but most likely **mpirun** just works. 

```
rxmd $ mpirun -np nprocessors ./rxmd
```

If you see following output, congratulations! You have everything working.

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

```

## License

This project is licensed under the GPU 3.0 license - see the [LICENSE.md](LICENSE.md) file for details


## Publications
* Mechanochemistry of shock-induced nanobubble collapse near silica in water
K. Nomura, R. K. Kalia, A. Nakano, and P. Vashishta
[Applied Physics Letters 101, 073108: 1-4  (2012)] (http://aip.scitation.org/doi/10.1063/1.4746270)

* Structure and dynamics of shock-induced nanobubble collapse in water
M. Vedadi, A. Choubey, K. Nomura, R. K. Kalia, A. Nakano, P. Vashishta, and A. C. T. van Duin 
[Physical Review Letters 105, 014503: 1-4  (2010)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.105.014503)

* Embrittlement of metal by solute segregation-induced amorphization
H. Chen,R. K. Kalia, E. Kaxiras, G. Lu, A. Nakano, K. Nomura, A. C. T. van Duin, P. Vashishta, and Z. Yuan
[Physical Review Letters 104, 155502: 1-4  (2010)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.104.155502)

* Metascalable molecular dynamics simulation of nano-mechano-chemistry
F. Shimojo, R. K. Kalia, A. Nakano, K. Nomura, and P. Vashishta
[Journal of Physics: Condensed Matter 20, 294204: 1-9  (2008)](http://iopscience.iop.org/article/10.1088/0953-8984/20/29/294204)

* A scalable parallel algorithm for large-scale reactive force-field molecular dynamics simulations
K. Nomura, R. K. Kalia, A. Nakano, and P. Vashishta
[Computer Physics Communications 178, 73-87  (2008)](http://www.sciencedirect.com/science/article/pii/S0010465507003748)