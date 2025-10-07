# rxmd : Linear Scalable Parallel ReaxFF Molecular Dynamics Simulator

**rxmd** has been developed to simulate large-scale Reactive Force Field molecular dynamics (MD) simulations on from commodity laptops to high-end supercomputing platforms. **rxmd** has been used in a various class of material studies, such as shock-induced chemical reactions, stress corrosion cracking, underwater bubble collapse, fracture of self-healing ceramics and oxidation of nanoparticles. 

## 0. Prerequisites

**rxmd** is designed to be simple, portable and minimally dependent on 3rd party library. You will need 1) a Fortran compiler that supports OpenMP, and 2) MPI (Message Passing Interface) library for parallel and distributed simulation. Modern Fortran compilers natively support OpenMP, and you can find many freely available MPI libraries online. Please refer to MPI library developer website about how to install their library. 3) PyTorch machine learning framework 4) C++ Compiler

**rxmd** has been tested on following environments.

### Fortan Compiler: 
```
GNU Fortran (GCC) 6.1.0
Intel Fortran (IFORT) 17.0.4
IBM XL Fortran V14.1
Intel Fortran Compiler 2025.1
```

### MPI library: 
```
OpenMPI 1.8.8
MPICH2
MVAPICH2 
Cray Mpich 7.6.0
Intel MPI Library 2025.1
```

### PyTorch
```
PyTorch 2.8
```

### C++ Compiler
```
Intel oneAPI DPC++/C++ Compiler 2025.1
```

## 1. Getting Started

To get started,  clone this repository to your computer. 
```
~$ git clone https://github.com/USCCACS/rxmd.git
```

## 2. How to build RXMD

### 2.1 Working Directory
Frist, change working directory to **rxmd/**
```
~$ cd rxmd-master
```
you will see following files and directories.

```
rxmd $ ls
.
├── build.sh
├── CMakeLists.txt
├── config
├── DAT
├── docs
├── examples
├── init
├── LICENSE.md
├── pot
├── README.md
├── regtests
├── rxmdtorch
├── src
├── unittests
└── util
```

Here, two directories, **src/** and **init/**, are especially important for you. **src/** contains all rxmd source codes and **init/** has a program and input files to generate an initial configurations for simulation. 

### 2.2 Build with CMake

Use the script `build.sh` uses CMake 3.16 to manage building the source code on your hardware.

This script will first build the `rxmdtorch` library and link this to the executable `rxmd` in the `bin` directory.

In order for CMake to find the libraries and packages at the time of build generation, ensure the adding directories to the PATH environment variable before using the build script.

#### 2.2.1 Example build on Aurora

Load a version of oneAPI and Python compiler with
```
module load cmake
module load frameworks
module use /soft/compilers/oneapi/2025.1.3/modulefiles/oneapi/public/
module load 2025.1.3
```

Install IPEX and PyTorch version which are suitable for the oneAPI version which is loaded. It is advisable to do this on your virtual environment by doing `source /path/to/new/venv/bin/activate`
```
python -m pip install torch==2.8.0 torchvision==0.23.0 torchaudio==2.8.0 --index-url https://download.pytorch.org/whl/xpu
python -m pip install intel-extension-for-pytorch==2.8.10+xpu oneccl_bind_pt==2.8.0+xpu --extra-index-url https://pytorch-extension.intel.com/release-whl/stable/xpu/us/
```
At any given time newer versins can be installed from the [IPEX page](https://github.com/intel/intel-extension-for-pytorch).

Finally, build by using the CMakeLists.txt in the project directory
```
cmake .. -DCMAKE_PREFIX_PATH=`python -c 'import torch; print(torch.utils.cmake_prefix_path)'`
make -j 16
```
or the build script
```
~/RXMD> sh build.sh
```

### 2.3 Prepare Initial Geometry

Next step is to generate initial MD geometry. This is done with the help of the `geninit.py` script which is located inside the `init` directory.

`geninit.py` reads a geometry file (input.xyz by default), and a replication parameter to replicate the geometry and save the entire initial MD geometry into `rxff.bin` file. The generated `rxff.bin` should then be placed in the `DAT/` directory for the `rxmd` executable to find.

sample usage
```
python geninit.py -i ../tobe.xyz -em ../rxmdnn.in
cp rxff.bin ../DAT/
```

### 2.4 Load Machine Learning Interatomic Potential Model

The first line of the `rxmdnn.in` file is written in the following format:
```
<model name> <path to MLIP model> <atomic number> [<species> <atomic number>]
```

One such model which is available to use can be downloaded from [here](https://zenodo.org/records/14915165/files/afm256_01_HL.pt)

With the `rxmd` executable, initial geomerty input `DAT/rxff.bin` in place and the MLIP potential specified in `rxmdnn.in` with the species, you are ready to start a simulation.

## 3. How to run

Default input parameters are set to run a single process job. In **rxmd.in**, the parameter **vprocs** defines how many MPI ranks in x, y, and z directions. Make sure you have **1 1 1** here. 

```
rxmd $ grep vprocs rxmd.in 
1 1 1                <vprocs>
```

To run single MPI rank job on a typical Linux computer, from the project directory you can simply type
```
rxmd $ ./bin/rxmd
```

How to run a multi process job depends on which MPI library you use, but most likely **mpirun** just works for you. 

```
rxmd $ mpirun -np nprocessors ./bin/rxmd
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

## 4. License

This project is licensed under the GPL v3 license - see the [LICENSE.md](https://github.com/USCCACS/rxmd/blob/master/LICENSE.md) file for details


## 5. Publications
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
