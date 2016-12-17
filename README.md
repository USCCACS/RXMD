# rxmd
ReaxFF MD code repository 

## Currently tested env
> source /usr/usc/intel/16.0/setup.sh

> source /usr/usc/openmpi/1.8.7/setup.sh.intel

## mk.sh and init.sh
There are several helper scripts in util directory. 

To compile the executable. Check Makefile in 'src' directory before compile. 
> ./util/mk.sh 

To prepare initial configuration. init.sh compiles geninit.F90, generate initial config binary 'rxff.bin' and copy it in 'DAT' directory following Makefile in 'init' directory. Again please check which compiler you use in the Makefile. 
> ./util/init.sh

## mulithreaded potential energy, bond-order calc, neighborlist, charge equilibration support
current version supports OpenMP mulithreading in all major functions. You can control the number of OpenMP threads by passing OMP_NUM_THREADS environment variable to mpirun.

Example) 32 MPI ranks on 2 node case, 64 threads in total. 

> qsub -I -d . -l nodes=2:ppn=16:priya_IB,walltime=2:00:00

> mpirun -x OMP_NUM_THREADS=2 --bind-to none -npernode 16 -np 32 ./rxmd

###sample output
```
Macintosh:rxmd knomura$ ./rxmd 
              rxmd has started
----------------------------------------------------------------
         req/alloc # of procs:        1  /        1
         req proc arrengement:        1        1        1
                parameter set:Reactive MD-force field: nitramines (RDX/HMX/TATB/PETN)               
                time step[fs]:    2.50E-01
 MDMODE CURRENTSTEP NTIMESTPE:  1         0       100
  isQEq,QEq_tol,NMAXQEq,qstep:     1   1.0E-06   500     1
                Lex_fqs,Lex_k:   1.000   2.000
            treq,vsfact,sstep:     300.000   1.000      100
                  fstep,pstep:   100    10
               NATOMS GNATOMS:                    4536                    4536
                         LBOX:       1.000       1.000       1.000
                  Hmatrix [A]:         39.540          0.000          0.000
                  Hmatrix [A]:          0.000         34.710          0.000
                  Hmatrix [A]:          0.000          0.000         32.130
               lata,latb,latc:      39.540      34.710      32.130
          lalpha,lbeta,lgamma:      90.000      90.000      90.000
               density [g/cc]:    1.8061
         # of linkedlist cell:    12    10    10
            maxrc, lcsize [A]:     3.160        3.29      3.47      3.21
    # of linkedlist cell (NB):    13    11    10
              lcsize [A] (NB):      3.04      3.16      3.21
     MAXNEIGHBS, MAXNEIGHBS10:    30   700
            NMINCELL, NBUFFER:     3    30000
    FFPath, DataDir, ParmPath:      ffield          DAT      rxmd.in
          # of atoms per type:         648 - 1        1296 - 2        1296 - 3        1296 - 4
----------------------------------------------------------------
nstep  TE  PE  KE: 1-Ebond 2-(Elnpr,Eover,Eunder) 3-(Eval,Epen,Ecoa) 4-(Etors,Econj) 5-Ehbond 6-(Evdw,EClmb,Echarge)
        0 -9.82458E+01 -9.82458E+01  0.00000E+00 -1.369E+02  1.287E+00 -1.361E+00  5.209E-01 -1.408E-03  3.821E+01     0.00    0.00   -0.00  37    0.21    0.46
       10 -9.82459E+01 -9.82461E+01  2.27213E-04 -1.369E+02  1.290E+00 -1.363E+00  5.215E-01 -1.407E-03  3.821E+01     0.08    0.00    0.00  30    0.21    1.40
       20 -9.82460E+01 -9.82465E+01  4.74545E-04 -1.369E+02  1.286E+00 -1.365E+00  5.202E-01 -1.419E-03  3.821E+01     0.16    0.00   -0.00   1    0.21    1.19
       30 -9.82462E+01 -9.82470E+01  7.78722E-04 -1.369E+02  1.283E+00 -1.365E+00  5.186E-01 -1.417E-03  3.821E+01     0.26    0.00   -0.00   1    0.21    1.05
       40 -9.82461E+01 -9.82474E+01  1.29512E-03 -1.369E+02  1.295E+00 -1.364E+00  5.202E-01 -1.425E-03  3.821E+01     0.43    0.00    0.00  24    0.21    1.31
       50 -9.82461E+01 -9.82479E+01  1.78842E-03 -1.369E+02  1.298E+00 -1.364E+00  5.202E-01 -1.418E-03  3.822E+01     0.60    0.00    0.00   1    0.21    1.18
       60 -9.82459E+01 -9.82483E+01  2.36084E-03 -1.369E+02  1.291E+00 -1.365E+00  5.180E-01 -1.290E-03  3.823E+01     0.79    0.00   -0.00   1    0.21    1.13
       70 -9.82460E+01 -9.82492E+01  3.13668E-03 -1.369E+02  1.292E+00 -1.368E+00  5.176E-01 -1.289E-03  3.822E+01     1.05    0.00   -0.00   1    0.21    1.16
       80 -9.82461E+01 -9.82500E+01  3.95321E-03 -1.369E+02  1.293E+00 -1.367E+00  5.162E-01 -1.290E-03  3.821E+01     1.33    0.00    0.00   1    0.21    1.26
       90 -9.82462E+01 -9.82510E+01  4.85092E-03 -1.369E+02  1.293E+00 -1.364E+00  5.147E-01 -1.297E-03  3.822E+01     1.63    0.00   -0.00  30    0.21    1.63
----------------------------------------------
        MAXNEIGHBS:           12
      MAXNEIGHBS10:          447
  MAXNBUFFER(MOVE):         4563
  MAXNBUFFER(COPY):        17226

               QEq:       4.6220         4.6220
    qeq_initialize:       1.6910         1.6910
           get_hsh:       1.6430         1.6430
      get_gradient:       0.8850         0.8850

        LINKEDLIST:       0.0870         0.0870
         COPYATOMS:       0.3540         0.3540
      NEIGHBORLIST:       0.3330         0.3330
     GetNBPairList:       0.8400         0.8400

            BOCALC:       0.3650         0.3650
            ENbond:       1.3500         1.3500
             Ebond:       0.0550         0.0550
             Elnpr:       0.3700         0.3700
               Ehb:       0.9450         0.9450
               E3b:       1.3320         1.3320
               E4b:       2.1300         2.1300
  ForceBondedTerms:       0.1010         0.1010

          WriteBND:       0.0000         0.0000
          WritePDB:       0.0000         0.0000
           ReadBIN:       0.0180         0.0180
          WriteBIN:       0.0000         0.0000
       Memory (GB):       0.2145         0.2145

       total (sec):      12.4920        12.4920
----------------------------------------------
    rxmd successfully finished
'''
