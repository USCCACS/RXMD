# rxmd
ReaxFF MD code repository 

## Currently tested env
source /usr/usc/intel/14.0/setup.sh

source /usr/usc/openmpi/1.8.7/setup.sh.intel


## threaded non-bonding pair support 
this version requires OpenMP enabled as default. If you are going to use two Infiniband nodes in priya queue to run a 16 MPIrank rxmd job,

> qsub -I -d . -l nodes=2:ppn=16:priya_IB,walltime=2:00:00

> mpirun -x OMP_NUM_THREADS=2 --bind-to none -npernode 8 -np 16 ./rxmd

###sample output

```
[knomura@hpc4018 rxmd]$ mpirun -x OMP_NUM_THREADS=2 --bind-to none -npernode 8 -np 16 ./rxmd
INFO: mdmode==0, setting isQEQ is 1. Atomic velocities are scaled to      300.000 [K] every    100 steps.
         # of atoms per type:   1           1536
         # of atoms per type:   2           3072
         # of atoms per type:   3           3072
         # of atoms per type:   4           3072
         # of atoms per type:   5              0
         # of atoms per type:   6              0
         # of atoms per type:   7              0
----------------------------------------------------------------
         req/alloc # of procs:       16  /       16
         req proc arrengement:        4        2        2
                parameter set:Reactive MD-force field: nitramines (RDX/HMX/TATB/PETN)               
                time step[fs]:    2.50E-01
                       MDMODE:    0
  isQEq,QEq_tol,NMAXQEq,qstep:     1   1.0E-06   500  1000
                Lex_fqs,Lex_k:   1.000   2.000
            treq,vsfact,sstep:     300.000   1.000      100
                  fstep,pstep:  1000   100
                  CURRENTSTEP:       2000
                    NTIMESTPE:         1000
               NATOMS GNATOMS:                     673                   10752
         NBUFFER_N, NBUFFER_P:     20000      2000
                         LBOX:       0.250       0.500       0.500
                  Hmatrix [A]:         52.720          0.000          0.000
                  Hmatrix [A]:          0.000         46.280          0.000
                  Hmatrix [A]:          0.000          0.000         42.840
               lata,latb,latc:      52.720      46.280      42.840
          lalpha,lbeta,lgamma:      90.000      90.000      90.000
               density [g/cc]:    1.8061
         # of linkedlist cell:     4     7     6
            maxrc, lcsize [A]:     3.160        3.29      3.31      3.57
            NMINCELL, NCELL10:     3     4
     MAXNEIGHBS, MAXNEIGHBS10:    30   600
         NBUFFER_P, NBUFFER_N:     2000    20000
  FFPath, DataPath, ParmPath: ffield DAT rxmd.in
----------------------------------------------------------------
nstep  TE  PE  KE: 1-Ebond 2-(Elnpr,Eover,Eunder) 3-(Eval,Epen,Ecoa) 4-(Etors,Econj) 5-Ehbond 6-(Evdw,EClmb,Echarge)
 
     2000 -9.65941E+01 -9.74858E+01  8.91767E-01 -1.367E+02  1.553E+00  8.188E-01  4.304E+01 -1.414E+01  9.438E+00   299.28    0.00   -0.00  56    0.52
     2100 -9.65915E+01 -9.74722E+01  8.80722E-01 -1.366E+02  1.556E+00  8.108E-01  4.302E+01 -1.412E+01  9.438E+00   295.57    0.00   -0.00   1    8.22
     2200 -9.65784E+01 -9.74624E+01  8.84019E-01 -1.367E+02  1.557E+00  8.115E-01  4.305E+01 -1.412E+01  9.438E+00   296.68    0.00   -0.00   1    7.37
     2300 -9.65686E+01 -9.74574E+01  8.88758E-01 -1.367E+02  1.552E+00  8.076E-01  4.305E+01 -1.411E+01  9.438E+00   298.27    0.00   -0.00   1    7.39
     2400 -9.65635E+01 -9.74572E+01  8.93672E-01 -1.367E+02  1.560E+00  8.079E-01  4.306E+01 -1.411E+01  9.438E+00   299.92    0.00   -0.00   1    7.39
     2500 -9.65626E+01 -9.74524E+01  8.89764E-01 -1.366E+02  1.543E+00  7.981E-01  4.304E+01 -1.411E+01  9.438E+00   298.60    0.00   -0.00   1    7.38
     2600 -9.65583E+01 -9.74492E+01  8.90916E-01 -1.366E+02  1.548E+00  8.130E-01  4.304E+01 -1.411E+01  9.438E+00   298.99    0.00   -0.00   1    7.36
     2700 -9.65552E+01 -9.74349E+01  8.79693E-01 -1.366E+02  1.556E+00  8.153E-01  4.305E+01 -1.411E+01  9.438E+00   295.22    0.00   -0.00   1    7.43
     2800 -9.65403E+01 -9.74311E+01  8.90748E-01 -1.367E+02  1.550E+00  8.145E-01  4.307E+01 -1.411E+01  9.438E+00   298.93    0.00   -0.00   1    7.44
     2900 -9.65372E+01 -9.74360E+01  8.98805E-01 -1.367E+02  1.554E+00  8.249E-01  4.305E+01 -1.411E+01  9.438E+00   301.64    0.00   -0.00   1    7.46
Max MAXNEIGHBS, Max MAXNEIGHBS10, Max NBUFFER_P, Max NBUFFER_N:           14         456         690        9384
               QEq:       0.4926         0.4718
          QEq_COPY:       0.1701         0.0800
        LINKEDLIST:       0.0219         0.0170
         COPYATOMS:       3.3327         2.4334
      NEIGHBORLIST:      11.8363         6.2388
GetNonbondingPairs:      21.4818        15.8759
            BOCALC:      18.9233        12.4252
            ENbond:      15.5581         7.0433
             Ebond:       0.7056         0.4654
             Elnpr:       1.7178         1.3502
               Ehb:       8.3507         5.0880
               E3b:       6.7780         6.3131
               E4b:       9.1922         8.5893
  ForceBondedTerms:       1.2148         0.8521
     COPYATOMS(-1):      19.8682         2.1551
       total (sec):      74.9829        74.9627
    program finished
```

