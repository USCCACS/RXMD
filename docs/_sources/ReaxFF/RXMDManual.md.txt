## 1. Initial System Preparation 

  1. Obtain the unit cell or the prepared system coordinates in xyz file format.  
  For example: If properties of a system such as RDX, FOX-7, TATB, CaDPA etc. is to be determined, obtain the unit cell for the desired system from literature. Or, if the simulation is to be performed for some other systems containing voids, vacuums etc. prepare the system accordingly and verify the coordinates in xyz format using some visualization software like VMD, Ovito, VESTA etc. to check if the system looks right or not.  
  
  2. After verification of the system by visualizing the coordinates in xyz file format with some molecular visualization software as mentioned above, the system is read by geninit.f90 code. This code will prepare initial binary files which will be read by rxff code to perform the simulation.
    
**Parameters and format required by geninit.f90 file:**

**2.1.Specific format read by the code**
    
    As mentioned above the coordinates of the system under consideration should be in xyz file format and they should be unscaled units. The format of xyz file read by geninit.f90 file requires some extra specifications as compared to a regular xyz file.  
    The format is as follows:  
    Line 1: Natoms (number of atoms)	 tag (some string to describe your system)  
    Line 2: System size (lata	latb	latc	alpha	beta	gamma)  
    Line 3 – Line Natoms + 2: Atomic coordinates (atype	x	y	z)  
    where atype is a string for describing the atom type such as C, H, O, N, S, Mo, Ni, Al and x, y, z specifies the Cartesian coordinates.  
    
    
**2.2.Parameters to be specified for geninit.f90 file**
    
    1.Number of Processors: The number of processors required in x, y and z format should be specified as arguments to vprocs array as:
    
    vprocs(3) = (/2,2,2/) if the simulation requires a total of 8 processors such that 2 processors are required in each of x, y and z directions.

    Important: To estimate the number of processor for correct computation please ensure that you divide the system length in x, y and z (in angstrom (A)) direction by 12 A and use the quotient as the arguments to vprocs(3).

    Note: If the number of processor  in direction  comes out to be odd then number of processors in direction  should be equal to . This is important for MPI communications performed by the code.

    For Example: 
    1.If your system is  then the total number of processors in x, y and z direction should be .

    2.If your system is  then the total number of processors in x, y and z direction should be . Since 3 and 5 processors cannot be used in x and y direction as they will violate the MPI communications performed in the code and the computation will experience a deadlock due to communication error.
    
    
**2.3.Multiplicity in x, y and z directions**  

    If the simulation needs to be performed on a scaled system then the integral multiples in x, y and z directions should be specified as arguments to array mc(3).  
  
    For example:
    1.If the system is  13x13x17 Aand you specify mc(3) = (/2,3,2/) then the final system size will be  and hence the number of processors  	should be specified accordingly as mentioned above in 2.2.  
    
    
**These steps will ensure that the binaries created are correct and hence will ensure appropriate simulation results, if the simulation is performed properly**


## 2. Input file description (rxmd.in)

  **The input file rxmd.in has following format:**
	
  <mode>
  <dt>	<time_step>
  <treq>  <vsfact>  <sstep>
  <fstep> <pstep>
  <vprocs>
  <isQEq> <NMAXQEq> <QEq_tol> <qstep>
  <isBinary> <isBondFile> <isPDB>
  <ftol>
  
 ***Here the keywords mean:***  
 
  
      dt  :Increment in time during a single molecular dynamics step (usually in fs).  
      
      time_step: Total number of time step the simulation needs to be run. So the total time of simulation will be (time_step*dt) fs.  
      
      treq: The temperature at which simulation is required to be done.  
      
      vsfact: The factor by which velocity is to be scaled. This factor will ensure quenching or heating of the system, as a value above 1 will ensure heating while a value below 1 will ensure cooling of the system.   
      
      sstep: This parameter dictates the rate of heating or cooling of the system. A higher 
             value will mean gradual heating or cooling while a lower value will mean rapid heating or cooling.  
             
      fstep: This parameter controls when the trajectory in various output file formats will be saved to a file.   
      
      pstep:  This parameter controls the print out information to the standard output or to a file related to various properties such as energy, temperature pressure etc.  
      
      vprocs: This parameter should be same as specified for geninit.f90 file so that the number of processors remain consistent.  
      
      Mode : mode for doing the simulation. The Rxff code provides following mode the perform simulations.
           1. Mode 0: This mode provide velocity to the atoms distributed according to a Gaussian distribution samples using box-muller algorithm and scaled to be at temperature = treq. Also, it initializes the charges on each atom to estimate electro static forces experienced on each atom.
           2. Mode 4: Scale velocities of the atoms every sstep by a factor of vsfact. 
           3. Mode 5: Scale velocities of the atoms according to the total kinetic energy of the system.
           4. Mode 6: This mode provide velocity to the atoms distributed according to a Gaussian distribution samples using box muller algorithm and scaled to be at temperature = treq every sstep.
           5. Mode 7: This mode scales the velocities of the atoms based on their individual kinetic energies every sstep to obtain an NVE ensemble.  
        
      isQEq: This parameter is true or false specifying whether the charge equilibration scheme is required or not.  
      
      QEq_tol: The required tolerance for electro-static term relative to its previous value  
      
      NMAXQEq: The maximum number of iterations to perform for charge equilibration step if electro-static terms do not approach to a value less than QEq_tol between any 	two consecutive iterations.  
      
      qstep: This parameter defines the frequency at which charge equilibration is to be done during the simulation.  
      
      isBinary: The parameter can be given a true value or a false value and defines
	              whether a binary file will be saved or not.
                
      isBondfile: This parameter can be given a true or a false value and defines whether	information related to the bond list of every atom will be saved to a file or not.  
      
      isPDB: This parameter can be given a true or a false value and defines whether information related to atomic trajectories is saved to the file in PDB format or not.


## 3.Running a parallel job:

1.The first step is to create binaries corresponding to your system which can be done 	by following the steps in 1 in init folder.  
2.After creating binaries ensure that you are at the root directory such that init 	directory is one level below you. Then create a folder called DAT and copy all the 
	binaries of the form of rxff* from init to DAT by doing 
  `cp init/rxff* DAT/.`  
3.Then run mkdir.sh in the unit folder as sh mkdir.sh from the root directory as 	mentioned in step 2. This will create as many folders as binary files for every 	processor to read and write data.  
4. Then prepare your input file by specifying reasonable values for various parameters 	as mentioned in 2.  
5. Then in the src folder compile the code using make command and prepare an 	executable. This executable rxmd will be used for doing the computation.  
6. From your root directory type this command for generating an output log file.  

  `mpirun –np ${nprocs} ./rxmd > log`   
mpirun is the compiler which was used for compiling the code and preparing the 	executable in step 5. You can use other compilers as well which are compatible to the 	requirements of our code.  
${nprocs} is the total number of processors required for performing the simuation.  
  log is the name of output file where the output results will be saved.




