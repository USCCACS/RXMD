!-------------------------------------------------------------------------------------------
module atoms
!-------------------------------------------------------------------------------------------
#ifdef INTEROP
use iso_c_binding 
#else 
#define c_double 8
#define c_int 4
#endif

include 'mpif.h'

!--- command arguments 
logical :: isFF=.false., isData=.false., isMDparm=.false.
character(64) :: FFPath="ffield", DataPath="DAT", ParmPath="rxmd.in"

!--- Size of I/O group for collective I/O 
integer,parameter :: nio=1

!--- For array size statistics
!  1-NATOMS, 2-nbrlist, 3-nbrlist for qeq, 4-NBUFFER_P, 5-not used, 
!  6-NBUFFER_P for qeq
integer,parameter :: nmaxas=5
integer,allocatable :: maxas(:,:)

!--- lattice parameters 
real(8) :: lata,latb,latc,lalpha,lbeta,lgamma

integer :: myid, nprocs, ierr, myparity(3), vID(3)

!<sbuffer> send buffer, <rbuffer> receive buffer
real(8),allocatable :: sbuffer(:), rbuffer(:)
   
!<ns> # of atoms to be sent, <nr> # of atoms to be received, <na> # of all of transfered atoms.
!<ne> # of elements for one atom.
!     Example) In case atom type, position and velocity to be sent,  ne = 1+3+3 = 7
integer :: ns, nr, na, ne

!<NE_COPY>,<NE_MOVE>,<NE_CPBK> :: Number of Elements to COPY, MOVE atoms and CoPy BacK force. 
!  integer,parameter :: NE_COPY=6, NE_MOVE=10, NE_CPBK=10
integer,parameter :: NE_COPY=13, NE_MOVE=12, NE_CPBK=10

!<MAXLAYERS> MAXimum # of linkedlist cell LAYERS.
integer,parameter :: MAXLAYERS=5
        
!<scl([xyz], start-end cell, direction flag, # of layers)> :: cell indices to be sent
!<scl1()> <scl2()> store opposite side layer indices.
integer :: scl1(3,2,6,0:MAXLAYERS), scl2(3,2,6,0:MAXLAYERS)
        
! <target_node> stores partner node ID in the 6-communications. 
! if targe_node(i)==-1, the node doesn't have a partner in i-direction.
integer :: target_node(6)

! For benchmarking, <vprocs> and <mc> will be read from vprocs.in
integer :: vprocs(3)

!<mc(3)> # of unit cells in each directions
integer :: mc(3)

real(8),allocatable:: rc(:), rc2(:)   !<RCUT>: cutoff length for sigma-bonding.
real(8),allocatable:: rcpi(:), rcpp(:)!      : cutoff length for other bonding.

real(8),parameter :: MINBOSIG = 1d-3      !<minBOsig>: criterion to decide <rc> 
real(8),parameter :: MINBO0 = 1d-4       !<minBO0>: cutoff bond order 
real(8),parameter :: cutof2_esub = 1d-4
!real(8),parameter :: cutof2_bo = 1.d-2
real(8),parameter :: cutof2_bo = 1.d-3
integer,parameter :: is_idEh = 1

!real(8),parameter :: MINBOSIG = 1d-4      !<minBOsig>: criterion to decide <rc> 
!real(8),parameter :: MINBO0 = 0.d0       !<minBO0>: cutoff bond order 
!real(8),parameter :: cutof2_esub = 0.d0
!real(8),parameter :: cutof2_bo = 1d-4
!integer,parameter :: is_idEh = 0

!integer :: NBUFFER_P=5000
!integer,parameter :: MAXNEIGHBS=50  !<MAXNEIGHBS>: Max # of Ngbs one atom may have. 
!integer,parameter :: MAXNEIGHBS10=200 !<MAXNEIGHBS>: Max # of Ngbs within 10[A]. 

integer :: NBUFFER_P=20000
integer,parameter :: MAXNEIGHBS=30  !<MAXNEIGHBS>: Max # of Ngbs one atom may have. 
integer,parameter :: MAXNEIGHBS10=600 !<MAXNEIGHBS>: Max # of Ngbs within 10[A]. 

integer,parameter :: NMINCELL=3  !<NMINCELL>: Nr of minimum linkedlist cell <-> minimum grain size.
integer :: NCELL10               !<NCELL10>: Nr of linkedlist cell to cover up 10[A].
real(8),parameter :: MAXANGLE= 0.999999999999d0 
real(8),parameter :: MINANGLE=-0.999999999999d0
real(8),parameter :: NSMALL = 1.d-10
real(8) :: maxrc                        !<maxRCUT>: Max cutoff length. used to decide lcsize.

real(8),parameter :: pi=3.14159265358979d0

! position, atom type, velocity, force & charge
real(kind=c_double),allocatable,target :: pos(:,:),atype(:)
real(8),allocatable :: v(:,:), f(:,:), q(:)

! atomic stress tensor
real(8),allocatable :: astr(:,:) 
real(8) :: pint(3,3)

!--- coefficient of bonding energy derivative 
real(8),allocatable :: ccbnd(:)

!--- weight of each component on stress calculation. These numbers depends on which system been computed.
integer :: ncmp
real(8) :: acmp(0:10)

real(8) :: HH(3,3,0:1), HHi(3,3), MDBOX, LBOX(0:3), OBOX(1:3) !MD box, local MD box, origin of box.
integer :: NATOMS         !local # of atoms
integer(8) :: GNATOMS     !global # of atoms

!<llist> Linked List
!<header> header atom of linkedlist cell.
!<nacell> Nr of atoms in a likedlist cell.
integer,allocatable :: llist(:), header(:,:,:), nacell(:,:,:)

!<nbllist> Linked List for non-bonding interaction
!<nbheader> header atom of linkedlist cell for non-bonding interaction
!<nbnacell> Nr of atoms in a likedlist cell for non-bonding interaction
integer,allocatable :: nbllist(:), nbheader(:,:,:), nbnacell(:,:,:)

!<nbrlist> neighbor list, <nbrindx> neighbor index
integer(kind=c_int),pointer :: nbrlist(:,:), nbrindx(:,:)

!<nbplist> neighbor list of nonbonding interaction, non-bonding pair list
integer,allocatable :: nbplist(:,:)


!<BO> Bond Order of atoms i-j (nearest neighb only) - (Eq 3a-3d)
real(kind=c_double),pointer :: BO(:,:,:) 

real(8),allocatable :: delta(:)

!--- Output variables from the BOp_CALC() subroutine:
real(8),allocatable :: deltap(:,:)

real(8),pointer :: dln_BOp(:,:,:)

real(8),pointer :: dBOp(:,:)


!--- For NEW DBO calc:
real(8),allocatable :: exp_delt1(:,:), exp_delt2(:,:)  ! exp( -pboc#(inxn) * deltap(i,1) ) - {inxn, i}   

!--- A[0123] coefficients for force calculation 
real(8),pointer:: A0(:,:),A1(:,:), A2(:,:), A3(:,:) 

!--- <switch> flag to omit pi and double pi bond in bond-order prime calculation.
real(8),allocatable :: switch(:,:) 

!--- Passed between Elnpr and E3body
real(8),allocatable :: nlp(:), dDlp(:) !Number of Lone Pairs, its derivatives.

! TE: Total Energy,  KE: Kinetic Energy,  PE :: Potential Energies
!  0-Esystem, 1-Ebond, 2-Elp, 3-Eover, 4-Eunder, 5-Eval, 6-Epen
!  7-Ecoa,  8-Etors, 9-Econj, 10-Ehbond, 11-Evdwaals, 12-Ecoulomb 13-Echarge
real(8) :: TE, KE, PE(0:13)
real(8) :: GTE, GKE, GPE(0:13)

!--- output file format 
logical :: isBinary, isBondFile, isPDB

!--- one vector charge equlibration
! g: gradient,  h: conjugate direction, hsh: hessian * h
!real(8),allocatable :: g(:), h(:), hsh(:)

! Two vectors electrostatic energy minimization 
real(8),allocatable :: qs(:),qt(:),gs(:), gt(:), hs(:), ht(:), hshs(:), hsht(:)

!--- variables for extended Lagrangian method ---
!<Lex_fqs> fraction between two QEq vectors
!<Lex_w> spring constant
real(8),allocatable :: qsfp(:),qsfv(:),qtfp(:),qtfv(:) 
real(8) :: Lex_fqs=1.0, Lex_w=1.d0, Lex_w2=1.d0, Lex_k=2.d0
 
integer :: ast ! Allocation STatus of allocatable variables.

! <lcsize> Linked list Cell SIZE. <cc> Nr of likedlist cell in local node.
real(8) :: lcsize(3), nblcsize(3)
integer :: cc(3), nbcc(3)

integer :: nmesh, nbnmesh
integer,allocatable :: mesh(:,:), nbmesh(:,:)

!--- Unit convetors. In original ReaxFF program, the units of length is [A]
!--- energy is [kcal/mol] and mass is [amu] respectively.
!--- Most of numbers written here are obtained below site.
!--- http://physics.nist.gov/cuu/Constants/Table/allascii.txt
!--- Bohr length
real(8),parameter :: Lbohr_a  = 0.5291772108d0   ! [A]
real(8),parameter :: Lbohr_m  = 0.5291772108d-10 ! [m]

!--- Electron rest mass
real(8),parameter :: Merest_amu = 5.48580152d-4  ! [amu]
real(8),parameter :: Merest_kg  = 9.10938d-31    ! [kg]
!--- Energy in Hartree
real(8),parameter :: Ehrtre_km = 6.2751d2        ! [kcal/mol]
real(8),parameter :: Ehrtre_ev  = 27.2113845d0   ! [eV]
real(8),parameter :: Ehrtre_j = 4.35974417d-18   ! [J] 
real(8),parameter :: Eev_kcal = 23.060538d0      ! [kcal/mol]

real(8),parameter :: Ekcal_j = 6.95016611d-21  ! [J]

!--- Boltzmann Constant
real(8),parameter :: BLTZMN = 1.3806503d-23  ! [m^2 kg s^-2 K^-1 ] 

real(8),parameter :: UTEMP0 = 503.398008d0    ! Ekcal_j/BLZMN [K]
real(8),parameter :: UTEMP = UTEMP0*2.d0/3.d0 ! [K]
real(8),parameter :: USTRS = 6.94728103d0     ! [GPa]
real(8),parameter :: UDENS = 1.66053886d0     ! [g/cc]
real(8),parameter :: UTIME = 1.d3/20.455d0    ! 1 = 1/20.445[ps] = 48.88780[fs]

!--- QEq variables. 
!<isQEq> flag to run QEq routine: 0-No QEq, 1-CG, 2-Extended Lagrangian
integer :: isQEq
!<NMAXQEq> Number of MAXimum iteration in QEq routine
integer :: NMAXQEq
!<QEq_thrsld> energy criterion in QEq routine
real(8) :: QEq_tol
!<nstep_qeq> counter of iteration
integer :: nstep_qeq, qstep

!-- variables for timing
integer,parameter :: Ntimer=20
integer :: it_timer(Ntimer)=0, it_timer_max(Ntimer)=0, it_timer_min(Ntimer)=0

!-- for timing in PRINTE
real(8) ::  wt0

!---
! <mdmode> determines MD mode
integer :: mdmode
! <nstep> current MD step, <ntime_step> Total # of time steps in one MD run.
! <current_step> will be used for subsequent runs.
integer :: nstep=0, ntime_step, current_step
!<vsfact> velocity scaling factor, <dt> one time step
real(8) :: treq, vsfact, dt, dmt
integer :: sstep

!--- output format flags, explained in 'rxmdopt.in'
integer :: fstep, pstep

!--- <frcindx> FoRCe INDeX. Index to return calculated force to original atoms.
integer,allocatable :: frcindx(:)
integer :: frcptr(0:6)

!--- parameter & structure file description
character(70) :: pfile

!--- stress components
real(8) :: xx,yy,zz,yz,zx,xy
!--- atomic stress index
integer :: ia,ja

!--- conjugate_gradient
REAL(8) :: ftol   ! tolerance of energy convergence 

!--- cutoff range calculation. 
integer(8),allocatable :: natoms_per_type(:)

!--- dthm=dt/(2*mass), hmas=mass/2
real(8),allocatable :: dthm(:), hmas(:)

!--- potential teble
integer,parameter :: NTABLE=5000
real(8),allocatable :: TBL_Eclmb(:,:,:), TBL_Evdw(:,:,:), TBL_Eclmb_QEq(:,:)
real(8) :: UDR, UDRi

integer,allocatable :: ibuf(:), ibuf1(:)
integer(8),allocatable :: ibuf8(:)

end module atoms

!-------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------
module parameters
!-------------------------------------------------------------------------------------------
!  This module stores the parameters input from the rxmda.in file.  This file is 
!    intentionally very similar to the input setup used by Adri van Duin at Caltech (to aid
!    in cross-checking and updates).  However, please note that there are a few changes made
!    aside from the obvious f77 -> f90 switch.  I have tried to clearly note these. 
!-------------------------------------------------------------------------------------------

!Independant Parameters

integer :: nso    !Number of different types of atoms
integer :: nboty  !Number of different bonds given

! Atom Dependant (ie where they appear in input file - not implementation in code)
character(2),allocatable :: atmname(:)      !holds the Chemical Abbrev for each atomtype
real(8),allocatable :: Val(:),Valboc(:)  !Valency of atomtype (norm, boc) 
real(8),allocatable :: mass(:)           !mass of atomtype

real(8),allocatable :: pbo1(:), pbo2(:), pbo3(:)   !Bond Order terms
real(8),allocatable :: pbo4(:), pbo5(:), pbo6(:)   !Bond Order terms
real(8),allocatable :: pboc1(:), pboc2(:), pboc3(:)  !Bond Order correction terms (f1-5)
real(8),allocatable :: pboc4(:), pboc5(:)            !Bond Order correction terms (f1-5)  
real(8),allocatable :: v13cor(:) !<kn>

real(8),allocatable :: rat(:),rapt(:),vnq(:)   !r0s/r0p/r0pp for like bonds 
real(8),allocatable :: ovc(:)  !a flag to apply fn4 and fn5 !<kn>

real(8),allocatable :: Desig(:),Depi(:),Depipi(:)  !Bond Energy parameters (eq. 6)
real(8),allocatable :: pbe1(:),pbe2(:)             !Bond Energy parameters (eq. 6) 

real(8),allocatable :: Vale(:)                      !Lone Pair Energy parameters (eq. 7)
real(8),allocatable :: plp1(:), nlpopt(:), plp2(:)  !Lone Pair Energy parameters (eq.8-10)  

real(8),allocatable :: povun1(:), povun2(:), povun3(:), povun4(:)   !Overcoordination Energy (eq. 11)
real(8),allocatable :: povun5(:), povun6(:), povun7(:), povun8(:)   !Undercoordination Energy (eq. 12)

real(8),allocatable :: pval1(:), pval2(:), pval3(:), pval4(:), pval5(:)   !Valency Angle Energy (eq. 13a-g)
real(8),allocatable :: pval6(:), pval7(:), pval8(:), pval9(:), pval10(:)   
real(8),allocatable :: Valangle(:), theta00(:)

real(8),allocatable :: ppen1(:), ppen2(:), ppen3(:), ppen4(:)   !Penalty Energy (eq. 14ab)

real(8),allocatable :: pcoa1(:), pcoa2(:), pcoa3(:), pcoa4(:)   !Conjugation (3 body) Energy (eq.15)
real(8),allocatable :: Valval(:)

real(8),allocatable :: ptor1(:), ptor2(:), ptor3(:), ptor4(:)   !Torsional Energy Terms (eq.16abc)
real(8),allocatable :: V1(:), V2(:), V3(:)

real(8),allocatable :: pcot1(:), pcot2(:)  !Conjugation (4body) Energy (eq. 17ab)

real(8),allocatable :: phb1(:), phb2(:), phb3(:), r0hb(:)   !Hydrogren Bond Energy (eq. 18)

real(8),allocatable :: Dij(:,:), alpij(:,:), rvdW(:,:), gamW(:,:)  !Van der Waals Energy (eq. 21ab)
real(8) :: pvdW1, pvdW1h, pvdW1inv

!Taper function 
real(8),parameter :: rctap0 = 10.d0 ![A]
real(8) :: rctap, rctap2, CTap(0:7)

! hydrogen bonding interaction cutoff
real(8),parameter :: rchb = 10.d0 ![A]
real(8),parameter :: rchb2 = rchb*rchb

!Coulomb Energy (eq. 22)  
real(8),parameter:: Cclmb0 = 332.0638d0 ! [kcal/mol/A] line 2481 in poten.f
real(8),parameter:: Cclmb0_qeq = 14.4d0 ! [ev]
real(8),parameter:: CEchrge = 23.02d0   ! [ev]
real(8) :: Cclmb = Cclmb0

real(8),allocatable :: gam(:), gamij(:,:)  

!Charge Equilibration part, <chi>  electronegativity   <eta> stiffness
real(8),allocatable :: chi(:), eta(:)

!End Lost Parameters Listing

!Not Understood Parameters: 
real(8),allocatable :: bom(:)

! 2-atom combo dependant: 
real(8),allocatable :: r0s(:,:)       !Bond Order terms
real(8),allocatable :: r0p(:,:)       !  "" 
real(8),allocatable :: r0pp(:,:)      !Bond Order terms 

! <inxn2> is type of 2-body interaction 1=C-C, 2=H-C, 3=H-H
integer,allocatable :: inxn2(:,:)

integer,allocatable :: inxn3(:,:,:)
integer,allocatable :: inxn3hb(:,:,:)
integer,allocatable :: inxn4(:,:,:,:)

!Saved calculations to prevent having to recalc lots of times.
real(8),allocatable :: cBOp1(:), cBOp3(:), cBOp5(:)
real(8),allocatable :: pbo2h(:), pbo4h(:), pbo6h(:)  

!NOTE: for debugging purpose variables
real(8)  :: vpar30,vpar1,vpar2

end module parameters 

!!------------------------------------------------------------------------------------------
!module ustruct
!!  unit cell variables 
!!------------------------------------------------------------------------------------------
!integer :: unitN                 !<unitN> # of atoms in a unit cell
!real(8),allocatable :: pos0(:,:) !<x0> atom coordinate array for unit cell
!integer,allocatable :: atype0(:) !<atype0> atom type array for unit cell
!                                 ! 1=C, 2=H, 3=O, 4=N, 5=S, 6=Si, 7=Al, X=8
!real(8),allocatable :: q0(:)     ! Initial charge distributin (not used) 
!real(8) :: ul(3,3)               !<ul> unit lattice vector
!real(8),parameter :: rshift = 0.0d0 !<rshift> shift atoms coordinate (not used)
!
!end module ustruct

#ifdef INTEROP
!------------------------------------------------------------------------------------------
module interop
!------------------------------------------------------------------------------------------
use iso_c_binding
implicit none

interface
    subroutine CppInterface(param_ptr, pos_ptr, atype_ptr, nbrlist_ptr, bo_ptr) &
               bind(c, name="CppInterface")
        use iso_c_binding
            type(c_ptr),value :: param_ptr
            type(c_ptr),value :: pos_ptr
            type(c_ptr),value :: atype_ptr
            type(c_ptr),value :: nbrlist_ptr
            type(c_ptr),value :: bo_ptr
    end subroutine
end interface

contains

!-------------------------------------------------------------------------
   subroutine print_atom()
!-------------------------------------------------------------------------
   use atoms
   use iso_c_binding
   implicit none
   
   integer,parameter :: NumInterOpParams=32
   
   type(c_ptr) :: param_ptr, pos_ptr, atype_ptr
   type(c_ptr) :: nbrlist_ptr, bo_ptr
   
   integer(c_int),allocatable,target :: ParamCppInterface(:)
   
   allocate(ParamCppInterface(NumInterOpParams))
   
   ParamCppInterface(1)=myid
   ParamCppInterface(2)=NATOMS
   ParamCppInterface(3)=NBUFFER_P
   ParamCppInterface(4)=NBUFFER_P
   ParamCppInterface(5)=MAXNEIGHBS
   
   param_ptr = c_loc(paramCppInterface(1))
   pos_ptr = c_loc(pos(1,1))
   atype_ptr = c_loc(atype(1))
   nbrlist_ptr = c_loc(nbrlist(1,0))
   bo_ptr = c_loc(bo(0,1,1))
   
   call CppInterface(param_ptr, pos_ptr, atype_ptr, nbrlist_ptr, bo_ptr)

   deallocate(ParamCppInterface)

   call MPI_FINALIZE(ierr)
   stop
   
   end subroutine

end module
#endif
