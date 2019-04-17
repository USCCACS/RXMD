!-------------------------------------------------------------------------------------------
module base
!-------------------------------------------------------------------------------------------
! position, atom type, velocity, force & charge
real(8),allocatable,dimension(:),target :: atype, q
real(8),allocatable,dimension(:,:),target :: pos, v, f

character(2),allocatable :: atmname(:)  
real(8),allocatable :: mass(:)          

end module

!-------------------------------------------------------------------------------------------
module atoms
!-------------------------------------------------------------------------------------------
use utils

#ifdef NOMPI
  use nompi
#else
  include 'mpif.h'
#endif

interface potentialtable_interface
  subroutine potentialtable_ctor()
  end subroutine
end interface

interface forcefield_param_interface
  subroutine forcefield_param(path)
    character(len=:),allocatable,intent(in) :: path
  end subroutine
end interface

interface charge_model_interface
  subroutine charge_model(atype, pos, q)
    real(8),intent(in out),allocatable :: atype(:), pos(:,:), q(:)
  end subroutine
end interface

interface force_model_interface
  subroutine force_model(atype, pos, f, q)
    real(8),intent(in out),allocatable :: atype(:), pos(:,:), q(:), f(:,:)
  end subroutine
end interface

interface velocity_scaling_interface
  subroutine velocity_scaling(atype, v, factor)
    real(8),intent(in),allocatable :: atype(:)
    real(8),intent(in out),allocatable :: v(:,:)
    real(8),intent(in),optional :: factor
  end subroutine
end interface

procedure(charge_model),pointer :: charge_model_func => null()
procedure(force_model),pointer :: force_model_func => null()
procedure(forcefield_param),pointer :: get_forcefield_param => null()
procedure(potentialtable_ctor),pointer :: set_potentialtable => null()
procedure(velocity_scaling),pointer :: velocity_scaling_func => null()

!--- supported force types : 'reaxff', 'nff'
character(MAXSTRLENGTH) :: force_type='reaxff'

!--- For array size statistics
!  1-NATOMS, 2-nbrlist, 3-nbrlist for qeq, 4-NBUFFER for move, 5-NBUFFER for copy
!  6-NBUFFER for qeq
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
integer,parameter :: MODE_COPY = 1, MODE_MOVE = 2, MODE_CPBK = 3
integer,parameter :: MODE_QCOPY1 = 4, MODE_QCOPY2 = 5

integer,parameter :: NE_CPBK = 4

!<MAXLAYERS> MAXimum # of linkedlist cell LAYERS.
integer,parameter :: MAXLAYERS=5
integer,parameter :: MAXLAYERS_NB=10
        
! <target_node> stores partner node ID in the 6-communications. 
! if targe_node(i)==-1, the node doesn't have a partner in i-direction.
integer :: target_node(6)

! For benchmarking, <vprocs> and <mc> will be read from vprocs.in
integer :: vprocs(3)

real(8),allocatable:: rc(:), rc2(:)   !<RCUT>: cutoff length for sigma-bonding.

!integer :: NBUFFER=5000
!integer,parameter :: MAXNEIGHBS=50  !<MAXNEIGHBS>: Max # of Ngbs one atom may have. 
!integer,parameter :: MAXNEIGHBS10=200 !<MAXNEIGHBS>: Max # of Ngbs within the taper function cutoff. 

integer :: NBUFFER=30000
integer,parameter :: MAXNEIGHBS=30  !<MAXNEIGHBS>: Max # of Ngbs one atom may have. 
integer,parameter :: MAXNEIGHBS10=1500 !<MAXNEIGHBS>: Max # of Ngbs within 10[A]. 

integer,parameter :: NMINCELL=4  !<NMINCELL>: Nr of minimum linkedlist cell <-> minimum grain size.
real(8),parameter :: MAXANGLE= 0.999999999999d0 
real(8),parameter :: MINANGLE=-0.999999999999d0
real(8),parameter :: NSMALL = 1.d-10
real(8) :: maxrc                        !<maxRCUT>: Max cutoff length. used to decide lcsize.

! stress tensor
real(8) :: astr(6) = 0.d0

real(8) :: HH(3,3,0:1), HHi(3,3), MDBOX, LBOX(0:3), OBOX(1:3) !MD box, local MD box, origin of box.
integer :: NATOMS         !local # of atoms
integer(8) :: GNATOMS     !global # of atoms

! TE: Total Energy,  KE: Kinetic Energy,  PE :: Potential Energies
!  0-Esystem, 1-Ebond, 2-Elp, 3-Eover, 4-Eunder, 5-Eval, 6-Epen
!  7-Ecoa,  8-Etors, 9-Econj, 10-Ehbond, 11-Evdwaals, 12-Ecoulomb 13-Echarge
real(8) :: TE, KE, PE(0:13)
real(8) :: GTE, GKE, GPE(0:13)

!--- output file format 
logical :: isBinary=.false., isBondFile=.false., isPDB=.false., isXYZ=.false.

!<llist> Linked List
!<header> header atom of linkedlist cell.
!<nacell> Nr of atoms in a likedlist cell.
integer,allocatable :: llist(:), header(:,:,:), nacell(:,:,:)

 
! <lcsize> Linked list Cell SIZE. <cc> Nr of likedlist cell in local node.
real(8) :: lcsize(3), nblcsize(3)
integer :: cc(3), nbcc(3)

integer :: nmesh, nbnmesh
integer,allocatable :: mesh(:,:), nbmesh(:,:)


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
integer,parameter :: Ntimer=30
integer :: it_timer(Ntimer)=0, it_timer_max(Ntimer)=0, it_timer_min(Ntimer)=0
real(8) :: wt0=0.d0

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

!--- conjugate_gradient tolerance of energy convergence 
REAL(8) :: ftol   

!--- <frcindx> FoRCe INDeX. Index to return calculated force to original atoms.
real(8),allocatable,target :: frcindx(:)
integer :: copyptr(0:6)


!--- cutoff range calculation. 
integer(8),allocatable :: natoms_per_type(:)

!--- dthm=dt/(2*mass), hmas=mass/2
real(8),allocatable :: dthm(:), hmas(:)

!--- potential teble
integer,parameter :: NTABLE=5000
real(8),allocatable :: TBL_Eclmb(:,:,:), TBL_Evdw(:,:,:), TBL_Eclmb_QEq(:,:)
real(8) :: UDR, UDRi

!--- from cmdline_args
logical :: isSpring=.false.
real(8) :: springConst=0.d0
logical :: hasSpringForce(16)
real(8),allocatable :: ipos(:,:)

logical :: isFF=.false., isData=.false., isMDparm=.false.
character(MAXSTRLENGTH) :: FFPath0='ffield', DataDir0='DAT', ParmPath0='rxmd.in'
character(len=:),allocatable :: FFPath, DataDir, ParmPath

logical :: isLG=.false.

logical :: isEfield=.false.


!<nbllist> Linked List for non-bonding interaction
!<nbheader> header atom of linkedlist cell for non-bonding interaction
!<nbnacell> Nr of atoms in a likedlist cell for non-bonding interaction
integer,allocatable :: nbllist(:), nbheader(:,:,:), nbnacell(:,:,:)

!<nbrlist> neighbor list, <nbrindx> neighbor index
integer,allocatable :: nbrlist(:,:), nbrindx(:,:)

!<nbplist> neighbor list of nonbonding interaction, non-bonding pair list
integer,allocatable :: nbplist(:,:)

!--- coefficient of bonding energy derivative 
real(8),allocatable :: ccbnd(:), cdbnd(:)

!<BO> Bond Order of atoms i-j (nearest neighb only) - (Eq 3a-3d)
real(8),allocatable :: BO(:,:,:) 
real(8),allocatable :: delta(:)

!--- Output variables from the BOp_CALC() subroutine:
real(8),allocatable :: deltap(:,:)
real(8),allocatable :: dln_BOp(:,:,:)
real(8),allocatable :: dBOp(:,:)

!--- A[0123] coefficients for force calculation 
real(8),allocatable :: A0(:,:),A1(:,:), A2(:,:), A3(:,:) 

!--- Passed between Elnpr and E3body
real(8),allocatable :: nlp(:), dDlp(:) !Number of Lone Pairs, its derivatives.
real(8),allocatable :: deltalp(:)

! Two vectors electrostatic energy minimization 
real(8),allocatable,target :: qs(:),qt(:),gs(:), gt(:), hs(:), ht(:), hshs(:), hsht(:)

!--- variables for extended Lagrangian method ---
!<Lex_fqs> fraction between two QEq vectors
!<Lex_w> spring constant
real(8),allocatable,target :: qsfp(:),qsfv(:),qtfp(:),qtfv(:) 
real(8),allocatable :: hessian(:,:)
real(8) :: Lex_fqs=1.0, Lex_w=1.d0, Lex_w2=1.d0, Lex_k=2.d0


!--- Taper function 
real(8),parameter :: rctap0 = 10.d0 ![A]
real(8),parameter :: rctap0_pqeq = 12.5d0 ![A]
real(8) :: rctap, rctap2, CTap(0:7)

!--- PQEQ variables
real(8),allocatable,target :: spos(:,:)
integer,allocatable :: inxnpqeq(:,:)

real(8),allocatable :: TBL_Eclmb_pcc(:,:,:),TBL_Eclmb_psc(:,:,:),TBL_Eclmb_pss(:,:,:)

integer :: ntype_pqeq, ntype_pqeq2
logical,allocatable :: isPolarizable(:)
character(len=:),allocatable :: Elempqeq(:)
real(8),allocatable :: X0pqeq(:),J0pqeq(:)

real(8),allocatable :: Zpqeq(:), Rcpqeq(:), Rspqeq(:), Kspqeq(:)
real(8),allocatable :: alphacc(:,:), alphasc(:,:), alphass(:,:)
real(8) :: lambda_pqeq = 0.462770d0

real(8) :: eFieldStrength=0.0
integer :: eFieldDir=1

logical :: isPQEq = .false.
character(MAXSTRLENGTH) :: PQEqParmPath

contains

!------------------------------------------------------------------------------------------
subroutine reaxff_mdvariables_allocator()
use memory_allocator_mod
!------------------------------------------------------------------------------------------
implicit none

!--- extra lists for ReaxFF 
call allocator(nbrindx,1,NBUFFER,1,MAXNEIGHBS)
call allocator(nbplist,0,MAXNEIGHBS10,1,NBUFFER)

!--- Bond Order Prime and deriv terms:
call allocator(dln_BOp,1,3,1,NBUFFER,1,MAXNEIGHBS)
call allocator(dBOp,1,NBUFFER,1,MAXNEIGHBS)
call allocator(deltap,1,NBUFFER,1,3)
call allocator(deltalp,1,NBUFFER)

!--- Bond Order terms
call allocator(BO,0,3,1,NBUFFER,1,MAXNEIGHBS)
call allocator(delta,1,NBUFFER)
call allocator(A0,1,NBUFFER,1,MAXNEIGHBS)
call allocator(A1,1,NBUFFER,1,MAXNEIGHBS)
call allocator(A2,1,NBUFFER,1,MAXNEIGHBS)
call allocator(A3,1,NBUFFER,1,MAXNEIGHBS)
call allocator(nlp,1,NBUFFER)
call allocator(dDlp,1,NBUFFER)
call allocator(ccbnd,1,NBUFFER)
call allocator(cdbnd,1,NBUFFER)

!--- 2 vector QEq varialbes
call allocator(qs,1,NBUFFER)
call allocator(gs,1,NBUFFER)
call allocator(qt,1,NBUFFER)
call allocator(gt,1,NBUFFER)
call allocator(hs,1,NBUFFER)
call allocator(hshs,1,NBUFFER)
call allocator(ht,1,NBUFFER)
call allocator(hsht,1,NBUFFER)
call allocator(hessian,1,MAXNEIGHBS10,1,NBUFFER)

!--- Varaiable for extended Lagrangian method
call allocator(qtfp,1,NBUFFER)
call allocator(qtfv,1,NBUFFER)

return
end subroutine


!-----------------------------------------------------------------------------------------------------------------------
character(len=256) function rankToString(irank)
!-----------------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(in) :: irank
write(rankToString,*) irank
rankToString = adjustl(rankToString)
end function rankToString

!-------------------------------------------------------------------------------------------
function GetFileNameBase(DataDir,nstep) result(fileNameBase)
!-------------------------------------------------------------------------------------------
integer,intent(in) :: nstep
character(len=:),allocatable :: DataDir,fileNameBase
character(9) :: a9

if(nstep>=0) then
  write(a9,'(i9.9)') nstep
  fileNameBase=trim(adjustl(DataDir))//"/"//a9
else
  fileNameBase=trim(adjustl(DataDir))//"/rxff"
endif

end function

end module atoms

