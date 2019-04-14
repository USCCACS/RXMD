!------------------------------------------------------------------------------------------
module init
!------------------------------------------------------------------------------------------

  use fileio

contains

!------------------------------------------------------------------------------------------
SUBROUTINE INITSYSTEM(atype, pos, v, f, q)
! This subroutine takes care of setting up initial system configuration.
! Unit conversion of parameters (energy, length & mass) are also done here.
!------------------------------------------------------------------------------------------
use qeq_mod
use pqeq_mod
use force_mod
use reaxff_param_mod
use velocity_modifiers_mod
use memory_allocator_mod

implicit none

real(8),allocatable,dimension(:) :: atype, q
real(8),allocatable,dimension(:,:) :: pos,v,f

integer :: i,j,k, ity, l(3), ist=0
real(8) :: mm, dns, mat(3,3)
integer(8) :: i8
real(8) :: rcsize(3), maxrcell

wt0 = MPI_WTIME()

!--- an error trap
if(vprocs(1)*vprocs(2)*vprocs(3) /= nprocs ) then
  if(myid==0) write(6,'(a60,3i3,i5)')  &
  "ERROR: requested/allocated # of procs are not consistent: ", vprocs(1:3), nprocs
  call MPI_FINALIZE(ierr)
  stop
endif

!--- setup the vector ID and parity for processes, in x, y and z order.
vID(1)=mod(myid,vprocs(1))
vID(2)=mod(myid/vprocs(1),vprocs(2))
vID(3)=myid/(vprocs(1)*vprocs(2))
myparity(1)=mod(vID(1),2)
myparity(2)=mod(vID(2),2)
myparity(3)=mod(vID(3),2)

k=0
do i=1,3
!--- Add (1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1) to vector ID applying PBC. 
   do j=1,-1,-2

!--- save vector ID.
      l(1:3) = vID(1:3)

!--- apply PBC.
      l(i)=mod( (vID(i) + j + vprocs(i)), vprocs(i) )
     
!--- get direction index.
      k=k+1 ![123456] 

!--- convert vector ID to sequential ID.
      target_node(k) = l(1)+l(2)*vprocs(1)+l(3)*vprocs(1)*vprocs(2) 
   enddo         
                   
enddo    

!--- read atom and MD box data
call allocator(atype,1,NBUFFER)
call allocator(q,1,NBUFFER)
call allocator(pos,1,NBUFFER,1,3)
call allocator(v,1,NBUFFER,1,3)
call allocator(f,1,NBUFFER,1,3)
!--- index array for returning reaction force
call allocator(frcindx,1,NBUFFER)


call ReadBIN(atype, pos, v, q, f, trim(DataDir)//"/rxff.bin")

!--- get global number of atoms
GNATOMS = NATOMS ! Convert 4 byte to 8 byte
call MPI_ALLREDUCE(MPI_IN_PLACE, GNATOMS, 1, MPI_INTEGER8, MPI_SUM,  MPI_COMM_WORLD, ierr)



!--- set force field parameters
get_forcefield_param => get_reaxff_param

!--- set force model
force_model_func => force_reaxff

!--- set potentialtable constructor
set_potentialtable => set_reaxff_potentialtables

!--- set charge model
charge_model_func => QEq
rctap = rctap0

if(isPQEq) then
  charge_model_func => PQEq
  rctap = rctap0_pqeq

  call allocator(spos,1,NBUFFER,1,3)

  call initialize_pqeq(chi,eta)
  if(isEfield) call initialize_eField(myid)

endif

!--- set taper function for the vdw and coulomb terms
rctap2 = rctap**2

CTap(0:7)=(/1.d0, 0.d0, 0.d0, 0.d0,   -35.d0/(rctap)**4, &
          84.d0/(rctap)**5, -70.d0/(rctap)**6, &
          20.d0/(rctap)**7 /)


!--- dt/2*mass, mass/2
call allocator(dthm, 1, size(mass))
call allocator(hmas, 1, size(mass))
do ity=1, size(mass)
   if(mass(ity) > 0.d0) then
      dthm(ity) = dt*0.5d0/mass(ity)
      hmas(ity) = 0.5d0*mass(ity)
   endif
enddo

!--- Varaiable for extended Lagrangian method
call allocator(qtfp,1,NBUFFER)
call allocator(qtfv,1,NBUFFER)

!--- get total number of atoms per type. This will be used to determine
!--- subroutine cutofflength() 
call allocator(natoms_per_type, 1, nso)
do i=1, NATOMS
   ity=nint(atype(i))
   natoms_per_type(ity)=natoms_per_type(ity)+1
enddo

call MPI_ALLREDUCE(MPI_IN_PLACE, natoms_per_type, size(natoms_per_type), &
                   MPI_INTEGER8, MPI_SUM,  MPI_COMM_WORLD, ierr)


!--- determine cutoff distances only for exsiting atom pairs. cutoff cleanup using natoms_per_type()
call get_bondorder_cutoff(rc, rc2, maxrc, natoms_per_type)

!--- setup potential table
call set_potentialtable()

!--- update box-related variables
call UpdateBoxParams()

!--- Linked List & Near Neighb Parameters
call allocator(nbrlist,1,NBUFFER,0,MAXNEIGHBS)
call allocator(llist,1,NBUFFER)
call allocator(header,-MAXLAYERS,cc(1)-1+MAXLAYERS,-MAXLAYERS,cc(2)-1+MAXLAYERS,-MAXLAYERS,cc(3)-1+MAXLAYERS)
call allocator(nacell,-MAXLAYERS,cc(1)-1+MAXLAYERS,-MAXLAYERS,cc(2)-1+MAXLAYERS,-MAXLAYERS,cc(3)-1+MAXLAYERS)



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

!--- setup 10[A] radius mesh to avoid visiting unecessary cells 
call GetNonbondingMesh()

!--- get density 
mm = 0.d0
do i=1, NATOMS
   ity = nint(atype(i))
   mm = mm + mass(ity)
enddo
call MPI_ALLREDUCE(MPI_IN_PLACE, mm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
dns=mm/MDBOX*UDENS

!--- allocate & initialize Array size Stat variables
call allocator(maxas, 1, ntime_step/pstep+1, 1, nmaxas)

!--- for spring force
if (isSpring) then
  call allocator(ipos,1,NBUFFER,1,3)
  ipos(1:NATOMS,1:3)=pos(1:NATOMS,1:3)
endif

!--- print out parameters and open data file
if(myid==0) then
   write(6,'(a)') "----------------------------------------------------------------"
   write(6,'(a30,i9,a3,i9)') "req/alloc # of procs:", vprocs(1)*vprocs(2)*vprocs(3), "  /",nprocs
   write(6,'(a30,3i9)')      "req proc arrengement:", vprocs(1),vprocs(2),vprocs(3)
   write(6,'(a30,es12.2)')   "time step[fs]:",dt*UTIME
   write(6,'(a30,i3, i10, i10)') "MDMODE CURRENTSTEP NTIMESTPE:", &
                                  mdmode, current_step, ntime_step
   write(6,'(a30,i6,es10.1,i6,i6)') "isQEq,QEq_tol,NMAXQEq,qstep:", &
                                     isQEq,QEq_tol,NMAXQEq,qstep
   write(6,'(a30,f8.3,f8.3)') 'Lex_fqs,Lex_k:',Lex_fqs,Lex_k
   write(6,'(a30,f12.3,f8.3,i9)') 'treq,vsfact,sstep:',treq*UTEMP0, vsfact, sstep
   write(6,'(a30,2i6)') 'fstep,pstep:', fstep,pstep
   write(6,'(a30,i24,i24)') "NATOMS GNATOMS:", NATOMS, GNATOMS
   write(6,'(a30,3f12.3)') "LBOX:",LBOX(1:3)
   write(6,'(a30,3f15.3)') "Hmatrix [A]:",HH(1:3,1,0)
   write(6,'(a30,3f15.3)') "Hmatrix [A]:",HH(1:3,2,0)
   write(6,'(a30,3f15.3)') "Hmatrix [A]:",HH(1:3,3,0)
   print'(a30,3f12.3)', 'lata,latb,latc:', lata,latb,latc
   print'(a30,3f12.3)', 'lalpha,lbeta,lgamma:', lalpha,lbeta,lgamma
   write(6,'(a30,3f10.4)') "density [g/cc]:",dns
   write(6,'(a30,3i6)')  '# of linkedlist cell:', cc(1:3)
   write(6,'(a30,f10.3,2x,3f10.2)') "maxrc, lcsize [A]:", &
   maxrc,lata/cc(1)/vprocs(1),latb/cc(2)/vprocs(2),latc/cc(3)/vprocs(3)
   write(6,'(a30,3i6)')  '# of linkedlist cell (NB):', nbcc(1:3)
   write(6,'(a30,3f10.2)') "lcsize [A] (NB):", &
   lata/nbcc(1)/vprocs(1),latb/nbcc(2)/vprocs(2),latc/nbcc(3)/vprocs(3)
   write(6,'(a30,2i6)') "MAXNEIGHBS, MAXNEIGHBS10:", MAXNEIGHBS,MAXNEIGHBS10
   write(6,'(a30,i6,i9)') "NMINCELL, NBUFFER:", NMINCELL, NBUFFER
   write(6,'(a30,3(a12,1x))') "FFPath, DataDir, ParmPath:", &
                          trim(FFPath), trim(DataDir), trim(ParmPath)

   print'(a30 $)','# of atoms per type:'
   do ity=1, nso
      if(natoms_per_type(ity)>0) print'(i12,a2,i2 $)',natoms_per_type(ity),' -',ity
   enddo
   print*

   if(isSpring) then
     print'(a,f8.2 $)','springConst [kcal/mol] ', springConst 
     do ity=1, size(hasSpringForce) 
        if(hasSpringForce(ity)) print'(i3 $)',ity
     enddo
     print*
   endif

   print'(a)', "----------------------------------------------------------------"
   write(6,'(a)')  &
   "nstep  TE  PE  KE: 1-Ebond 2-(Elnpr,Eover,Eunder) 3-(Eval,Epen,Ecoa) 4-(Etors,Econj) 5-Ehbond 6-(Evdw,EClmb,Echarge)"

endif

END SUBROUTINE
end module

!----------------------------------------------------------------
subroutine GetNonbondingMesh()
use atoms
use reaxff_param_mod
use memory_allocator_mod
! setup 10[A] radius mesh to avoid visiting unecessary cells 
!----------------------------------------------------------------
implicit none

integer :: i,j,k

real(8) :: latticePerNode(3), rr(3), dr2
real(8) :: maxrcell
integer :: imesh(3), maximesh, ii(3), i1

!--- initial estimate of LL cell dims
nblcsize(1:3)=3d0

!--- get mesh resolution which is close to the initial value of rlc.
latticePerNode(1)=lata/vprocs(1)
latticePerNode(2)=latb/vprocs(2)
latticePerNode(3)=latc/vprocs(3)
nbcc(1:3)=int(latticePerNode(1:3)/nblcsize(1:3))
nblcsize(1:3)=latticePerNode(1:3)/nbcc(1:3)
maxrcell = maxval(nblcsize(1:3))

!--- get # of linked list cell to cover up the non-bonding cutoff length
imesh(1:3)  = int(rctap/nblcsize(1:3)) + 1
maximesh = maxval(imesh(1:3))

!--- List up only cell indices within the cutoff range.
!--- pre-compute nmesh to get exact array size.
nbnmesh=0
do i=-imesh(1), imesh(1)
do j=-imesh(2), imesh(2)
do k=-imesh(3), imesh(3)
   ii(1:3) = [i,j,k]
   do i1 = 1, 3
      if(ii(i1)>0) then
         ii(i1)=ii(i1)-1
      else if(ii(i1)<0) then
         ii(i1)=ii(i1)+1
      endif
   enddo
   rr(1:3) = ii(1:3)*nblcsize(1:3)
   dr2 = sum(rr(1:3)*rr(1:3))
   if(dr2 <= rctap**2) nbnmesh = nbnmesh + 1
enddo; enddo; enddo

call allocator(nbmesh,1,3,1,nbnmesh)

nbmesh(:,:)=0
nbnmesh=0
do i=-imesh(1), imesh(1)
do j=-imesh(2), imesh(2)
do k=-imesh(3), imesh(3)
   ii(1:3) = [i,j,k]
   do i1 = 1, 3
      if(ii(i1)>0) then
         ii(i1)=ii(i1)-1
      else if(ii(i1)<0) then
         ii(i1)=ii(i1)+1
      endif
   enddo
   rr(1:3) = ii(1:3)*nblcsize(1:3)
   dr2 = sum(rr(1:3)*rr(1:3))
   if(dr2 <= rctap**2) then
      nbnmesh = nbnmesh + 1
      nbmesh(1:3,nbnmesh) = (/i, j, k/)
   endif
enddo; enddo; enddo

call allocator(nbllist,1,NBUFFER)
call allocator(nbheader, &
                -MAXLAYERS_NB,nbcc(1)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,nbcc(2)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,nbcc(3)-1+MAXLAYERS_NB)
call allocator(nbnacell, &
                -MAXLAYERS_NB,nbcc(1)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,nbcc(2)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB,nbcc(3)-1+MAXLAYERS_NB)

!--- normalize nblcsize, like lcsize.
nblcsize(1:3)=nblcsize(1:3)/(/lata,latb,latc/)

end subroutine

!----------------------------------------------------------------
subroutine GetBoxParams(H,la,lb,lc,angle1,angle2,angle3)
!----------------------------------------------------------------
implicit none
real(8),intent(out) :: H(3,3)
real(8),intent(in) :: la,lb,lc, angle1,angle2,angle3
real(8) :: hh1, hh2 , lal, lbe, lga
real(8) :: pi=atan(1.d0)*4.d0

!--- convet unit for angles
lal=angle1*pi/180.d0
lbe=angle2*pi/180.d0
lga=angle3*pi/180.d0

!--- construct H-matrix
hh1=lc*(cos(lal)-cos(lbe)*cos(lga))/sin(lga)
hh2=lc*sqrt( 1.d0-cos(lal)**2-cos(lbe)**2-cos(lga)**2 + &
             2*cos(lal)*cos(lbe)*cos(lga) )/sin(lga)

H(1,1)=la;          H(2,1)=0.d0;        H(3,1)=0.d0
H(1,2)=lb*cos(lga); H(2,2)=lb*sin(lga); H(3,2)=0.d0
H(1,3)=lc*cos(lbe); H(2,3)=hh1;         H(3,3)=hh2

return
end subroutine

!----------------------------------------------------------------
subroutine UpdateBoxParams()
use atoms
!----------------------------------------------------------------
implicit none

!--- get volume 
MDBOX = &
HH(1,1,0)*(HH(2,2,0)*HH(3,3,0) - HH(3,2,0)*HH(2,3,0)) + &
HH(2,1,0)*(HH(3,2,0)*HH(1,3,0) - HH(1,2,0)*HH(3,3,0)) + &
HH(3,1,0)*(HH(1,2,0)*HH(2,3,0) - HH(2,2,0)*HH(1,3,0))

!--- get inverse of H-matrix
call matinv(HH,HHi)

!--- local box dimensions (a temporary use of lbox)
LBOX(1)=lata/vprocs(1)
LBOX(2)=latb/vprocs(2)
LBOX(3)=latc/vprocs(3)

!--- get the number of linkedlist cell per domain
cc(1:3)=int(LBOX(1:3)/maxrc)

!--- local system size in the unscaled coordinate.
LBOX(1:3) = 1.d0/vprocs(1:3)

!--- get the linkedlist cell dimensions (normalized)
lcsize(1:3) = LBOX(1:3)/cc(1:3)

!--- get origin of local MD box in the scaled coordiate.
OBOX(1:3) = LBOX(1:3)*vID(1:3)

return
end
