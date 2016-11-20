!------------------------------------------------------------------------------------------
SUBROUTINE INITSYSTEM(NBUFFER, atype, pos, v, f, q)
! This subroutine takes care of setting up initial system configuration.
! Unit conversion of parameters (energy, length & mass) are also done here.
!------------------------------------------------------------------------------------------
!use parameters; use atoms; use ustruct
use parameters; use atoms
implicit none

integer,intent(in) :: NBUFFER
real(8),allocatable,dimension(:) :: atype, q
real(8),allocatable,dimension(:,:) :: pos,v,f

integer :: i,j,k, ix,iy,iz,ii, n, m, m3, p,s, inxn, ity, jty, l(3), sID, ist=0
real(8) :: dr, dr2, rr(3), mm, gmm, dns, mat(3,3)
integer(8) :: i8
real(8) :: rcsize(3), maxrcell, rcmesh2
integer :: imesh(3), maximesh

character(8) :: fname0
character(2) :: ctype
character(6) :: a6

character(64) :: argv

!--- read FF file, output dir, MD parameter file paths from command line
do i=1, command_argument_count()
   call get_command_argument(i,argv)
   select case(adjustl(argv))
     case("--help","-h")
       if(myid==0) print'(a)', "--ffield ffield --outDir DAT --rxmdin rxmd.in"
       stop
     case("--ffield", "-ff")
       call get_command_argument(i+1,argv)
       FFPath=adjustl(argv)
     case("--outDir", "-o")
       call get_command_argument(i+1,argv)
       DataDir=adjustl(argv)
     case("--rxmdin", "-in")
       call get_command_argument(i+1,argv)
       ParmPath=adjustl(argv)
     case default
   end select

enddo

!--- read MD control parameters
open(1, file=trim(ParmPath), status="old")
read(1,*) mdmode
read(1,*) dt, ntime_step
read(1,*) treq, vsfact, sstep
read(1,*) fstep, pstep
read(1,*) vprocs(1:3)
read(1,*) isQEq, NMAXQEq, QEq_tol, qstep
read(1,*) Lex_fqs, Lex_k
read(1,*) isBinary, isBondFile, isPDB
read(1,*) ftol
close(1)

!--- an error trap
if(vprocs(1)*vprocs(2)*vprocs(3) /= nprocs ) then
  if(myid==0) write(6,'(a60,3i3,i5)')  &
  "ERROR: requested/allocated # of procs are not consistent: ", vprocs(1:3), nprocs
  call MPI_FINALIZE(ierr)
  stop
endif

!--- initialize charge with QEq
if(mdmode==0) then
  if(myid==0) then
    print'(a,f12.3,a,i6,a)', &
         'INFO: mdmode==0, setting isQEQ is 1. Atomic velocities are scaled to ', &
          treq, ' [K] every ', sstep, ' steps.'
  endif
  isQEq=1
endif

!--- time unit conversion from [fs] -> time unit
dt = dt/UTIME

!--- square the spring const in the extended Lagrangian method 
Lex_w2=2.d0*Lex_k/dt/dt

!--- get reduced temperature from [K]
treq=treq/UTEMP0

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

!--- dt/2*mass, mass/2
allocate(dthm(nso),hmas(nso),stat=ast); ist=ist+ast
do ity=1, nso
   dthm(ity) = dt*0.5d0/mass(ity)
   hmas(ity) = 0.5d0*mass(ity)
enddo

!--- Particle Parameters
allocate(pos(3,NBUFFER),f(3,NBUFFER),v(3,NBUFFER),stat=ast); ist=ist+ast
allocate(atype(NBUFFER),q(NBUFFER), stat=ast); ist=ist+ast
f(:,:)=0.d0

!--- OpenMP private data array for ENbond
allocate(fnb(3,NBUFFER),stat=ast); ist=ist+ast
fnb(:,:)=0.d0

!--- Varaiable for extended Lagrangian method
allocate(qsfp(NBUFFER), qsfv(NBUFFER), stat=ast); ist=ist+ast
allocate(qtfp(NBUFFER), qtfv(NBUFFER), stat=ast); ist=ist+ast
qsfp(:)=0.d0; qsfv(:)=0.d0; qtfp(:)=0.d0; qtfv(:)=0.d0

call ReadBIN(NBUFFER, atype, pos, v, f, q)

call getbox(mat,lata,latb,latc,lalpha,lbeta,lgamma)
do i=1, 3
do j=1, 3
   HH(i,j,0)=mat(i,j)
enddo; enddo

!--- get total number of atoms per type. This will be used to determine
!--- subroutine cutofflength() 
allocate(natoms_per_type(nso),ibuf8(nso))
natoms_per_type(:)=0
do i=1, NATOMS
   ity=nint(atype(i))
   natoms_per_type(ity)=natoms_per_type(ity)+1
enddo

call MPI_ALLREDUCE(natoms_per_type, ibuf8, nso, MPI_INTEGER8, MPI_SUM,  MPI_COMM_WORLD, ierr)
natoms_per_type(1:nso)=ibuf8(1:nso)
deallocate(ibuf8)

!--- determine cutoff distances only for exsiting atom pairs
call CUTOFFLENGTH()

!--- update box-related variables
call updatebox()

!--- scaled to unscaled coordinates
call xs2xu(NBUFFER, pos)

!--- get global number of atoms
i8=NATOMS ! Convert 4 byte to 8 byte
call MPI_ALLREDUCE(i8, GNATOMS, 1, MPI_INTEGER8, MPI_SUM,  MPI_COMM_WORLD, ierr)

#ifdef STRESS
!--- stress variables
allocate(astr(6,NBUFFER),stat=ast); ist=ist+ast
astr(:,:)=0.d0; 
#endif

!--- Linked List & Near Neighb Parameters
allocate(nbrlist(NBUFFER,0:MAXNEIGHBS), nbrindx(NBUFFER, MAXNEIGHBS),stat=ast); ist=ist+ast

allocate(nbplist(NBUFFER,0:MAXNEIGHBS10),stat=ast); ist=ist+ast

allocate(llist(NBUFFER),stat=ast);ist=ist+ast 
allocate(header(-MAXLAYERS:cc(1)-1+MAXLAYERS, -MAXLAYERS:cc(2)-1+MAXLAYERS, -MAXLAYERS:cc(3)-1+MAXLAYERS), stat=ast); ist=ist+ast
allocate(nacell(-MAXLAYERS:cc(1)-1+MAXLAYERS, -MAXLAYERS:cc(2)-1+MAXLAYERS, -MAXLAYERS:cc(3)-1+MAXLAYERS), stat=ast); ist=ist+ast

!--- Bond Order Prime and deriv terms:
allocate(dln_BOp(3,NBUFFER, MAXNEIGHBS), dBOp(NBUFFER,MAXNEIGHBS), stat=ast); ist=ist+ast

allocate(deltap(NBUFFER, 3), stat=ast); ist=ist+ast

!--- Bond Order terms
allocate(BO(0:3,NBUFFER,MAXNEIGHBS), delta(NBUFFER), stat=ast); ist=ist+ast
allocate(A0(NBUFFER, MAXNEIGHBS), stat=ast); ist=ist+ast
allocate(A1(NBUFFER, MAXNEIGHBS), stat=ast); ist=ist+ast 
allocate(A2(NBUFFER, MAXNEIGHBS), stat=ast); ist=ist+ast 
allocate(A3(NBUFFER, MAXNEIGHBS), stat=ast); ist=ist+ast 

allocate(nlp(NBUFFER), dDlp(NBUFFER), stat=ast); ist=ist+ast

allocate(ccbnd(NBUFFER), stat=ast); ist=ist+ast
ccbnd(:)=0.d0

!--- 2 vector QEq varialbes
allocate(qs(NBUFFER), gs(NBUFFER), stat=ast); ist=ist+ast
allocate(qt(NBUFFER), gt(NBUFFER), stat=ast); ist=ist+ast
allocate(hs(NBUFFER), hshs(NBUFFER), stat=ast); ist=ist+ast
allocate(ht(NBUFFER), hsht(NBUFFER), stat=ast); ist=ist+ast
qs(:)=0.d0; qt(:)=0.d0; gs(:)=0.d0; gt(:)=0.d0; hs(:)=0.d0; ht(:)=0.d0; hshs(:)=0.d0; hsht(:)=0.d0

!--- returning force index array 
allocate(frcindx(NBUFFER), stat=ast); ist=ist+ast 

!--- Calculate constants used in loop (to avoid slow down later)
allocate(cBOp1(nboty), cBOp3(nboty), cBOp5(nboty), stat=ast); ist=ist+ast
allocate(pbo2h(nboty), pbo4h(nboty), pbo6h(nboty), stat=ast); ist=ist+ast

!--- setup potential table
call POTENTIALTABLE()

!--- <switch> flag to omit pi and double pi bond.
allocate(switch(1:3,nboty), stat=ast); ist=ist+ast

switch(:,:)=0
do i=1,nso
   do j=1,nso
   inxn = inxn2(i,j)

   if(inxn/=0) then

!!--- In BOp calculation, <switch> will be multiplied to <BOp> to remove
!!--- BOpi and BOpipi for bonding interaction of atoms with a hydrogen.
       if((rat(i)>0.d0)  .and. rat(j)>0.d0 )  switch(1,inxn)=1
       if((rapt(i)>0.d0) .and. rapt(j)>0.d0 ) switch(2,inxn)=1
       if((vnq(i)>0.d0)  .and. vnq(j)>0.d0 )  switch(3,inxn)=1

      if(r0s(i,j)<=0.d0) then 
         cBOp1(inxn) = 0.d0
      else
         cBOp1(inxn) = pbo1(inxn)/(r0s(i,j)**pbo2(inxn))
      endif
      if(r0p(i,j)<=0.d0) then
         cBOp3(inxn) = 0.d0
      else
         cBOp3(inxn) = pbo3(inxn)/(r0p(i,j)**pbo4(inxn))
      endif

      if(r0pp(i,j)<=0.d0) then
         cBOp5(inxn) = 0.d0
      else
         cBOp5(inxn) = pbo5(inxn)/(r0pp(i,j)**pbo6(inxn))
      endif
 
      pbo2h(inxn) = 0.5d0*pbo2(inxn)
      pbo4h(inxn) = 0.5d0*pbo4(inxn)
      pbo6h(inxn) = 0.5d0*pbo6(inxn)
   endif
   enddo
enddo

!--- get real size of linked list cell
rcsize(1) = lata/vprocs(1)/cc(1)
rcsize(2) = latb/vprocs(2)/cc(2)
rcsize(3) = latc/vprocs(3)/cc(3)
maxrcell = maxval(rcsize(1:3))

!--- setup 10[A] radius mesh to avoid visiting unecessary cells 
call GetNonbondingMesh(NBUFFER,pos)

!--- get density 
mm = 0.d0
do i=1, NATOMS
   ity = atype(i)
   mm = mm + mass(ity)
enddo
call MPI_ALLREDUCE (mm, gmm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
dns=gmm/MDBOX*UDENS

!--- allocate & initialize Array size Stat variables
i=ntime_step/pstep+1
allocate(maxas(i,nmaxas))
maxas(:,:)=0

!--- print out parameters and open data file
if(myid==0) then
   write(6,'(a)') "----------------------------------------------------------------"
   write(6,'(a30,i9,a3,i9)') "req/alloc # of procs:", vprocs(1)*vprocs(2)*vprocs(3), "  /",nprocs
   write(6,'(a30,3i9)')      "req proc arrengement:", vprocs(1),vprocs(2),vprocs(3)
   write(6,'(a30,a70)')      "parameter set:", pfile
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
   write(6,'(a30,2i6)') "MAXNEIGHBS, MAXNEIGHBS10:", MAXNEIGHBS,MAXNEIGHBS10
   write(6,'(a30,i6,i9)') "NMINCELL, NBUFFER:", NMINCELL, NBUFFER
   write(6,'(a30,3(a12,1x))') "FFPath, DataDir, ParmPath:", &
                          trim(FFPath), trim(DataDir), trim(ParmPath)

   print'(a30 $)','# of atoms per type:'
   do ity=1, nso
      if(natoms_per_type(ity)>0) print'(i12,a2,i2 $)',natoms_per_type(ity),' -',ity
   enddo
   print*

   print'(a)', "----------------------------------------------------------------"
   write(6,'(a)')  &
   "nstep  TE  PE  KE: 1-Ebond 2-(Elnpr,Eover,Eunder) 3-(Eval,Epen,Ecoa) 4-(Etors,Econj) 5-Ehbond 6-(Evdw,EClmb,Echarge)"

endif

!--- set initial time
wt0 = MPI_WTIME()

END SUBROUTINE

!------------------------------------------------------------------------------------------
SUBROUTINE INITVELOCITY(NBUFFER, atype, v)
use parameters; use atoms
! Generate gaussian distributed velocity as an initial value  using Box-Muller algorithm
!------------------------------------------------------------------------------------------
implicit none

integer,intent(in) :: NBUFFER
real(8) :: atype(NBUFFER)
real(8) :: v(3,NBUFFER)

integer :: i,j,k, ity
real(8) :: vv(2), vsqr, vsl, rndm(2)
real(8) :: vCM(3), GvCM(3), mm, Gmm
real(8) :: vfactor

!--- assign velocity to two atoms together with BM algoritm. 
!--- If <NATOMS> is odd, the <NATOMS> + 1 element will be the ignored in later calculations.

do i=1, NATOMS, 2

  do k=1,3 ! three directions
     !--- generate gaussian distributed velocity
     vsqr=0.d0
     do while ( (vsqr >= 1.d0) .or. (vsqr==0.d0) ) 
        call random_number(rndm)
        vv(1) = 2.d0 * rndm(1) - 1.d0
        vv(2) = 2.d0 * rndm(2) - 1.d0
        vsqr = vv(1)**2 + vv(2)**2
     enddo

     vsl = sqrt(-2.d0 * log(vsqr)/vsqr)
     v(k,i)   = vv(1)*vsl
     v(k,i+1) = vv(2)*vsl
  enddo
  
enddo

!--- get the local momentum and mass.
vCM(:)=0.d0;  mm = 0.d0
do i=1, NATOMS
   ity = atype(i)
   vCM(1:3)=vCM(1:3) + mass(ity)*v(1:3,i)
   mm = mm + mass(ity)
enddo
 
call MPI_ALLREDUCE (vCM, GvCM, size(vCM), MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE (mm, Gmm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)

!--- get the global momentum
GvCM(:)=GvCM(:)/Gmm

!--- set the total momentum to be zero and get the current kinetic energy. 
KE = 0.d0
do i=1, NATOMS
   v(1:3,i) = v(1:3,i) - GvCM(1:3)

   ity = atype(i)
   KE = KE + hmas(ity)*sum( v(1:3,i)*v(1:3,i) )
enddo

call MPI_ALLREDUCE (KE, GKE, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
GKE = GKE/GNATOMS

!--- scale the obtained velocity to get the initial temperature.
vfactor = sqrt(1.5d0*treq/GKE)
v(:,:) = vfactor * v(:,:)

end subroutine

!------------------------------------------------------------------------------------------
subroutine CUTOFFLENGTH()
use atoms; use parameters
!------------------------------------------------------------------------------------------
implicit none
integer :: i,j,nn,ity,jty,inxn
real(8) :: dr,BOsig

!--- get the cutoff length based on sigma bonding interaction.

! --- Remark --- 
! sigma bond before correction is the longest, namely longer than than pi and double-pi bonds
! thus check only sigma bond convergence.
allocate(rc(nboty), rc2(nboty), stat=ast)

rc(:)=0.d0; rc2(:)=0.d0
do ity=1,nso
do jty=ity, nso

   inxn=inxn2(ity,jty)
   if(inxn==0) cycle

   dr = 1.0d0
   BOsig=1.d0
   do while (BOsig > MINBOSIG) 
      dr = dr + 0.01d0
      BOsig = exp( pbo1(inxn)*(dr/r0s(ity,jty))**pbo2(inxn) ) !<- sigma bond prime
   enddo

   rc(inxn)  = dr
   rc2(inxn) = dr*dr
enddo
enddo

!----------------------------------------------------------------------------
! In some cases, an atom that do not exist in simulation gives
! the longest bond-order cutoff length. the check below is to ignore
! such atoms to keep the linkedlist cell dimensions as small as possible.
!----------------------------------------------------------------------------
do ity=1, nso
   if(natoms_per_type(ity)==0) then
      do jty=1, nso
         inxn=inxn2(ity,jty)
         if(inxn/=0) rc(inxn)=0.d0
         inxn=inxn2(jty,ity)
         if(inxn/=0) rc(inxn)=0.d0
      enddo
   endif
enddo

!--- get the max cutoff length 
maxrc=maxval(rc(:))

end subroutine

!------------------------------------------------------------------------------------------
subroutine POTENTIALTABLE()
use atoms; use parameters
!------------------------------------------------------------------------------------------
implicit none
integer :: i, n, ity,jty,inxn
real(8) :: dr1, dr2, dr3, dr4, dr5, dr6, dr7

real(8) :: exp1, exp2, fsum(3)
real(8) :: gamwinvp, gamWij, alphaij, Dij0, rvdW0
real(8) :: Tap, dTap, fn13, dfn13, dr3gamij, CEvdw, CEclmb, rij_vd1
real(8) :: rres

!--- first element in table 0: potential
!---                        1: derivative of potential
allocate(TBL_EClmb(0:1,NTABLE,nboty), TBL_Evdw(0:1,NTABLE, nboty), stat=ast)
allocate(TBL_EClmb_QEq(NTABLE,nboty), stat=ast)

!--- unit distance in r^2 scale
UDR = rctap2/NTABLE
UDRi = 1.d0/UDR

do ity=1, nso
do jty=ity, nso

   inxn = inxn2(ity,jty)
   if(inxn/=0) then
      do i=1, NTABLE

         dr2 = UDR*i
         dr1 = sqrt(dr2)

!--- Interaction Parameters:
         gamWij = gamW(ity,jty)
         alphaij = alpij(ity,jty)
         Dij0 = Dij(ity,jty)
         rvdW0 = rvdW(ity,jty) 
         gamwinvp = (1.d0/gamWij)**pvdW1

         dr3 = dr1*dr2
         dr4 = dr2*dr2
         dr5 = dr1*dr2*dr2
         dr6 = dr2*dr2*dr2
         dr7 = dr1*dr2*dr2*dr2 

         rij_vd1 = dr2**pvdW1h
         Tap = CTap(7)*dr7 + CTap(6)*dr6 + &
               CTap(5)*dr5 + CTap(4)*dr4 + CTap(0)
         fn13 = (rij_vd1 + gamwinvp)**pvdW1inv
         exp1 = exp( alphaij*(1.d0 - fn13 / rvdW0) )
         exp2 = sqrt(exp1)

         dr3gamij = ( dr3 + gamij(ity,jty) )**( -1.d0/3.d0 )

!!--- Energy Calculation:
!      PEvd = Tap*Dij0*(exp1 - 2d0*exp2)      
!      PEclmb = Tap*Cclmb*q(i)*q(j)*dr3gamij
!if(myid==0) print*,i, Tap*Dij0*(exp1 - 2d0*exp2), Tap*Cclmb*dr3gamij

          TBL_Evdw(0,i,inxn) = Tap*Dij0*(exp1 - 2d0*exp2)      
          TBL_Eclmb(0,i,inxn) = Tap*Cclmb*dr3gamij
          TBL_Eclmb_QEq(i,inxn) = Tap*Cclmb0_qeq*dr3gamij

!if(inxn==1.and.myid==0) print*,i,TBL_Evdw(i,inxn,0), TBL_Eclmb(i,inxn,0)

!--- Force Calculation:
         dTap = 7d0*CTap(7)*dr5 + 6d0*CTap(6)*dr4 + &
                5d0*CTap(5)*dr3 + 4d0*CTap(4)*dr2

         dfn13 = ((rij_vd1 + gamwinvp)**(pvdW1inv-1.d0)) * (dr2**(pvdW1h-1.d0)) 

!      CEvdw = Dij0*( dTap*(exp1 - 2.d0*exp2)  &
!           - Tap*(alphaij/rvdW0)*(exp1 - exp2)*dfn13 )
!      CEclmb = Cclmb*q(i)*q(j)*dr3gamij*( dTap - (dr3gamij**3)*Tap*dr(0) ) 

         TBL_Evdw(1,i,inxn) = Dij0*( dTap*(exp1 - 2.d0*exp2)  &
                            - Tap*(alphaij/rvdW0)*(exp1 - exp2)*dfn13 )
         TBL_Eclmb(1,i,inxn) = Cclmb*dr3gamij*( dTap - (dr3gamij**3)*Tap*dr1 ) 

      enddo
   endif

enddo
enddo

end subroutine

!----------------------------------------------------------------
subroutine GetNonbondingMesh(NBUFFER,pos)
use atoms; use parameters
! setup 10[A] radius mesh to avoid visiting unecessary cells 
!----------------------------------------------------------------
implicit none

integer,intent(in) :: NBUFFER
real(8) :: pos(3,NBUFFER)

integer :: i,j,k

real(8) :: latticePerNode(3), rr(3), dr2
real(8) :: rcsize(3),maxrcell,rcmesh2
integer :: ilc(3), imesh(3), maximesh 

!--- initial estimate of LL cell dims
nblcsize(1:3)=2.d0

!--- get mesh resolution which is close to the initial value of rlc.
latticePerNode(1)=lata/vprocs(1)
latticePerNode(2)=latb/vprocs(2)
latticePerNode(3)=latc/vprocs(3)
nbcc(1:3)=latticePerNode(1:3)/nblcsize(1:3)
nblcsize(1:3)=latticePerNode(1:3)/nbcc(1:3)
maxrcell = maxval(nblcsize(1:3))

!--- get # of linked list cell to cover up the non-bonding (10[A]) cutoff length
imesh(1:3)  = int(rctap/nblcsize(1:3)) + 1
maximesh = maxval(imesh(1:3))

!--- List up only cell indices within the cutoff range.
!--- pre-compute nmesh to get exact array size.
nbnmesh=0
do i=-imesh(1), imesh(1)
do j=-imesh(2), imesh(2)
do k=-imesh(3), imesh(3)
   rr(1:3) = (/i,j,k/)*nblcsize(1:3)
   dr2 = sum(rr(1:3)*rr(1:3))
   if(dr2 <= (rctap*1.2d0)**2) nbnmesh = nbnmesh + 1
   !if(dr2 <= rctap**2) nbnmesh = nbnmesh + 1
enddo; enddo; enddo

allocate(nbmesh(3,nbnmesh),stat=ast)

nbmesh(:,:)=0
nbnmesh=0
do i=-imesh(1), imesh(1)
do j=-imesh(2), imesh(2)
do k=-imesh(3), imesh(3)
   rr(1:3) = (/i,j,k/)*nblcsize(1:3)
   dr2 = sum(rr(1:3)*rr(1:3))
   if(dr2 <= (rctap*1.2d0)**2) then
   !if(dr2 <= rctap**2) then
      nbnmesh = nbnmesh + 1
      nbmesh(1:3,nbnmesh) = (/i, j, k/)
   endif
enddo; enddo; enddo

allocate(nbllist(NBUFFER),stat=ast)
allocate(nbheader( &
                -MAXLAYERS_NB:nbcc(1)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB:nbcc(2)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB:nbcc(3)-1+MAXLAYERS_NB), stat=ast)
allocate(nbnacell( &
                -MAXLAYERS_NB:nbcc(1)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB:nbcc(2)-1+MAXLAYERS_NB, &
                -MAXLAYERS_NB:nbcc(3)-1+MAXLAYERS_NB), stat=ast)

!--- normalize nblcsize, like lcsize.
nblcsize(1:3)=nblcsize(1:3)/(/lata,latb,latc/)

end subroutine

!----------------------------------------------------------------
subroutine getbox(H,la,lb,lc,angle1,angle2,angle3)
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
subroutine updatebox()
use atoms; use parameters
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
cc(1:3)=LBOX(1:3)/maxrc

!--- local system size in the unscaled coordinate.
LBOX(1:3) = 1.d0/vprocs(1:3)

!--- get the linkedlist cell dimensions (normalized)
lcsize(1:3) = LBOX(1:3)/cc(1:3)

!--- get origin of local MD box in the scaled coordiate.
OBOX(1:3) = LBOX(1:3)*vID(1:3)

return
end
