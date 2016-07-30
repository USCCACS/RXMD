!------------------------------------------------------------------------------
subroutine QEq()
use atoms; use parameters
! Two vector electronegativity equilization routine
!
! The linkedlist cell size is determined by the cutoff length of bonding 
! interaction <rc> = 3A. Since the non-bonding interaction cutoff <Rcut> = 10A,
! need to take enough layers to calculate non-bonding interactoins.
!
!<Gnew>, <Gold> :: NEW and OLD squre norm of Gradient vector.
!<Est> :: ElectroSTatic energy
!-------------------------------------------------------------------------------
implicit none
integer :: i1,j1,k1, nmax
real(8) :: Gnew(2), Gold(2) 
real(8) :: Est, GEst1, GEst2,lmin(2), g_h(2), h_hsh(2)
real(8) :: buf(4), Gbuf(4)
real(8) :: ssum, tsum, Gssum, Gtsum, qsum, gqsum, mu
real(8) :: qwtime
real(8) :: QCopyDr(3)

call system_clock(i1,k1)

QCopyDr(1:3)=10d0/(/lata,latb,latc/)

!--- Initialize <s> vector with current charge and <t> vector with zero.
!--- isQEq==1 Normal QEq, isQEq==2 Extended Lagrangian method, DEFAULT skip QEq 
select case(isQEq)

!=== original QEq ===!
  case (1) 
!--- In the original QEq, fictitious charges are initialized with real charges
!--- and set zero.
    qsfp(1:NATOMS)=q(1:NATOMS)
    qsfv(1:NATOMS)=0.d0
!--- Initialization of the two vector QEq 
    qs(1:NATOMS)=q(1:NATOMS)
    qt(1:NATOMS)=0.d0
    nmax=NMAXQEq

!=== Extended Lagrangian method ===!
  case(2)
!--- charge mixing.
    qs(1:NATOMS)=Lex_fqs*qsfp(1:NATOMS)+(1.d0-Lex_fqs)*q(1:NATOMS)
!--- the same as the original QEq, set t vector zero
    qt(1:NATOMS)=0.d0
!--- just run one step
    nmax=1

!=== else, just return ===!
  case default
     return

end select

!--- copy atomic coords and types from neighbors, used in qeq_initialize()
call xu2xs()
call COPYATOMS(MODE_COPY, QCopyDr)
call LINKEDLIST()
call NBLINKEDLIST()
call xs2xu()

call qeq_initialize()

!--- after the initialization, only the normalized coords are necessary for COPYATOMS()
!--- The atomic coords are converted back to real at the end of this function.
call xu2xs()
call COPYATOMS(MODE_QCOPY1,QCopyDr)
call get_gradient(Gnew)

!--- Let the initial CG direction be the initial gradient direction
hs(1:NATOMS) = gs(1:NATOMS)
ht(1:NATOMS) = gt(1:NATOMS)

call COPYATOMS(MODE_QCOPY2,QCopyDr)

GEst2=1.d99
!do nstep_qeq=0, NMAXQEq-1
do nstep_qeq=0, nmax-1

!  qsum = sum(q(1:NATOMS))
!  call MPI_ALLREDUCE(qsum, gqsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)

  call get_hsh(Est)

  call MPI_ALLREDUCE(Est, GEst1, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
  !if(myid==0) print'(i5,4es25.15)',nstep_qeq, 0.5d0*log(Gnew(1:2)/NATOMS), GEst1, GEst2

  if( ( 0.5d0*(abs(GEst2)+abs(GEst1))<QEq_tol) .or. (abs(GEst1/GEst2-1.d0) < QEq_tol) ) exit
  GEst2 = GEst1

!--- line minimization factor of <s> vector
  g_h(1) = dot_product(gs(1:NATOMS), hs(1:NATOMS))
  h_hsh(1) = dot_product(hs(1:NATOMS), hshs(1:NATOMS))

!--- line minimization factor of <t> vector
  g_h(2) = dot_product(gt(1:NATOMS), ht(1:NATOMS))
  h_hsh(2) = dot_product(ht(1:NATOMS), hsht(1:NATOMS))

  buf(1)=g_h(1);   buf(2)=g_h(2)
  buf(3)=h_hsh(1); buf(4)=h_hsh(2)
  Gbuf(:)=0.d0
  call MPI_ALLREDUCE(buf, Gbuf, 4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
  g_h(1) = Gbuf(1);   g_h(2) = Gbuf(2)
  h_hsh(1) = Gbuf(3); h_hsh(2) = Gbuf(4)

  lmin(1:2) = g_h(1:2)/h_hsh(1:2)

!--- line minimization for each vector
  qs(1:NATOMS) = qs(1:NATOMS) + lmin(1)*hs(1:NATOMS)
  qt(1:NATOMS) = qt(1:NATOMS) + lmin(2)*ht(1:NATOMS)

!--- get a current electronegativity <mu>
  ssum = sum(qs(1:NATOMS))
  tsum = sum(qt(1:NATOMS))
  buf(1) = ssum; buf(2) = tsum

  Gbuf(:)=0.d0
  call MPI_ALLREDUCE(buf, Gbuf, 2, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
  ssum=Gbuf(1); tsum=Gbuf(2)

  mu = ssum/tsum

!--- update atom charges
  q(1:NATOMS) = qs(1:NATOMS) - mu*qt(1:NATOMS)

!--- update new charges of buffered atoms.
  call COPYATOMS(MODE_QCOPY1,QCopyDr)

!--- save old residues.  
  Gold(:) = Gnew(:)
  call get_gradient(Gnew)

!--- get new conjugate direction
  hs(1:NATOMS) = gs(1:NATOMS) + (Gnew(1)/Gold(1))*hs(1:NATOMS)
  ht(1:NATOMS) = gt(1:NATOMS) + (Gnew(2)/Gold(2))*ht(1:NATOMS)

!--- update new conjugate direction for buffered atoms.
  call COPYATOMS(MODE_QCOPY2,QCopyDr)

enddo

call qeq_finalize()
call xs2xu()

call system_clock(j1,k1)
it_timer(1)=it_timer(1)+(j1-i1)

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------
subroutine qeq_initialize()
use atoms; use parameters
! This subroutine create a neighbor list with cutoff length = 10[A] and save the hessian into <A0>.  
! <nbrlist> and <A0> will be used for different purpose later.
!-----------------------------------------------------------------------------------------------------------------------
implicit none
integer :: i,j,k, ity, jty, n, m, mn, nn,nn1, ist=0
integer :: c1,c2,c3, c4,c5,c6, ic(3)
real(8) :: dr(3), dr1,dr2,dr3,dr4, Tap, CEst, eta_ity, hsan
real(8) :: CTdr4, CTdr5, CTdr6, CTdr7, dr3gt
real(8) :: third = -1.d0/3.d0
real(8) :: drtb, drtb1
integer :: itb, itb1, inxn, l2g

deallocate(BO,A0,A1,A2,A3,nbrlist,nbrindx,stat=ast); ist=ist+ast
deallocate(dln_BOp,dBOp,stat=ast); ist=ist+ast
allocate(nbrlist(0:MAXNEIGHBS10,NATOMS),stat=ast); ist=ist+ast
allocate(A0(MAXNEIGHBS10,NATOMS),stat=ast); ist=ist+ast

if(ist/=0) then
   print*,'Error @ qeq_initialize: ', myid
   call MPI_FINALIZE(ierr)
   stop
endif

nn=0
nbrlist(0,:) = 0
do c1=0, nbcc(1)-1
do c2=0, nbcc(2)-1
do c3=0, nbcc(3)-1

   i = nbheader(c1,c2,c3)
   do m = 1, nbnacell(c1,c2,c3)

   ity=atype(i)

   do mn = 1, nbnmesh
      c4 = c1 + nbmesh(1,mn)
      c5 = c2 + nbmesh(2,mn)
      c6 = c3 + nbmesh(3,mn)

      j = nbheader(c4,c5,c6)
      do n=1, nbnacell(c4,c5,c6)

         if(i/=j) then
            dr(1:3) = pos(1:3,i) - pos(1:3,j)
            dr2 =  sum(dr(1:3)*dr(1:3))

            if(dr2 < rctap2) then

               jty = atype(j)

!--- make a neighbor list with cutoff length = 10[A]
               nbrlist(0,i) = nbrlist(0,i) + 1
               nbrlist(nbrlist(0,i),i) = j


!--- get table index and residual value
               itb = int(dr2*UDRi)
               drtb = dr2 - itb*UDR
               drtb = drtb*UDRi

!--- save hessian into A0
               inxn = inxn2(ity, jty)
#ifdef DEBUG
if(inxn==0) print'(a,4i9)','inxn==0,myid,inxn,ity,jty: ',myid,inxn,ity,jty
if(itb==0) print'(5i,6f10.5)',myid,l2g(atype(i)),l2g(atype(j)),i,j,pos(1:3,i), pos(1:3,j)
#endif 
               hsan = (1.d0-drtb)*TBL_Eclmb_QEq(itb,inxn) + drtb*TBL_Eclmb_QEq(itb+1,inxn)
               A0(nbrlist(0,i),i) = hsan
            endif
         endif

         j=nbllist(j)
      enddo
   enddo !   do mn = 1, nbnmesh

   i=nbllist(i)
   enddo
enddo; enddo; enddo

!--- for array size stat
if(mod(nstep,pstep)==0) then
  nn=maxval(nbrlist(0,1:NATOMS))
  i=nstep/pstep+1
  maxas(i,3)=nn
endif

end subroutine 

!-----------------------------------------------------------------------------------------------------------------------
subroutine qeq_finalize()
use atoms
integer :: l2g
integer :: iast
!-----------------------------------------------------------------------------------------------------------------------
deallocate(A0, nbrlist,stat=ast)

iast=0
allocate(BO(0:3,NBUFFER_P, MAXNEIGHBS), stat=ast); iast=iast+ast
allocate(A0(NBUFFER_P, MAXNEIGHBS), stat=ast); iast=iast+ast
allocate(A1(NBUFFER_P, MAXNEIGHBS), stat=ast); iast=iast+ast
allocate(A2(NBUFFER_P, MAXNEIGHBS), stat=ast); iast=iast+ast
allocate(A3(NBUFFER_P, MAXNEIGHBS), stat=ast); iast=iast+ast
allocate(dln_BOp(3,NBUFFER_P, MAXNEIGHBS), dBOp(NBUFFER_P,MAXNEIGHBS), stat=ast); iast=iast+ast

allocate(nbrlist(NBUFFER_P,-1:MAXNEIGHBS), stat=ast); iast=iast+ast
allocate(nbrindx(NBUFFER_P,-1:MAXNEIGHBS), stat=ast); iast=iast+ast

if (iast/=0) then
   if (myid==0) print*, 'ERROR: qeq_finalize', iast
   call MPI_FINALIZE(ierr)
   stop 
endif

end subroutine

!-----------------------------------------------------------------------------------------------------------------------
subroutine get_hsh(Est)
use atoms; use parameters
! This subroutine updates hessian*cg array <hsh> and the electrostatic energy <Est>.  
!-----------------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(OUT) :: Est
integer :: i,j,j1, ity
real(8) :: eta_ity, Est1

Est = 0.d0
do i=1,NATOMS
   ity = atype(i)
   eta_ity = eta(ity)

   hshs(i) = eta_ity*hs(i)
   hsht(i) = eta_ity*ht(i)

   Est = Est + chi(ity)*q(i) + 0.5d0*eta_ity*q(i)*q(i)

   do j1 = 1, nbrlist(0,i)
      j = nbrlist(j1,i)
      hshs(i) = hshs(i) + A0(j1,i)*hs(j)
      hsht(i) = hsht(i) + A0(j1,i)*ht(j)
!--- get half of potential energy, then sum it up if atoms are resident.
      Est1 = 0.5d0*A0(j1,i)*q(i)*q(j)
      Est = Est + Est1
      if(j>0) Est = Est + Est1
   enddo

enddo

end subroutine 

!-----------------------------------------------------------------------------------------------------------------------
subroutine get_gradient(Gnew)
! Update gradient vector <g> and new residue <Gnew>
!-----------------------------------------------------------------------------------------------------------------------
use atoms; use parameters
implicit none
real(8),intent(OUT) :: Gnew(2)
real(8) :: eta_ity, ggnew(2)
integer :: i,j,j1, ity

do i=1,NATOMS
   ity = atype(i)
   eta_ity = eta(ity)

!--- Initialize a gradient vector
   gs(i) = - chi(ity) - eta_ity*qs(i)
   gt(i) = - 1.d0     - eta_ity*qt(i)

   do j1=1, nbrlist(0,i) 
      j = nbrlist(j1,i)
      gs(i) = gs(i) - A0(j1,i)*qs(j)
      gt(i) = gt(i) - A0(j1,i)*qt(j)
   enddo

enddo 

ggnew(1) = dot_product(gs(1:NATOMS), gs(1:NATOMS))
ggnew(2) = dot_product(gt(1:NATOMS), gt(1:NATOMS))
call MPI_ALLREDUCE(ggnew, Gnew, size(ggnew), MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)

end subroutine

end subroutine QEq
