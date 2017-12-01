!------------------------------------------------------------------------------
subroutine PQEq(atype, pos, q)
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

real(8),intent(in) :: atype(NBUFFER), pos(NBUFFER,3)
real(8),intent(out) :: q(NBUFFER)

real(8) :: fpqeq(NBUFFER)
real(8) :: vdummy(1,1), fdummy(1,1)

integer :: i,j,l2g
integer :: i1,j1,k1, nmax
real(8) :: Gnew(2), Gold(2) 
real(8) :: Est, GEst1, GEst2, g_h(2), h_hsh(2)
real(4) :: lmin(2)
real(8) :: buf(4), Gbuf(4)
real(8) :: ssum, tsum, mu
real(8) :: qsum, gqsum
real(8) :: QCopyDr(3)

call system_clock(i1,k1)

QCopyDr(1:3)=rctap/(/lata,latb,latc/)

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
    qs(:)=0.d0
    qt(:)=0.d0
    qs(1:NATOMS)=q(1:NATOMS)
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

#ifdef QEQDUMP 
open(91,file="qeqdump"//trim(rankToString(myid))//".txt")
#endif

!--- copy atomic coords and types from neighbors, used in qeq_initialize()
call COPYATOMS(MODE_COPY, QCopyDr, atype, pos, vdummy, fdummy, q)
call LINKEDLIST(atype, pos, nblcsize, nbheader, nbllist, nbnacell, nbcc, MAXLAYERS_NB)

call qeq_initialize()

#ifdef QEQDUMP 
do i=1, NATOMS
   do j1=1,nbplist(i,0)
      j = nbplist(i,j1)
      write(91,'(4i6,4es25.15)') -1, l2g(atype(i)),nint(atype(i)),l2g(atype(j)),hessian(j1,i)
   enddo
enddo
#endif

!--- after the initialization, only the normalized coords are necessary for COPYATOMS()
!--- The atomic coords are converted back to real at the end of this function.
call COPYATOMS(MODE_QCOPY1,QCopyDr, atype, pos, vdummy, fdummy, q)
call get_gradient(Gnew)

!--- Let the initial CG direction be the initial gradient direction
hs(1:NATOMS) = gs(1:NATOMS)
ht(1:NATOMS) = gt(1:NATOMS)

call COPYATOMS(MODE_QCOPY2,QCopyDr, atype, pos, vdummy, fdummy, q)

GEst2=1.d99
do nstep_qeq=0, nmax-1


#ifdef QEQDUMP 
  qsum = sum(q(1:NATOMS))
  call MPI_ALLREDUCE(qsum, gqsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
#endif

  call get_hsh(Est)

  call MPI_ALLREDUCE(Est, GEst1, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)

#ifdef QEQDUMP 
  if(myid==0) print'(i5,5es25.15)', nstep_qeq, 0.5d0*log(Gnew(1:2)/GNATOMS), GEst1, GEst2, gqsum
#endif

  if( ( 0.5d0*( abs(GEst2) + abs(GEst1) ) < QEq_tol) ) exit 
  if( abs(GEst2) > 0.d0 .and. (abs(GEst1/GEst2-1.d0) < QEq_tol) ) exit
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
  call COPYATOMS(MODE_QCOPY1,QCopyDr, atype, pos, vdummy, fdummy, q)

!--- save old residues.  
  Gold(:) = Gnew(:)
  call get_gradient(Gnew)

!--- get new conjugate direction
  hs(1:NATOMS) = gs(1:NATOMS) + (Gnew(1)/Gold(1))*hs(1:NATOMS)
  ht(1:NATOMS) = gt(1:NATOMS) + (Gnew(2)/Gold(2))*ht(1:NATOMS)

!--- update new conjugate direction for buffered atoms.
  call COPYATOMS(MODE_QCOPY2,QCopyDr, atype, pos, vdummy, fdummy, q)

enddo

!--- for PQEq
call update_shell_positions()

call system_clock(j1,k1)
it_timer(1)=it_timer(1)+(j1-i1)

! save # of QEq iteration 
it_timer(24)=it_timer(24)+nstep_qeq

#ifdef QEQDUMP 
close(91)
#endif

return 

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------
subroutine update_shell_positions()
implicit none
!-----------------------------------------------------------------------------------------------------------------------
integer :: i,ity,j,jty,j1,inxn
real(8) :: shelli(3),shellj(3), qjc, clmb, dclmb, dr2
real(8) :: sforce(NATOMS,3), sf(3)

do i=1, NATOMS

   sforce(i,1:3)=0.d0
   ity = nint(atype(i))

   ! if i-atom is not polarizable, no force acting on i-shell. 
   if( isPolarizable(ity) ) cycle 

   sforce(i,1:3) = sforce(i,1:3) - Kspqeq(ity)*spos(i,1:3) ! Eq. (37)

   shelli(1:3) = pos(i,1:3) + spos(i,1:3)

   do j1 = 1, nbplist(i,0)

      j = nbplist(i,j1)
      jty = nint(atype(j))

      inxn = inxn2(ity, jty)

      qjc = q(j) + Zpqeq(jty)
      shellj(1:3) = pos(j,1:3) + spos(j,1:3)

      ! j-atom can be either polarizable or non-polarizable. In either case,
      ! there will be force on i-shell from j-core.  qjc takes care of the difference.  Eq. (38)
      !call get_coulomb_and_dcoulomb(shelli(1:3)-pos(j,1:3),inxn,sf)
      sforce(i,1:3)=sforce(i,1:3)-sf(1:3)*qjc*Zpqeq(ity)

      ! if j-atom is polarizable, there will be force on i-shell from j-shell. Eq. (38)
      if( isPolarizable(jty) ) then 
         !call get_coulomb_and_dcoulomb(shelli(1:3)-shellj(1:3),inxn,sf)
         sforce(i,1:3)=sforce(i,1:3)-sf(1:3)*Zpqeq(ity)*Zpqeq(jty)
      endif

   enddo

enddo

!--- update shell positions after finishing the shell-force calculation.  Eq. (39)
do i=1, NATOMS
   ity = nint(atype(i))
   if( isPolarizable(ity) ) spos(i,1:3) = spos(i,1:3) + sforce(i,1:3)/Kspqeq(ity)
enddo

end subroutine

!-----------------------------------------------------------------------------------------------------------------------
subroutine qeq_initialize()
use atoms; use parameters; use MemoryAllocator
! This subroutine create a neighbor list with cutoff length = 10[A] and save the hessian into <hessian>.  
! <nbrlist> and <hessian> will be used for different purpose later.
!-----------------------------------------------------------------------------------------------------------------------
implicit none
integer :: i,j, ity, jty, n, m, mn, nn
integer :: c1,c2,c3, c4,c5,c6
real(4) :: dr2
real(8) :: dr(3), drtb
real(8) :: alphaij, pqeqc, pqeqs, ff(3)
integer :: itb, inxn

integer :: ti,tj,tk

call system_clock(ti,tk)

nbplist(:,0) = 0

!$omp parallel do schedule(runtime), default(shared), &
!$omp private(i,j,ity,jty,n,m,mn,nn,c1,c2,c3,c4,c5,c6,dr,dr2,drtb,itb,inxn)
do c1=0, nbcc(1)-1
do c2=0, nbcc(2)-1
do c3=0, nbcc(3)-1

   i = nbheader(c1,c2,c3)
   do m = 1, nbnacell(c1,c2,c3)

   ity=nint(atype(i))

   fpqeq(i)=0.d0

   do mn = 1, nbnmesh
      c4 = c1 + nbmesh(1,mn)
      c5 = c2 + nbmesh(2,mn)
      c6 = c3 + nbmesh(3,mn)

      j = nbheader(c4,c5,c6)
      do n=1, nbnacell(c4,c5,c6)

         if(i/=j) then
            dr(1:3) = pos(i,1:3) - pos(j,1:3)
            dr2 =  sum(dr(1:3)*dr(1:3))

            if(dr2 < rctap2) then

               jty = nint(atype(j))

!--- make neighbor-list upto the taper function cutoff
!$omp atomic
               nbplist(i,0) = nbplist(i,0) + 1
               nbplist(i,nbplist(i,0)) = j

!--- get table index and residual value
               itb = int(dr2*UDRi)
               drtb = dr2 - itb*UDR
               drtb = drtb*UDRi

               inxn = inxn2(ity, jty)

!--- PEQq : 
!print'(a,4i6,i9,3es15.5)','i,j,ity,jty,itb,pqeqc',i,j,ity,jty,itb,sqrt(dr2),pqeqc,fpqeq(i)
               ! contribution from core(i)-core(j)
               call get_coulomb_and_dcoulomb_pqeq(dr,alphacc(ity,jty),pqeqc,ff)

               hessian(nbplist(i,0),i) = pqeqc

               ! contribution from C(r_icjc) and C(r_icjs) if j-atom is polarizable
               if( isPolarizable(jty) ) then ! pos(i,1:3)-(pos(j,1:3)+spos(j,1:3))  
                  call get_coulomb_and_dcoulomb_pqeq(dr-spos(j,1:3),alphasc(jty,ity),pqeqs,ff)

                  fpqeq(i) = fpqeq(i) + (pqeqc - pqeqs) * Zpqeq(jty) ! Eq. 30
               endif


            endif
         endif

         j=nbllist(j)
      enddo
   enddo !   do mn = 1, nbnmesh

   i=nbllist(i)
   enddo
enddo; enddo; enddo
!$omp end parallel do

!--- for array size stat
if(mod(nstep,pstep)==0) then
  nn=maxval(nbplist(1:NATOMS,0))
  i=nstep/pstep+1
  maxas(i,3)=nn
endif

call system_clock(tj,tk)
it_timer(16)=it_timer(16)+(tj-ti)

end subroutine 

!-----------------------------------------------------------------------------------------------------------------------
subroutine get_hsh(Est)
use atoms; use parameters
! This subroutine updates hessian*cg array <hsh> and the electrostatic energy <Est>.  
!-----------------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(OUT) :: Est
integer :: i,j,j1, ity, jty, inxn
real(8) :: eta_ity, Est1, dr2

real(8) :: Ccicj,Ccisj,Csicj,Csisj,shelli(3),shellj(3),qic,qjc,ff(3)

integer :: ti,tj,tk
call system_clock(ti,tk)

Est = 0.d0
!$omp parallel do default(shared), schedule(runtime), private(i,j,j1,ity,eta_ity,Est1),reduction(+:Est)
do i=1, NATOMS
   ity = nint(atype(i))
   eta_ity = eta(ity)

   hshs(i) = eta_ity*hs(i)
   hsht(i) = eta_ity*ht(i)

!--- for PQEq
   qic = q(i) + Zpqeq(ity)
   shelli(1:3) = pos(i,1:3) + spos(i,1:3)

   dr2 = sum(spos(i,1:3)*spos(i,1:3)) ! distance between core-and-shell for i-atom

   Est = Est + chi(ity)*q(i) + 0.5d0*eta_ity*q(i)*q(i) + 0.5d0*Kspqeq(ity)*dr2

   do j1 = 1, nbplist(i,0)
      j = nbplist(i,j1)
      jty = nint(atype(j))

!--- for PQEq
      qjc = q(j) + Zpqeq(jty)
      shellj(1:3) = pos(j,1:3) + spos(j,1:3)

      Ccicj = hessian(j1,i)
      call get_coulomb_and_dcoulomb_pqeq(shelli(1:3)-pos(j,1:3),alphasc(ity,jty),Csicj,ff)
      call get_coulomb_and_dcoulomb_pqeq(pos(i,1:3)-shellj(1:3),alphasc(jty,ity),Ccisj,ff)
      call get_coulomb_and_dcoulomb_pqeq(shelli(1:3)-shellj(1:3),alphass(ity,jty),Csisj,ff)

      hshs(i) = hshs(i) + hessian(j1,i)*hs(j)
      hsht(i) = hsht(i) + hessian(j1,i)*ht(j)

!--- get half of potential energy, then sum it up if atoms are resident.
      Est1 = Ccicj*qic*qjc - Ccisj*qic*Zpqeq(jty) - Csicj*Zpqeq(ity)*qjc + Csisj*Zpqeq(ity)*Zpqeq(jty)

      Est = Est + 0.5d0*Est1
      if(j<=NATOMS) Est = Est + 0.5d0*Est1
   enddo

enddo
!$omp end parallel do

call system_clock(tj,tk)
it_timer(18)=it_timer(18)+(tj-ti)

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

real(8) :: gssum, gtsum

integer :: ti,tj,tk
call system_clock(ti,tk)

!$omp parallel do default(shared), schedule(runtime), private(gssum, gtsum, eta_ity,i,j,j1,ity)
do i=1,NATOMS

   gssum=0.d0
   gtsum=0.d0
   do j1=1, nbplist(i,0) 
      j = nbplist(i,j1)
      gssum = gssum + hessian(j1,i)*qs(j)
      gtsum = gtsum + hessian(j1,i)*qt(j)
   enddo

   ity = nint(atype(i))
   eta_ity = eta(ity)

   gs(i) = - chi(ity) - eta_ity*qs(i) - gssum - fpqeq(i)
   gt(i) = - 1.d0     - eta_ity*qt(i) - gtsum

enddo 
!$omp end parallel do

ggnew(1) = dot_product(gs(1:NATOMS), gs(1:NATOMS))
ggnew(2) = dot_product(gt(1:NATOMS), gt(1:NATOMS))
call MPI_ALLREDUCE(ggnew, Gnew, size(ggnew), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

call system_clock(tj,tk)
it_timer(19)=it_timer(19)+(tj-ti)

end subroutine

end subroutine PQEq
