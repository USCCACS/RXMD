module pqeq_mod

  use atoms
  use reaxff_param_mod
  use memory_allocator_mod
  use lists_mod

contains
!------------------------------------------------------------------------------
subroutine PQEq(atype, pos, q)
  use communication_mod
!use atoms
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

real(8),allocatable,intent(in out) :: atype(:), pos(:,:)
real(8),allocatable,intent(in out) :: q(:)

real(8) :: fpqeq(NBUFFER)
real(8),allocatable :: vdummy(:,:), fdummy(:,:) 

integer :: i,j
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
call LINKEDLIST(atype, pos, nblcsize, nbheader, nbllist, nbnacell)

call qeq_initialize()

#ifdef QEQDUMP 
do i=1, NATOMS
   do j1=1,nbplist(0,i)
      j = nbplist(j1,i)
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
  call MPI_ALLREDUCE(MPI_IN_PLACE, qsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
  gqsum = qsum
#endif

  call get_hsh(Est)

  call MPI_ALLREDUCE(MPI_IN_PLACE, Est, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
  GEst1 = Est

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
  call MPI_ALLREDUCE(MPI_IN_PLACE, buf, 4, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
  g_h(1) = buf(1);   g_h(2) = buf(2)
  h_hsh(1) = buf(3); h_hsh(2) = buf(4)

  lmin(1:2) = g_h(1:2)/h_hsh(1:2)

!--- line minimization for each vector
  qs(1:NATOMS) = qs(1:NATOMS) + lmin(1)*hs(1:NATOMS)
  qt(1:NATOMS) = qt(1:NATOMS) + lmin(2)*ht(1:NATOMS)

!--- get a current electronegativity <mu>
  ssum = sum(qs(1:NATOMS))
  tsum = sum(qt(1:NATOMS))
  buf(1) = ssum; buf(2) = tsum

  call MPI_ALLREDUCE(MPI_IN_PLACE, buf, 2, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, ierr)
  ssum=buf(1); tsum=buf(2)

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
real(8),parameter :: MAX_SHELL_DISPLACEMENT=1d-3

integer :: i,ity,j,jty,j1,inxn
real(8) :: shelli(3),shellj(3), qjc, clmb, dclmb, ddr
real(8) :: sforce(NATOMS,3), sf(3), Esc, Ess
real(8) :: ff(3), dr(3)

sforce(1:NATOMS,1:3)=0.d0
do i=1, NATOMS

   ity = nint(atype(i))

   ! if i-atom is not polarizable, no force acting on i-shell. 
   if( .not. isPolarizable(ity) ) cycle 

   if(isEfield) sforce(i,eFieldDir) = sforce(i,eFieldDir) - Zpqeq(ity)*eFieldStrength*Eev_kcal

   sforce(i,1:3) = sforce(i,1:3) - Kspqeq(ity)*spos(i,1:3) ! Eq. (37)
   shelli(1:3) = pos(i,1:3) + spos(i,1:3)

   do j1 = 1, nbplist(0,i)

      j = nbplist(j1,i)
      jty = nint(atype(j))

      qjc = q(j) + Zpqeq(jty)
      shellj(1:3) = pos(j,1:3) + spos(j,1:3)

      ! j-atom can be either polarizable or non-polarizable. In either case,
      ! there will be force on i-shell from j-core.  qjc takes care of the difference.  Eq. (38)
      dr(1:3)=shelli(1:3)-pos(j,1:3)
      call get_coulomb_and_dcoulomb_pqeq(dr,alphasc(ity,jty),Esc, inxnpqeq(ity, jty), TBL_Eclmb_psc,sf)

      ff(1:3)=-Cclmb0*sf(1:3)*qjc*Zpqeq(ity)
      sforce(i,1:3)=sforce(i,1:3)-ff(1:3)

      ! if j-atom is polarizable, there will be force on i-shell from j-shell. Eq. (38)
      if( isPolarizable(jty) ) then 
         dr(1:3)=shelli(1:3)-shellj(1:3)
         call get_coulomb_and_dcoulomb_pqeq(dr,alphass(ity,jty),Ess, inxnpqeq(ity, jty), TBL_Eclmb_pss,sf)

         ff(1:3)=Cclmb0*sf(1:3)*Zpqeq(ity)*Zpqeq(jty)
         sforce(i,1:3)=sforce(i,1:3)-ff(1:3)

      endif

   enddo

enddo

!--- update shell positions after finishing the shell-force calculation.  Eq. (39)
do i=1, NATOMS

   ity = nint(atype(i))

   dr(1:3)=sforce(i,1:3)/Kspqeq(ity)
   ddr = sqrt(sum(dr(1:3)*dr(1:3)))

   ! check the shell displacement per MD step for stability
   if(ddr>MAX_SHELL_DISPLACEMENT) then
      !print'(a,i6,i9,f12.6)', &
      !     '[WARNING] large shell displacement found : myid,i,ddr : ', myid, i, ddr 
      dr(1:3) = dr(1:3)/ddr*MAX_SHELL_DISPLACEMENT
   endif

   if( isPolarizable(ity) ) spos(i,1:3) = spos(i,1:3) + dr(1:3)
enddo


end subroutine

!-----------------------------------------------------------------------------------------------------------------------
subroutine qeq_initialize()
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

nbplist(0,:) = 0

!$omp parallel do schedule(runtime), default(shared), &
!$omp private(i,j,ity,jty,n,m,mn,nn,c1,c2,c3,c4,c5,c6,dr,dr2,drtb,itb,inxn,pqeqc,pqeqs,ff)
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
               nbplist(0,i) = nbplist(0,i) + 1
               nbplist(nbplist(0,i),i) = j

!--- get table index and residual value
               itb = int(dr2*UDRi)
               drtb = dr2 - itb*UDR
               drtb = drtb*UDRi

!--- PEQq : 
               ! contribution from core(i)-core(j)
               call get_coulomb_and_dcoulomb_pqeq(dr,alphacc(ity,jty),pqeqc,inxnpqeq(ity, jty),TBL_Eclmb_pcc,ff)

               hessian(nbplist(0,i),i) = Cclmb0_qeq * pqeqc

               fpqeq(i) = fpqeq(i) + Cclmb0_qeq * pqeqc * Zpqeq(jty) ! Eq. 30

               ! contribution from C(r_icjc) and C(r_icjs) if j-atom is polarizable
               if( isPolarizable(jty) ) then 
                  dr(1:3)=pos(i,1:3) - pos(j,1:3) - spos(j,1:3) ! pos(i,1:3)-(pos(j,1:3)+spos(j,1:3))  
                  call get_coulomb_and_dcoulomb_pqeq(dr,alphasc(jty,ity),pqeqs,inxnpqeq(jty, ity),TBL_Eclmb_psc,ff)

                  fpqeq(i) = fpqeq(i) - Cclmb0_qeq * pqeqs * Zpqeq(jty) ! Eq. 30
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
  nn=maxval(nbplist(0,1:NATOMS))
  i=nstep/pstep+1
  maxas(i,3)=nn
endif

call system_clock(tj,tk)
it_timer(16)=it_timer(16)+(tj-ti)

end subroutine 

!-----------------------------------------------------------------------------------------------------------------------
subroutine get_hsh(Est)
! This subroutine updates hessian*cg array <hsh> and the electrostatic energy <Est>.  
!-----------------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(OUT) :: Est
integer :: i,j,j1, ity, jty, inxn
real(8) :: eta_ity, Est1, dr2, dr(3)

real(8) :: Ccicj,Csicj,Csisj,shelli(3),shellj(3),qic,qjc,ff(3)
real(8) :: Eshell

integer :: ti,tj,tk
call system_clock(ti,tk)

Est = 0.d0
!$omp parallel do default(shared), reduction(+:Est) &
!$omp private(i,j,j1,ity,jty,eta_ity,Est1,Eshell,Ccicj,Csicj,Csisj,shelli,shellj,qic,qjc,ff,dr,dr2)
do i=1, NATOMS
   ity = nint(atype(i))
   eta_ity = eta(ity)

   hshs(i) = eta_ity*hs(i)
   hsht(i) = eta_ity*ht(i)

!--- for PQEq
   qic = q(i) + Zpqeq(ity)
   shelli(1:3) = pos(i,1:3) + spos(i,1:3)

   dr2 = sum(spos(i,1:3)*spos(i,1:3)) ! distance between core-and-shell for i-atom

   Est = Est + chi(ity)*q(i) + 0.5d0*eta_ity*q(i)*q(i)

   do j1 = 1, nbplist(0,i)
      j = nbplist(j1,i)
      jty = nint(atype(j))

!--- for PQEq
      qjc = q(j) + Zpqeq(jty)
      shellj(1:3) = pos(j,1:3) + spos(j,1:3)

      Ccicj = 0.d0; Csicj=0.d0; Csisj=0.d0

      Ccicj = hessian(j1,i)*qic*qjc ! hessian() is in [eV]

      if(isPolarizable(ity)) then
         dr(1:3)=shelli(1:3)-pos(j,1:3)
         call get_coulomb_and_dcoulomb_pqeq(dr,alphasc(ity,jty),Csicj,inxnpqeq(ity,jty),TBL_Eclmb_psc,ff)
         Csicj=-Cclmb0_qeq*Csicj*qjc*Zpqeq(ity)

         if(isPolarizable(jty)) then
             dr(1:3)=shelli(1:3)-shellj(1:3)
             call get_coulomb_and_dcoulomb_pqeq(dr,alphass(ity,jty),Csisj,inxnpqeq(ity,jty),TBL_Eclmb_pss,ff)
             Csisj=Cclmb0_qeq*Csisj*Zpqeq(ity)*Zpqeq(jty)
         endif
      endif

      hshs(i) = hshs(i) + hessian(j1,i)*hs(j)
      hsht(i) = hsht(i) + hessian(j1,i)*ht(j)

!--- get half of potential energy, then sum it up if atoms are resident.
      Est1 = 0.5d0*(Ccicj + Csisj)

      Est = Est + Est1 + Csicj
!--- nbplist does not distinguish i,j pairs from intra-/inter-node atoms.
      !if(j<=NATOMS) Est = Est + Est1 
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
implicit none
real(8),intent(OUT) :: Gnew(2)
real(8) :: eta_ity
integer :: i,j,j1, ity

real(8) :: gssum, gtsum

integer :: ti,tj,tk
call system_clock(ti,tk)

!$omp parallel do default(shared), schedule(runtime), private(gssum, gtsum, eta_ity,i,j,j1,ity)
do i=1,NATOMS

   gssum=0.d0
   gtsum=0.d0
   do j1=1, nbplist(0,i) 
      j = nbplist(j1,i)
      gssum = gssum + hessian(j1,i)*qs(j)
      gtsum = gtsum + hessian(j1,i)*qt(j)
   enddo

   ity = nint(atype(i))
   eta_ity = eta(ity)

   gs(i) = - chi(ity) - eta_ity*qs(i) - gssum - fpqeq(i)
   gt(i) = - 1.d0     - eta_ity*qt(i) - gtsum

enddo 
!$omp end parallel do

gnew(1) = dot_product(gs(1:NATOMS), gs(1:NATOMS))
gnew(2) = dot_product(gt(1:NATOMS), gt(1:NATOMS))
call MPI_ALLREDUCE(MPI_IN_PLACE, Gnew, size(Gnew), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

call system_clock(tj,tk)
it_timer(19)=it_timer(19)+(tj-ti)

end subroutine

end subroutine PQEq


!-------------------------------------------------------------------------------------------
subroutine initialize_eField(myid)
implicit none
!-------------------------------------------------------------------------------------------
integer,intent(in) :: myid

if(myid==0) then
   print'(a)','-----------------------------------------------------------'
   print'(a,f12.6,i6)','eField [V/A], eFiled direction : ', eFieldStrength, eFieldDir
   print'(a)','-----------------------------------------------------------'
endif

return
end subroutine


!-------------------------------------------------------------------------------------------
subroutine EEfield(Etotal,NATOMS,pos,q,f,atype,Eev_kcal)
implicit none 
!-------------------------------------------------------------------------------------------
integer,intent(in) :: NATOMS
real(8),intent(in) :: pos(NATOMS,3),q(NATOMS),atype(NATOMS),Eev_kcal

real(8) :: f(NATOMS,3), Etotal

integer :: i, ity
real(8) :: Eenergy, Eforce, qic

! NOTE: the energy function from an eField is to be determined and not computed here. 
do i=1, NATOMS

   ity = nint(atype(i))
   qic = q(i) + Zpqeq(ity)

   Eforce = -qic*eFieldStrength*Eev_kcal

   f(i,eFieldDir) = f(i,eFieldDir) + Eforce

enddo

return
end subroutine

!-------------------------------------------------------------------------------------------
subroutine get_coulomb_and_dcoulomb_pqeq(rr,alpha,Eclmb,inxn,TBL_Eclmb_PQEq,ff)
implicit none
!-------------------------------------------------------------------------------------------
real(8),intent(in) :: rr(3),alpha,TBL_Eclmb_PQEq(ntype_pqeq2,NTABLE,0:1)
integer,intent(in) :: inxn
real(8) :: Ecc,Esc,Ess,dEcc,dEsc,dEss,Eclmb,dEclmb,ff(3)

integer :: i,itb, itb1
real(8) :: drtb, drtb1, dr1, dr2, dr1i 
real(8) :: Tap, dTap, clmb, dclmb, screen, dscreen 

real(8) :: dr3,dr4,dr5,dr6,dr7

dr2 = sum(rr(1:3)*rr(1:3))

! note that shell-shell distance could be greater than cutoff though core-core is within the cutoff.
if(dr2>rctap2) return 

dr1 = sqrt(dr2)

itb = int(dr2*UDRi)
itb1 = itb+1
drtb = dr2 - itb*UDR
drtb = drtb*UDRi
drtb1= 1.d0-drtb

Eclmb = drtb1*TBL_Eclmb_PQEq(inxn,itb,0)  + drtb*TBL_Eclmb_PQEq(inxn,itb1,0)
dEclmb = drtb1*TBL_Eclmb_PQEq(inxn,itb,1)  + drtb*TBL_Eclmb_PQEq(inxn,itb1,1)

ff(1:3)=dEclmb*rr(1:3)

return

dr1i = 1.d0/dr1
clmb = dr1i
!TODO 1. Tabulate the coulomb term, 2. change the derivatives to our convention, i.e. u(r)'/r^2*dr(1:3)
!dclmb = -dr1i*dr1i*dr1i 
dclmb = -dr1i*dr1i

screen = erf(alpha*dr1)
!dscreen = 2.d0*alpha*sqrtpi_inv*exp(-alpha*alpha*dr2)*dr1i
dscreen = 2.d0*alpha*sqrtpi_inv*exp(-alpha*alpha*dr2)

!--- core-core distance is withing the taper cutoff, but core-shell & shell-shell distance can be beyond the cutoff.
!--- Directly computing the taper function for now. 
dr3 = dr1*dr2
dr4 = dr2*dr2
dr5 = dr1*dr2*dr2
dr6 = dr2*dr2*dr2
dr7 = dr1*dr2*dr2*dr2
Tap = CTap(7)*dr7 + CTap(6)*dr6 + CTap(5)*dr5 + CTap(4)*dr4 + CTap(0)
dTap = 7d0*CTap(7)*dr6 + 6d0*CTap(6)*dr5 + 5d0*CTap(5)*dr4 + 4d0*CTap(4)*dr3

Eclmb = clmb*screen*Tap
dEclmb = dclmb*screen*Tap + clmb*dscreen*Tap + clmb*screen*dTap

ff(1:3)=dEclmb*rr(1:3)/dr1

return
end subroutine


!-------------------------------------------------------------------------------------------
subroutine set_alphaij_pqeq()
implicit none
!-------------------------------------------------------------------------------------------
integer :: ity,jty

real(8) :: alpha_ci,alpha_cj,alpha_si,alpha_sj 

alphacc(:,:)=0.d0
alphasc(:,:)=0.d0
alphass(:,:)=0.d0

do ity=1,ntype_pqeq

   alpha_ci=0.5d0*lambda_pqeq/Rcpqeq(ity)**2
   alpha_si=0.5d0*lambda_pqeq/Rspqeq(ity)**2

   do jty=1,ntype_pqeq

      alpha_cj=0.5d0*lambda_pqeq/Rcpqeq(jty)**2
      alpha_sj=0.5d0*lambda_pqeq/Rspqeq(jty)**2

      ! core(i)-core(j) term.
      alphacc(ity,jty)=sqrt( (alpha_ci*alpha_cj)/(alpha_ci + alpha_cj) )

      ! shell(i)-shell(j) term.
      if( isPolarizable(ity) .and. isPolarizable(jty) ) &
          alphass(ity,jty)=sqrt( (alpha_si*alpha_sj)/(alpha_si + alpha_sj) )

      ! shell(i)-core(j) term. the first index must be polarizable atom type. 
      ! C(r_si_cj) is fine but use alphasc(jty,ity) when refer to alphasc for C(r_ci_sj).
      if( isPolarizable(ity) ) &
          alphasc(ity,jty)=sqrt( (alpha_si*alpha_cj)/(alpha_si + alpha_cj) )

   enddo

enddo

end subroutine


!-------------------------------------------------------------------------------------------
subroutine initialize_pqeq(chi,eta)
implicit none
!-------------------------------------------------------------------------------------------
integer :: ity,jty,icounter
real(8) :: chi(ntype_pqeq),eta(ntype_pqeq)

integer :: i,inxn
real(8) :: Acc,Asc,Ass
real(8) :: dr1,dr2,dr3,dr4,dr5,dr6,dr7,dr1i,UDR,UDRi
real(8) :: clmb,dclmb,screen,dscreen,Tap,dTap,Eclmb,dEclmb

call set_alphaij_pqeq()

!--- for PQEq
do ity = 1, ntype_pqeq
  if( .not. isPolarizable(ity) ) then
    if(myid==0) &
       print'(a,i3,a)','atom type ', ity, ' is not polarizable. Setting Z & K to zero.'
    Zpqeq(ity)=0.d0
    Kspqeq(ity)=0.d0
  else
    if(myid==0) then
       print'(a $)','updating chi and eta '

       print'(a,i3,1x,a2,1x,4f12.6)', &
         'name,chi,X0pqeq,eta,J0pqeq :', &
          ity,Elempqeq(ity), chi(ity),X0pqeq(ity), eta(ity),J0pqeq(ity)
    endif

    chi(ity)=X0pqeq(ity)
    eta(ity)=J0pqeq(ity)
  endif
enddo

!--- 2x factor in param.F90
eta(:)=2.d0*eta(:)

allocate(inxnpqeq(ntype_pqeq,ntype_pqeq))

icounter=0
do ity=1, ntype_pqeq
do jty=ity, ntype_pqeq

   icounter = icounter + 1

   inxnpqeq(ity,jty)=icounter
   inxnpqeq(jty,ity)=inxnpqeq(ity,jty)
enddo; enddo

allocate(TBL_Eclmb_pcc(ntype_pqeq2,NTABLE,0:1))
allocate(TBL_Eclmb_psc(ntype_pqeq2,NTABLE,0:1))
allocate(TBL_Eclmb_pss(ntype_pqeq2,NTABLE,0:1))

!--- unit distance in r^2 scale
UDR = rctap2/NTABLE
UDRi = 1.d0/UDR

do ity=1, ntype_pqeq
do jty=ity, ntype_pqeq

   Acc = alphacc(ity,jty)
   Asc = alphasc(ity,jty)
   Ass = alphass(ity,jty)

   inxn = inxnpqeq(ity,jty)

   do i=1,NTABLE

         dr2 = UDR*i
         dr1 = sqrt(dr2)

!--- Interaction Parameters:
         dr3 = dr1*dr2
         dr4 = dr2*dr2
         dr5 = dr1*dr2*dr2
         dr6 = dr2*dr2*dr2
         dr7 = dr1*dr2*dr2*dr2 

         Tap = CTap(7)*dr7 + CTap(6)*dr6 + &
               CTap(5)*dr5 + CTap(4)*dr4 + CTap(0)

!--- Force Calculation:
         dTap = 7d0*CTap(7)*dr5 + 6d0*CTap(6)*dr4 + &
                5d0*CTap(5)*dr3 + 4d0*CTap(4)*dr2

         dr1i = 1.d0/dr1

         clmb = dr1i
         dclmb = -dr1i*dr1i*dr1i 

         !--- core-core interaction
         screen = erf(Acc*dr1)
         dscreen = 2.d0*Acc*sqrtpi_inv*exp(-Acc*Acc*dr2)*dr1i

         Eclmb = clmb*screen*Tap
         dEclmb = dclmb*screen*Tap + clmb*dscreen*Tap + clmb*screen*dTap

         TBL_Eclmb_pcc(inxn,i,0) = Eclmb
         TBL_Eclmb_pcc(inxn,i,1) = dEclmb

         !--- core-shell interaction
         screen = erf(Asc*dr1)
         dscreen = 2.d0*Asc*sqrtpi_inv*exp(-Asc*Asc*dr2)*dr1i

         Eclmb = clmb*screen*Tap
         dEclmb = dclmb*screen*Tap + clmb*dscreen*Tap + clmb*screen*dTap

         TBL_Eclmb_psc(inxn,i,0) = Eclmb
         TBL_Eclmb_psc(inxn,i,1) = dEclmb

         !--- shell-shell interaction
         screen = erf(Ass*dr1)
         dscreen = 2.d0*Ass*sqrtpi_inv*exp(-Ass*Ass*dr2)*dr1i

         Eclmb = clmb*screen*Tap
         dEclmb = dclmb*screen*Tap + clmb*dscreen*Tap + clmb*screen*dTap

         TBL_Eclmb_pss(inxn,i,0) = Eclmb
         TBL_Eclmb_pss(inxn,i,1) = dEclmb

   enddo

enddo
enddo

end subroutine


end module
