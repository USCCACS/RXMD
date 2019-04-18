!------------------------------------------------------------------------------
program rxmd
!------------------------------------------------------------------------------
use base
use init
use mpi_mod
use reaxff_param_mod
use cmdline_args
use velocity_modifiers_mod
use communication_mod
!use CG
!------------------------------------------------------------------------------
implicit none
integer :: i,ity,it1,it2,irt

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

if(myid==0)  print'(a30)', 'rxmd has started'

!--- process command line arguments
call get_cmdline_args(myid, eFieldDir, eFieldStrength)

!--- get forcefield-independent parameters 
call get_rxmd_parms(ParmPath)

!--- initialize the MD system
call mdcontext_base(atype, pos, v, f, q)

!--- Enter Main MD loop 
call system_clock(it1,irt)
call mddriver_func(ntime_step)
call system_clock(it2,irt)
it_timer(Ntimer)=(it2-it1)

call FinalizeMD(irt)

call MPI_FINALIZE(ierr)
end PROGRAM

!------------------------------------------------------------------------------
subroutine FinalizeMD(irt)
use base, only : myid, ierr
use mpi_mod
use atoms
use memory_allocator_mod
!------------------------------------------------------------------------------
implicit none
integer :: i
integer,intent(in) :: irt ! time resolution
integer,allocatable :: ibuf(:),ibuf1(:)

allocate(ibuf(nmaxas),ibuf1(nmaxas))
ibuf(:)=0
do i=1,nmaxas
   ibuf(i)=maxval(maxas(:,i))
enddo
call MPI_ALLREDUCE(ibuf, ibuf1, nmaxas, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

call MPI_ALLREDUCE(it_timer, it_timer_max, Ntimer, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(it_timer, it_timer_min, Ntimer, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)

if(myid==0) then
   print'(a)','----------------------------------------------'
   print'(a20,i12)', 'MAXNEIGHBS: ', ibuf1(2)
   print'(a20,i12)', 'MAXNEIGHBS10: ', ibuf1(3)
   print'(a20,i12)', 'MAXNBUFFER(MOVE): ', ibuf1(1)+ibuf1(4)
   print'(a20,i12)', 'MAXNBUFFER(COPY): ', ibuf1(1)+ibuf1(5)
   print'(a20,i12)', 'QEq Iterations: ', it_timer_max(24)
   print'(a20,f12.2)','Memory (MB): ', GetTotalMemory()*1d-6
   print*

   print'(a20,f12.4,3x,f12.4)','QEq: ',  dble(it_timer_max(1))/irt, dble(it_timer_min(1))/irt
   print'(a20,f12.4,3x,f12.4)','qeq_initialize: ',  dble(it_timer_max(16))/irt, dble(it_timer_min(16))/irt
   print'(a20,f12.4,3x,f12.4)','get_hsh: ',  dble(it_timer_max(18))/irt, dble(it_timer_min(18))/irt
   print'(a20,f12.4,3x,f12.4)','get_gradient: ',  dble(it_timer_max(19))/irt, dble(it_timer_min(19))/irt
   print*

   print'(a20,f12.4,3x,f12.4)','LINKEDLIST: ',  dble(it_timer_max(3))/irt, dble(it_timer_min(3))/irt
   print'(a20,f12.4,3x,f12.4)','COPYATOMS: ',    dble(it_timer_max(4))/irt, dble(it_timer_min(4))/irt
   print'(a20,f12.4,3x,f12.4)','send_rec: ', dble(it_timer_max(25))/irt, dble(it_timer_min(25))/irt
   print'(a20,f12.4,3x,f12.4)','store_atoms: ', dble(it_timer_max(26))/irt, dble(it_timer_min(26))/irt
   print'(a20,f12.4,3x,f12.4)','append_atoms: ', dble(it_timer_max(27))/irt, dble(it_timer_min(27))/irt

   print'(a20,f12.4,3x,f12.4)','NEIGHBORLIST: ', dble(it_timer_max(5))/irt, dble(it_timer_min(5))/irt
   print'(a20,f12.4,3x,f12.4)','GetNBPairList: ', dble(it_timer_max(15))/irt, dble(it_timer_min(15))/irt
   print*

   print'(a20,f12.4,3x,f12.4)','BOCALC: ', dble(it_timer_max(6))/irt, dble(it_timer_min(6))/irt
   print'(a20,f12.4,3x,f12.4)','ENbond: ', dble(it_timer_max(7))/irt, dble(it_timer_min(7))/irt
   print'(a20,f12.4,3x,f12.4)','Ebond: ', dble(it_timer_max(8))/irt, dble(it_timer_min(8))/irt
   print'(a20,f12.4,3x,f12.4)','Elnpr: ', dble(it_timer_max(9))/irt, dble(it_timer_min(9))/irt
   print'(a20,f12.4,3x,f12.4)','Ehb: ', dble(it_timer_max(10))/irt, dble(it_timer_min(10))/irt
   print'(a20,f12.4,3x,f12.4)','E3b: ', dble(it_timer_max(11))/irt, dble(it_timer_min(11))/irt
   print'(a20,f12.4,3x,f12.4)','E4b: ', dble(it_timer_max(12))/irt, dble(it_timer_min(12))/irt
   print'(a20,f12.4,3x,f12.4)','ForceBondedTerms: ', dble(it_timer_max(13))/irt, dble(it_timer_min(13))/irt
   print*

   print'(a20,f12.4,3x,f12.4)','WriteBND: ', dble(it_timer_max(20))/irt, dble(it_timer_min(20))/irt
   print'(a20,f12.4,3x,f12.4)','WritePDB: ', dble(it_timer_max(21))/irt, dble(it_timer_min(21))/irt
   print'(a20,f12.4,3x,f12.4)','ReadBIN: ', dble(it_timer_max(22))/irt, dble(it_timer_min(22))/irt
   print'(a20,f12.4,3x,f12.4)','WriteBIN: ', dble(it_timer_max(23))/irt, dble(it_timer_min(23))/irt
   print*

   print'(a20,f12.4,3x,f12.4)','total (sec): ',dble(it_timer_max(Ntimer))/irt, dble(it_timer_min(Ntimer))/irt

   print'(a)','----------------------------------------------'

   print'(a30)', 'rxmd successfully finished'
endif

deallocate(ibuf,ibuf1)

end subroutine

!----------------------------------------------------------------------------------------
subroutine PRINTE(atype, v, q)
use mpi_mod
use base, only : hh, hhi, natoms, gnatoms, mdbox, myid, ierr
use atoms
use reaxff_param_mod
use memory_allocator_mod
! calculate the kinetic energy and sum up all of potential energies, then print them.
!----------------------------------------------------------------------------------------
implicit none

real(8),intent(in) :: atype(NBUFFER), q(NBUFFER)
real(8),intent(in) :: v(NBUFFER,3)

integer :: i,ity,cstep
real(8) :: qq=0.d0,tt=0.d0,ss=0.d0,buf(0:23)

i=nstep/pstep+1
maxas(i,1)=NATOMS

KE=0.d0
do i=1, NATOMS
   ity=nint(atype(i))
   KE = KE + hmas(ity)*sum(v(i,1:3)*v(i,1:3))
enddo
qq=sum(q(1:NATOMS))

!--- pressure 
ss=sum(astr(1:3))/3.d0

!--- potential energy 
PE(0)=sum(PE(1:13))

!--- copy data into buffer
buf(0:13) = PE(0:13)
buf(14) = KE; buf(15) = ss; buf(16) = qq
buf(17:22)=astr(1:6)
call MPI_ALLREDUCE (MPI_IN_PLACE, buf, size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!--- copy data from buffer
GPE(0:13) = buf(0:13)
GKE = buf(14); ss = buf(15); qq = buf(16)
astr(1:6)=buf(17:22)

!--- compute properties
GPE(:)=GPE(:)/GNATOMS
GKE=GKE/GNATOMS
tt=GKE*UTEMP
astr(1:6)=astr(1:6)/MDBOX*USTRS/pstep
ss=ss/MDBOX*USTRS/pstep

!--- total energy
GTE = GKE + GPE(0)
if(myid==0) then
   
   cstep = nstep + current_step 

   write(6,'(a,i9,3es13.5,6es11.3,1x,3f8.2,i4,f8.2,f8.2)') 'MDstep: ',cstep,GTE,GPE(0),GKE, &
   GPE(1),sum(GPE(2:4)),sum(GPE(5:7)),sum(GPE(8:9)),GPE(10),sum(GPE(11:13)), &
   tt, ss, qq, nstep_qeq, GetTotalMemory()*1e-9, MPI_WTIME()-wt0

   !write(6,'(a,i9,6f12.6)') 'stress : ',cstep,astr(1:6)

endif

!--- reset stress tensor accumulator
astr(1:6)=0.d0

!--- save current time
wt0 = MPI_WTIME()
end subroutine
