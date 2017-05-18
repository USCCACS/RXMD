!------------------------------------------------------------------------------
program rxmd
use base; use atoms; use parameters; use CG
!------------------------------------------------------------------------------
implicit none
integer :: i,it1,it2,irt,provided
real(8) :: ctmp, dr(3)

!call MPI_INIT(ierr)
call MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,provided,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

if(myid==0)  print'(a30)', 'rxmd has started'

!--- read ffield file
CALL GETPARAMS(FFPath,FFDescript)

!--- initialize the MD system
CALL INITSYSTEM(atype, pos, v, f, q)

if(mdmode==10) call ConjugateGradient(atype,pos)

call QEq(atype, pos, q)
call FORCE(atype, pos, f, q)

!--- Enter Main MD loop 
call system_clock(it1,irt)

do nstep=0, ntime_step-1

   if(mod(nstep,pstep)==0) then
       call PRINTE(atype, v, q)
       if(saveRunProfile) call SaveRunProfileData(RunProfileFD, nstep)
   endif
   if(mod(nstep,fstep)==0) &
        call OUTPUT(atype, pos, v, f, q, GetFileNameBase(current_step+nstep))

   if(mod(nstep,sstep)==0.and.mdmode==4) &
      v(1:3,1:NATOMS)=vsfact*v(1:3,1:NATOMS)

   if(mod(nstep,sstep)==0.and.mdmode==5) then
      ctmp = (treq*UTEMP0)/( GKE*UTEMP )
      v(1:3,1:NATOMS)=sqrt(ctmp)*v(1:3,1:NATOMS)
   endif

   if(mod(nstep,sstep)==0.and.(mdmode==0.or.mdmode==6)) &
      call INITVELOCITY(atype, v)

!--- correct the c.o.m motion
   if(mod(nstep,sstep)==0.and.mdmode==7) &
      call ScaleTemperature(atype, v)

!--- update velocity
   call vkick(1.d0, atype, v, f) 

!--- update coordinates
   qsfv(1:NATOMS)=qsfv(1:NATOMS)+0.5d0*dt*Lex_w2*(q(1:NATOMS)-qsfp(1:NATOMS))
   qsfp(1:NATOMS)=qsfp(1:NATOMS)+dt*qsfv(1:NATOMS)

   pos(1:3,1:NATOMS)=pos(1:3,1:NATOMS)+dt*v(1:3,1:NATOMS)

!--- migrate atoms after positions are updated
   call COPYATOMS(MODE_MOVE,[0.d0, 0.d0, 0.d0],atype, pos, v, f, q)
   
   if(mod(nstep,qstep)==0) call QEq(atype, pos, q)
   call FORCE(atype, pos, f, q)

!--- update velocity
   call vkick(1.d0, atype, v, f) 
   qsfv(1:NATOMS)=qsfv(1:NATOMS)+0.5d0*dt*Lex_w2*(q(1:NATOMS)-qsfp(1:NATOMS))

enddo

!--- save the final configurations
call OUTPUT(atype, pos, v, f, q,  GetFileNameBase(current_step+nstep))

!--- update rxff.bin in working directory for continuation run
call WriteBIN(atype, pos, v, q, GetFileNameBase(-1))

call system_clock(it2,irt)
it_timer(Ntimer)=(it2-it1)

call FinalizeMD(irt)

call MPI_FINALIZE(ierr)
end PROGRAM

!------------------------------------------------------------------------------
subroutine SaveRunProfileData(fd, MDstep)
use base; use atoms
!------------------------------------------------------------------------------
implicit none
integer :: i, fd, MDstep

write(fd,'(i9,a3)') MDstep," : " 

do i=1,3
   write(fd,'(3es16.8)'), HH(1:3,i,0)
enddo

write(fd,'(14es16.8)') GPE(0:13)

!do i=1, NATOMS
!   write(fd,'(es25.13,9f15.5)') atype(i),pos(1:3,i),v(1:3,i),f(1:3,i)
!enddo

write(fd,*) 
write(fd,*) 

end subroutine

!------------------------------------------------------------------------------
subroutine FinalizeMD(irt)
use atoms; use MemoryAllocator
!------------------------------------------------------------------------------
implicit none
integer :: i
integer,intent(in) :: irt ! time resolution
integer,allocatable :: ibuf(:),ibuf1(:)

!--- close summary file
if(saveRunProfile) close(RunProfileFd)

allocate(ibuf(nmaxas),ibuf1(nmaxas))
ibuf(:)=0
do i=1,nmaxas
   ibuf(i)=maxval(maxas(:,i))
enddo
call MPI_ALLREDUCE (ibuf, ibuf1, nmaxas, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

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

!------------------------------------------------------------------------------
subroutine vkick(dtf, atype, v, f)
use atoms
!------------------------------------------------------------------------------
implicit none

real(8) :: atype(NBUFFER),v(3,NBUFFER),f(3,NBUFFER)

integer :: i, ity
real(8) :: dtf

do i=1,NATOMS
   ity = nint(atype(i))
   v(1:3,i) = v(1:3,i) + dtf*dthm(ity)*f(1:3,i)
enddo

end subroutine

!----------------------------------------------------------------------------------------
subroutine PRINTE(atype, v, q)
use atoms; use parameters; use MemoryAllocator
! calculate the kinetic energy and sum up all of potential energies, then print them.
!----------------------------------------------------------------------------------------
implicit none

real(8),intent(in) :: atype(NBUFFER), q(NBUFFER)
real(8),intent(in) :: v(3,NBUFFER)

integer :: i,ity,cstep
real(8),save :: wt0
real(8) :: qq=0.d0,tt=0.d0,ss=0.d0,buf(0:20),Gbuf(0:20)

i=nstep/pstep+1
maxas(i,1)=NATOMS

KE=0.d0
do i=1, NATOMS
   ity=nint(atype(i))
   KE = KE + hmas(ity)*sum(v(1:3,i)*v(1:3,i))
enddo
qq=sum(q(1:NATOMS))

#ifdef STRESS
!--- pressure 
ss=sum(astr(1:3,1:NATOMS))
#endif

!--- potential energy 
PE(0)=sum(PE(1:13))

!--- copy data into buffer
buf(0:13) = PE(0:13)
buf(14) = KE; buf(15) = ss; buf(16) = qq
call MPI_ALLREDUCE (buf, Gbuf, size(buf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

!--- copy data from buffer
GPE(0:13) = Gbuf(0:13)
GKE = Gbuf(14); ss = Gbuf(15); qq = Gbuf(16)

!--- compute properties
GPE(:)=GPE(:)/GNATOMS
GKE=GKE/GNATOMS
tt=GKE*UTEMP
#ifdef STRESS
ss=ss/3.d0/MDBOX*USTRS
#endif 

!--- total energy
GTE = GKE + GPE(0)
if(myid==0) then
   
   cstep = nstep + current_step 

#ifdef PROFILE
   write(6,'(i9,3es13.5,13es11.3,1x,3f8.2,i4,f8.2,f8.2)') cstep,GTE,GPE(0),GKE, &
   GPE(1:13), tt, ss, qq, nstep_qeq, GetTotalMemory()*1e-9, MPI_WTIME()-wt0 
#else
   write(6,'(i9,3es13.5,6es11.3,1x,3f8.2,i4,f8.2,f8.2)') cstep,GTE,GPE(0),GKE, &
   GPE(1),sum(GPE(2:4)),sum(GPE(5:7)),sum(GPE(8:9)),GPE(10),sum(GPE(11:13)), &
   tt, ss, qq, nstep_qeq, GetTotalMemory()*1e-9, MPI_WTIME()-wt0 
#endif

#ifdef STRESS
   write(6,'(6es13.5)') pint(1,1)*USTRS, pint(2,2)*USTRS, pint(3,3)*USTRS, &
                        pint(2,3)*USTRS, pint(3,1)*USTRS, pint(1,2)*USTRS
#endif

endif

!--- save current time
wt0 = MPI_WTIME()
end subroutine

!----------------------------------------------------------------------------------------
subroutine LINKEDLIST(atype, rreal, cellDims, headAtom, atomList, NatomPerCell, Ncells, NLAYERS)
use atoms
! partitions the volume into linked-list cells <lcsize>
!----------------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: atype(NBUFFER), rreal(3,NBUFFER), cellDims(3)

integer,intent(in) :: Ncells(3), NLAYERS
integer,intent(out) :: atomList(NBUFFER)
integer,intent(out) :: NatomPerCell(-NLAYERS:Ncells(1)-1+NLAYERS, &
                                    -NLAYERS:Ncells(2)-1+NLAYERS, &
                                    -NLAYERS:Ncells(3)-1+NLAYERS) 
integer,intent(out) :: headAtom(-NLAYERS:Ncells(1)-1+NLAYERS, & 
                                -NLAYERS:Ncells(2)-1+NLAYERS, &
                                -NLAYERS:Ncells(3)-1+NLAYERS) 

real(8) :: rnorm(3,NBUFFER)
integer :: n, l(3), j

integer :: ti,tj,tk
call system_clock(ti,tk)

call xu2xs(rreal,rnorm)

headAtom(:,:,:) = -1; atomList(:) = 0; NatomPerCell(:,:,:)=0

!--- copyptr(6) stores the last atom index copied in COPYATOMS.
do n=1, copyptr(6) 

   if(nint(atype(n))==0) cycle

   l(1:3) = floor(rnorm(1:3,n)/cellDims(1:3))

   atomList(n) = headAtom(l(1), l(2), l(3))
   headAtom(l(1), l(2), l(3)) = n
   NatomPerCell(l(1), l(2), l(3)) = NatomPerCell(l(1), l(2), l(3)) + 1
enddo

call system_clock(tj,tk)
it_timer(3)=it_timer(3)+(tj-ti)

end subroutine 

!----------------------------------------------------------------------
subroutine NEIGHBORLIST(nlayer, atype, pos)
use atoms; use parameters
! calculate neighbor list for atoms witin cc(1:3, -nlayer:nlayer) cells.
!----------------------------------------------------------------------
implicit none
integer,intent(in) :: nlayer
real(8),intent(in) :: atype(NBUFFER), pos(3,NBUFFER)

integer :: c1,c2,c3, ic(3), c4, c5, c6
integer :: n, n1, m, m1, nty, mty, inxn
real(8) :: dr(3), dr2

integer :: i,j,i1,j1
logical :: isFound

integer :: ti,tj,tk
call system_clock(ti,tk)

nbrlist(:,0) = 0

!$omp parallel do default(shared) collapse(3) & 
!$omp private(c1,c2,c3,ic,c4,c5,c6,n,n1,m,m1,nty,mty,inxn,dr,dr2) 
DO c1=-nlayer, cc(1)-1+nlayer
DO c2=-nlayer, cc(2)-1+nlayer
DO c3=-nlayer, cc(3)-1+nlayer

  m = header(c1, c2, c3)
  do m1=1, nacell(c1, c2, c3)
     mty = nint(atype(m))

     do c4 = -1, 1
     do c5 = -1, 1
     do c6 = -1, 1
        ic(1:3) = [c1, c2, c3] + [c4, c5, c6]

        n = header(ic(1),ic(2),ic(3))
        do n1=1, nacell(ic(1), ic(2), ic(3))

           if(n/=m) then
             nty = nint(atype(n))
             inxn = inxn2(mty, nty)

             dr(1:3) = pos(1:3,n) - pos(1:3,m) 
             dr2 = sum(dr(1:3)*dr(1:3))

             if(dr2<rc2(inxn)) then 
                nbrlist(m, 0) = nbrlist(m, 0) + 1
                nbrlist(m, nbrlist(m, 0)) = n
             endif 
           endif

           n=llist(n) 
        enddo
     enddo; enddo; enddo

     m = llist(m)
  enddo
enddo; enddo; enddo
!$omp end parallel do 

!--- to get the reverse information (i.e. from i,j1&j to i1), store <i1> into <nbrindx>.

!$omp parallel do default(shared) private(i,i1,j,j1,isFound)
do i=1, copyptr(6)
   do i1 = 1, nbrlist(i,0)
      j = nbrlist(i,i1)
      isFound=.false.
      do j1 = 1, nbrlist(j,0)
         if(i == nbrlist(j,j1)) then
            nbrindx(i,i1)=j1
            isFound=.true.
         endif
      enddo
      if(.not.isFound) &
      print'(a,i6,30i4)','ERROR: inconsistency between nbrlist and nbrindx found', &
           myid, i,nbrlist(i,0:nbrlist(i,0)), j, nbrlist(j,0:nbrlist(j,0))
   enddo
enddo
!$omp end parallel do

!--- error trap
n=maxval(nbrlist(1:NATOMS,0))
if(n > MAXNEIGHBS) then
   write(6,'(a45,2i5)') "ERROR: overflow of max # in neighbor list, ", myid, n
   call MPI_FINALIZE(ierr)
   stop
endif

!--- for array size stat
if(mod(nstep,pstep)==0) then
  maxas(nstep/pstep+1,2)=maxval(nbrlist(1:NATOMS,0))
endif

call system_clock(tj,tk)
it_timer(5)=it_timer(5)+(tj-ti)

end subroutine

!----------------------------------------------------------------------
subroutine GetNonbondingPairList(pos)
use atoms; use parameters
!----------------------------------------------------------------------
implicit none

real(8),intent(in) :: pos(3,NBUFFER)

integer :: c1,c2,c3,c4,c5,c6,i,j,m,n,mn,iid,jid
integer :: l2g
real(8) :: dr(3), dr2

integer :: ti,tj,tk
call system_clock(ti,tk)

! reset non-bonding pair list
nbplist(:,0)=0

!$omp parallel do default(shared),private(c1,c2,c3,c4,c5,c6,i,j,m,n,mn,iid,jid,dr,dr2)
do c1=0, nbcc(1)-1
do c2=0, nbcc(2)-1
do c3=0, nbcc(3)-1

   i = nbheader(c1,c2,c3)
   do m = 1, nbnacell(c1,c2,c3)

      do mn = 1, nbnmesh
         c4 = c1 + nbmesh(1,mn)
         c5 = c2 + nbmesh(2,mn)
         c6 = c3 + nbmesh(3,mn)

         j = nbheader(c4,c5,c6)
         do n=1, nbnacell(c4,c5,c6)

            !if(i<j .or. NATOMS<j) then
            if(i/=j) then
               dr(1:3) = pos(1:3,i) - pos(1:3,j)
               dr2 = sum(dr(1:3)*dr(1:3))

               if(dr2<=rctap2) then
                 nbplist(i,0)=nbplist(i,0)+1
                 nbplist(i,nbplist(i,0))=j
               endif

            endif

            j=nbllist(j)
         enddo
       enddo

      i=nbllist(i)
   enddo
enddo; enddo; enddo
!$omp end parallel do

call system_clock(tj,tk)
it_timer(15)=it_timer(15)+(tj-ti)

end subroutine

!----------------------------------------------------------------------
subroutine angular_momentum(atype, pos, v)
use atoms; use parameters
!----------------------------------------------------------------------
implicit none

real(8) :: atype(NBUFFER), pos(3,NBUFFER),v(3,NBUFFER)

integer :: i,ity
real(8) :: com(3), Gcom(3), intsr(3,3), Gintsr(3,3), intsr_i(3,3), angm(3), Gangm(3), angv(3), mm, Gmm
real(8) :: dr(3), dv(3)

!--- get center of mass
com(:)=0.d0;     Gcom(:)=0.d0
mm=0.d0; Gmm=0.d0

do i=1, NATOMS
   ity = nint(atype(i))
   mm = mm + mass(ity)
   com(1:3) = mass(ity)*pos(1:3,i)
enddo

call MPI_ALLREDUCE(mm, Gmm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(com, Gcom, 3, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
Gcom(1:3) = Gcom(1:3)/Gmm


!--- get the angular momentum and inertia tensor from the com
angm(:)=0.d0;    Gangm(:)=0.d0
intsr(:,:)=0.d0; Gintsr(:,:)=0.d0

do i=1, NATOMS
   dr(1:3) = pos(1:3,i) - Gcom(1:3)
   
   angm(1) = mass(ity)*( dr(2)*v(3,i)-dr(3)*v(2,i) )
   angm(2) = mass(ity)*( dr(3)*v(1,i)-dr(1)*v(3,i) )
   angm(3) = mass(ity)*( dr(1)*v(2,i)-dr(2)*v(1,i) )

   intsr(1,1) = mass(ity)*( dr(2)**2+dr(3)**2 )
   intsr(2,2) = mass(ity)*( dr(3)**2+dr(1)**2 )
   intsr(3,3) = mass(ity)*( dr(1)**2+dr(2)**2 )

   intsr(1,2) =-mass(ity)*( dr(1)*dr(2) )
   intsr(1,3) =-mass(ity)*( dr(1)*dr(3) )
   intsr(2,3) =-mass(ity)*( dr(2)*dr(3) )

   intsr(2,1) = intsr(1,2)
   intsr(3,1) = intsr(1,3)
   intsr(3,2) = intsr(2,3)
enddo

call MPI_ALLREDUCE(angm, Gangm, 3, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE(intsr, Gintsr, 9, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)

!--- get angular velocity
call matinv(Gintsr, intsr_i)

angv(1) = sum(intsr_i(1,1:3)*angm(1:3))
angv(2) = sum(intsr_i(2,1:3)*angm(1:3))
angv(3) = sum(intsr_i(3,1:3)*angm(1:3))


!--- correct rotational motion wrt CoM.
do i=1,NATOMS
   dr(1:3) = pos(1:3,i) - Gcom(1:3)
   dv(1) = angv(2)*dr(3) - angv(3)*dr(2)
   dv(2) = angv(3)*dr(1) - angv(1)*dr(3)
   dv(3) = angv(1)*dr(2) - angv(2)*dr(1)

   v(1:3,i) = v(1:3,i) - dv(1:3)
enddo

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine matinv(m1,m2)
! get inverse of m1 and save to m2
!--------------------------------------------------------------------------------------------------------------
implicit none
real(8) :: m1(3,3), m2(3,3), detm

m2(1,1) = m1(2,2)*m1(3,3)-m1(2,3)*m1(3,2)
m2(1,2) = m1(1,3)*m1(3,2)-m1(1,2)*m1(3,3)
m2(1,3) = m1(1,2)*m1(2,3)-m1(1,3)*m1(2,2)
m2(2,1) = m1(2,3)*m1(3,1)-m1(2,1)*m1(3,3)
m2(2,2) = m1(1,1)*m1(3,3)-m1(1,3)*m1(3,1)
m2(2,3) = m1(1,3)*m1(2,1)-m1(1,1)*m1(2,3)
m2(3,1) = m1(2,1)*m1(3,2)-m1(2,2)*m1(3,1)
m2(3,2) = m1(1,2)*m1(3,1)-m1(1,1)*m1(3,2)
m2(3,3) = m1(1,1)*m1(2,2)-m1(1,2)*m1(2,1)

detm = m1(1,1)*m1(2,2)*m1(3,3) + m1(1,2)*m1(2,3)*m1(3,1) &
     + m1(1,3)*m1(2,1)*m1(3,2) - m1(1,3)*m1(2,2)*m1(3,1) &
     - m1(1,2)*m1(2,1)*m1(3,3) - m1(1,1)*m1(2,3)*m1(3,2) 

m2(:,:) = m2(:,:)/detm

end subroutine

!--------------------------------------------------------------------------------------------------------------
function l2g(atype)
implicit none
!convert Local ID to Global ID 
!--------------------------------------------------------------------------------------------------------------
real(8),intent(IN) :: atype
integer :: l2g,ity

ity = nint(atype)
l2g = nint((atype-ity)*1d13)

return
end function

!--------------------------------------------------------------------------------------------------------------
subroutine xu2xs(rreal, rnorm)
! update normalized coordinate from real coordinate. Subtract obox to make them local. 
use atoms
real(8),intent(in) :: rreal(3,NBUFFER)
real(8),intent(out) :: rnorm(3,NBUFFER)

!--------------------------------------------------------------------------------------------------------------
real(8) :: rr(3)

do i=1, NBUFFER
   rr(1:3) = rreal(1:3,i)
   rnorm(1,i)=sum(HHi(1,1:3)*rr(1:3))
   rnorm(2,i)=sum(HHi(2,1:3)*rr(1:3))
   rnorm(3,i)=sum(HHi(3,1:3)*rr(1:3))
   rnorm(1:3,i) = rnorm(1:3,i) - OBOX(1:3)
enddo

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine xs2xu(rnorm,rreal)
! update real coordinate from normalized coordinate
use atoms
!--------------------------------------------------------------------------------------------------------------
real(8),intent(in) :: rnorm(3,NBUFFER)
real(8),intent(out) :: rreal(3,NBUFFER)

real(8) :: rr(3)

do i=1, NBUFFER
   rr(1:3) = rnorm(1:3,i) + OBOX(1:3)
   rreal(1,i)=sum(HH(1,1:3,0)*rr(1:3))
   rreal(2,i)=sum(HH(2,1:3,0)*rr(1:3))
   rreal(3,i)=sum(HH(3,1:3,0)*rr(1:3))
enddo

end subroutine

!-----------------------------------------------------------------------
subroutine ScaleTemperature(atype, v)
use atoms; use parameters
!-----------------------------------------------------------------------
implicit none
real(8) :: atype(NBUFFER), v(3,NBUFFER)

integer :: i,ity
real(8) :: Ekinetic, ctmp

do i=1, NATOMS
   ity=nint(atype(i))
   Ekinetic=0.5d0*mass(ity)*sum(v(1:3,i)*v(1:3,i))
   ctmp = (treq*UTEMP0)/( Ekinetic*UTEMP )
   v(1:3,i)=sqrt(ctmp)*v(1:3,i)
enddo

call LinearMomentum(atype, v)

return
end

!-----------------------------------------------------------------------
subroutine LinearMomentum(atype, v)
use atoms; use parameters
!-----------------------------------------------------------------------
implicit none

real(8) :: atype(NBUFFER), v(3,NBUFFER)

integer :: i,ity
real(8) :: mm,vCM(3),sbuf(4),rbuf(4)

!--- get the local momentum and mass.
vCM(:)=0.d0;  mm = 0.d0
do i=1, NATOMS
   ity = nint(atype(i))
   vCM(1:3)=vCM(1:3) + mass(ity)*v(1:3,i)
   mm = mm + mass(ity)
enddo

sbuf(1)=mm; sbuf(2:4)=vCM(1:3)
call MPI_ALLREDUCE (sbuf, rbuf, size(sbuf), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
mm=rbuf(1); vCM(1:3)=rbuf(2:4)

!--- get the global momentum
vCM(:)=vCM(:)/mm

!--- set the total momentum to be zero 
do i=1, NATOMS
   v(1:3,i) = v(1:3,i) - vCM(1:3)
enddo

return
end
