!------------------------------------------------------------------------------
program rxmd
use atoms; use parameters
#ifdef INTEROP
use interop 
#endif
!------------------------------------------------------------------------------
implicit none
integer :: i,i1, j,j1, k, n,ity,jty,it1,it2,irt
real(8) :: ctmp
integer :: l2g, igd

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

!--- read ffield file
CALL GETPARAMS()

!--- initialize the MD system
CALL INITSYSTEM()
call OUTPUT(0)

if(mdmode==10) call conjugate_gradient()

call QEq(NCELL10)
call FORCE()

#ifdef INTEROP
!--- A test for Fortran03/C++ interoperability
call print_atom()
#endif

!--- Enter Main MD loop 
call system_clock(it1,irt)
do nstep=0, ntime_step-1

   if(mod(nstep,pstep)==0) call PRINTE()
   if(mod(nstep,fstep)==0) call OUTPUT(2)

   if(mod(nstep,sstep)==0.and.mdmode==4) &
      v(1:3,1:NATOMS)=vsfact*v(1:3,1:NATOMS)

   if(mod(nstep,sstep)==0.and.mdmode==5) then
      ctmp = (treq*UTEMP0)/( GKE*UTEMP )
      v(1:3,1:NATOMS)=sqrt(ctmp)*v(1:3,1:NATOMS)
   endif

   if(mod(nstep,sstep)==0.and.(mdmode==0.or.mdmode==6)) &
      call INITVELOCITY()

!--- correct the c.o.m motion
   if(mod(nstep,sstep)==0.and.mdmode==7) &
      call ScaleTemperature()

!--- update velocity
   call vkick(1.d0) 

!--- update coordinates
   qsfv(1:NATOMS)=qsfv(1:NATOMS)+0.5d0*dt*Lex_w2*(q(1:NATOMS)-qsfp(1:NATOMS))
   qsfp(1:NATOMS)=qsfp(1:NATOMS)+dt*qsfv(1:NATOMS)

   pos(1:3,1:NATOMS)=pos(1:3,1:NATOMS)+dt*v(1:3,1:NATOMS)

!--- migrate atoms after positions are updated
   call xu2xs()
   call LINKEDLIST()
   call COPYATOMS(0)
   call xs2xu()
   
   if(mod(nstep,qstep)==0) call QEq(NCELL10)
   call FORCE()


!--- update velocity
   call vkick(1.d0) 
   qsfv(1:NATOMS)=qsfv(1:NATOMS)+0.5d0*dt*Lex_w2*(q(1:NATOMS)-qsfp(1:NATOMS))

enddo

!--- save the final configurations
call OUTPUT(-1)

call system_clock(it2,irt)
it_timer(Ntimer)=(it2-it1)
call MPI_ALLREDUCE (it_timer, it_timer_max, Ntimer, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
call MPI_ALLREDUCE (it_timer, it_timer_min, Ntimer, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)

!--- simulaiton finished w/o problem. Printing some information, array sizes, timings. 
allocate(ibuf(nmaxas),ibuf1(nmaxas))
ibuf(:)=0
do i=1,nmaxas
   ibuf(i)=maxval(maxas(:,i))
   !if(myid==0) print'(a,20i)','i,maxas: ', i,maxas(:,i)
enddo
call MPI_ALLREDUCE (ibuf, ibuf1, nmaxas, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
if(myid==0) then
   !print'(a,10i12)', '1. max array size: ', ibuf1(1:nmaxas)
   print'(a,10i12)', 'Max MAXNEIGHBS, Max MAXNEIGHBS10, Max NBUFFER_P, Max NBUFFER_N: ', &
                      ibuf1(2), ibuf1(3), ibuf1(1)+ibuf1(4), ibuf1(5)
endif

!do i=1,nmaxas
!   ibuf(i)=minval(maxas(:,i))
!enddo
!call MPI_ALLREDUCE (ibuf, ibuf1, nmaxas, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
!if(myid==0) then
!   !print'(a,10i12)', '2. min array size: ', ibuf1(1:nmaxas)
!   print'(a,10i12)', 'Min MAXNEIGHBS, Min MAXNEIGHBS10, Min NBUFFER_P, Min NBUFFER_N: ', &
!                      ibuf1(2), ibuf1(3), ibuf1(1)+ibuf1(4), ibuf1(5)
!endif

deallocate(ibuf,ibuf1)

if(myid==0) then
   print'(a20,f12.4,3x,f12.4)','QEq: ',  dble(it_timer_max(1))/irt, dble(it_timer_min(1))/irt
   print'(a20,f12.4,3x,f12.4)','QEq_COPY: ',  dble(it_timer_max(2))/irt, dble(it_timer_min(2))/irt
   print'(a20,f12.4,3x,f12.4)','LINKEDLIST: ',  dble(it_timer_max(3))/irt, dble(it_timer_min(3))/irt
   print'(a20,f12.4,3x,f12.4)','COPYATOMS: ',    dble(it_timer_max(4))/irt, dble(it_timer_min(4))/irt
   print'(a20,f12.4,3x,f12.4)','NEIGHBORLIST: ', dble(it_timer_max(5))/irt, dble(it_timer_min(5))/irt
   print'(a20,f12.4,3x,f12.4)','BOCALC: ', dble(it_timer_max(6))/irt, dble(it_timer_min(6))/irt
   print'(a20,f12.4,3x,f12.4)','ENbond: ', dble(it_timer_max(7))/irt, dble(it_timer_min(7))/irt
   print'(a20,f12.4,3x,f12.4)','Ebond: ', dble(it_timer_max(8))/irt, dble(it_timer_min(8))/irt
   print'(a20,f12.4,3x,f12.4)','Elnpr: ', dble(it_timer_max(9))/irt, dble(it_timer_min(9))/irt
   print'(a20,f12.4,3x,f12.4)','Ehb: ', dble(it_timer_max(10))/irt, dble(it_timer_min(10))/irt
   print'(a20,f12.4,3x,f12.4)','E3b: ', dble(it_timer_max(11))/irt, dble(it_timer_min(11))/irt
   print'(a20,f12.4,3x,f12.4)','E4b: ', dble(it_timer_max(12))/irt, dble(it_timer_min(12))/irt
   print'(a20,f12.4,3x,f12.4)','ForceBondedTerms: ', dble(it_timer_max(13))/irt, dble(it_timer_min(13))/irt
   print'(a20,f12.4,3x,f12.4)','COPYATOMS(-1): ', dble(it_timer_max(14))/irt, dble(it_timer_min(14))/irt
   print'(a20,f12.4,3x,f12.4)','total (sec): ',dble(it_timer_max(Ntimer))/irt, dble(it_timer_min(Ntimer))/irt
   print'(a20)', 'program finished'
endif

call MPI_FINALIZE(ierr)
end PROGRAM

!------------------------------------------------------------------------------
subroutine vkick(dtf)
use atoms
!------------------------------------------------------------------------------
implicit none
integer :: i, ity
real(8) :: dtf

do i=1,NATOMS
   ity = atype(i)
   v(1:3,i) = v(1:3,i) + dtf*dthm(ity)*f(1:3,i)
enddo

end subroutine


!----------------------------------------------------------------------------------------
subroutine PRINTE()
use atoms; use parameters
! calculate the kinetic energy and sum up all of potential energies, then print them.
!----------------------------------------------------------------------------------------
implicit none
integer :: i,j,ity, cstep
real(8) :: qq=0.d0,tt=0.d0,ss=0.d0,buf(0:20),Gbuf(0:20)

i=nstep/pstep+1
maxas(i,1)=NATOMS

KE=0.d0
do i=1, NATOMS
   ity=atype(i)
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

   write(6,'(i9,3es13.5,6es11.3,1x,3f8.2,i4,f8.2)') cstep,GTE,GPE(0),GKE, &
   GPE(1), sum(GPE(2:3)), sum(GPE(3:10)), GPE(11:13), &
   tt, ss, qq, nstep_qeq, MPI_WTIME()-wt0 

#ifdef STRESS
   write(6,'(6es13.5)') pint(1,1)*USTRS, pint(2,2)*USTRS, pint(3,3)*USTRS, &
                        pint(2,3)*USTRS, pint(3,1)*USTRS, pint(1,2)*USTRS
#endif

endif

!--- save current time
wt0 = MPI_WTIME()

end subroutine


!----------------------------------------------------------------------------------------
subroutine OUTPUT(imode)
use atoms; use parameters
! imode == 0: finalize
!
!----------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: imode 
integer :: i, ity, j, j1, jty, m, n,n3, cstep
integer :: l2g
real(8) :: ri(3), rj(3), bndordr(MAXNEIGHBS), tt=0.d0, ss=0.d0
integer :: igd,jgd,bndlist(0:MAXNEIGHBS)
character(8) :: fname0
character(6) :: a6
character(9) :: a9
character(128) :: FileName

integer (kind=MPI_OFFSET_KIND) :: offsetIO
integer :: localDataSize
integer :: fh ! file handler

integer,parameter :: PDBLineSize=67
character(PDBLineSize) :: PDBOneLine

integer :: BNDLineSize
integer,parameter :: MaxBNDLineSize=512
character(MaxBNDLineSize) :: BNDOneLine
real(8),parameter :: BNDcutoff=0.3d0

write(a6(1:6),'(i6.6)') myid
write(a9(1:9),'(i9.9)') nstep + current_step
FileName=trim(DataPath)//"/"//a6//"/rxff"//a6//"-"//a9

!--- binary ------------------------------------------------------------------
if(isBinary) then
  call xu2xs()
  call coio_write(imode)
  call xs2xu()
endif
!------------------------------------------------------------------ binary ---


!--- BondFile -------------------------------------------------------------
if(isBondFile) then

    ! precompute the total # of neighbors
    m=0
    do i=1, NATOMS
        do j1 = 1, nbrlist(i,0)
!--- don't count if BO is less than BNDcutoff.
            if(BO(0,i,j1) > BNDcutoff) then 
                m=m+1
            endif
        enddo
    enddo

200 format(i9,1x,3f9.3,1x,2i3,20(i9,f6.3)) 

    ! get local datasize based on above format and the total # of neighbors
    localDataSize=NATOMS*(9+1+3*9+1+2*3 +1)+m*(9+6)

    ! offsetIO will point the end of local write after the scan
    call MPI_Scan(localDataSize,offsetIO,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

    ! set offsetIO at the beginning of the local write
    offsetIO=offsetIO-localDataSize

    call MPI_File_Open(MPI_COMM_WORLD,trim(DataPath)//"/"//a9//".bnd", &
        MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

   !open(10,file=trim(FileName)//".bnd")

   do i=1, NATOMS
      ity = atype(i)
!--- get global ID for i-atom
      igd = l2g(atype(i))

!--- count the number bonds to be shown.
      bndlist(0)=0
      do j1 = 1, nbrlist(i,0)
         j = nbrlist(i,j1)
         jty = atype(j)

!--- get global ID for j-atom
         jgd = l2g(atype(j))

!--- if bond order is less than 0.3, ignore the bond.
         if( BO(0,i,j1) < 0.3d0) cycle

         bndlist(0) = bndlist(0) + 1
         bndlist(bndlist(0)) = jgd
         bndordr(bndlist(0)) = BO(0,i,j1)
      enddo

        BNDOneLine=""
        write(BNDOneLine,200) igd, pos(1:3,i),int(atype(i)),bndlist(0), &
            (bndlist(j1),bndordr(j1),j1=1,bndlist(0))
        ! remove space and add new_line
        BNDOneLine=trim(BNDOneLine)//NEW_LINE('A')
        BNDLineSize=len(trim(BNDOneLine))

        call MPI_File_Seek(fh,offsetIO,MPI_SEEK_SET,ierr)
        call MPI_File_Write(fh,BNDOneLine,BNDLineSize, &
            MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)

        offsetIO=offsetIO+BNDLineSize

   enddo
   !close(10)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_File_Close(fh,ierr)

endif
!------------------------------------------------------------ BondFile ----

!--- PDB ---------------------------------------------------------------------
if(isPDB) then

    ! get local datasize
    localDataSize=NATOMS*PDBLineSize

    ! offsetIO will point the end of local write after the scan
    call MPI_Scan(localDataSize,offsetIO,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

    if(myid==nprocs-1) print*,'offsetIO',offsetIO

    ! set offsetIO at the beginning of the local write
    offsetIO=offsetIO-localDataSize

    call MPI_File_Open(MPI_COMM_WORLD,trim(DataPath)//"/"//a9//".pdb", &
        MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

    do i=1, NATOMS

      ity = atype(i)
!--- calculate atomic temperature 
      tt = hmas(ity)*sum(v(1:3,i)*v(1:3,i))
      tt = tt*UTEMP*1d-2 !scale down to use two decimals in PDB format 

!--- sum up diagonal atomic stress components 
#ifdef STRESS
      ss = sum(astr(1:3,i))/3.d0
#endif
      ss = ss*USTRS

      ss = q(i)*10 ! 10x atomic charge

      igd = l2g(atype(i))
      select case(ity)
        case(1) 
          write(PDBOneLine,100)'ATOM  ',0, 'C', igd, pos(1:3,i), tt, ss
        case(2) 
          write(PDBOneLine,100)'ATOM  ',0, 'H', igd, pos(1:3,i), tt, ss
        case(3) 
          write(PDBOneLine,100)'ATOM  ',0, 'O', igd, pos(1:3,i), tt, ss
        case(4) 
          write(PDBOneLine,100)'ATOM  ',0, 'N', igd, pos(1:3,i), tt, ss
        case(5) 
          write(PDBOneLine,100)'ATOM  ',0, 'S', igd, pos(1:3,i), tt, ss
        case(6) 
          write(PDBOneLine,100)'ATOM  ',0,'Si', igd, pos(1:3,i), tt, ss
        case(7) 
          write(PDBOneLine,100)'ATOM  ',0,'Al', igd, pos(1:3,i), tt, ss
      end select
        PDBOneLine(PDBLineSize:PDBLineSize)=NEW_LINE('A')

        call MPI_File_Seek(fh,offsetIO,MPI_SEEK_SET,ierr)
        call MPI_File_Write(fh,PDBOneLine,PDBLineSize, &
            MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)

        offsetIO=offsetIO+PDBLineSize
    enddo

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_File_Close(fh,ierr)
endif
100 format(A6,I5,1x,A2,i12,4x,3f8.3,f6.2,f6.2)
!-------------------------------------------------------------------- PDB ----

end subroutine

!----------------------------------------------------------------------------------------
subroutine LINKEDLIST()
use atoms
! partitions the volume into linked-list cells with minimum cless size <lcsize>
!----------------------------------------------------------------------------------------
implicit none
integer :: n,n0,l(3), c1,c2,c3, j

header(:,:,:) = -1; llist(:) = 0; nacell(:,:,:)=0

do n=1, NATOMS
   l(1:3) = pos(1:3,n)/lcsize(1:3)
   do j=1,3
      if(pos(j,n)<0.d0) l(j) = l(j) - 1
   enddo
   llist(n) = header(l(1), l(2), l(3))
   header(l(1), l(2), l(3)) = n
   nacell(l(1), l(2), l(3)) = nacell(l(1), l(2), l(3)) + 1
enddo

end subroutine 


!----------------------------------------------------------------------
subroutine NEIGHBORLIST(nlayer)
! calculate neighbor list for atoms witin cc(1:3, -nlayer:nlayer) cells.
use atoms; use parameters
!----------------------------------------------------------------------
implicit none
integer,intent(IN) :: nlayer
integer :: c1,c2,c3, ic(3), c4, c5, c6
integer :: n, n1, m, m1, nty, mty, inxn
real(8) :: dr(3), dr2
integer :: l2g

nbrlist(:,:) = 0
nbrindx(:,:) = 0

DO c1=-nlayer, cc(1)-1+nlayer
DO c2=-nlayer, cc(2)-1+nlayer
DO c3=-nlayer, cc(3)-1+nlayer

  m = header(c1, c2, c3)
  do m1=1, nacell(c1, c2, c3)
     mty = atype(m)

     do c4 = -1, 1
     do c5 = -1, 1
     do c6 = -1, 1
        ic(1:3) = (/c1, c2, c3/) + (/c4, c5, c6/)

        n = header(ic(1),ic(2),ic(3))
        do n1=1, nacell(ic(1), ic(2), ic(3))

           if(n<m) then
             nty = atype(n)
             inxn = inxn2(mty, nty)

             dr(1:3) = pos(1:3,n) - pos(1:3,m) 
             dr2 = sum(dr(1:3)*dr(1:3))

             if(dr2<rc2(inxn)) then 
                nbrlist(m, 0) = nbrlist(m, 0) + 1
                nbrlist(n, 0) = nbrlist(n, 0) + 1
                nbrlist(m, nbrlist(m, 0)) = n
                nbrlist(n, nbrlist(n, 0)) = m
!--- to get the reverse information (i.e. from i,j1&j to i1), store <i1> into <nbrindx>.
                nbrindx(m, nbrlist(m, 0)) = nbrlist(n, 0)
                nbrindx(n, nbrlist(n, 0)) = nbrlist(m, 0)
             endif 
           endif
           n=llist(n) 
        enddo
     enddo; enddo; enddo

     m = llist(m)
  enddo
enddo; enddo; enddo

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

!do n=1,NATOMS
!  print*,'------------------------------------------'
!  print*,l2g(atype(n)),nbrlist(n,0),  nbrlist(n,-1), atype(n)
!  print*,(mod(l2g(atype(nbrlist(n,m)))-1,168)+1,m=1,nbrlist(n,0) )
!  print*,int(atype(nbrlist(n,1:nbrlist(n,0))))
!  print*,nbrindx(n,1:nbrlist(n,0))
!enddo
!stop 'neighbor'

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine angular_momentum()
use atoms; use parameters
!--------------------------------------------------------------------------------------------------------------
implicit none
integer :: i,ity
real(8) :: com(3), Gcom(3), intsr(3,3), Gintsr(3,3), intsr_i(3,3), angm(3), Gangm(3), angv(3), mm, Gmm
real(8) :: dr(3), dv(3)

!--- get center of mass
com(:)=0.d0;     Gcom(:)=0.d0
mm=0.d0; Gmm=0.d0

do i=1, NATOMS
   ity = atype(i) 
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
subroutine COPYATOMS(nlayer)
use atoms
! In this subroutine, boundary atoms are copied to neighbour nodes. Three subroutines, 
! <store_atoms()>, <send_recv()> and  <append_atoms()> work together, shareing variables 
! <sbuffer>, <rbuffer>, <ne>, <na>, <ns> and <nr>. 
!
!--- Variables --- 
!<dflag>:  Direction FLAG specifing the direction of communication.
!          1 +x, 2 -x, 3 +y, 4 -y, 5 +z, 6 -z. 
!<parity>:  PARITY of a node in <dflag> direction. 
!<nlayer>:  Nr of LAYERs to be copied. COPYATOMS() features two different modes, migration mode and copy mode.
!           <nlayer> will be used as a flag to distinguish them. 
!
!--- Shared Variables ---
!<ne>: # of elements one atom has.
!   Example) in copy mode, three positions[xyz] + charge + atomtype = 5 elements. 
!<na>, <ns> & <nr>: Nr of All of transfered elemetns, Number of elements to be Sent, and Recieved one respectivly.
!        
!--- Subroutines ---
!<store_atoms()>: store boundary atoms information.
!<send_recv()):   send & receive. 
!<append_atoms()>: append received infromation into array.  deallocate <sbuffer> & <rbuffer>
!
!--------------------------------------------------------------------------------------------------------------
implicit none

integer,intent(IN) :: nlayer
integer :: i,tn, dflag, parity
integer :: ni, ity
integer :: i1,j1,k1

integer :: c1,c2,c3,n

!--- clear total # of copied atoms, sent atoms, recieved atoms
na=0;ns=0;nr=0

!--- set Nr of elements during this communication. 
if(nlayer==0) ne=NE_MOVE
if(nlayer> 0) ne=NE_COPY
if(nlayer< 0) ne=NE_CPBK

do dflag=1, 6
   tn = target_node(dflag)
   i=(dflag+1)/2  !<- [123]
   parity=myparity(i)

   if(nlayer<0) then            ! communicate with neighbors in reversed order
      tn = target_node(7-dflag) ! <-[654321] 
      i=(6-dflag)/2 + 1         ! <-[321]
      parity=myparity(i)
   endif

   call store_atoms(tn, dflag, parity, nlayer)
   call send_recv(tn, dflag, parity)
   call append_atoms(dflag, parity, nlayer)
enddo

if(nlayer==0) then
!--- remove atoms which are transfered to neighbor nodes.
   ni=0
   do i=1, NATOMS + na/ne
      ity = atype(i)
!--- if atype is smaller than zero (this is done in store_atoms), ignore the atom.
      if(ity>0) then
        ni=ni+1
        pos(1:3,ni) = pos(1:3,i)
        v(1:3,ni) = v(1:3,i)
        atype(ni) = atype(i)
        q(ni) = q(i)
        qs(ni) = qs(i)
        qt(ni) = qt(i)
        qsfp(ni) = qsfp(i)
        qsfv(ni) = qsfv(i)
      endif
   enddo 

!--- Update Nr of local atom
   NATOMS=ni
endif

!--- for array size stat
if(mod(nstep,pstep)==0) then
  ni=nstep/pstep+1
  if(nlayer==0) maxas(ni,4)=na/ne
  !if(nlayer>0) maxas(ni,5)=na/ne
endif

end subroutine COPYATOMS

!--------------------------------------------------------------------------------------------------------------
subroutine send_recv(tn, dflag, parity)
use atoms
! shared variables::  <ns>, <nr>, <na>, <sbuffer()>, <rbuffer()>
! This subroutine only takes care of communication part. won't be affected by wether atom migration or atom 
! copy mode. 
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) ::tn, dflag, parity
integer :: recv_stat(MPI_STATUS_SIZE)

!--- if the traget node is the node itself, atoms informations are already copied 
!--- to <rbuffer> in <store_atoms()>. Nothing to do here. Just return from this subroutine.
if(myid==tn) return 

if (parity == 0) then

     call MPI_SEND(ns, 1, MPI_INTEGER, tn, 10, MPI_COMM_WORLD, ierr)
     if (ns > 0) &
       call MPI_SEND(sbuffer, ns, MPI_DOUBLE_PRECISION, tn, 11, MPI_COMM_WORLD, ierr)

     call MPI_RECV(nr, 1, MPI_INTEGER, tn, 12, MPI_COMM_WORLD, recv_stat, ierr)
     if (nr > 0) then
       allocate(rbuffer(nr), stat=ast); if(ast/=0) print*, "ERROR: rbuffer@send_recv: 1"
       call MPI_RECV(rbuffer, nr, MPI_DOUBLE_PRECISION, tn, &
                     13, MPI_COMM_WORLD, recv_stat, ierr)
     endif

elseif (parity == 1) then

       call MPI_RECV(nr, 1, MPI_INTEGER, tn, 10, MPI_COMM_WORLD, recv_stat, ierr)
       if (nr > 0) then
         allocate(rbuffer(nr), stat=ast); if(ast/=0) print*,"ERROR: rbuffer@send_recv: 2"
         call MPI_RECV(rbuffer, nr, MPI_DOUBLE_PRECISION, tn, &
                      11, MPI_COMM_WORLD, recv_stat, ierr)
       endif

       call MPI_SEND(ns, 1, MPI_INTEGER, tn, 12, MPI_COMM_WORLD, ierr)
       if (ns > 0) &
       call MPI_SEND(sbuffer, ns, MPI_DOUBLE_PRECISION, tn, 13, MPI_COMM_WORLD, ierr)

endif

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine store_atoms(tn, dflag, parity, nlayer)
use atoms
! <nlayer> will be used as a flag to change the behavior of this subroutine. 
!    <nlayer>==0 migration mode
!            > 0 copy mode
! shared variables::  <ns>, <nr>, <na>, <ne>, <sbuffer()>, <rbuffer()>
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: tn, dflag, parity, nlayer
integer :: i,j,k,m,n,n1,ni,is,l(3,2)
real(8) :: sft, xshift

if(nlayer>=0) then
!--- if node parity is even, the node starts communicating in positive direction first. use <scl1>.
!--- if it's odd, negative direction first. use <scl2>. 
   if(parity==0) l(1:3, 1:2)=scl1(1:3, 1:2, dflag, nlayer)
   if(parity==1) l(1:3, 1:2)=scl2(1:3, 1:2, dflag, nlayer)

!--- get the Number of atoms to be Sent in the subjected layer. 
   ns=sum(nacell(l(1,1):l(1,2), l(2,1):l(2,2), l(3,1):l(3,2)))

!--- get the Number of Elements to be sent.
   ns=ns*ne

!--- <sbuffer> will deallocated in store_atoms.
   if(ns>0) allocate(sbuffer(ns),stat=ast)

!--- get the coordinate Index to be Shifted.
   is = int((dflag-1)/2) !<- [012]

!--- When atom moves to neighbor nodes, their coordinates must be shifted. 
!--- function <xshift> returns the edge length of one node <sft>, assuming all of node size
!--- is same.
   sft=xshift(dflag, parity)
endif

!--- Store subjected atoms infromation into <sbuffer>. Layer indices <l(:,:)> is given 
!--- by <scl1> or <scl2>, depending on the node parity. Even parity nodes send 
!--- positive direction layer fast.

!--- Initialize index of array <sbuffer>
ni=0

!====== MIGRATION MODE =================================================================
if(nlayer==0) then 
   do i=l(1,1), l(1,2)
   do j=l(2,1), l(2,2)
   do k=l(3,1), l(3,2)

      n=header(i,j,k)
      do n1=1, nacell(i,j,k)
         sbuffer(ni+1:ni+3) = pos(1:3,n)
         sbuffer(ni+1+is) = sbuffer(ni+1+is) + sft
         sbuffer(ni+4:ni+6) = v(1:3,n)
         sbuffer(ni+7) = atype(n)
         sbuffer(ni+8) = q(n)
         sbuffer(ni+9) = qs(n)
         sbuffer(ni+10) = qt(n)
         sbuffer(ni+11) = qsfp(n)
         sbuffer(ni+12) = qsfv(n)

!--- In append_atoms subroutine, atoms with <atype>==-1 will be removed
         atype(n) = -1.d0 

!--- chenge index to point next atom.
         ni=ni+ne

         n=llist(n)
      enddo
   enddo; enddo; enddo

endif
!================================================================= MIGRATION MODE ====
         
!===== COPY MODE =====================================================================
if(nlayer>0) then 
   do i=l(1,1), l(1,2)
   do j=l(2,1), l(2,2)
   do k=l(3,1), l(3,2)

      n=header(i,j,k)
      do n1=1, nacell(i,j,k)

         sbuffer(ni+1:ni+3) = pos(1:3,n)
         sbuffer(ni+1+is) = sbuffer(ni+1+is) + sft
         sbuffer(ni+4) = atype(n)
         sbuffer(ni+5) = q(n)
         sbuffer(ni+6) = dble(n)
         sbuffer(ni+7) = qs(n)
         sbuffer(ni+8) = qt(n)
         sbuffer(ni+9) = hs(n)
         sbuffer(ni+10) = ht(n)

!--- chenge index to point next atom
         ni=ni+ne

         n=llist(n)
      enddo
   enddo; enddo; enddo

endif
!====================================================================== COPY MODE ===

!===== FORCE COPYBACK MODE ==========================================================
if(nlayer<0) then

   is = 7 - dflag !<- [654321] reversed order direction flag
   allocate(sbuffer(size(frcindx)*ne),stat=ast)

   do n=frcindx(is-1)-1, frcindx(is), -1
      sbuffer(ni+1) = dble(frcindx(n))
      sbuffer(ni+2:ni+4) = f(1:3,n)
#ifdef STRESS
      sbuffer(ni+5:ni+10) = astr(1:6,n)
#endif
!--- chenge index to point next atom.
      ni=ni+ne
   enddo

!--- update Nr of atoms to be sent.
   ns = ni

endif
!=========================================================== FORCE COPYBACK MODE ====


!--- if myid is the same of target-node ID, don't use MPI call.
!--- Just copy <sbuffer> to <rbuffer>. Because <send_recv()> will not be used,
!--- <nr> has to be updated here for <append_atoms()>.
   if(myid==tn) then
     if(ns>0) then
        nr=ns
        allocate(rbuffer(nr), stat=ast)
        rbuffer(:) = sbuffer(:)
     else
        nr=0
     endif
   endif

end subroutine store_atoms

!--------------------------------------------------------------------------------------------------------------
subroutine  append_atoms(dflag, parity, nlayer)
use atoms
! <append_atoms> append copied information into arrays
! shared variables::  <ns>, <nr>, <na>, <ne>, <sbuffer()>, <rbuffer()>
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: dflag, parity, nlayer
integer :: m, i, ine, j, l(3)

if( (na+nr)/ne > NBUFFER_N) then
    print'(a26,i4,5i8)', "over capacity @append_atoms",myid,na,nr,ne,(na+nr)/ne, NBUFFER_N
    call MPI_FINALIZE(ierr)
    stop
endif

!====== MIGRATION MODE =================================================================   
if(nlayer==0) then  

     do i=0, nr/ne-1
!--- get current index <ine> in <rbuffer(1:nr)>.
        ine=i*ne

!--- Transferred atoms will be appended to the positive direction after <NATOMS>. 
        m = NATOMS + 1 + (na/ne+i)

        pos(1:3,m) = rbuffer(ine+1:ine+3)
        v(1:3,m) = rbuffer(ine+4:ine+6)
        atype(m) = rbuffer(ine+7)
        q(m)  = rbuffer(ine+8)
        qs(m) = rbuffer(ine+9)
        qt(m) = rbuffer(ine+10)
        qsfp(m) = rbuffer(ine+11)
        qsfv(m) = rbuffer(ine+12)

        l(1:3)=int(pos(1:3,m)/lcsize(1:3))

!--- if atoms are stored in negative indices cells, these indices must be substracted by 1
!--- to distinguish int(-1:0) and int[0:1).
        do j=1,3 
          if(pos(j,m)<0.d0) l(j) = l(j) - 1
        enddo

        llist(m) = header(l(1), l(2), l(3))
        header(l(1), l(2), l(3)) = m
        nacell(l(1), l(2), l(3)) = nacell(l(1), l(2), l(3)) + 1
     enddo

endif
!=================================================================== MIGRATION MODE ===

!===== COPY MODE ======================================================================   
if(nlayer>0) then 

      do i=0, nr/ne-1
!--- get current index <ine> in <rbuffer(1:nr)>.
         ine=i*ne
!--- Transferred atoms will be appended to the negative direction from <m>=-1. 
         m = -(na/ne+i) - 1

         pos(1:3,m) = rbuffer(ine+1:ine+3)
         atype(m) = rbuffer(ine+4)
         q(m)  = rbuffer(ine+5)
         frcindx(m) = anint(rbuffer(ine+6))
         qs(m) = rbuffer(ine+7)
         qt(m) = rbuffer(ine+8)
         hs(m) = rbuffer(ine+9)
         ht(m) = rbuffer(ine+10)

         l(1:3)=int(pos(1:3,m)/lcsize(1:3))

!--- if atoms are stored in negative indices cells, these indices must be substracted by 1
!--- to distinguish int(-1:0) and int[0:1).
         do j=1,3 
           if(pos(j,m)<0) l(j) = l(j) - 1
         enddo

         llist(m) = header(l(1), l(2), l(3))
         header(l(1), l(2), l(3)) = m
         nacell(l(1), l(2), l(3)) = nacell(l(1), l(2), l(3)) + 1
      enddo

!--- Save the last transfered atom index to <dflag> direction.
      if(nr==0) m = frcindx(dflag-1)
      frcindx(dflag) = m

endif
!========================================================================  COPY MODE  ===


!===== FORCE COPYBACK MODE =============================================================
if(nlayer<0) then

  do i=0, nr/ne-1
!--- get current index <ine> in <rbuffer(1:nr)>.
     ine=i*ne
!--- Append the transferred forces into the original position of force array.
     m = anint(rbuffer(ine+1))
     f(1:3,m) = f(1:3,m) + rbuffer(ine+2:ine+4)
#ifdef STRESS
     astr(1:6,m) = astr(1:6,m) + rbuffer(ine+5:ine+10)
#endif
  enddo

endif   
!============================================================== FORCE COPYBACK MODE  ===

!--- update the total # of transfered elements.
na=na+nr

!--- update the total Nr of transfered atoms.
llist(0) = na/ne

!--- In case node doesn't have a partner node for certain directions, <sbuffer> and <rbuffer> will not be allocated. 
!--- check the allocation status of <sbuffer> and <rbuffer>.
if(allocated(sbuffer)) deallocate(sbuffer)
if(allocated(rbuffer)) deallocate(rbuffer)

end subroutine append_atoms

!--------------------------------------------------------------------------------------------------------------
real(8) function xshift(dflag, parity)
use atoms
! This subroutine returns a correction value in coordinate for transfered atoms.
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: dflag, parity
integer :: i, j
  
  i = mod(dflag,2)  !<- i=1 positive, i=0 negative for even parity nodes
  j = int( (dflag+1)/2 )  !<- [xyz] = [123]
  if(parity==0) then
    if(i==1) xshift =-LBOX(j)
    if(i==0) xshift = LBOX(j)
  else if(parity==1) then
    if(i==1) xshift = LBOX(j)
    if(i==0) xshift =-LBOX(j)
  endif
   
!print'(a,i,3f10.3)',"shift: ",myid, xshift
 
end function

!--------------------------------------------------------------------------------------------------------------
function l2g(atype)
implicit none
!convert Local ID to Global ID 
!--------------------------------------------------------------------------------------------------------------
real(8),intent(IN) :: atype
integer :: l2g,ity

ity = nint(atype)
l2g= anint((atype-ity)*1d13)

return
end function

!--------------------------------------------------------------------------------------------------------------
subroutine xu2xs()
! unscaled coordinate to scaled coordinate. shift coordinate to scaled-local.
use atoms
!--------------------------------------------------------------------------------------------------------------
real(8) :: rr(3)
real(8) :: rx,ry,rz

do i=-NBUFFER_N, NBUFFER_P
   rr(1:3) = pos(1:3,i)
   pos(1,i)=sum(HHi(1,1:3)*rr(1:3))
   pos(2,i)=sum(HHi(2,1:3)*rr(1:3))
   pos(3,i)=sum(HHi(3,1:3)*rr(1:3))
   pos(1:3,i) = pos(1:3,i) - OBOX(1:3)
enddo

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine xs2xu()
! scaled coordinate to unscaled coordinate. shift coordinate to unscaled-global.
use atoms
!--------------------------------------------------------------------------------------------------------------
real(8) :: rr(3)

do i=-NBUFFER_N, NBUFFER_P
   rr(1:3) = pos(1:3,i) + OBOX(1:3)
   pos(1,i)=sum(HH(1,1:3,0)*rr(1:3))
   pos(2,i)=sum(HH(2,1:3,0)*rr(1:3))
   pos(3,i)=sum(HH(3,1:3,0)*rr(1:3))
enddo

end subroutine

!-----------------------------------------------------------------------
subroutine ScaleTemperature()
use atoms; use parameters
!-----------------------------------------------------------------------
implicit none
integer :: i,ity
real(8) :: Ekinetic, ctmp

do i=1, NATOMS
   ity=atype(i)
   Ekinetic=0.5d0*mass(ity)*sum(v(1:3,i)*v(1:3,i))
   ctmp = (treq*UTEMP0)/( Ekinetic*UTEMP )
   v(1:3,i)=sqrt(ctmp)*v(1:3,i)
enddo

call LinearMomentum()

return
end

!-----------------------------------------------------------------------
subroutine LinearMomentum()
use atoms; use parameters
!-----------------------------------------------------------------------
implicit none
integer :: i,ity
real(8) :: mm,vCM(3),sbuf(4),rbuf(4)

!--- get the local momentum and mass.
vCM(:)=0.d0;  mm = 0.d0
do i=1, NATOMS
   ity = atype(i)
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
