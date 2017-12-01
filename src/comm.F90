!--------------------------------------------------------------------------------------------------------------
subroutine COPYATOMS(imode, dr, atype, rreal, v, f, q)
use atoms; use pqeq_vars
!
! TODO: update notes here
!
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

integer,intent(IN) :: imode 
real(8),intent(IN) :: dr(3)
real(8) :: atype(NBUFFER), q(NBUFFER)
real(8) :: rreal(NBUFFER,3),v(NBUFFER,3),f(NBUFFER,3)

real(8) :: pos(NBUFFER,3) ! <- normalized coordinate

integer :: i,tn1,tn2, dflag
integer :: ni, ity

integer :: ti,tj,tk,tti,ttj

integer,parameter :: dinv(6)=(/2,1,4,3,6,5/)

call system_clock(tti,tk)

!--- Normalized local coordinate will be used through this function. 
!--- How many atom coords need to be normalized depends on who calls this.
!--- For example function calls during QEq (MODE_QCOPY1 & MODE_QCOPY2) assume
!--- atom information of the extended domain regions. 
call xu2xs(rreal,pos,max(NATOMS,copyptr(6)))

!--- clear total # of copied atoms, sent atoms, recieved atoms
na=0;ns=0;nr=0

!--- REMARK: note that MODE_CPBK depends on copyptr() generated during MODE_COPY.
!--- Since cached atoms are stored after the resident atoms (i.e. i > NATOMS), 
!--- initialize the cache atom pointer 0th elem with NATOMS.
copyptr(0)=NATOMS

!--- set the number of data per atom 
select case(imode)
   case(MODE_COPY)
      ne = NE_COPY
   case(MODE_MOVE)
      ne = NE_MOVE
   case(MODE_CPBK)
      ne = NE_CPBK
   case(MODE_QCOPY1)
      ne = NE_QCOPY1
   case(MODE_QCOPY2)
      ne = NE_QCOPY2
   case(MODE_STRESSCALC)
      ne = NE_STRESSCALC
   case default
      print'(a,i3)', "ERROR: imode doesn't match in COPYATOMS: ", imode
end select

do dflag=1, 6
   tn1 = target_node(dflag)
   tn2 = target_node(dinv(dflag))
   i = (dflag+1)/2  !<- [123]

   if(imode==MODE_CPBK) then  ! communicate with neighbors in reversed order
      tn1 = target_node(7-dinv(dflag)) ! <-[563412] 
      tn2 = target_node(7-dflag) ! <-[654321] 
      i = (6-dflag)/2 + 1         ! <-[321]
   endif

   call store_atoms(tn1, dflag, imode, dr)
   call send_recv(tn1, tn2, myparity(i))
   call append_atoms(dflag, imode)

enddo


if(imode==MODE_MOVE) then
!--- remove atoms which are transfered to neighbor nodes.
   ni=0
   do i=1, NATOMS + na/ne
      ity = nint(atype(i))
!--- if atype is smaller than zero (this is done in store_atoms), ignore the atom.
      if(ity>0) then
        ni=ni+1
        pos(ni,1:3) = pos(i,1:3)
        v(ni,1:3) = v(i,1:3)
        atype(ni) = atype(i)
        q(ni) = q(i)
        qs(ni) = qs(i)
        qt(ni) = qt(i)
        qsfp(ni) = qsfp(i)
        qsfv(ni) = qsfv(i)
!--- For PQEq
        spos(ni,1:3) = spos(i,1:3)
      endif
   enddo 

!--- update the number of resident atoms
   NATOMS=ni

endif

!--- by here, we got new atom positions in the normalized coordinate, need to update real coordinates.
if(imode== MODE_COPY .or. imode == MODE_MOVE) call xs2xu(pos,rreal,copyptr(6))

!--- for array size stat
if(mod(nstep,pstep)==0) then
  ni=nstep/pstep+1
  if(imode==MODE_MOVE) maxas(ni,4)=na/ne
  if(imode==MODE_COPY) maxas(ni,5)=na/ne
endif

call system_clock(ttj,tk)
it_timer(4)=it_timer(4)+(ttj-tti)

return
CONTAINS 

!--------------------------------------------------------------------------------------------------------------
subroutine send_recv(tn1, tn2, mypar)
use atoms
! shared variables::  <ns>, <nr>, <na>, <sbuffer()>, <rbuffer()>
! This subroutine only takes care of communication part. won't be affected by wether atom migration or atom 
! copy mode. 
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) ::tn1, tn2, mypar 
integer :: recv_stat(MPI_STATUS_SIZE)
real(8) :: recv_size

!--- if the traget node is the node itself, atoms informations are already copied 
!--- to <rbuffer> in <store_atoms()>. Nothing to do here. Just return from this subroutine.

if(myid==tn1) return 

call system_clock(ti,tk)

if (mypar == 0) then

     ! the number of elements per data packet has to be greater than 1, for example NE_COPY = 10.
     ! if ns == 0, send one double to tell remote rank that there will be no atom data to be sent. 
     if (ns > 0) then 
       call MPI_SEND(sbuffer, ns, MPI_DOUBLE_PRECISION, tn1, 10, MPI_COMM_WORLD, ierr)
     else
       call MPI_SEND(1, 1, MPI_DOUBLE_PRECISION, tn1, 10, MPI_COMM_WORLD, ierr)
     endif

     call MPI_Probe(tn2, 11, MPI_COMM_WORLD, recv_stat, ierr)
     call MPI_Get_count(recv_stat, MPI_DOUBLE_PRECISION, nr, ierr)

     call CheckSizeThenReallocate(rbuffer,nr)

     call MPI_RECV(rbuffer, nr, MPI_DOUBLE_PRECISION, tn2, 11, MPI_COMM_WORLD, recv_stat, ierr)

     ! the number of elements per data packet has to be greater than 1, for example NE_COPY = 10.
     ! nr == 1 means no atom data to be received. 
     if(nr==1) nr=0 

elseif (mypar == 1) then

     call MPI_Probe(tn2, 10, MPI_COMM_WORLD, recv_stat, ierr)
     call MPI_Get_count(recv_stat, MPI_DOUBLE_PRECISION, nr, ierr)

     call CheckSizeThenReallocate(rbuffer,nr)

     call MPI_RECV(rbuffer, nr, MPI_DOUBLE_PRECISION, tn2, 10, MPI_COMM_WORLD, recv_stat, ierr)

     ! the number of elements per data packet has to be greater than 1, for example NE_COPY = 10.
     ! nr == 1 means no atom data to be received. 
     if(nr==1) nr=0 

     if (ns > 0) then
        call MPI_SEND(sbuffer, ns, MPI_DOUBLE_PRECISION, tn1, 11, MPI_COMM_WORLD, ierr)
     else
        call MPI_SEND(1, 1, MPI_DOUBLE_PRECISION, tn1, 11, MPI_COMM_WORLD, ierr)
     endif

endif

call system_clock(tj,tk)
it_timer(25)=it_timer(25)+(tj-ti)

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine store_atoms(tn, dflag, imode, dr)
use atoms
! <nlayer> will be used as a flag to change the behavior of this subroutine. 
!    <nlayer>==0 migration mode
!            > 0 copy mode
! shared variables::  <ns>, <nr>, <na>, <ne>, <sbuffer()>, <rbuffer()>
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: tn, dflag, imode 
real(8),intent(IN) :: dr(3)

integer :: n,ni,is
real(8) :: sft 
integer :: cptridx

call system_clock(ti,tk)

!--- reset the number of atoms to be sent
ns=0

cptridx=((dflag-1)/2)*2 ! <- [002244]

if(imode/=MODE_CPBK) then

!--- # of elements to be sent. should be more than enough. 
   ni = copyptr(cptridx)*ne

!--- <sbuffer> will deallocated in store_atoms.
   call CheckSizeThenReallocate(sbuffer,ni)

!--- get the coordinate Index to be Shifted.
   is = int((dflag-1)/2) !<- [012] means [xyz]

!--- When atom moves to neighbor nodes, their coordinates must be shifted. 
!--- xshift() returns the edge length of one node assuming all of node size is same.
   sft=xshift(dflag)

!--- start buffering data depending on modes. all copy&move modes use buffer size, dr, to select atoms.
   do n=1, copyptr(cptridx)

      if(inBuffer(dflag,n,dr)) then

        select case(imode)
        case(MODE_MOVE)
           sbuffer(ns+1:ns+3) = pos(n,1:3)
           sbuffer(ns+1+is) = sbuffer(ns+1+is) + sft
           sbuffer(ns+4:ns+6) = v(n,1:3)
           sbuffer(ns+7) = atype(n)
           sbuffer(ns+8) = q(n)
           sbuffer(ns+9) = qs(n)
           sbuffer(ns+10) = qt(n)
           sbuffer(ns+11) = qsfp(n)
           sbuffer(ns+12) = qsfv(n)
           sbuffer(ns+13:ns+15) = spos(n,1:3)
  
!--- In append_atoms subroutine, atoms with <atype>==-1 will be removed
           atype(n) = -1.d0 

        case(MODE_COPY)
           sbuffer(ns+1:ns+3) = pos(n,1:3)
           sbuffer(ns+1+is) = sbuffer(ns+1+is) + sft
           sbuffer(ns+4) = atype(n)
           sbuffer(ns+5) = q(n)
           sbuffer(ns+6) = dble(n)
           sbuffer(ns+7) = qs(n)
           sbuffer(ns+8) = qt(n)
           sbuffer(ns+9) = hs(n)
           sbuffer(ns+10) = ht(n)
           sbuffer(ns+11:ns+13) = spos(n,1:3)

        case(MODE_QCOPY1)
           sbuffer(ns+1) = qs(n)
           sbuffer(ns+2) = qt(n)

        case(MODE_QCOPY2)
           sbuffer(ns+1) = hs(n)
           sbuffer(ns+2) = ht(n)
           sbuffer(ns+3) = q(n)

        case(MODE_STRESSCALC)
           sbuffer(ns+1:ns+6) = astr(1:6,n)

        end select 

!--- increment the number of atoms to be sent 
        ns = ns + ne
     endif

   enddo

!===== FORCE COPYBACK MODE ===========================================================
else if(imode==MODE_CPBK) then

   is = 7 - dflag !<- [654321] reversed order direction flag

   n = copyptr(is) - copyptr(is-1) + 1
   call CheckSizeThenReallocate(sbuffer,n*ne)

   do n=copyptr(is-1)+1, copyptr(is)
      sbuffer(ns+1) = dble(frcindx(n))
      sbuffer(ns+2:ns+4) = f(n,1:3)
#ifdef STRESS
      sbuffer(ns+5:ns+10) = astr(1:6,n)
#endif
!--- chenge index to point next atom.
      ns=ns+ne
   enddo
!=========================================================== FORCE COPYBACK MODE ====
endif

!--- if myid is the same of target-node ID, don't use MPI call.
!--- Just copy <sbuffer> to <rbuffer>. Because <send_recv()> will not be used,
!--- <nr> has to be updated here for <append_atoms()>.
if(myid==tn) then
   if(ns>0) then
      nr=ns
      call CheckSizeThenReallocate(rbuffer,nr)
      rbuffer(1:ns) = sbuffer(1:ns)
   else
      nr=0
   endif
endif

call system_clock(tj,tk)
it_timer(26)=it_timer(26)+(tj-ti)

end subroutine store_atoms

!--------------------------------------------------------------------------------------------------------------
subroutine append_atoms(dflag, imode)
use atoms
! <append_atoms> append copied information into arrays
! shared variables::  <ns>, <nr>, <na>, <ne>, <sbuffer()>, <rbuffer()>
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: dflag, imode 
integer :: m, i, ine

call system_clock(ti,tk)

if( (na+nr)/ne > NBUFFER) then
    print'(a,i4,5i8)', "ERROR: over capacity in append_atoms; myid,na,nr,ne,(na+nr)/ne,NBUFFER: ", &
         myid,na,nr,ne,(na+nr)/ne, NBUFFER
    call MPI_FINALIZE(ierr)
    stop
endif

if(imode /= MODE_CPBK) then  

!--- go over the buffered atom
   do i=0, nr/ne-1

!--- get current index <ine> in <rbuffer(1:nr)>.
      ine=i*ne

!--- current atom index; resident + 1 + stored atoms so far + atoms in buffer
      m = NATOMS + 1 + na/ne + i

      select case(imode)
         case(MODE_MOVE)
              pos(m,1:3) = rbuffer(ine+1:ine+3)
              v(m,1:3) = rbuffer(ine+4:ine+6)
              atype(m) = rbuffer(ine+7)
              q(m)  = rbuffer(ine+8)
              qs(m) = rbuffer(ine+9)
              qt(m) = rbuffer(ine+10)
              qsfp(m) = rbuffer(ine+11)
              qsfv(m) = rbuffer(ine+12)
              spos(m,1:3) = rbuffer(ine+13:ine+15)
      
         case(MODE_COPY)
              pos(m,1:3) = rbuffer(ine+1:ine+3)
              atype(m) = rbuffer(ine+4)
              q(m)  = rbuffer(ine+5)
              frcindx(m) = nint(rbuffer(ine+6))
              qs(m) = rbuffer(ine+7)
              qt(m) = rbuffer(ine+8)
              hs(m) = rbuffer(ine+9)
              ht(m) = rbuffer(ine+10)
              spos(m,1:3) = rbuffer(ine+11:ine+13)
      
           case(MODE_QCOPY1)
              qs(m) = rbuffer(ine+1)
              qt(m) = rbuffer(ine+2)
      
           case(MODE_QCOPY2)
              hs(m) = rbuffer(ine+1)
              ht(m) = rbuffer(ine+2)
              q(m)  = rbuffer(ine+3)

           case(MODE_STRESSCALC)
              astr(1:6,m) = rbuffer(ine+1:ine+6)
      
      end select

     enddo

!===== FORCE COPYBACK MODE =============================================================
else if(imode == MODE_CPBK) then

   do i=0, nr/ne-1
!--- get current index <ine> in <rbuffer(1:nr)>.
      ine=i*ne
!--- Append the transferred forces into the original position of force array.
      m = nint(rbuffer(ine+1))
      f(m,1:3) = f(m,1:3) + rbuffer(ine+2:ine+4)
#ifdef STRESS
      astr(1:6,m) = astr(1:6,m) + rbuffer(ine+5:ine+10)
#endif
   enddo

endif   
!============================================================== FORCE COPYBACK MODE  ===

!--- store the last transfered atom index to <dflag> direction.
!--- if no data have been received, use the index of previous direction.
if(imode /= MODE_CPBK) then
  if(nr==0) m = copyptr(dflag-1)
  copyptr(dflag) = m
endif

!--- update the total # of transfered elements.
na=na+nr

call system_clock(tj,tk)
it_timer(27)=it_timer(27)+(tj-ti)

end subroutine append_atoms

!--------------------------------------------------------------------------------------------------------------
real(8) function xshift(dflag)
use atoms
! This subroutine returns a correction value in coordinate for transfered atoms.
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: dflag 
integer :: i, j
  
  i = mod(dflag,2)  !<- i=1 means positive, i=0 means negative direction 
  j = int((dflag+1)/2)  !<- [123] means [xyz]

  if(i==1) then
      xshift =-lbox(j)
  else if(i==0) then
      xshift = lbox(j)
  endif
   
!print'(a,i,3f10.3)',"shift: ",myid, xshift
end function

!--------------------------------------------------------------------------------------------------------------
function inBuffer(dflag, idx, dr) result(isInside)
use atoms
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: dflag, idx
real(8),intent(IN) :: dr(3)
integer :: i
logical :: isInside

i = mod(dflag,2)  !<- i=1 means positive, i=0 means negative direction 
select case(dflag)
   case(1) 
      isInside = lbox(1) - dr(1) < pos(idx,1)
   case(2) 
      isInside = pos(idx,1) <= dr(1)
   case(3) 
      isInside = lbox(2) - dr(2) < pos(idx,2)
   case(4) 
      isInside = pos(idx,2) <= dr(2)
   case(5) 
      isInside = lbox(3) - dr(3) < pos(idx,3)
   case(6) 
      isInside = pos(idx,3) <= dr(3)
   case default
      write(6,*) "ERROR: no matching directional flag in isInside: ", dflag
end select

end function

!--------------------------------------------------------------------------------------------------------------
subroutine CheckSizeThenReallocate(buffer,nsize)
implicit none
!--------------------------------------------------------------------------------------------------------------
real(8),allocatable :: buffer(:)
integer,intent(IN) :: nsize
integer :: ast

if(allocated(buffer)) then
   if(nsize > size(buffer)) then
       deallocate(buffer)
       allocate(buffer(2*nsize), stat=ast)
   endif
else
   allocate(buffer(nsize), stat=ast)
endif

end subroutine

end subroutine COPYATOMS
