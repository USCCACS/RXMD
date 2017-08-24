!--------------------------------------------------------------------------------------------------------------
subroutine stress()
use atoms; use parameters
! calculate stress tensor components of kinetic energy part. 
!--------------------------------------------------------------------------------------------------------------
implicit none
integer :: i,ity, j,jty,icmp, m,n
integer :: c1,c2,c3,c4,c5,c6
real(8) :: dr(3),dr2
real(8) :: mh, cstr(0:6,10)
real(8) :: dbuf(6), Gdbuf(6)
!--- cutoff length of stress calculation range
real(8) :: rcstr2 = 5.d0**2

!--- update buffer atom's stress value
call COPYATOMS_str(4)

!--- get the potential contribution of internal pressure 
do i=1, NATOMS
  ity = atype(i)
  mh = mass(ity)
  xx = mh*v(1,i)*v(1,i)
  yy = mh*v(2,i)*v(2,i)
  zz = mh*v(3,i)*v(3,i)
  yz = mh*v(2,i)*v(3,i)
  zx = mh*v(1,i)*v(3,i)
  xy = mh*v(1,i)*v(2,i)

!--- one step values.
  astr(1,i) = astr(1,i) + xx
  astr(2,i) = astr(2,i) + yy
  astr(3,i) = astr(3,i) + zz
  astr(4,i) = astr(4,i) + yz
  astr(5,i) = astr(5,i) + zx
  astr(6,i) = astr(6,i) + xy
enddo

dbuf(1:6)=0.d0
do i=1, NATOMS
   dbuf(1:6)=dbuf(1:6)+astr(1:6,i)
enddo

call MPI_ALLREDUCE(dbuf,Gdbuf,size(dbuf),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

dbuf(1:6)=Gdbuf(1:6)/MDBOX
pint(1,1)=dbuf(1) !xx
pint(2,2)=dbuf(2) !yy
pint(3,3)=dbuf(3) !zz
pint(3,1)=dbuf(4) !zx
pint(1,3)=dbuf(4) !xz
pint(2,3)=dbuf(5) !yz
pint(3,2)=dbuf(5) !zy
pint(1,2)=dbuf(6) !xy
pint(2,1)=dbuf(6) !yx

return
end subroutine
!--------------------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------------------
subroutine COPYATOMS_str(nlayer)
! subset of <COPYATOMS>. Just let buffer atoms have stress value. 
use atoms 
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: nlayer
integer :: i,tn, dflag, parity

na=0;ns=0;nr=0

!--- set Nr of elements during this communication. 
ne=6
do dflag=1, 6
   tn = target_node(dflag)
   i=(dflag+1)/2  !<- [123]
   parity=myparity(i)

   call store_atoms_str(tn, dflag, parity, nlayer)
   call send_recv(tn, dflag, parity)
   call append_atoms_str(dflag, parity, nlayer)
enddo

end subroutine 

!--------------------------------------------------------------------------------------------------------------
subroutine store_atoms_str(tn, dflag, parity, nlayer)
use atoms 
! shared variables::  <ns>, <nr>, <na>, <ne>, <sbuffer()>, <rbuffer()>
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: tn, dflag, parity, nlayer
integer :: i,j,k,m,n,n1,ni,is,l(3,2)

if(parity==0) l(1:3, 1:2)=scl1(1:3, 1:2, dflag, nlayer)
if(parity==1) l(1:3, 1:2)=scl2(1:3, 1:2, dflag, nlayer)
ns=sum(nacell(l(1,1):l(1,2), l(2,1):l(2,2), l(3,1):l(3,2)))
ns=ns*ne
if(ns>0) allocate(sbuffer(ns),stat=ast)

ni=0
do i=l(1,1), l(1,2)
do j=l(2,1), l(2,2)
do k=l(3,1), l(3,2)

   n=header(i,j,k)
   do n1=1, nacell(i,j,k)
      sbuffer(ni+1:ni+6) = astr(1:6,n)
      ni=ni+ne
      n=llist(n)
   enddo
enddo; enddo; enddo

if(myid==tn) then
  if(ns>0) then
     nr=ns
     allocate(rbuffer(nr), stat=ast)
     rbuffer(:) = sbuffer(:)
  else
     nr=0
  endif
endif

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine  append_atoms_str(dflag, parity, nlayer)
use atoms
!--------------------------------------------------------------------------------------------------------------
implicit none
integer,intent(IN) :: dflag, parity, nlayer
integer :: m, i, ine

do i=0, nr/ne-1
   ine=i*ne
   m = -(na/ne+i) - 1
   astr(1:6,m) = rbuffer(ine+1:ine+6)
enddo

na=na+nr
llist(0) = na/ne

if(allocated(sbuffer)) deallocate(sbuffer)
if(allocated(rbuffer)) deallocate(rbuffer)

end subroutine 
