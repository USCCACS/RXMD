!--------------------------------------------------------------------------------------------------------------
subroutine stress(atype, pos, v, f)
use atoms; use parameters
! calculate stress tensor components of kinetic energy part. 
!--------------------------------------------------------------------------------------------------------------
implicit none

real(8) :: atype(NBUFFER),pos(3,NBUFFER),v(3,NBUFFER),f(3,NBUFFER)

integer :: i,ity, j,jty,icmp, m,n
integer :: c1,c2,c3,c4,c5,c6
real(8) :: dr(3),dr2
real(8) :: mh, cstr(0:6,10)
real(8) :: dbuf(6), Gdbuf(6)
!--- cutoff length of stress calculation range
real(8) :: rcstr2 = 5.d0**2

!--- update buffer atom's stress value
dr=4*lcsize(1:3)
call COPYATOMS(MODE_STRESSCALC, dr)

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
