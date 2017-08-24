!--------------------------------------------------------------------------------------------------------------
subroutine stress(atype, pos, v, f)
use atoms; use parameters
! calculate stress tensor components of kinetic energy part. 
!--------------------------------------------------------------------------------------------------------------
implicit none

real(8) :: atype(NBUFFER),pos(NBUFFER,3),v(NBUFFER,3),f(NBUFFER,3)

integer :: i,ity, j,jty,icmp, m,n
integer :: c1,c2,c3,c4,c5,c6
real(8) :: dr(3),dr2
real(8) :: mh, cstr(0:6,10)
real(8) :: dbuf(6), Gdbuf(6)
!--- cutoff length of stress calculation range
real(8) :: rcstr2 = 5.d0**2
real(8) :: qdummy(1)

!--- update buffer atom's stress value
dr=4*lcsize(1:3)
call COPYATOMS(atype, pos, v, f, qdummy, MODE_STRESSCALC, dr)

!--- get the potential contribution of internal pressure 
do i=1, NATOMS
  ity = atype(i)
  mh = mass(ity)
  xx = mh*v(i,1)*v(i,1)
  yy = mh*v(i,2)*v(i,2)
  zz = mh*v(i,3)*v(i,3)
  yz = mh*v(i,2)*v(i,3)
  zx = mh*v(i,1)*v(i,3)
  xy = mh*v(i,1)*v(i,2)

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
