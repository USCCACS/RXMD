module params
implicit none
integer,parameter :: vprocs(3)=(/8,8,8/)
integer,parameter :: mc(3)=(/1,1,1/)

integer,parameter :: nprocs=vprocs(1)*vprocs(2)*vprocs(3)
integer,parameter :: mctot=mc(1)*mc(2)*mc(3)
integer :: lnatoms(0:nprocs-1)  ! local # of atoms
integer :: lnatoms1(0:nprocs-1) ! prefix sum of above
integer :: lnatoms2(0:nprocs-1) ! temp counter
real(8) :: L1,L2,L3,Lalpha,Lbeta,Lgamma
real(8) :: lbox(3),obox(3), dtype
real(8) :: H(3,3), Hi(3,3), rr(3),rr1(3),vv(3),qq,rmin(3), rmax(3)
real(8) :: qfsp0=0.d0, qfsv0=0.d0
integer :: natoms,  ntot, ix,iy,iz
real(8),allocatable :: pos0(:,:), pos1(:,:)
character(3),allocatable :: ctype0(:), ctype1(:)
integer,allocatable :: itype0(:)
real(8),allocatable :: itype1(:)
character(256) :: fname 
character(256) :: fnote

contains

!----------------------------------------------------------------
subroutine getbox(la,lb,lc,angle1,angle2,angle3)
!----------------------------------------------------------------
implicit none
real(8),intent(in) :: la,lb,lc, angle1,angle2,angle3
real(8) :: hh1, hh2 , lal, lbe, lga
real(8) :: pi=datan(1.d0)*4.d0

!--- convet unit for angles
lal=angle1*pi/180.d0
lbe=angle2*pi/180.d0
lga=angle3*pi/180.d0

!--- construct H-matrix
hh1=lc*(cos(lal)-cos(lbe)*cos(lga))/sin(lga)
hh2=lc*sqrt( 1.d0-cos(lal)**2-cos(lbe)**2-cos(lga)**2 + &
             2*cos(lal)*cos(lbe)*cos(lga) )/sin(lga)

H(1,1)=la
H(2,1)=0.d0
H(3,1)=0.d0
H(1,2)=lb*cos(lga)
H(2,2)=lb*sin(lga)
H(3,2)=0.d0
H(1,3)=lc*cos(lbe)
H(2,3)=hh1
H(3,3)=hh2

call matinv(H,Hi)

return 
end subroutine


!----------------------------------------------------------------
subroutine matinv(m1,m2)
! get inverse of m1 and save to m2
!----------------------------------------------------------------
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

return
end subroutine

end module


!================================================================
program geninit
use params
implicit none
!================================================================
integer :: i, j, k, n, ia, myid, sID
integer(8) :: ii=0
character(6) :: a6

!--- get structure file name
call getarg(1,fname)

open(1,file=fname,form="formatted")
write(6,'(2a,1x,a,i9,3i6)') 'input file: ', trim(fname),' mctot,mc',mctot, mc(1:3)

!--- read # of atoms and file description
read(1,*) natoms, fnote
print'(i9,3x,a)',natoms, trim(fnote)

!--- read lattice parameters
read(1,*) L1, L2, L3, Lalpha, Lbeta, Lgamma
!print'(a,3x,6f12.3)','1, L2, L3, Lalpha, Lbeta, Lgamma: ', & 
!           L1, L2, L3, Lalpha, Lbeta, Lgamma

call getbox(L1,L2,L3,lalpha,lbeta,lgamma)

!--- allocate arrays for atom position and type
allocate(ctype0(natoms),pos0(3,natoms),itype0(natoms))
allocate(ctype1(natoms*mctot),pos1(3,natoms*mctot),itype1(natoms*mctot))
do i=1, natoms
   read(1,*) ctype0(i),pos0(1:3,i)
   ctype0(i)=adjustl(ctype0(i))
   select case(ctype0(i)) 
     case("C  ") 
       itype0(i)=1
     case("H  ") 
       itype0(i)=2
     case("O  ") 
       itype0(i)=3
     case("S  ") 
       itype0(i)=4
     case("Mo ") 
       itype0(i)=5
     case("Ni ") 
       itype0(i)=6
     case("Al ")
       itype0(i)=7   
   end select
enddo
close(1)

!--- repeat the unit cell
ntot=0
do ix=0, mc(1)-1
do iy=0, mc(2)-1
do iz=0, mc(3)-1
do i=1, natoms
   ntot=ntot+1
   !rr(1)=sum(Hi(1,1:3)*pos0(1:3,i))
   !rr(2)=sum(Hi(2,1:3)*pos0(1:3,i))
   !rr(3)=sum(Hi(3,1:3)*pos0(1:3,i))
   rr(1:3)=pos0(1:3,i) !!! THIS IS FOR NORMALIZED COORD XYZ !!!
   rr(1)=(rr(1)+ix)/mc(1)
   rr(2)=(rr(2)+iy)/mc(2)
   rr(3)=(rr(3)+iz)/mc(3)
   pos1(1:3,ntot)=rr(1:3)
   ctype1(ntot)=ctype0(i)
   itype1(ntot)=itype0(i)+ntot*1d-13+1d-14
   !print'(a3,1x,i3,3f8.3,3f)',ctype0(i),itype0(i),pos0(1:3,i), rr(1:3)
enddo 
enddo; enddo; enddo

!--- shift coordinates, then shift a bit to avoid zero coordinates
do i=1, 3
   rmin(i)=minval(pos1(i,:))
   pos1(i,:)=pos1(i,:)-rmin(i)
enddo

!--- wrap back atom coordinates if necessary
do i=1, ntot
   do j=1, 3
      pos1(j,i)=mod(pos1(j,i),1.d0)
   enddo
   pos1(1:3,i)=pos1(1:3,i)+1d-9 ! shift by a small value to avoid coordinates at 0
enddo

!--- check coordinates
do i=1, 3
   rmin(i)=minval(pos1(i,:))
   rmax(i)=maxval(pos1(i,:))
enddo
print'(a,3es15.5,3x,3es15.5)','rmin(1:3),rmax(1:3): ', rmin(1:3), rmax(1:3)

!--- update natoms & H-matrix 
natoms=natoms*mctot
L1=L1*mc(1); L2=L2*mc(2); L3=L3*mc(3)
call getbox(L1,L2,L3,lalpha,lbeta,lgamma)

!--- count how many atoms per MPI domain, get lnatoms()
lnatoms(:)=0
do n=1,natoms 
   i=pos1(1,n)*vprocs(1)
   j=pos1(2,n)*vprocs(2)
   k=pos1(3,n)*vprocs(3)
   sID = i + j*vprocs(1) + k*vprocs(1)*vprocs(2)
   lnatoms(sID)=lnatoms(sID)+1
enddo

!--- prefix sum
lnatoms1(0)=0
do i=1, nprocs-1
   lnatoms1(i)=lnatoms1(i-1)+lnatoms(i-1)
enddo

!--- To avoid opening too many files, sort atom data using disk
lbox(1:3)=1.d0/vprocs(1:3)
lnatoms2(:)=0 ! for atom counting & result check
open(1,file="all.bin",form="unformatted",access="stream")
do n=1,natoms 
   i=pos1(1,n)*vprocs(1)
   j=pos1(2,n)*vprocs(2)
   k=pos1(3,n)*vprocs(3)
   sID = i + j*vprocs(1) + k*vprocs(1)*vprocs(2)
   obox(1:3)=lbox(1:3)*(/i,j,k/)
   pos1(1:3,n)=pos1(1:3,n)-obox(1:3)

!--- total number of atoms before n-th atom, 
!--- each atom has 32 (=3*8 + 8) bytes data
   ii=lnatoms1(sID)+lnatoms2(sID) 
   write(1,pos=ii*32+1) pos1(1:3,n),itype1(n)

   lnatoms2(sID)=lnatoms2(sID)+1
enddo
close(1)

!--- error check
do i=0, nprocs-1
   if(lnatoms(i)/=lnatoms2(i)) then
      print'(a,3i9)', 'Error,i,lnatoms(i),lnatoms2(i)',i,lnatoms(i),lnatoms2(i)
      stop
   endif
enddo

qq=0.d0; vv(:)=0.d0; qfsp0=0.d0; qfsv0=0.d0
open(1,file="all.bin",form="unformatted",access="stream")
open(20,file="geninit.xyz")
open(30,file="rxff.bin",form="unformatted",access="stream")
write(30) nprocs, vprocs(1:3)
write(30) lnatoms(0:nprocs-1)
write(30) 0
write(30) L1, L2, L3, Lalpha, Lbeta, Lgamma

write(20,'(i12)') sum(lnatoms(:))
write(20,'(a)') trim(fname)
do myid=0,nprocs-1
   write(a6(1:6),'(i6.6)') myid

   i=mod(myid,vprocs(1))
   j=mod(myid/vprocs(1),vprocs(2))
   k=myid/(vprocs(1)*vprocs(2))
   obox(1:3)=lbox(1:3)*(/i,j,k/)

   do n=1,lnatoms(myid)
      read(1) rr(1:3),dtype

      write(30)rr(1:3)
      write(30)vv(1:3)
      write(30)qq
      write(30)dtype
      write(30)qfsp0
      write(30)qfsv0

      rr1(1:3)=rr(1:3)+obox(1:3)
      rr(1)=sum(H(1,1:3)*rr1(1:3))
      rr(2)=sum(H(2,1:3)*rr1(1:3))
      rr(3)=sum(H(3,1:3)*rr1(1:3))
      select case(nint(dtype)) 
        case(1) 
          write(20,'(a3, $)') "C  "
        case(2) 
          write(20,'(a3, $)') "H  "
        case(3) 
          write(20,'(a3, $)') "O  "
        case(4) 
          write(20,'(a3, $)') "N  "
        case(6) 
          write(20,'(a3, $)') "Si "
      end select
      write(20,'(3f12.5)') rr(1:3)
   enddo

enddo
close(1)
close(20)
close(30)

print*,'sum(lnatoms), lnatoms: ',sum(lnatoms),lnatoms(:)
print'(a,3x,6f12.3)','L1, L2, L3, Lalpha, Lbeta, Lgamma: ', & 
           L1, L2, L3, Lalpha, Lbeta, Lgamma

end

