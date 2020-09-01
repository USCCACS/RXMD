!-------------------------------------------------------------------------------------------
module utils
!-------------------------------------------------------------------------------------------

  integer,parameter :: MAXSTRLENGTH=256

  real(8),parameter :: pi=3.14159265358979d0
  real(8),parameter :: sqrtpi_inv=1.d0/sqrt(pi)

!--- Unit convetors. In original ReaxFF program, the units of length is [A]
!--- energy is [kcal/mol] and mass is [amu] respectively.
!--- Most of numbers written here are obtained below site.
!--- http://physics.nist.gov/cuu/Constants/Table/allascii.txt
!--- Bohr length
  real(8),parameter :: Lbohr_a  = 0.5291772108d0   ! [A]
  real(8),parameter :: Lbohr_m  = 0.5291772108d-10 ! [m]

!--- Electron rest mass
  real(8),parameter :: Merest_amu = 5.48580152d-4  ! [amu]
  real(8),parameter :: Merest_kg  = 9.10938d-31    ! [kg]
!--- Energy in Hartree
  real(8),parameter :: Ehrtre_km = 6.2751d2        ! [kcal/mol]
  real(8),parameter :: Ehrtre_ev  = 27.2113845d0   ! [eV]
  real(8),parameter :: Ehrtre_j = 4.35974417d-18   ! [J] 
  real(8),parameter :: Eev_kcal = 23.060538d0      ! [kcal/mol]

  real(8),parameter :: Ekcal_j = 6.95016611d-21  ! [J]

!--- Boltzmann Constant
  real(8),parameter :: BLTZMN = 1.3806503d-23  ! [m^2 kg s^-2 K^-1 ] 

  real(8),parameter :: UTEMP0 = 503.398008d0    ! Ekcal_j/BLZMN [K]
  real(8),parameter :: UTEMP = UTEMP0*2.d0/3.d0 ! [K]
  real(8),parameter :: USTRS = 6.94728103d0     ! [GPa]
  real(8),parameter :: UDENS = 1.66053886d0     ! [g/cc]
  real(8),parameter :: UTIME = 1.d3/20.455d0    ! 1 = 1/20.445[ps] = 48.88780[fs]

  type string_array
    character(len=:),allocatable :: str
  end type

contains

!-------------------------------------------------------------------------------------------
subroutine matinv(m1,m2)
! get inverse of m1 and save to m2
!-------------------------------------------------------------------------------------------
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

!-------------------------------------------------------------------------------------------
function str_gen(c) result(str)
!-------------------------------------------------------------------------------------------
implicit none
character(*) :: c
character(len=:),allocatable :: str

str = trim(adjustl(c))

return
end function

!-------------------------------------------------------------------------------------------
function int_to_str(i) result(str)
!-------------------------------------------------------------------------------------------
  integer,intent(in) :: i
  character(len=:),allocatable :: str
  character(len=32) :: str_buf

  write(str_buf,'(i32)') i
  str = trim(adjustl(str_buf))

  return
end function

!-------------------------------------------------------------------------------------------
integer function l2g(atype) result(global_id)
!convert Local ID to Global ID 
!-------------------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: atype
integer :: ity

ity = nint(atype)
global_id = nint((atype-ity)*1d13)

return
end function

!-------------------------------------------------------------------------------------------
function rank_to_string(irank) result(str)
!-------------------------------------------------------------------------------------------
implicit none
integer,intent(in) :: irank
character(len=:),allocatable :: str

write(str,*) irank
str = trim(adjustl(str))

end function 

!-------------------------------------------------------------------------------------------
function GetFileNameBase(DataDir,nstep) result(fileNameBase)
!-------------------------------------------------------------------------------------------
integer,intent(in) :: nstep
character(len=:),allocatable :: DataDir,fileNameBase
character(9) :: a9

if(nstep>=0) then
  write(a9,'(i9.9)') nstep
  fileNameBase=trim(adjustl(DataDir))//"/"//a9
else
  fileNameBase=trim(adjustl(DataDir))//"/rxff"
endif

end function


!-------------------------------------------------------------------------------------------
subroutine xu2xs(hhinv,box_origin,nmax,rreal,rnorm)
! update normalized coordinate from real coordinate. Subtract obox to make them local. 
!-------------------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: hhinv(3,3), box_origin(3)
integer,intent(in) :: nmax
real(8),intent(in) :: rreal(:,:)
real(8),intent(in out) :: rnorm(:,:)

real(8) :: rr(3)
integer :: i

do i=1,nmax
   rr(1:3) = rreal(i,1:3)
   rnorm(i,1)=sum(hhinv(1,1:3)*rr(1:3))
   rnorm(i,2)=sum(hhinv(2,1:3)*rr(1:3))
   rnorm(i,3)=sum(hhinv(3,1:3)*rr(1:3))
   rnorm(i,1:3) = rnorm(i,1:3) - box_origin(1:3)
enddo

end subroutine

!-------------------------------------------------------------------------------------------
subroutine xu2xs_inplace(hhinv,box_origin,nmax,rreal)
! update normalized coordinate from real coordinate. Subtract obox to make them local. 
!-------------------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: hhinv(3,3), box_origin(3)
integer,intent(in) :: nmax
real(8),intent(in out) :: rreal(:,:)
real(8) :: rr(3)
integer :: i

do i=1,nmax
   rr(1:3) = rreal(i,1:3)
   rreal(i,1)=sum(hhinv(1,1:3)*rr(1:3))
   rreal(i,2)=sum(hhinv(2,1:3)*rr(1:3))
   rreal(i,3)=sum(hhinv(3,1:3)*rr(1:3))
   rreal(i,1:3) = rreal(i,1:3) - box_origin(1:3)
enddo

end subroutine

!-------------------------------------------------------------------------------------------
subroutine xs2xu(hh,box_origin,nmax,rnorm,rreal)
! update real coordinate from normalized coordinate
!-------------------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: hh(3,3,0:1),box_origin(3)
integer,intent(in) :: nmax
real(8),intent(in) :: rnorm(:,:)
real(8),intent(in out) :: rreal(:,:)

real(8) :: rr(3)
integer :: i

do i=1,nmax 
   rr(1:3) = rnorm(i,1:3) + box_origin(1:3)
   rreal(i,1)=sum(hh(1,1:3,0)*rr(1:3))
   rreal(i,2)=sum(hh(2,1:3,0)*rr(1:3))
   rreal(i,3)=sum(hh(3,1:3,0)*rr(1:3))
enddo

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine xs2xu_inplace(hh,box_origin,nmax,rnorm)
! update real coordinate from normalized coordinate
!--------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: hh(3,3,0:1),box_origin(3)
integer,intent(in) :: nmax
real(8),intent(in out) :: rnorm(:,:)

real(8) :: rr(3)
integer :: i

do i=1,nmax 
   rr(1:3) = rnorm(i,1:3) + box_origin(1:3)
   rnorm(i,1)=sum(hh(1,1:3,0)*rr(1:3))
   rnorm(i,2)=sum(hh(2,1:3,0)*rr(1:3))
   rnorm(i,3)=sum(hh(3,1:3,0)*rr(1:3))
enddo

end subroutine

!--------------------------------------------------------------------------------------------------------------
subroutine get_boxparameters(H,la,lb,lc,angle1,angle2,angle3)
!--------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(in out) :: H(3,3)
real(8),intent(in) :: la,lb,lc, angle1,angle2,angle3
real(8) :: hh1, hh2 , lal, lbe, lga
real(8) :: pi=atan(1.d0)*4.d0

!--- convet unit for angles
lal=angle1*pi/180.d0
lbe=angle2*pi/180.d0
lga=angle3*pi/180.d0

!--- construct H-matrix
hh1=lc*(cos(lal)-cos(lbe)*cos(lga))/sin(lga)
hh2=lc*sqrt( 1.d0-cos(lal)**2-cos(lbe)**2-cos(lga)**2 + &
             2*cos(lal)*cos(lbe)*cos(lga) )/sin(lga)

H(1,1)=la;          H(2,1)=0.d0;        H(3,1)=0.d0
H(1,2)=lb*cos(lga); H(2,2)=lb*sin(lga); H(3,2)=0.d0
H(1,3)=lc*cos(lbe); H(2,3)=hh1;         H(3,3)=hh2

return
end subroutine

!----------------------------------------------------------------
subroutine update_box_params(vprocs, vid, hh, lata, latb, latc, maxrc, &
                             cc, lcsize, hhi, mdbox, lbox, obox)
!----------------------------------------------------------------
implicit none
real(8),intent(in) :: hh(3,3,0:1), lata, latb, latc, maxrc
integer,intent(in) :: vprocs(3), vid(3)
real(8),intent(in out) :: hhi(3,3), mdbox, lbox(3), obox(3), lcsize(3)
integer,intent(in out) :: cc(3)

!--- get volume 
MDBOX = &
HH(1,1,0)*(HH(2,2,0)*HH(3,3,0) - HH(3,2,0)*HH(2,3,0)) + &
HH(2,1,0)*(HH(3,2,0)*HH(1,3,0) - HH(1,2,0)*HH(3,3,0)) + &
HH(3,1,0)*(HH(1,2,0)*HH(2,3,0) - HH(2,2,0)*HH(1,3,0))

!--- get inverse of H-matrix
call matinv(HH,HHi)

!--- local box dimensions (a temporary use of lbox)
LBOX(1)=lata/vprocs(1)
LBOX(2)=latb/vprocs(2)
LBOX(3)=latc/vprocs(3)

!--- get the number of linkedlist cell per domain
cc(1:3)=int(LBOX(1:3)/maxrc)

!--- local system size in the unscaled coordinate.
LBOX(1:3) = 1.d0/vprocs(1:3)

!--- get the linkedlist cell dimensions (normalized)
lcsize(1:3) = LBOX(1:3)/cc(1:3)

!--- get origin of local MD box in the scaled coordiate.
OBOX(1:3) = LBOX(1:3)*vID(1:3)

return
end

!-------------------------------------------------------------------------------------------
integer function getstr(linein,lineout)
implicit none
!-------------------------------------------------------------------------------------------

character(len=:),allocatable,intent(in out) :: linein,lineout
integer :: pos1

!--- remove whitespace 
linein = adjustl(linein)

!--- return if black line
if(len(linein)==0) then
  getstr=-2
  return
endif

!--- return if it's a comment line or entirely whitespace
if(linein(1:1)=='#' .or. linein == repeat(' ', len(linein)) ) then
   getstr=-1
   return
endif

! find position in linein to get a token. if whitespace is not found, take entire line
pos1=index(linein,' ')-1
if(pos1==-1) pos1=len(linein)

lineout=linein(1:pos1)
linein=linein(pos1+1:)
getstr=len(lineout)

return
end

!-------------------------------------------------------------------------------------------
function get_command_argument_str(idx) result(string)
implicit none
!-------------------------------------------------------------------------------------------
integer,intent(in) :: idx
character(len=:),allocatable :: string
character(MAXSTRLENGTH) :: argv

call get_command_argument(idx,argv)
string = trim(adjustl(argv))

return
end function

!-------------------------------------------------------------------------------------------
logical function find_cmdline_argc(key,idx)
implicit none
!-------------------------------------------------------------------------------------------
integer,intent(inout) :: idx
character(*) :: key

integer :: i
character(MAXSTRLENGTH) :: argv

do i=1, command_argument_count()
   call get_command_argument(i,argv)
   if(index(argv,trim(adjustl(key))//' ')/=0) then ! trailing zero to distinguish '-foo ' and '-fooo'
      idx=i
      find_cmdline_argc=.true.
      return
   endif
enddo

idx=-1
find_cmdline_argc=.false.

return
end function

!-------------------------------------------------------------------------------------------
subroutine assert(condition, message, myrank, val)
!-------------------------------------------------------------------------------------------
logical,intent(in) :: condition
character(*),intent(in) :: message
integer,optional,intent(in) :: myrank
real(8),optional,intent(in) :: val

if(condition) return

print'(a)',repeat('!',60)

if(present(myrank)) write(6,fmt='(a)', advance='no') 'myrank '//int_to_str(myrank)//' :'
if(present(val)) write(6,fmt='(a,e15.5,a)', advance='no') ' value',val,' : '
print'(a)',message
print'(a)',repeat('!',60)
stop

return
end subroutine

end module
