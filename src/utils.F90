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

  character(len=:),allocatable :: token
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
if(present(val)) write(6,fmt='(a,es15.5,a)', advance='no') ' value',val,' : '
print'(a)',message
print'(a)',repeat('!',60)
stop

return
end subroutine

!-------------------------------------------------------------------------------------------
subroutine put_rng_seed(seed,myrank)
!-------------------------------------------------------------------------------------------
implicit none
integer,intent(in) :: seed, myrank

integer,allocatable :: seeds0(:), seeds(:)
integer :: i,seed_size 

 
call random_seed(size=seed_size)
allocate(seeds(seed_size),seeds0(seed_size))
seeds0=0; seeds=0

call random_seed(get=seeds0)
seeds = seeds0
seeds(1) = seed

if(myrank==0) then
  print'(a)','INFO :resetting RNG seeds'
  do i=1, seed_size
     print'(a,i4,2i12)','seeds,seeds0 ', i, seeds0(i), seeds(i)
  enddo
endif

call random_seed(put=seeds)

deallocate(seeds)

end subroutine

end module

!-------------------------------------------------------------------------------------------
module element_name_manager
!-------------------------------------------------------------------------------------------

  type single_element_type
     character(len=:),allocatable :: name
     integer :: atomic
     logical :: set
  endtype

  type elements_type
     type(single_element_type),allocatable :: e(:)
     character(len=:),allocatable :: str 
  end type

  type(elements_type) :: elements0, elements

  integer,parameter :: BUF_LEN=5
  character(len=BUF_LEN) :: a5

contains

  subroutine set_initial_lookup_table()
     integer :: i

     if(.not.allocated(elements0%e)) allocate(elements0%e(0))

     elements0%e = [elements0%e, single_element_type(name="H", atomic=0,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="He", atomic=1,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Li", atomic=2,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Be", atomic=3,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="B", atomic=4,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="C", atomic=5,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="N", atomic=6,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="O", atomic=7,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="F", atomic=8,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Ne", atomic=9,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Na", atomic=10,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Mg", atomic=11,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Al", atomic=12,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Si", atomic=13,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="P", atomic=14,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="S", atomic=15,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Cl", atomic=16,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Ar", atomic=17,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="K", atomic=18,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Ca", atomic=19,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Sc", atomic=20,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Ti", atomic=21,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="V", atomic=22,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Cr", atomic=23,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Mn", atomic=24,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Fe", atomic=25,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Co", atomic=26,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Ni", atomic=27,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Cu", atomic=28,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Zn", atomic=29,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Ga", atomic=30,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Ge", atomic=31,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="As", atomic=32,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Se", atomic=33,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Br", atomic=34,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Kr", atomic=35,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Rb", atomic=36,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Sr", atomic=37,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Y", atomic=38,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Zr", atomic=39,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Nb", atomic=40,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Mo", atomic=41,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Tc", atomic=42,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Ru", atomic=43,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Rh", atomic=44,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Pd", atomic=45,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Ag", atomic=46,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Cd", atomic=47,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="In", atomic=48,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Sn", atomic=49,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Sb", atomic=50,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Te", atomic=51,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="I", atomic=52,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Xe", atomic=53,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Cs", atomic=54,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Ba", atomic=55,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="La", atomic=56,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Ce", atomic=57,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Pr", atomic=58,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Nd", atomic=59,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Pm", atomic=60,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Sm", atomic=61,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Eu", atomic=62,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Gd", atomic=63,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Tb", atomic=64,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Dy", atomic=65,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Ho", atomic=66,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Er", atomic=67,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Tm", atomic=68,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Yb", atomic=69,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Lu", atomic=70,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Hf", atomic=71,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Ta", atomic=72,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="W", atomic=73,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Re", atomic=74,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Os", atomic=75,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Ir", atomic=76,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Pt", atomic=77,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Au", atomic=78,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Hg", atomic=79,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Tl", atomic=80,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Pb", atomic=81,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Bi", atomic=82,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Ac", atomic=83,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Th", atomic=84,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Pa", atomic=85,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="U", atomic=86,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Np", atomic=87,set=.false.)]
     elements0%e = [elements0%e, single_element_type(name="Pu", atomic=88,set=.false.)]

     do i=1, size(elements0%e)
        write(a5,'(a5)') elements0%e(i)%name
        elements0%str=elements0%str//adjustl(a5)
     enddo

  end subroutine

  subroutine add_elem(elems, elem)
     character(len=*),intent(in) :: elem
     type(elements_type),intent(in out) :: elems
     integer :: idx

     if(.not.allocated(elems%e)) allocate(elems%e(0))

     idx=get_loc(elements0%str, elem) 
     elements0%e(idx)%set = .true. 

     elems%e = [elems%e, elements0%e(idx)]

     write(a5,'(a5)') elements0%e(idx)%name
     elems%str= elems%str//adjustl(a5)
  end subroutine

  function get_loc(element_buffer, elem) result(loc)
     character(len=:),allocatable :: element_buffer
     character(len=*) :: elem
     integer :: loc
     loc = index(element_buffer, elem)
     if (loc > 0) then
        loc = (loc-1)/BUF_LEN + 1
     else
        print'(2a)','Undefined element found. Please add it to set_initial_lookup_table(): ', elem
        stop 
     endif
  end function

end module 
