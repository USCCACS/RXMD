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

contains

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
integer function l2g(atype) result(global_id)
implicit none
!convert Local ID to Global ID 
!--------------------------------------------------------------------------------------------------------------
real(8),intent(IN) :: atype
integer :: ity

ity = nint(atype)
global_id = nint((atype-ity)*1d13)

return
end function

!--------------------------------------------------------------------------------------------------------------
subroutine xu2xs(hhinv,box_origin,nmax,rreal,rnorm)
! update normalized coordinate from real coordinate. Subtract obox to make them local. 
!--------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: hhinv(3,3), box_origin(3)
integer,intent(in) :: nmax
real(8),intent(in) :: rreal(:,:)
real(8),intent(out) :: rnorm(:,:)

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

!--------------------------------------------------------------------------------------------------------------
subroutine xu2xs_inplace(hhinv,box_origin,nmax,rreal)
! update normalized coordinate from real coordinate. Subtract obox to make them local. 
!--------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: hhinv(3,3), box_origin(3)
integer,intent(in) :: nmax
real(8),intent(inout) :: rreal(:,:)
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


!--------------------------------------------------------------------------------------------------------------
subroutine xs2xu(hh,box_origin,nmax,rnorm,rreal)
! update real coordinate from normalized coordinate
!--------------------------------------------------------------------------------------------------------------
implicit none
real(8),intent(in) :: hh(3,3,0:1),box_origin(3)
integer,intent(in) :: nmax
real(8),intent(in) :: rnorm(:,:)
real(8),intent(out) :: rreal(:,:)

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

end module
