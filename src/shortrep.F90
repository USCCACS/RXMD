module mod_short_repulsion
use base
use utils

implicit none

type short_repulsion_type
  logical :: has_short_repulsion=.false.

  !real(8),parameter :: alpha=0.1d0, GeGe_x=10d0, SeSe_x=1.11d0
  !real(8),parameter :: A_GeSe = 3.5d0-20, Ge_sig=0.73d0, Se_sig=2.d0 ! [J] & [A]
  !integer,parameter :: GeGe_eta = 11, GeSe_eta = 9, SeGe_eta=9, SeSe_eta = 7

  integer :: num_types
  real(8) :: alpha
  real(8),allocatable :: sigma(:)
  integer,allocatable :: eta(:,:)
  real(8),allocatable :: multiplier(:,:)
  real(8),allocatable :: Aij(:,:)
  real(8),allocatable :: coef(:,:)
end type 

type(short_repulsion_type) short_rep

contains

!-----------------------------------------------------------------------------
subroutine initialize_short_repulsion(sr)
!-----------------------------------------------------------------------------
type(short_repulsion_type),intent(in out) :: sr

integer :: idx, ity, jty, funit
real(8) :: sigma_sum

character(len=:),allocatable :: filename
logical :: fexists


! num_types = 2
! multiplier = [[GeGe_x,0.d0], [0.d0, SeSe_x]]
! Aij = A_GeSe/Ekcal_j
! sigma = [Ge_sig,Se_sig]
! eta = [[GeGe_eta, GeSe_eta],[SeGe_eta, SeSe_eta]]

if(find_cmdline_argc('--short_rep',idx)) then
   filename = 'shortrep.in'
   inquire(file=filename, exist=sr%has_short_repulsion)

   if(.not. sr%has_short_repulsion) then

     filename = get_command_argument_str(idx+1)
     inquire(file=filename, exist=sr%has_short_repulsion)

     if(.not. sr%has_short_repulsion) then
       if(myid==0) then
          print'(a)', 'ERROR: short_rep flag is given but could not find shortrep.in file. Existing rxmd.'
          call MPI_FINALIZE(ierr)
          stop
       endif
     endif

   endif
endif

! the short_repulsion is not applied if input file is not found by here.
if(.not. sr%has_short_repulsion) return

open(newunit=funit, file=filename, form='formatted')
read(funit,*) sr%num_types

allocate(sr%Aij(sr%num_types,sr%num_types), sr%multiplier(sr%num_types,sr%num_types), sr%coef(sr%num_types,sr%num_types))
allocate(sr%eta(sr%num_types,sr%num_types), sr%sigma(sr%num_types))

read(funit,*) sr%alpha
read(funit,*) (sr%sigma(ity), ity=1, sr%num_types)
do ity = 1, sr%num_types
   read(funit,*) (sr%eta(ity,jty), jty=1, sr%num_types)
enddo
do ity = 1, sr%num_types
   read(funit,*) (sr%Aij(ity,jty), jty=1, sr%num_types)
enddo
do ity = 1, sr%num_types
   read(funit,*) (sr%multiplier(ity,jty), jty=1, sr%num_types)
enddo

do ity=1,sr%num_types
do jty=1,sr%num_types
   sigma_sum = sr%sigma(ity)+sr%sigma(jty)
   sr%coef(ity,jty) = sr%alpha*(sr%Aij(ity,jty)/Ekcal_j)*sigma_sum**sr%eta(ity,jty)

   sr%coef(ity,jty) = sr%multiplier(ity,jty)*sr%coef(ity,jty)
enddo; enddo

call print(sr)

end subroutine

!-----------------------------------------------------------------------------
subroutine print(sr)
!-----------------------------------------------------------------------------
implicit none
type(short_repulsion_type),intent(in) :: sr
integer :: ity,jty

if(myid==0) then
   print'(a)', repeat('-',60)
   print'(a)', 'parameters for short repulsion: '
   print'(a)', repeat('-',60)
   print'(a,i9)', 'num_types : ', sr%num_types

   print*

   print'(9(a4,i1,a1,1x,f10.3,3x))', ('sig(',ity,')',sr%sigma(ity), ity=1, sr%num_types)
   do ity=1, sr%num_types
      print'(9(a4,i1,a1,i1,a1,1x,i6,3x))', ('eta(',ity,',',jty,')', sr%eta(ity,jty), jty=1, sr%num_types)
   enddo
   do ity=1, sr%num_types
      print'(9(a4,i1,a1,i1,a1,1x,es13.2,3x))', ('Aij(',ity,',',jty,')', sr%Aij(ity,jty), jty=1, sr%num_types)
   enddo

   print*

   print'(a,f10.3)', 'alpha : ', sr%alpha
   do ity=1, sr%num_types
      print'(9(a4,i1,a1,i1,a1,1x,f8.3,2x))', ('mul(',ity,',',jty,')', sr%multiplier(ity,jty), jty=1, sr%num_types)
   enddo

   print*

   do ity=1, sr%num_types
      print'(9(a5,i1,a1,i1,a1,1x,f10.3,3x))', ('coef(',ity,',',jty,')',sr%coef(ity,jty), jty=1, sr%num_types)
   enddo
   print'(a)', repeat('-',60)
endif

end subroutine

!-----------------------------------------------------------------------------
subroutine short_repulsion(sr)
!-----------------------------------------------------------------------------
type(short_repulsion_type),intent(in) :: sr

real(8) :: coef(2,2), fcoef, ff0(3), dr(3), dr1, dr2, dri, dri_eta
integer :: i, j, j1, ity, jty

if(.not. sr%has_short_repulsion) return

do i=1, NATOMS

   ity = nint(atype(i))

   do j1 = 1, nbrlist(i,0)

      j = nbrlist(i,j1)
      jty = nint(atype(j))

      dr(1:3) = pos(i,1:3) - pos(j,1:3)
      dr2 = sum(dr(1:3)*dr(1:3))
      dr1 = sqrt(dr2)
      dri = 1d0/dr1
      dri_eta = dri**(sr%eta(ity,jty)+1)
      fcoef = sr%eta(ity,jty)*sr%coef(ity,jty)*dri_eta*dri
      !if(l2g(atype(i))==1.and.dr1<3.d0)
      !print'(2i6,2f12.5)',l2g(atype(i)),l2g(atype(j)),dr1,fcoef

      f(i,1:3) = f(i,1:3) + fcoef*dr(1:3)

   enddo
enddo

end subroutine

end module
