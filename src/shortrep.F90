module mod_short_repulsion
use base
use utils

implicit none

type short_repulsion_type1_params

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

  logical :: is_initialized = .false.

  contains
    procedure :: read => read_type1_params
    procedure :: print => print_type1_params
    procedure :: apply => apply_type1_short_repulsion
end type

type short_repulsion_type2_params

  real(8) :: alpha
  real(8) :: rc_oh

  real(8) :: Kr, Kq
  real(8) :: r0, q0

  integer :: htype, otype

  character(2),allocatable :: atom_name(:)

  logical :: is_initialized = .false.

  contains 
    procedure :: read => read_type2_params
    procedure :: print => print_type2_params
    procedure :: apply => apply_type2_short_repulsion
end type

type short_repulsion_type
  logical :: has_short_repulsion=.false.

  integer :: potential_type = 2

  type(short_repulsion_type1_params) :: p1
  type(short_repulsion_type2_params) :: p2

end type 

type(short_repulsion_type) short_rep

contains


!-----------------------------------------------------------------------------
subroutine initialize_short_repulsion(sr, atom_name)
!-----------------------------------------------------------------------------
type(short_repulsion_type),intent(in out) :: sr
character(2),allocatable,intent(in) :: atom_name(:)

integer :: idx, ity, jty, funit
real(8) :: sigma_sum

character(len=:),allocatable :: filename
logical :: fexists

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

read(funit,*) sr%potential_type

if(sr%potential_type==1) then
   call sr%p1%read(funit)
   call sr%p1%print()
else if(sr%potential_type==2) then
   call sr%p2%read(funit, atom_name)
   call sr%p2%print()
else
   print'(a)', 'ERROR: unsupported potential type. Exiting the code', sr%potential_type
   stop
endif

close(funit)

end subroutine

!-----------------------------------------------------------------------------
subroutine read_type2_params(this, funit, atom_name)
!-----------------------------------------------------------------------------
  implicit none
  class(short_repulsion_type2_params) :: this
  character(2),allocatable,intent(in) :: atom_name(:)

  integer,intent(in) :: funit 
  integer :: ity, jty
  real(8) :: sigma_sum

  do ity=1, size(atom_name)
     if(index(atom_name(ity),"H") /=0) this%htype=ity
     if(index(atom_name(ity),"O") /=0) this%otype=ity
  enddo

  read(funit,*) this%alpha, this%rc_oh
  read(funit,*) this%Kr, this%Kq
  read(funit,*) this%r0, this%q0

  this%is_initialized = .true.
end subroutine

!-----------------------------------------------------------------------------
subroutine print_type2_params(this)
!-----------------------------------------------------------------------------
  implicit none
  class(short_repulsion_type2_params) :: this

  if(myid==0) then
    print'(a)', repeat('-',60)
    print'(a30,a7,i3,1x,a7,i3)', "atom types: ", "H -", this%htype, "O -", this%otype
    print'(a30,2f10.3)', 'alpha, rc_oh: ', this%alpha, this%rc_oh
    print'(a30,2f10.3)', 'Kr, Kq: ', this%Kr, this%Kq
    print'(a30,2f10.3)', 'r0, q0: ', this%r0, this%q0
    print'(a)', repeat('-',60)
  endif

end subroutine


!-----------------------------------------------------------------------------
subroutine read_type1_params(this, funit)
!-----------------------------------------------------------------------------
  implicit none
  class(short_repulsion_type1_params) :: this
  integer,intent(in) :: funit 
  integer :: ity, jty
  real(8) :: sigma_sum

  read(funit,*) this%num_types
  
  allocate(this%Aij(this%num_types,this%num_types), &
           this%multiplier(this%num_types,this%num_types), &
           this%coef(this%num_types,this%num_types))
  allocate(this%eta(this%num_types,this%num_types), this%sigma(this%num_types))
  
  read(funit,*) this%alpha
  read(funit,*) (this%sigma(ity), ity=1, this%num_types)
  do ity = 1, this%num_types
     read(funit,*) (this%eta(ity,jty), jty=1, this%num_types)
  enddo
  do ity = 1, this%num_types
     read(funit,*) (this%Aij(ity,jty), jty=1, this%num_types)
  enddo
  do ity = 1, this%num_types
     read(funit,*) (this%multiplier(ity,jty), jty=1, this%num_types)
  enddo
  
  do ity=1,this%num_types
  do jty=1,this%num_types
     sigma_sum = this%sigma(ity)+this%sigma(jty)
     this%coef(ity,jty) = this%alpha*(this%Aij(ity,jty)/Ekcal_j)*sigma_sum**this%eta(ity,jty)
  
     this%coef(ity,jty) = this%multiplier(ity,jty)*this%coef(ity,jty)
  enddo; enddo

  this%is_initialized = .true.
end subroutine


!-----------------------------------------------------------------------------
subroutine print_type1_params(this)
!-----------------------------------------------------------------------------
implicit none
class (short_repulsion_type1_params),intent(in) :: this
integer :: ity,jty

if(myid==0) then
   print'(a)', repeat('-',60)
   print'(a)', 'parameters for short repulsion: '
   print'(a)', repeat('-',60)
   print'(a,i9)', 'num_types : ', this%num_types

   print*

   print'(9(a4,i1,a1,1x,f10.3,3x))', ('sig(',ity,')',this%sigma(ity), ity=1, this%num_types)
   do ity=1, this%num_types
      print'(9(a4,i1,a1,i1,a1,1x,i6,3x))', ('eta(',ity,',',jty,')', this%eta(ity,jty), jty=1, this%num_types)
   enddo
   do ity=1, this%num_types
      print'(9(a4,i1,a1,i1,a1,1x,es13.2,3x))', ('Aij(',ity,',',jty,')', this%Aij(ity,jty), jty=1, this%num_types)
   enddo

   print*

   print'(a,f10.3)', 'alpha : ', this%alpha
   do ity=1, this%num_types
      print'(9(a4,i1,a1,i1,a1,1x,f8.3,2x))', ('mul(',ity,',',jty,')', this%multiplier(ity,jty), jty=1, this%num_types)
   enddo

   print*

   do ity=1, this%num_types
      print'(9(a5,i1,a1,i1,a1,1x,f10.3,3x))', ('coef(',ity,',',jty,')',this%coef(ity,jty), jty=1, this%num_types)
   enddo
   print'(a)', repeat('-',60)
endif

end subroutine

!-----------------------------------------------------------------------------
subroutine apply_type2_short_repulsion(this)
!-----------------------------------------------------------------------------
class (short_repulsion_type2_params),intent(in) :: this

real(8) :: coef(2,2), fcoef, ff0(3), dr(3), dr1, dr2, dri, dri_eta
integer :: i, j, k, j1, k1, ity, jty

integer :: nbr_count, nbr_index(10)
real(8) :: nbr_dist(10), dist_1st, dist_2nd, rij(0:3), rik(0:3)
real(8) :: sine, cosine, theta, Krcoef, Kqcoef, ff(3)

real(8),parameter :: MAX_DIST = 1e9, pi=3.14159265358979d0

do i=1, NATOMS

   ity = nint(atype(i))

   if(ity/=this%otype) cycle

   !if(myid==0) print'(a,2i6)','ity, this%otype:',ity,this%otype

   nbr_count = 0; nbr_dist = MAX_DIST
   do j1 = 1, nbrlist(i,0)

      j = nbrlist(i,j1)
      jty = nint(atype(j))

      if(jty/=this%htype) cycle

      dr(1:3) = pos(i,1:3) - pos(j,1:3)
      dr2 = sum(dr(1:3)*dr(1:3))
      dr1 = sqrt(dr2)
      dri = 1d0/dr1

      if(dr1>this%rc_oh) cycle

      nbr_count = nbr_count + 1
      nbr_index(nbr_count) = j
      nbr_dist(nbr_count) = dr1

      !if(myid==0) print'(a,i6,2i6,f8.3)','jty, this%htype:',i,jty,this%htype, dr1

      !f(i,1:3) = f(i,1:3) + fcoef*dr(1:3)
   enddo
   !print*,nbr_index(1:nbr_count)

   j1 = minloc(nbr_dist,dim=1)
   dist_1st = nbr_dist(j1)
   nbr_dist(j1) = MAX_DIST
   j = nbr_index(j1) 

   k1 = minloc(nbr_dist,dim=1)
   dist_2nd = nbr_dist(k1)
   nbr_dist(k1) = MAX_DIST
   k = nbr_index(k1) 

   if(nbr_count > 2) then
      print'(a)',repeat('-=',30)
      print'(a,2i6)', 'ERROR: more than two neighbor hydrogen found : ', myid, nbr_count
      print'(a)',repeat('-=',30)
   endif

   rij(1:3) = pos(i,1:3)-pos(j,1:3)
   rij(0) = sqrt(sum(rij(1:3)*rij(1:3)))
   rij(1:3) = rij(1:3)/rij(0)

   rik(1:3) = pos(i,1:3)-pos(k,1:3)
   rik(0) = sqrt(sum(rik(1:3)*rik(1:3)))
   rik(1:3) = rik(1:3)/rik(0)

   cosine = sum(rij(1:3)*rik(1:3))
   theta = acos(cosine)
   sine = sin(theta)
  

   Krcoef = this%alpha * this%Kr*(rij(0)-this%r0) ! ij
   ff(1:3) = Krcoef*rij(1:3)
   f(i,1:3) = f(i,1:3) - ff(1:3)
   f(j,1:3) = f(j,1:3) + ff(1:3)
   !write(unit=6,fmt='(a,4f12.3)') 'Krcoef, ff(1:3): ', Krcoef, ff(1:3)

   Krcoef = this%alpha * this%Kr*(rik(0)-this%r0) ! ik
   ff(1:3) = Krcoef*rik(1:3)
   f(i,1:3) = f(i,1:3) - ff(1:3)
   f(k,1:3) = f(k,1:3) + ff(1:3)
   !write(unit=6,fmt='(a,4f12.3)') 'Krcoef, ff(1:3): ', Krcoef, ff(1:3)

   Kqcoef = this%alpha * this%Kq*(theta*180d0/pi-this%q0)*(-1.d0/sine)

   rij(1:3)=-rij(1:3)
   !print*,Kqcoef, theta*180d0/pi, this%q0
   call ForceA3(Kqcoef,j,i,k,rij,rik)

   !if(myid==0) print'(a,3i6,3f8.3)','i,j,k,dist_1st,dist_2nd, theta: ', &
   !    i, j, k, dist_1st, dist_2nd, theta


enddo

contains

!-----------------------------------------------------------------------
subroutine cross_product(dr1, dr2, crs)
! Calculate a cross product <dr1(1:3)> x <dr2(1:3)> = <crs(1:3)>
! <dr1> and <dr2> must have thier norm in 0th element.
!-----------------------------------------------------------------------
   implicit none
   real(8),INTENT(IN) :: dr1(0:3), dr2(0:3)
   real(8),INTENT(OUT) :: crs(0:3)
   real(8) :: ndr1(1:3), ndr2(1:3)
   real(8),parameter :: NSMALL = 1d-6

   ndr1(1:3) = dr1(1:3)/dr1(0)
   ndr2(1:3) = dr2(1:3)/dr2(0)

   crs(1) = ndr1(2)*ndr2(3) - ndr1(3)*ndr2(2)
   crs(2) = ndr1(3)*ndr2(1) - ndr1(1)*ndr2(3)
   crs(3) = ndr1(1)*ndr2(2) - ndr1(2)*ndr2(1)
   crs(0) = sqrt( sum(crs(1:3)*crs(1:3)) )
   if(crs(0)<NSMALL) crs(0) = NSMALL

end subroutine

!-----------------------------------------------------------------------
subroutine ForceA3(coeff,i,j,k,da0, da1)
! derivative of <cos_ijk>
!-----------------------------------------------------------------------
implicit none
! Caa(0,0)=Caa, Caa(0,-1)=Caa-1, Caa(-1,0)=Ca-1a, Caa(-1,-1)=Ca-1a-1

real(8),INTENT(IN) :: coeff, da0(0:3), da1(0:3)
integer,INTENT(IN) :: i,j,k
real(8) :: Caa(-2:0,-2:0), Ci(3), Ck(3)
real(8) :: fij(3), fjk(3), fijjk(3), rij(3), rjk(3)
real(8) :: CCisqr, coCC

Caa( 0,0) = da0(0)**2; Caa( 0,-1) = sum(da0(1:3)*da1(1:3))
Caa(-1,0) = Caa(0,-1); Caa(-1,-1) = da1(0)**2

CCisqr = 1.d0/( da0(0)*da1(0) )
coCC = coeff*CCisqr

!--- Some of calculations are unnecessary due to the action-reaction relation.
Ci(1) = -( Caa(0,-1)/Caa(0,0) )
Ci(2) =  1.d0
!Cj(1) =  Caa( 0,-1)/Caa( 0, 0) + 1.d0 
!Cj(2) = -( Caa(0,-1)/Caa(-1,-1) + 1.d0 ) 
Ck(1) = -1.d0
Ck(2) =  Caa(0,-1)/Caa(-1,-1)

rij(1:3) = da0(1:3)
rjk(1:3) = da1(1:3)

fij(1:3) = coCC*(Ci(1)*rij(1:3) + Ci(2)*rjk(1:3))
fjk(1:3) =-coCC*(Ck(1)*rij(1:3) + Ck(2)*rjk(1:3))
fijjk(1:3) =  -fij(1:3) + fjk(1:3) 

!$omp atomic
f(i,1) = f(i,1) + fij(1)
!$omp atomic
f(i,2) = f(i,2) + fij(2)
!$omp atomic
f(i,3) = f(i,3) + fij(3)

!$omp atomic
f(j,1) = f(j,1) + fijjk(1)
!$omp atomic
f(j,2) = f(j,2) + fijjk(2)
!$omp atomic
f(j,3) = f(j,3) + fijjk(3)

!$omp atomic
f(k,1) = f(k,1) - fjk(1)
!$omp atomic
f(k,2) = f(k,2) - fjk(2)
!$omp atomic
f(k,3) = f(k,3) - fjk(3)

!--- Check N3rd ---
!print'(a,6f20.13)','N3rd: ', &
!Ci(1)-Ci(2), -Cj(1), Ci(2), -Ck(1), Cj(2), Ck(1)-Ck(2)

end subroutine

end subroutine

!-----------------------------------------------------------------------------
subroutine apply_type1_short_repulsion(this)
!-----------------------------------------------------------------------------
class (short_repulsion_type1_params),intent(in) :: this

real(8) :: coef(2,2), fcoef, ff0(3), dr(3), dr1, dr2, dri, dri_eta
integer :: i, j, j1, ity, jty

do i=1, NATOMS

   ity = nint(atype(i))

   do j1 = 1, nbrlist(i,0)

      j = nbrlist(i,j1)
      jty = nint(atype(j))

      dr(1:3) = pos(i,1:3) - pos(j,1:3)
      dr2 = sum(dr(1:3)*dr(1:3))
      dr1 = sqrt(dr2)
      dri = 1d0/dr1
      dri_eta = dri**(this%eta(ity,jty)+1)
      fcoef = this%eta(ity,jty)*this%coef(ity,jty)*dri_eta*dri
      !if(l2g(atype(i))==1.and.dr1<3.d0)
      !print'(2i6,2f12.5)',l2g(atype(i)),l2g(atype(j)),dr1,fcoef

      f(i,1:3) = f(i,1:3) + fcoef*dr(1:3)

   enddo
enddo

end subroutine

!-----------------------------------------------------------------------------
subroutine short_repulsion(sr)
!-----------------------------------------------------------------------------
type(short_repulsion_type),intent(in) :: sr

if(.not. sr%has_short_repulsion) return

if (sr%potential_type==1) then
  call sr%p1%apply()
else if (sr%potential_type==2) then
  call sr%p2%apply()
else
  print*,'To be implemented'

endif


end subroutine

end module
