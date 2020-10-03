module mod_short_repulsion

use base
use fileio, only : output
use utils, only : l2g, find_cmdline_argc, get_command_argument_str, Ekcal_j, assert, int_to_str

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

  real(8) :: alpha_bond, alpha_angle, vscale
  real(8) :: rc_oh

  real(8) :: Kr, Kq
  real(8) :: r0, q0

  real(8) :: rc_freeze

  real(8) :: beta_1, beta_2
  real(8) :: beta_s1, beta_l1, beta_s2, beta_l2

  integer :: htype, otype

  real(8) :: fcut_o, fcut_h, ffactor

  character(2),allocatable :: atom_name(:)

  logical :: is_initialized = .false.

  logical :: does_freeze_atoms = .false.
  logical :: does_flip_velocity = .false.

  logical,allocatable :: frozen_atom_flag(:)
  real(8),allocatable :: r_oh(:,:)

  contains 
    procedure :: read => read_type2_params
    procedure :: print => print_type2_params
    procedure :: apply2 => apply_type2_short_repulsion
    procedure :: apply3or4 => apply_type3or4_short_repulsion
    procedure :: flipv => flip_frozenh_velocity
    procedure :: freezex => freeze_h_in_overstreched_bond
end type

type short_repulsion_type

  logical :: has_short_repulsion=.false.

  integer :: potential_type = 0 

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
          print'(a)', 'ERROR: short_rep flag is given but could not find shortrep.in file. Exiting rxmd.'
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

if(myid==0) then
  print'(a)', repeat('-',60)
  print'(a)', 'parameters for short repulsion: '
  print'(a,i3)', 'potential type: ', sr%potential_type
  print'(a)', repeat('-',60)
endif

if(sr%potential_type==1) then
   call sr%p1%read(funit)
   call sr%p1%print()
else if(sr%potential_type==2.or.sr%potential_type==3.or.sr%potential_type==4.or.sr%potential_type==5) then
   call sr%p2%read(funit, atom_name, sr%potential_type) 
   call sr%p2%print(sr%potential_type)
else
   print'(a)', 'ERROR: unsupported potential type. Exiting the code', sr%potential_type
   stop
endif

if(myid==0) print'(a)', repeat('-',60)

close(funit)

end subroutine

!-----------------------------------------------------------------------------
subroutine read_type2_params(this, funit, atom_name, potential_type)
!-----------------------------------------------------------------------------
  implicit none
  class(short_repulsion_type2_params) :: this
  character(2),allocatable,intent(in) :: atom_name(:)
  integer,intent(in) :: potential_type

  integer,intent(in) :: funit 
  integer :: ity, jty
  real(8) :: sigma_sum

  call assert(2<=potential_type .and. potential_type<=5, &
     'unsupported potential type in read_type2_params(). Exiting rxmd.', myid)

  this%htype=-1; this%otype=-1; 
  do ity=1, size(atom_name)
     if(index(atom_name(ity),"H") /=0) this%htype=ity
     if(index(atom_name(ity),"O") /=0) this%otype=ity
  enddo
  call assert(this%htype==-1, 'could not find H type', myid)
  call assert(this%htype==-1, 'could not find O type', myid)

  read(funit,*) this%alpha_bond, this%alpha_angle, this%rc_oh
  read(funit,*) this%Kr, this%Kq
  read(funit,*) this%r0, this%q0

  if(potential_type == 3) then
     read(funit,*) this%does_freeze_atoms, this%rc_freeze
     read(funit,*) this%does_flip_velocity, this%vscale

     allocate(this%frozen_atom_flag(NBUFFER))
     allocate(this%r_oh(NBUFFER,0:3))
  else if(potential_type == 4.or.potential_type == 5) then
     read(funit,*) this%beta_1, this%beta_s1, this%beta_l1
     read(funit,*) this%beta_2, this%beta_s2, this%beta_l2
  endif
  if(potential_type==5) then
     read(funit,*) this%fcut_o, this%fcut_h, this%ffactor
  endif

  this%is_initialized = .true.

end subroutine

!-----------------------------------------------------------------------------
subroutine print_type2_params(this, potential_type)
!-----------------------------------------------------------------------------
  implicit none
  class(short_repulsion_type2_params) :: this
  integer,intent(in) :: potential_type

  if(myid==0) then
    print'(a30,a7,i3,1x,a7,i3)', "atom types: ", "H -", this%htype, "O -", this%otype
    print'(a30,3f10.3)', 'alpha_bond, alpha_angle, rc_oh: ', this%alpha_bond, this%alpha_angle, this%rc_oh
    print'(a30,2f10.3)', 'Kr, Kq: ', this%Kr, this%Kq
    print'(a30,2f10.3)', 'r0, q0: ', this%r0, this%q0
    if(potential_type==3) then 
       print'(a30,l6,f10.3)',  'rc_freeze: ', this%does_freeze_atoms, this%rc_freeze
       print'(a30,l6,f10.3)',  'vscale: ', this%does_flip_velocity, this%vscale
    endif
    if(potential_type==4.or.potential_type==5) then
       print'(a30,3f10.3)','beta_1, beta_s1, beta_l1: ', this%beta_1, this%beta_s1, this%beta_l1
       print'(a30,3f10.3)','beta_2, beta_s2, beta_l2: ', this%beta_2, this%beta_s2, this%beta_l2
    endif
    if(potential_type==5) then
       print'(a30,3f10.3)','fcut_o, fcut_h, ffactor: ', this%fcut_o, this%fcut_h, this%ffactor
    endif
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

!-----------------------------------------------------------------------------
subroutine apply_type3or4_short_repulsion(this, potential_type)
!-----------------------------------------------------------------------------
class (short_repulsion_type2_params),intent(in) :: this
integer,intent(in) :: potential_type

real(8) :: coef(2,2), fcoef, ff0(3), dr(3), dr1, dr2, dri, dri_eta
integer :: i, j, k, l, l1, ity, jty, kty, lty, igid, jgid, kgid

real(8) :: rij(0:3), rik(0:3)
real(8) :: sine, cosine, theta, Krcoef, Kqcoef, ff(3)

real(8),parameter :: MAX_DIST = 1e9, pi=3.14159265358979d0

logical :: is_found_k, is_found_j
real(8) :: dr_k, dr_j

real(8) :: bond_coeff

real(8),parameter :: rc_inner=0.97d0, rc_outer=1.05d0

character(len=:),allocatable :: info_ij, info_ik
character(len=:),allocatable :: filebase

!call set_force_zero(rc_inner,rc_outer)

do i=1, NATOMS                           ! O

   ity = nint(atype(i))

   if(ity/=this%otype) cycle

   !if(myid==0) print'(a,2i6)','ity, this%otype:',ity,this%otype

   igid = l2g(atype(i))

   is_found_j=.false.; jgid=-1; dr_j=0d0; jty=-1
   is_found_k=.false.; kgid=-1; dr_k=0d0; kty=-1

   do l1 = 1, nbrlist(i,0)

      l = nbrlist(i,l1)
      lty = nint(atype(l))

      if(lty/=this%htype) cycle

      if(l2g(atype(l))==igid+1) then
        j = l   !H1
        is_found_j = .true.
        jgid = l2g(atype(j))
        jty = nint(atype(j))
      endif
      if(l2g(atype(l))==igid+2) then
        k = l   !H2
        is_found_k = .true.
        kgid = l2g(atype(k))
        kty = nint(atype(k))
      endif
   enddo

   info_ij = int_to_str(myid)//' '//int_to_str(igid)//' '//int_to_str(jgid)//' '//int_to_str(ity)//' '//int_to_str(jty)
   info_ik = int_to_str(myid)//' '//int_to_str(igid)//' '//int_to_str(kgid)//' '//int_to_str(ity)//' '//int_to_str(kty)

   call assert(is_found_j, 'ERROR: wrong global ID for j '//info_ij)
   call assert(is_found_k, 'ERROR: wrong global ID for k '//info_ik)

   call assert(jty==this%htype, 'ERROR: j-atom type is not H '//info_ij)
   call assert(kty==this%htype, 'ERROR: k-atom type is not H '//info_ik)

   rij(1:3) = pos(i,1:3)-pos(j,1:3)
   rij(0) = sqrt(sum(rij(1:3)*rij(1:3)))
   rij(1:3) = rij(1:3)/rij(0)

   call assert(rij(0)>0.8d0, 'ERROR: j-atom is too close '//info_ij, val=rij(0))
   call assert(rij(0)<1.25d0, 'ERROR: j-atom is too far '//info_ik, val=rij(0))

   rik(1:3) = pos(i,1:3)-pos(k,1:3)
   rik(0) = sqrt(sum(rik(1:3)*rik(1:3)))
   rik(1:3) = rik(1:3)/rik(0)

   call assert(rik(0)>0.8d0, 'ERROR: k-atom is too close '//info_ij, val=rik(0))
   call assert(rik(0)<1.25d0, 'ERROR: k-atom is too far '//info_ik, val=rik(0))

   cosine = sum(rij(1:3)*rik(1:3))
   theta = acos(cosine)
   sine = sin(theta)

   call assert(abs(sine)>1d-9, 'ERROR: too small H-O-H angle ', myid)

   bond_coeff = this%alpha_bond
   if(potential_type==4) then
     if(rij(0)<this%beta_s1.or.this%beta_l1<rij(0)) bond_coeff = bond_coeff + this%beta_1
     if(rij(0)<this%beta_s2.or.this%beta_l2<rij(0)) bond_coeff = bond_coeff + this%beta_2
   endif

   !Krcoef = bond_coeff * this%Kr*(rij(0)-this%r0) ! ij
   Krcoef = bond_coeff * get_smooth_spring_coeff(rij(0),this%Kr,this%r0) ! ij
   if(potential_type==4) then
     if(rij(0)<this%beta_s1.or.this%beta_l1<rij(0)) &
       print'(a,f8.3,f8.3,e15.5)', 'beta1: myid,i,j,ity,jty: '//info_ij,bond_coeff,rij(0),Krcoef
     if(rij(0)<this%beta_s2.or.this%beta_l2<rij(0)) &
       print'(a,f8.3,f8.3,e15.5)', 'beta2: myid,i,j,ity,jty: '//info_ij,bond_coeff,rij(0),Krcoef
   endif

   ff(1:3) = Krcoef*rij(1:3)
   f(i,1:3) = f(i,1:3) - ff(1:3)
   f(j,1:3) = f(j,1:3) + ff(1:3)
   !write(unit=6,fmt='(a,4f12.3)') 'Krcoef, ff(1:3): ', Krcoef, ff(1:3)

   bond_coeff = this%alpha_bond
   if(potential_type==4) then
     if(rik(0)<this%beta_s1.or.this%beta_l1<rik(0)) bond_coeff = bond_coeff + this%beta_1
     if(rik(0)<this%beta_s2.or.this%beta_l2<rik(0)) bond_coeff = bond_coeff + this%beta_2
   endif

   !Krcoef = bond_coeff * this%Kr*(rik(0)-this%r0) ! ik
   Krcoef = bond_coeff * get_smooth_spring_coeff(rik(0),this%Kr,this%r0) ! ik
   if(potential_type==4) then
     if(rik(0)<this%beta_s1.or.this%beta_l1<rik(0)) &
        print'(a,f8.3,f8.3,e15.5)','beta1: myid,i,k,ity,kty: '//info_ij,bond_coeff,rik(0),Krcoef
     if(rik(0)<this%beta_s2.or.this%beta_l2<rik(0)) &
        print'(a,f8.3,f8.3,e15.5)','beta2: myid,i,k,ity,kty: '//info_ik,bond_coeff,rik(0),Krcoef
   endif

   ff(1:3) = Krcoef*rik(1:3)
   f(i,1:3) = f(i,1:3) - ff(1:3)
   f(k,1:3) = f(k,1:3) + ff(1:3)
   !write(unit=6,fmt='(a,4f12.3)') 'Krcoef, ff(1:3): ', Krcoef, ff(1:3)

   Kqcoef = this%alpha_angle * this%Kq*(theta*180d0/pi-this%q0)*(-1.d0/sine)

   rij(1:3)=-rij(1:3)
   !print*,Kqcoef, theta*180d0/pi, this%q0
   call ForceA3(Kqcoef,j,i,k,rij,rik)

enddo

contains 

!-----------------------------------------------------------------------------
subroutine set_force_zero(rc_inner, rc_outer)
!-----------------------------------------------------------------------------
real(8),intent(in):: rc_inner, rc_outer
character(len=:),allocatable :: info_ij, info_ik

do i=1, NATOMS                           ! O

   ity = nint(atype(i))

   if(ity/=this%otype) cycle

   !if(myid==0) print'(a,2i6)','ity, this%otype:',ity,this%otype

   igid = l2g(atype(i))

   is_found_j=.false.; jgid=-1; dr_j=0d0; jty=-1
   is_found_k=.false.; kgid=-1; dr_k=0d0; kty=-1

   do l1 = 1, nbrlist(i,0)

      l = nbrlist(i,l1)
      lty = nint(atype(l))

      if(lty/=this%htype) cycle

      if(l2g(atype(l))==igid+1) then
        j = l   !H1
        is_found_j = .true.
        jgid = l2g(atype(j))
        jty = nint(atype(j))
      endif
      if(l2g(atype(l))==igid+2) then
        k = l   !H2
        is_found_k = .true.
        kgid = l2g(atype(k))
        kty = nint(atype(k))
      endif
   enddo

   info_ij = int_to_str(myid)//' '//int_to_str(igid)//' '//int_to_str(jgid)//' '//int_to_str(ity)//' '//int_to_str(jty)
   info_ik = int_to_str(myid)//' '//int_to_str(igid)//' '//int_to_str(kgid)//' '//int_to_str(ity)//' '//int_to_str(kty)

   call assert(is_found_j, 'ERROR(sfz): sfz wrong global ID for j '//info_ij)
   call assert(is_found_k, 'ERROR(sfz): wrong global ID for k '//info_ij)

   call assert(jty==this%htype, 'ERROR(sfz): j-atom type is not H '//info_ij)
   call assert(kty==this%htype, 'ERROR(sfz): k-atom type is not H '//info_ij)

   rij(1:3) = pos(i,1:3)-pos(j,1:3)
   rij(0) = sqrt(sum(rij(1:3)*rij(1:3)))

   if(rij(0)<rc_inner.or.rc_outer<rij(0)) then
     print'(a,f6.3,7e11.3)','fzero: myid,i,k,ity,kty,rij,fi,fj: '//info_ij,rij(0),f(i,1:3),f(j,1:3)
     !f(i,1:3)=0.d0
     f(j,1:3)=0.d0
   endif

   rik(1:3) = pos(i,1:3)-pos(k,1:3)
   rik(0) = sqrt(sum(rik(1:3)*rik(1:3)))

   if(rik(0)<rc_inner.or.rc_outer<rik(0)) then
     print'(a,f6.3,7e11.3)','fzero: myid,i,k,ity,kty,rik,fi,fk: '//info_ij,rik(0),f(i,1:3),f(k,1:3)
     !f(i,1:3)=0.d0
     f(k,1:3)=0.d0
   endif

   return
enddo


end subroutine


!-----------------------------------------------------------------------------
function get_smooth_spring_coeff(r, Kr, r0_oh) result(NewSp)
!-----------------------------------------------------------------------------
!real(8),parameter :: A=9580.16d0, B=0.0286d0 ! r0_oh = 0.971d0
!real(8),parameter :: A=24520.6d0, B=0.049d0  ! r0_oh = 0.971d0
!real(8),parameter :: A=92294.2d0, B=-0.02d0    ! r0_oh = 0.970d0
!real(8),parameter :: A=92294.2d0, B=0.05d0    ! r0_oh = 0.970d0 !<- (9/11)
!real(8),parameter :: A=92294.2d0, B=0.02d0    ! r0_oh = 0.970d0 !<- (9/12)
!real(8),parameter :: A=92294.2d0, B=0.05d0    ! r0_oh = 0.970d0 !<- (9/13)
!real(8),parameter :: A=92294.2d0, B=0.03d0    ! r0_oh = 0.970d0 !<- (9/21)
real(8),parameter :: A=92294.2d0, B=0.05d0    ! r0_oh = 0.970d0 !<- (9/22)
!real(8),parameter :: A=277.78d0, B=0.0d0    ! r0_oh = 0.970d0
real(8),intent(in) :: r, Kr, r0_oh

real(8) :: S, Sp, F, Fp, NewSp

!Energy = S*F

F = A*(r - r0_oh - B)**4
Fp = 4d0*A*(r - r0_oh - B)**3
!F = A*(r - r0_oh)**2
!Fp = 2d0*A*(r - r0_oh)

S = 0.5d0*Kr*(r - r0_oh)**2
Sp = Kr*(r - r0_oh)

NewSp = S*Fp + Sp*F

return

end function

end subroutine

!-----------------------------------------------------------------------------
subroutine flip_frozenh_velocity(this)
!-----------------------------------------------------------------------------
class (short_repulsion_type2_params),intent(in out) :: this

integer :: i, ity
real(8) :: vv

if(.not. this%does_flip_velocity) return

call update_frozen_atom_flag(this)

do i=1, NATOMS
   if(this%frozen_atom_flag(i)) then
     ity = nint(atype(i))
     call assert(ity==this%htype, 'wrong atom type to flip velocity '//int_to_str(ity)//' '//int_to_str(i), myid)
     vv = sqrt(sum(v(i,1:3)*v(i,1:3)))
     !v(i,1:3) = vv*this%r_oh(i,1:3)
     v(i,1:3) = this%vscale*vv*this%r_oh(i,1:3)
     print'(a,3i6,4f8.3,1x,3(3f8.3,1x))','V: myid,i,ity,this%r_oh(0:3): ', &
        myid,l2g(atype(i)),ity,this%r_oh(i,0:3), pos(i,1:3),v(i,1:3),f(i,1:3)
   endif
enddo

end subroutine

!-----------------------------------------------------------------------------
subroutine freeze_h_in_overstreched_bond(this)
!-----------------------------------------------------------------------------
class (short_repulsion_type2_params),intent(in out) :: this

integer :: i, ity

if(.not. this%does_freeze_atoms) return

call update_frozen_atom_flag(this)

do i=1, NATOMS
   if(this%frozen_atom_flag(i)) then

     ity = nint(atype(i))
     call assert(ity==this%htype, 'wrong atom type to freeze position '//int_to_str(ity)//' '//int_to_str(i), myid)

     pos(i,1:3) = pos(i,1:3)-dt*v(i,1:3)

     print'(a,3i6,4f8.3,1x,3(3f8.3,1x))','X: myid,i,ity,this%r_oh(0:3): ', &
        myid,l2g(atype(i)),ity,this%r_oh(i,0:3),pos(i,1:3),v(i,1:3),f(i,1:3)
   endif
enddo

end subroutine

!-----------------------------------------------------------------------------
subroutine update_frozen_atom_flag(this)
!-----------------------------------------------------------------------------
class (short_repulsion_type2_params),intent(in out) :: this

real(8) :: ff0(0:3), fr, rij(0:3), rik(0:3)
integer :: i, j, k, l, ity, lty, l1, igid, jgid

this%frozen_atom_flag = .false.

do i=1, NATOMS                           ! O

   ity = nint(atype(i))

   igid = l2g(atype(i))

   if(ity==this%otype) then
   
      do l1 = 1, nbrlist(i,0)
   
         l = nbrlist(i,l1)
         lty = nint(atype(l))
   
         if(lty/=this%htype) cycle
   
         if(l2g(atype(l))==igid+1) j = l   !H1
         if(l2g(atype(l))==igid+2) k = l   !H2
      enddo
   
      rij(1:3) = pos(i,1:3)-pos(j,1:3)
      rij(0) = sqrt(sum(rij(1:3)*rij(1:3)))
      rij(1:3) = rij(1:3)/rij(0)
      if(rij(0) > this%rc_freeze) then
         this%frozen_atom_flag(j) = .true.
         this%r_oh(j,0:3) = rij(0:3)
         print'(a,3i6,4f8.3,2(3f8.3,1x))','o->h1:myid,i,j,rij(0:3),pos(i,1:3),pos(j,1:3) ', &
            myid,l2g(atype(i)),l2g(atype(j)),rij(0:3),pos(i,1:3),pos(j,1:3)
      endif
   
      rik(1:3) = pos(i,1:3)-pos(k,1:3)
      rik(0) = sqrt(sum(rik(1:3)*rik(1:3)))
      rik(1:3) = rik(1:3)/rik(0)
      if(rik(0) > this%rc_freeze) then
         this%frozen_atom_flag(k) = .true.
         this%r_oh(k,0:3) = rik(0:3)
         print'(a,3i6,4f8.3,1x,2(3f8.3,1x))','o->h2:myid,i,k,rik(0:3),pos(i,1:3),pos(k,1:3) ', &
            myid,l2g(atype(i)),l2g(atype(k)),rik(0:3),pos(i,1:3),pos(k,1:3)
      endif
   
   else if(ity==this%htype) then
   
      do l1 = 1, nbrlist(i,0)
         l = nbrlist(i,l1)
         lty = nint(atype(l))
   
         if(lty/=this%otype) cycle
         jgid = l2g(atype(l))
   
         if(mod(igid,3)==2 .and. jgid == igid - 1) j = l ! O for H1
         if(mod(igid,3)==0 .and. jgid == igid - 2) j = l ! O for H2
      enddo
   
      rij(1:3) = pos(j,1:3)-pos(i,1:3)        ! r_O - r_H
      rij(0) = sqrt(sum(rij(1:3)*rij(1:3)))
      rij(1:3) = rij(1:3)/rij(0)

      if(rij(0) > this%rc_freeze) then
         this%frozen_atom_flag(i) = .true.    ! freeze H
         this%r_oh(i,0:3) = rij(0:3)          ! r_OH
         print'(a,3i6,4f8.3,2(3f8.3,1x))','h->o:myid,i,j,rij(0:3),pos(i,1:3),pos(j,1:3) ', &
            myid,l2g(atype(i)),l2g(atype(j)),rij(0:3),pos(j,1:3),pos(i,1:3)
      endif
   
   endif

enddo

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
   if(nbr_count < 2) then
      print'(a)',repeat('-=',30)
      print'(a,2i6)', 'ERROR: less than two neighbor hydrogen found : ', myid, nbr_count
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

   Krcoef = this%alpha_bond * this%Kr*(rij(0)-this%r0) ! ij
   ff(1:3) = Krcoef*rij(1:3)
   f(i,1:3) = f(i,1:3) - ff(1:3)
   f(j,1:3) = f(j,1:3) + ff(1:3)
   !write(unit=6,fmt='(a,4f12.3)') 'Krcoef, ff(1:3): ', Krcoef, ff(1:3)

   Krcoef = this%alpha_bond * this%Kr*(rik(0)-this%r0) ! ik
   ff(1:3) = Krcoef*rik(1:3)
   f(i,1:3) = f(i,1:3) - ff(1:3)
   f(k,1:3) = f(k,1:3) + ff(1:3)
   !write(unit=6,fmt='(a,4f12.3)') 'Krcoef, ff(1:3): ', Krcoef, ff(1:3)

   call assert(abs(sine)>1d-9, 'ERROR: too small H-O-H angle', myid)

   Kqcoef = this%alpha_angle * this%Kq*(theta*180d0/pi-this%q0)*(-1.d0/sine)

   rij(1:3)=-rij(1:3)
   !print*,Kqcoef, theta*180d0/pi, this%q0
   call ForceA3(Kqcoef,j,i,k,rij,rik)

   !if(myid==0) print'(a,3i6,3f8.3)','i,j,k,dist_1st,dist_2nd, theta: ', &
   !    i, j, k, dist_1st, dist_2nd, theta


enddo

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
character(len=:),allocatable :: filebase
real(8),allocatable :: ftmp(:,:)

if(.not. sr%has_short_repulsion) return

if (sr%potential_type==1) then
  call sr%p1%apply()
else if (sr%potential_type==2) then
  call sr%p2%apply2()
else if (sr%potential_type==3.or.sr%potential_type==4.or.sr%potential_type==5) then
  call sr%p2%apply3or4(sr%potential_type)
else
  print*,'To be implemented'
endif

end subroutine

subroutine force_cutoff(num_atoms, atype, f, myrank, sr)
   implicit none
   integer,intent(in) :: num_atoms, myrank
   real(8),allocatable,intent(in) ::  atype(:)
   real(8),allocatable,intent(in out) ::  f(:,:)
   type(short_repulsion_type),intent(in) :: sr

   integer :: i,ia,ity
   real(8) :: fcut, fset

   do i = 1, num_atoms
      ity = nint(atype(i))

      if(ity==sr%p2%htype) then
          fcut = sr%p2%fcut_h
          fset = fcut*sr%p2%ffactor
      else if(ity==sr%p2%otype) then
          fcut = sr%p2%fcut_o
          fset = fcut*sr%p2%ffactor
      else
          print*,'ERROR: atomtype was not found.', myid, ity, i, l2g(atype(i))
          stop
      endif

      do ia=1,3
         if(f(i,ia) > fcut) then
           print'(a,5i9,es15.5,2f8.2)','max force cutoff applied.', myid, ity, i, l2g(atype(i)), ia, f(i,ia), fcut, fset
           f(i,ia) = fset
         endif
         if(f(i,ia) < -fcut) then
           print'(a,5i9,es15.5,2f8.2)','min force cutoff applied.', myid, ity, i, l2g(atype(i)), ia, f(i,ia), fcut, fset
           f(i,ia) = -fset
         endif
      enddo

   enddo
end subroutine


end module
