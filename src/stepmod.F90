module step_modifier

use utils, only : getstr, get_boxparameters, update_box_params, xs2xu_inplace, xu2xs_inplace
use base

type strain_type
  real(8) :: lattice(6)
  contains
    procedure :: apply => apply_strain
end type

type efield_type
  real(8) :: e(3)
  contains
    procedure :: apply => apply_efield
end type

type step_type
  integer :: nsteps = 0
  real(8),allocatable :: charge(:)
  type(strain_type),pointer :: strain
  type(efield_type),pointer :: efield
  contains 
    procedure :: show => show_step_params
end type

type(step_type),allocatable :: step_params(:) 

character(len=:),allocatable,private :: token

contains 

!------------------------------------------------------------------------------
subroutine apply_efield(this, num_atoms, q, f)
!------------------------------------------------------------------------------
  class(efield_type),intent(in) :: this
  integer,intent(in) :: num_atoms
  real(8),intent(in out),allocatable :: f(:,:)
  real(8),intent(in),allocatable :: q(:)

  integer :: ia,ib

  do ia=1, num_atoms
     do ib=1, 3 
        f(ia,ib) = f(ia,ib) + q(ib)*this%e(ib)
     enddo
  enddo
end subroutine

!------------------------------------------------------------------------------
subroutine apply_strain(this, num_atoms, pos)
!------------------------------------------------------------------------------
  class(strain_type),intent(in) :: this
  integer,intent(in) :: num_atoms
  real(8),intent(in out),allocatable :: pos(:,:)

  integer :: ia, ib
  real(8) :: mat(3,3)

  call xu2xs_inplace(hhi,obox,num_atoms,pos)

  lata = lata*(1d0+this%lattice(1))
  latb = latb*(1d0+this%lattice(2))
  latc = latc*(1d0+this%lattice(3))
  lalpha = lalpha*(1d0+this%lattice(4))
  lbeta = lbeta*(1d0+this%lattice(5))
  lgamma = lgamma*(1d0+this%lattice(6))

  call get_boxparameters(mat,lata,latb,latc,lalpha,lbeta,lgamma)

  do ia = 1, 3
  do ib = 1, 3
     HH(ia,ib,0)=mat(ia,ib)
  enddo; enddo

  !print'(a,6es15.6)','lata,latb,latc,lalpha,lbeta,lgamma', lata,latb,latc,lalpha,lbeta,lgamma
  !print'(a,3f15.3)','', HH(1:3,1,0)
  !print'(a,3f15.3)','', HH(1:3,2,0)
  !print'(a,3f15.3)','', HH(1:3,3,0)

  call update_box_params(vprocs, vid, hh, lata, latb, latc, maxrc, &
                       cc, lcsize, hhi, mdbox, lbox, obox)

  call xs2xu_inplace(hh,obox,num_atoms,pos)

end subroutine

!------------------------------------------------------------------------------
subroutine show_step_params(this)
!------------------------------------------------------------------------------
  class(step_type),intent(in) :: this
  if(associated(this%strain)) &
    print'(a,6es10.2)', 'strain(lattice) : ', this%strain%lattice
  if(associated(this%efield)) &
    print'(a,3es10.2)', 'efiled(e) : ', this%efield%e(1:3)
end subroutine

!------------------------------------------------------------------------------
subroutine get_tokens_and_append_step_params(linein, steps)
!------------------------------------------------------------------------------
  character(len=:),allocatable,intent(in out) :: linein
  type(step_type),allocatable,intent(in out) :: steps(:)

  type(step_type) :: c

  if (getstr(linein, token) < 0) stop 'error while reading nsteps'
  read(token, *) c%nsteps

  do while (getstr(linein, token) > 0) 
     select case(token)
       case("strain")
         if(.not.associated(c%strain)) allocate(c%strain)

         if (getstr(linein, token) < 0) stop 'error while reading strain lattice1'
         read(token, *) c%strain%lattice(1)
         if (getstr(linein, token) < 0) stop 'error while reading strain lattice2'
         read(token, *) c%strain%lattice(2)
         if (getstr(linein, token) < 0) stop 'error while reading strain lattice3'
         read(token, *) c%strain%lattice(3)
         if (getstr(linein, token) < 0) stop 'error while reading strain lattice4'
         read(token, *) c%strain%lattice(4)
         if (getstr(linein, token) < 0) stop 'error while reading strain lattice5'
         read(token, *) c%strain%lattice(5)
         if (getstr(linein, token) < 0) stop 'error while reading strain lattice6'
         read(token, *) c%strain%lattice(6)

       case("efield")
         if(.not.associated(c%efield)) allocate(c%efield)

         if (getstr(linein, token) < 0) stop 'error while reading efield x'
         read(token, *) c%efield%e(1)
         if (getstr(linein, token) < 0) stop 'error while reading efield y'
         read(token, *) c%efield%e(2)
         if (getstr(linein, token) < 0) stop 'error while reading efield z'
         read(token, *) c%efield%e(3)

     end select
  enddo

  steps = [steps, c]

end subroutine

!-------------------------------------------------------------------------------------------
function step_params_ctor(filename) result(steps)
!-------------------------------------------------------------------------------------------
!character(len=:),allocatable,intent(in) :: filename
character(*) :: filename
character(256) :: linein0
character(len=:),allocatable :: linein

integer :: iunit

type(step_type),allocatable :: steps(:) 

allocate(steps(0))

open(newunit=iunit, file=filename, form='formatted', status='old')

do while (.true.)
  read(iunit,'(a)',end=10) linein0
  linein = trim(adjustl(linein0))
  if (getstr(linein, token) > 0) then
     if(token=='step') call get_tokens_and_append_step_params(linein, steps)
  endif
end do
10 close(iunit)

if(myid==0) then
  do ia=1, size(steps)
    print'(a,i3,a)',"=== step ", ia, " ==="
    call steps(ia)%show()
  enddo
endif

end function

end module
