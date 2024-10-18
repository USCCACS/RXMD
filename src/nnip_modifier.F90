module nnip_modifier
  use mpi_mod
  implicit none

  type :: ceiling_params_type
     logical :: flag = .false.
     real(8) :: max_disp, max_velocity
  contains
     procedure :: show => show_ceiling_params
  end type

  type :: shortrep_params_type
     logical :: flag = .false.
     integer :: num_types = 0
     real(8) :: A(3,3),sigma(3,3),rc(3,3)
     integer :: n(3,3)
  contains
     procedure :: show => show_shortrep_params
  end type

  type :: springpot_params_type
     logical :: flag = .false.
     integer :: num_types = 0
     real(8) :: A(3,3),r0(3,3),rc_in(3,3),rc_out(3,3)
  contains
     procedure :: show => show_springpot_params
  end type

  type(ceiling_params_type) :: ceparams
  type(shortrep_params_type) :: shparams
  type(springpot_params_type) :: spparams

contains

!==========================================================================================
  subroutine show_springpot_params(this)
!==========================================================================================
     class(springpot_params_type),intent(in) :: this
     integer :: i,j

     print'(2a)', '--- springpot parameter ', repeat('-',60)

     print'(a,i6)', '  num_types', this%num_types
     do i=1, this%num_types
     do j=1, this%num_types
        print'(a,2i3,4f10.3)','  ity,jty,A,r0,rc_in,rc_out: ', &
             i, j, this%A(i,j), this%r0(i,j), this%rc_in(i,j), this%rc_out(i,j)
     enddo; enddo

     print'(a)', repeat('-',80)
  end subroutine

!------------------------------------------------------------------------------------------
  function springpot_params_ctor(filename) result(p)
!------------------------------------------------------------------------------------------
     character(len=:),allocatable,intent(in) :: filename
     type(springpot_params_type) :: p
     integer :: iunit, i, j

     open(newunit=iunit, file=filename)
     read(iunit,*) p%num_types

     p%A=-1d0; p%r0=-1d0; p%rc_in=-1d0; p%rc_out=-1d0

     do 
       read(iunit,*,end=99) i,j,p%A(i,j), p%r0(i,j), p%rc_in(i,j), p%rc_out(i,j)
     enddo 
     99 close(iunit)

     do i=1, p%num_types-1
     do j=i+1, p%num_types
        p%A(j,i) = p%A(i,j)
        p%r0(j,i) = p%r0(i,j)
        p%rc_in(j,i) = p%rc_in(i,j)
        p%rc_out(j,i) = p%rc_out(i,j)
     enddo; enddo

     p%A = p%A*23.0609d0 ! eV to kcal/mol conversion

     p%flag = .true.

  end function

!------------------------------------------------------------------------------------------
  subroutine springpot(myrank, num_atoms, atype, pos, nbrlist, f, p, NBUFFER, MAXNEIGHBS) 
!------------------------------------------------------------------------------------------
  integer,intent(in) :: myrank, num_atoms, NBUFFER, MAXNEIGHBS
  real(8),intent(in),allocatable :: atype(:), pos(:,:)
  real(8),intent(in out),allocatable :: f(:,:)
  integer,intent(in) :: nbrlist(1:NBUFFER,0:MAXNEIGHBS)
  type(springpot_params_type),intent(in) :: p

  integer :: i,j,j1,ity,jty
  real(8) :: rr(0:3), ff(3), E, coef

  do i=1, num_atoms 
     ity = atype(i)
     do j1 = 1, nbrlist(i,0)
        j = nbrlist(i,j1)
        jty = atype(j)

        if (p%A(ity,jty)<0d0) cycle

        rr(1:3) = pos(i,1:3)-pos(j,1:3)
        rr(0) = sqrt(sum(rr(1:3)*rr(1:3)))

        if (p%rc_in(ity,jty) < rr(0) .and. rr(0) < p%rc_out(ity,jty)) then

        coef = -2d0*p%A(ity,jty)*(rr(0) - p%r0(ity,jty))/rr(0)
        ff(1:3) = coef*rr(1:3)

        !print'(a,2i6,es15.5,2f8.3,1x,8f8.3)', 
        !     'ity,jty,coef,r0,rc_in,rc_out,rr(0:3),ff(1:3) : ', &
        !      ity,jty,coef,p%r0(ity,jty),p%rc_in(ity,jty),p%rc_out(ity,jty),rr(0:3),ff(1:3)

        f(i,1:3) = f(i,1:3) + coef*rr(1:3)

        endif
     enddo
  enddo

  end subroutine

!==========================================================================================
  subroutine show_shortrep_params(this)
!==========================================================================================
     class(shortrep_params_type),intent(in) :: this
     integer :: i,j
     print'(2a)', '--- shortrep parameter ', repeat('-',60)
     print'(a,i6)', '  num_types', this%num_types
     do i=1, this%num_types
     do j=1, this%num_types
        print'(a,2i3,2f10.3,i4,f10.3)','  ity,jty,A,sigma,n,rc: ', &
             i, j, this%A(i,j), this%sigma(i,j), this%n(i,j), this%rc(i,j)
     enddo; enddo
     print'(a)', repeat('-',80)
  end subroutine
!------------------------------------------------------------------------------------------
  function shortrep_params_ctor(filename) result(p)
!------------------------------------------------------------------------------------------
     character(len=:),allocatable,intent(in) :: filename
     type(shortrep_params_type) :: p
     integer :: iunit, i, j
     open(newunit=iunit, file=filename)
     read(iunit,*) p%num_types
     p%A=1d99; p%n=1; p%sigma=1d0
     do
       read(iunit,*,end=99) i,j,p%A(i,j), p%sigma(i,j), p%n(i,j), p%rc(i,j)
     enddo
     99 close(iunit)
     do i=1, p%num_types-1
     do j=i+1, p%num_types
        p%A(j,i) = p%A(i,j)
        p%n(j,i) = p%n(i,j)
        p%sigma(j,i) = p%sigma(i,j)
        p%rc(j,i) = p%rc(i,j)
     enddo; enddo
     p%A = p%A*23.0609d0 ! eV to kcal/mol conversion
     p%flag = .true.
  end function

!------------------------------------------------------------------------------------------
  subroutine shortrep(myrank, num_atoms, atype, pos, nbrlist, f, p, NBUFFER, MAXNEIGHBS)
!------------------------------------------------------------------------------------------
  integer,intent(in) :: myrank, num_atoms, NBUFFER, MAXNEIGHBS
  real(8),intent(in),allocatable :: atype(:), pos(:,:)
  real(8),intent(in out),allocatable :: f(:,:)
  integer,intent(in) :: nbrlist(1:NBUFFER,0:MAXNEIGHBS)
  type(shortrep_params_type),intent(in) :: p
  integer :: i,j,j1,ity,jty
  real(8) :: rr(0:3), ff(3), E, f1, f2, arg, coef
  real(8) :: hpi=1.5707963267948966d0
  do i=1, num_atoms
     ity = atype(i)
     do j1 = 1, nbrlist(i,0)
        j = nbrlist(i,j1)
        jty = atype(j)
        rr(1:3) = pos(i,1:3)-pos(j,1:3)
        rr(0) = sqrt(sum(rr(1:3)*rr(1:3)))
        if (p%rc(ity,jty) < rr(0)) cycle
        f1 = p%A(ity,jty)*(p%sigma(ity,jty)/rr(0))**p%n(ity,jty)
        arg = hpi/p%rc(ity,jty)
        f2 = cos(arg*rr(0))
        E = f1*f2
        coef = (1d0/rr(0))*( f1*f2*p%n(ity,jty) + arg*sin(arg*rr(0)) )
        ff(1:3) = coef*rr(1:3)
        !print'(a,2i6,es15.5,2f8.3,1x,8f8.3)', 'ity,jty,coef,f1,f2,rc,rr(0:3),ff(1:3) : ', &
        !                              ity,jty,coef,f1,f2,p%rc(ity,jty),rr(0:3),ff(1:3)
        f(i,1:3) = f(i,1:3) + coef*rr(1:3)
     enddo
  enddo
  end subroutine

!==========================================================================================
  subroutine show_ceiling_params(this)
!==========================================================================================
     class(ceiling_params_type),intent(in) :: this
     print'(a)', repeat('-',80)
     print'(2a)', '--- ceiling parameter ', repeat('-',60)
     print'(a,2f10.3)', '  max_disp, max_velocity: ', this%max_disp, this%max_velocity
     print'(a)', repeat('-',80)
  end subroutine

!------------------------------------------------------------------------------------------
  function ceiling_params_ctor(filename) result(p)
!------------------------------------------------------------------------------------------
     character(len=:),allocatable,intent(in) :: filename
     type(ceiling_params_type) :: p
     integer :: iunit

     open(newunit=iunit, file=filename)
     read(iunit,*) p%max_disp, p%max_velocity
     close(iunit)

     p%flag = .true.

  end function

!------------------------------------------------------------------------------------------
  subroutine ceiling_array(myrank, num_atoms, array, multiplier, threshold, atype, dthm)
!------------------------------------------------------------------------------------------
     integer,intent(in) :: myrank, num_atoms
     real(8),intent(in) :: threshold, multiplier
     real(8),intent(in out),allocatable :: array(:,:)
     real(8),intent(in),allocatable,optional:: atype(:),dthm(:) 

     integer :: i, num_samples, ity, ierr
     real(8) :: ary(0:3), multiplier_
     logical :: is_force_ceiling = .false.

     ! must use dthm for force ceiling
     is_force_ceiling = present(dthm)

     num_samples = 0
     do i=1, num_atoms
        ary(1:3) = array(i,1:3)
        ary(0) = sqrt(sum(ary(1:3)*ary(1:3)))

        ! for force ceiling. dthm = 0.5*dt/mass
        multiplier_ = multiplier
        if (is_force_ceiling) then
           ity = nint(atype(i))
           multiplier_ = dthm(ity)
        endif

        if (ary(0)*multiplier > threshold) then
          array(i,1:3) = ary(1:3)/ary(0)*threshold
          num_samples = num_samples + 1
        endif
     enddo

     call MPI_ALLREDUCE(MPI_IN_PLACE, num_samples, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

     if(myrank==0) then
        if (is_force_ceiling) then
          print'(a,es10.2,i6)','INFO:max velocity, num_samples:     ', threshold, num_samples
        else
          print'(a,es10.2,i6)','INFO:max displacement, num_samples: ', threshold, num_samples
        endif
     endif

  end subroutine

end module
