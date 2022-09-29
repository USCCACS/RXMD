!------------------------------------------------------------------------------
module harmonic_potential_mod
!------------------------------------------------------------------------------

  use utils, only : find_cmdline_argc, getstr, l2g, token

  implicit none

  type single_atom
    integer :: id,k_type
    real(8) :: pos(3)
  end type

  type harmonic_potential_type

    logical :: apply = .false.
    character(len=:),allocatable :: filename 

    real(8),allocatable :: K(:) ! spring constants

    type(single_atom),allocatable :: atoms(:)

  end type

  type(harmonic_potential_type) :: harmo_pot

contains

!------------------------------------------------------------------------------
  subroutine show_harmonic_potential_params(hp)
!------------------------------------------------------------------------------
    type(harmonic_potential_type),intent(in) :: hp
    integer :: idx

    print'(a)',repeat('-',80)
    print'(a,l)','apply: ', hp%apply
    print'(a,a)','filename: ', hp%filename
    print'(a,i3,10f10.5)','num_K, K: ', size(hp%K), hp%K

    print'(a,i3)','num(atoms): ', size(hp%atoms)
    do idx = 1, size(hp%atoms)
       print'(a,i9,i3,3f10.5)','id,k_type,pos: ', hp%atoms(idx)%id,hp%atoms(idx)%k_type, hp%atoms(idx)%pos
    enddo
    print'(a)',repeat('-',80)

  end subroutine

  function harmonic_potential_ctor() result(harmo_pot)
    integer :: idx, iunit
    character(256) :: argv, linein0
    character(len=:),allocatable :: linein 

    type(harmonic_potential_type) :: harmo_pot

    if(find_cmdline_argc('--harmo_pot',idx).or.find_cmdline_argc('-harmo',idx)) then

      call get_command_argument(idx+1,argv)
      harmo_pot%filename = trim(adjustl(argv))

      open(newunit=iunit, file=harmo_pot%filename, form='formatted', status='old')

      do while (.true.)
        read(iunit,'(a)',end=10) linein0
        linein = trim(adjustl(linein0))

        if (getstr(linein, token) > 0) then
          if(token=='K') call get_K_params(linein, harmo_pot%K)
          if(token=='atom') call get_atom_params(linein, harmo_pot%atoms)
        endif 
      enddo
      10 close(iunit)

      if(size(harmo_pot%atoms) > 0) harmo_pot%apply = .true.
    endif
  end function

!------------------------------------------------------------------------------
  subroutine get_atom_params(linein, atoms)
!------------------------------------------------------------------------------
    character(len=:),allocatable,intent(in out) :: linein
    type(single_atom),allocatable,intent(in out) :: atoms(:)
    type(single_atom) :: a
    real(8) :: dummy

    if (getstr(linein, token) < 0) stop 'error while reading id in get_atom_params'
    read(token, *) a%id
    if (getstr(linein, token) < 0) stop 'error while reading k_type in get_atom_params'
    read(token, *) a%k_type
    if (getstr(linein, token) < 0) stop 'error while reading pos(1) in get_atom_params'
    read(token, *) dummy; a%pos(1) = dummy
    if (getstr(linein, token) < 0) stop 'error while reading pos(2) in get_atom_params'
    read(token, *) dummy; a%pos(2) = dummy
    if (getstr(linein, token) < 0) stop 'error while reading pos(3) in get_atom_params'
    read(token, *) dummy; a%pos(3) = dummy

    ! allocate zero-sized array
    if(.not.allocated(atoms)) allocate(atoms(0)) 

    atoms = [atoms, a]

  end subroutine

!------------------------------------------------------------------------------
  subroutine get_K_params(linein, K)
!------------------------------------------------------------------------------
    character(len=:),allocatable,intent(in out) :: linein
    real(8),allocatable :: K(:) ! spring constants
    integer :: idx, num_K

    if (getstr(linein, token) < 0) stop 'error while reading num_K'
    read(token, *) num_K

    allocate(K(num_K))

    do idx = 1, num_K
       if (getstr(linein, token) < 0) stop 'error while reading Kval'
       read(token, *) K(idx)
    enddo

  end subroutine

!------------------------------------------------------------------------------
  subroutine apply_harmonic_potential(hp, natoms, atype, pos, f)
!------------------------------------------------------------------------------
    type(harmonic_potential_type) :: hp
    integer,intent(in) :: natoms
    real(8),allocatable,intent(in) :: atype(:), pos(:,:)
    real(8),allocatable,intent(in out) :: f(:,:)

    integer :: i, idx, igd
    real(8) :: dr(3), ff(3)

    do i=1, natoms
       igd = l2g(atype(i))
       do idx = 1, size(hp%atoms)
          if (igd == hp%atoms(idx)%id) then
             dr(1:3) = pos(i,1:3) - hp%atoms(idx)%pos(1:3)
             ff(1:3) = hp%K(hp%atoms(idx)%k_type)*dr(1:3)
             f(i,1:3) = f(i,1:3) - ff(1:3)
             !print'(2i,3f)',i,igd,ff(1:3)
          endif
       enddo
    enddo
  end subroutine
  
end module

!program main
!  use harmonic_potential 
!  harmo_pot = harmonic_potential_ctor() 
!  call show_harmonic_potential_params(harmo_pot)
!end program 
