module msd_mod
  use utils
  use mpi_mod
 
  implicit none

  type msd_type

    logical :: is_msd
    real(8) :: unit_time

    real(8),allocatable :: dat(:,:)
    integer,allocatable :: num_samples(:,:)
    type(string_array),allocatable :: msd_atom_name(:)
  

    contains 
      procedure :: measure => measure_msd
      procedure :: save => save_msd

  end type

  type(msd_type) :: msd_data

contains

  function msd_type_ctor(atom_name, msd_size, unit_time) result(m)
     character(2),allocatable,intent(in) :: atom_name(:)
     integer,intent(in) :: msd_size
     real(8),intent(in) :: unit_time
     type(msd_type) :: m

     integer :: i,idx

     if(find_cmdline_argc('--msd',idx).or.find_cmdline_argc('-msd',idx)) then 
       m%is_msd = .true.
     else
       m%is_msd = .false.
       return
     endif

     if(msd_size<1) then
       write(6,fmt='(a,i9)') 'ERROR: msd_size <1', msd_size
       stop
     endif

     m%unit_time = unit_time
     
     allocate(m%msd_atom_name(size(atom_name)))
     allocate(m%dat(size(atom_name),msd_size), m%num_samples(size(atom_name),msd_size))

     do i=1, size(atom_name)
        m%msd_atom_name(i)%str = atom_name(i)
     enddo

     m%dat = 0.d0; m%num_samples = 0  

  end function

  subroutine measure_msd(this, num_atoms, atype, pos, pos0)
     class(msd_type),intent(in out) :: this 
     integer,intent(in) :: num_atoms 
     real(8),allocatable,intent(in) :: atype(:), pos(:,:), pos0(:,:) 

     integer,save :: counter = 0

     integer :: i, ity
     real(8) :: dr(3), dr2

     counter = counter + 1

     do i=1, num_atoms
        dr(1:3) = pos(i,1:3) - pos0(i,1:3)
        ity = nint(atype(i)) 
        dr2 = sum(dr(1:3)*dr(1:3))
        this%dat(ity,counter) = this%dat(ity,counter) + dr2
        this%num_samples(ity,counter) = this%num_samples(ity,counter) + 1
     enddo

  end subroutine

  subroutine save_msd(this)
    class(msd_type) :: this
    integer :: i,t,iunit,ierr, num_msdsize, num_type

    real(8),allocatable :: a(:)
    integer,allocatable :: b(:)

    open(newunit=iunit,file='msd.dat',form='formatted')
    
    num_type = size(this%dat,dim=1)
    num_msdsize = size(this%dat,dim=2)

    a = reshape(this%dat, [num_type*num_msdsize]) 
    b = reshape(this%num_samples, [num_type*num_msdsize])

    call MPI_ALLREDUCE(MPI_IN_PLACE, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, b, size(b), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    this%dat = reshape(a, [num_type,num_msdsize])
    this%num_samples = reshape(b, [num_type,num_msdsize])

    write(iunit,fmt='(a)',advance='no') ' time(fs) '
    do i=1, size(this%msd_atom_name)
       write(iunit,fmt='(a15)',advance='no') this%msd_atom_name(i)%str
    enddo
    do i=1, size(this%msd_atom_name)
       write(iunit,fmt='(a9)',advance='no') 'Num_'//this%msd_atom_name(i)%str
    enddo
    write(iunit,fmt=*)

    do t=1, num_msdsize-1
       write(iunit,fmt='(f10.2)',advance='no') t*this%unit_time
       do i=1, num_type
          if(this%num_samples(i,t)>0) then
             write(iunit,fmt='(es18.5)',advance='no') &
                   this%dat(i,t)/this%num_samples(i,t)
          else
             write(iunit,fmt='(es18.5)',advance='no') 0.d0
          endif
       enddo
       do i=1, num_type
          if(this%num_samples(i,t)>0) then
             write(iunit,fmt='(i12)',advance='no') this%num_samples(i,t)
          else
             write(iunit,fmt='(i12)',advance='no') 0
          endif
       enddo
       write(iunit,fmt=*)
    enddo
    close(iunit)

  end subroutine 

end module
