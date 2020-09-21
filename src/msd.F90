module msd_mod
  use base, only : myid
  use utils
  use mpi_mod
 
  implicit none

  type msd_type

    logical :: is_msd = .false.
    real(8) :: length_fs, interval_fs, freq_fs
    integer :: length_step, interval_step, freq_step

    real(8) :: onestep_fs = 0.d0

    integer :: num_init_pos = 0
    integer :: num_data_points= 0
    logical,allocatable :: has_init_pos(:)

    real(8),allocatable :: dat(:,:,:)
    integer,allocatable :: num_samples(:,:,:)
    type(string_array),allocatable :: msd_atom_name(:)

  end type

  type(msd_type) :: msd_data

contains

  subroutine msd_initialize(m, atom_name, total_steps, onestep_fs)

     type(msd_type),intent(in out) :: m
     character(2),allocatable,intent(in) :: atom_name(:)
     integer,intent(in) :: total_steps
     real(8),intent(in) :: onestep_fs 

     integer :: i,idx, num_atoms
     character(256) :: argv

     if(find_cmdline_argc('--msd',idx).or.find_cmdline_argc('-msd',idx)) then 
       m%is_msd = .true.
       call get_command_argument(idx+1,argv); read(argv,*) m%length_fs
       call get_command_argument(idx+2,argv); read(argv,*) m%interval_fs
       call get_command_argument(idx+3,argv); read(argv,*) m%freq_fs
     endif

     if (.not. m%is_msd) return

     m%onestep_fs= onestep_fs
     m%length_step = m%length_fs/m%onestep_fs
     m%interval_step = m%interval_fs/m%onestep_fs
     m%freq_step = m%freq_fs/m%onestep_fs

     call assert(total_steps > m%length_step,'ERROR: total_steps is less than MSD correlation length')
     call assert(m%length_step>0,  'ERROR: MSD correlation length <= 0')
     call assert(m%interval_step>0,'ERROR: MSD correlation interval <= 0')
     call assert(m%freq_step>0,    'ERROR: MSD measurement frequency <= 0')
     call assert(mod(m%interval_step,m%freq_step) == 0, 'ERROR: interval must be a multiple of frequency')

     ! if total_steps > m%length_step is true, at least one initial position is there
     ! regardless of interval and measurement frequency.
     m%num_init_pos = max(int((total_steps - m%length_step)/m%interval_step), 1)

     allocate(m%has_init_pos(m%num_init_pos))
     m%has_init_pos = .false.

     num_atoms = size(atom_name)

     allocate(m%msd_atom_name(num_atoms))
     do i=1, num_atoms
        m%msd_atom_name(i)%str = atom_name(i)
     enddo

     m%num_data_points = int(m%length_step/m%freq_step)
     call assert(m%num_data_points>0, 'ERROR: MSD num_data_points <= 0')

     allocate(m%dat(m%num_init_pos, num_atoms, m%num_data_points))
     allocate(m%num_samples(m%num_init_pos, num_atoms, m%num_data_points))

     m%dat = 0.d0; m%num_samples = 0

     call msd_print(m)

  end subroutine

  subroutine msd_print(m)

     type(msd_type),intent(in) :: m 

     if(myid==0) then
       print'(a)', repeat('-',60)
       if(m%is_msd) write(6,fmt='(a)') 'start MSD measurement'
       print'(a,l3)', 'is_msd: ', m%is_msd
       print'(a,f9.3,a,i9,a)', 'length   : ', m%length_fs,  ' [fs] ', m%length_step,   ' [MDstep]'
       print'(a,f9.3,a,i9,a)', 'interval : ', m%interval_fs,' [fs] ', m%interval_step, ' [MDstep]'
       print'(a,f9.3,a,i9,a)', 'frequency: ', m%freq_fs,    ' [fs] ', m%freq_step,     ' [MDstep]'
       print'(a,i6)', '# of initial positions: ', m%num_init_pos
       print'(a,i9)', '# of data per initial positions: ', m%num_data_points
       print'(a)', repeat('-',60)
     endif

  end subroutine

  subroutine msd_add_initial_pos(m, current_step, num_atoms, pos, pos0)
     type(msd_type),intent(in out) :: m 
     integer,intent(in) :: num_atoms, current_step ! current_step is zero-indexed
     real(8),allocatable,intent(in) :: pos(:,:)
     real(8),allocatable,intent(in out) :: pos0(:,:,:)

     integer :: ia

     if(.not. m%is_msd) return
     if(mod(current_step, m%interval_step) /= 0) return

     ia = current_step/m%interval_step + 1

     if(ia > m%num_init_pos) return

     m%has_init_pos(ia)=.true.
     pos0(ia,1:num_atoms,1:3) = pos(1:num_atoms,1:3)

     if(myid==0) then
       print'(a,i4,a1,i4,a,i9)', & 
           'INFO: A new initial position has been added. ia & current_step = ', &
           ia, '/', m%num_init_pos, ' current_step = ', current_step 
     endif

  end subroutine

  subroutine msd_measure(m, current_step, num_atoms, atype, pos, pos0)
     type(msd_type),intent(in out) :: m 
     integer,intent(in) :: current_step ! current_step is zero-indexed
     integer,intent(in) :: num_atoms 
     real(8),allocatable,intent(in) :: atype(:), pos(:,:), pos0(:,:,:) 

     integer :: i, ia, ity, idx 
     real(8) :: dr(3), dr2

     if (.not. m%is_msd) return

     if(mod(current_step, m%freq_step) /= 0) return

     do ia = 1, m%num_init_pos

        if(.not. m%has_init_pos(ia)) cycle

        idx = (current_step - (ia-1)*m%interval_step)/m%freq_step + 1

        if(0<idx .and. idx <= m%num_data_points) then 

           do i=1, num_atoms
              ity = nint(atype(i)) 

              dr(1:3) = pos(i,1:3) - pos0(ia,i,1:3)
              dr2 = sum(dr(1:3)*dr(1:3))

              m%dat(ia, ity, idx) = m%dat(ia, ity, idx) + dr2
              m%num_samples(ia, ity, idx) = m%num_samples(ia, ity, idx) + 1
           enddo

        endif

     enddo

  end subroutine

  subroutine msd_save(m)
     type(msd_type),intent(in out) :: m
     integer :: i,j,t,iunit,ierr, num_msdsize, num_types
 
     real(8),allocatable :: a(:)
     integer,allocatable :: b(:)
 
     if (.not. m%is_msd) return
 
     open(newunit=iunit,file='msd.dat',form='formatted')
     
     num_types = size(m%dat,dim=2)
     num_msdsize = size(m%dat,dim=3)
 
     a = reshape(m%dat, [m%num_init_pos*num_types*num_msdsize]) 
     b = reshape(m%num_samples, [m%num_init_pos*num_types*num_msdsize])
 
     call MPI_ALLREDUCE(MPI_IN_PLACE, a, size(a), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE, b, size(b), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
 
     m%dat = reshape(a, [m%num_init_pos, num_types, num_msdsize])
     m%num_samples = reshape(b, [m%num_init_pos, num_types, num_msdsize])
 
     ! sum up all init positions and store in the first element of dat&num_samples
     do i=1, num_types
     do j=1, num_msdsize
        m%dat(1,i,j) = sum(m%dat(1:m%num_init_pos,i,j))
        m%num_samples(1,i,j) = sum(m%num_samples(1:m%num_init_pos,i,j))
     enddo; enddo
 
     if(myid==0) then

        write(iunit,fmt='(a)',advance='no') ' time(fs) '
        do i=1, size(m%msd_atom_name)
           write(iunit,fmt='(a15)',advance='no') m%msd_atom_name(i)%str
        enddo
        do i=1, size(m%msd_atom_name)
           write(iunit,fmt='(a9)',advance='no') 'Num_'//m%msd_atom_name(i)%str
        enddo
        write(iunit,fmt=*)
    
        do t=1, num_msdsize
   
           write(iunit,fmt='(f10.2)',advance='no') (t-1)*m%freq_fs
   
           do i=1, num_types
              if(m%num_samples(1,i,t)>0) then
                 write(iunit,fmt='(es18.5)',advance='no') &
                       m%dat(1,i,t)/m%num_samples(1,i,t)
              else
                 write(iunit,fmt='(es18.5)',advance='no') 0.d0
              endif
           enddo
   
           do i=1, num_types
              if(m%num_samples(1,i,t)>0) then
                 write(iunit,fmt='(i12)',advance='no') m%num_samples(1,i,t)
              else
                 write(iunit,fmt='(i12)',advance='no') 0
              endif
           enddo
   
           write(iunit,fmt=*)
   
        enddo

        close(iunit)

     endif


  end subroutine 

end module
