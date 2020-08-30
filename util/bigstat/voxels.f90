module voxels_mod

  use stat_mod

  use mpi

  implicit none

  !include 'mpif.h'

  integer :: myrank=0, num_ranks=1, ierr=0

  type voxel
    character(len=:),allocatable :: filename
    type(string_array),allocatable :: elems(:)
    integer :: num_atoms
    real(8),allocatable :: x(:),y(:),z(:),q(:)
    integer,allocatable :: id(:), itype(:)
  end type

!--- global vars to keep negative indices of type(voxels)
  integer,parameter :: OFFSET=1
  integer :: DIMS0(6), DIMS(6)

contains

!-----------------------------------------------------------------------------------------
  subroutine get_voxel_params_from_mdframe(oneframe, num_cells, cell_size, vox_rc)
!-----------------------------------------------------------------------------------------

    type(mdframe),intent(in) :: oneframe
    integer,intent(in out) :: num_cells(3)
    real(8),intent(in out) :: cell_size(3)
    real(8),intent(in) :: vox_rc

    num_cells(1:3)=int(oneframe%lattice(1:3)/vox_rc)
    cell_size(1:3)=oneframe%lattice(1:3)/num_cells(1:3)
    DIMS0(1:3)=0;      DIMS0(4:6)=num_cells(1:3)-1
    DIMS(1:3)=-OFFSET; DIMS(4:6)=DIMS0(4:6)+OFFSET
  
    print'(a)', repeat('-',60)
    print'(a,3i6)', 'num_cells: ', num_cells
    print'(a,3f8.3)', 'cell_size: ', cell_size
    print'(a,3f8.3)', 'lattice: ', oneframe%lattice(1:3)
    print'(a,6i6)', 'DIMS0: ', DIMS0
    print'(a,6i6)', 'DIMS:  ', DIMS
    print'(a)', repeat('-',60)

  end subroutine

!-----------------------------------------------------------------------------------------
  subroutine save_voxels_into_xyz(v, lookup_name, filename)
!-----------------------------------------------------------------------------------------
    type(voxel),intent(in out) :: v(DIMS(1):DIMS(4),DIMS(2):DIMS(5),DIMS(3):DIMS(6))
    type(string_array),allocatable :: lookup_name(:)
    character(*),intent(in) :: filename
    integer :: iunit, i,j,k,n,ity
   
    open(newunit=iunit,file=trim(adjustl(filename)), form='formatted')
    write(iunit,fmt=*) sum(v(:,:,:)%num_atoms)
    write(iunit,fmt=*)

    do i=DIMS(1),DIMS(4)
    do j=DIMS(2),DIMS(5)
    do k=DIMS(3),DIMS(6)
       do n=1, v(i,j,k)%num_atoms 
          ity = v(i,j,k)%itype(n)
          !print'(4i6,3x,a)',i,j,k,ity,lookup_name(ity)%str
          write(iunit,fmt=*) lookup_name(ity)%str, v(i,j,k)%x(n), v(i,j,k)%y(n), v(i,j,k)%z(n)
       enddo
    enddo; enddo; enddo

    close(iunit) 

  end subroutine

!-----------------------------------------------------------------------------------------
  subroutine cache_boundary_voxels(v, lattice)
!-----------------------------------------------------------------------------------------
    type(voxel),intent(in out) :: v(DIMS(1):DIMS(4),DIMS(2):DIMS(5),DIMS(3):DIMS(6))
    real(8),intent(in) :: lattice(6)
    integer :: i,j,k

    associate(d1=>DIMS0(1),d2=>DIMS0(2),d3=>DIMS0(3),d4=>DIMS0(4),d5=>DIMS0(5),d6=>DIMS0(6))

      do j=d2,d5
      do k=d3,d6
         v(DIMS(4),j,k)=v(0,j,k)
         v(DIMS(1),j,k)=v(d4,j,k)

         v(DIMS(4),j,k)%x=v(0,j,k)%x + lattice(1)
         v(DIMS(1),j,k)%x=v(d4,j,k)%x - lattice(1)
      enddo; enddo
  
      do i=DIMS(1),DIMS(4)
      do k=d3,d6
         v(i,DIMS(5),k)=v(i,0,k)
         v(i,DIMS(2),k)=v(i,d5,k)

         v(i,DIMS(5),k)%y=v(i,0,k)%y + lattice(2)
         v(i,DIMS(2),k)%y=v(i,d5,k)%y - lattice(2)
      enddo; enddo
  
      do i=DIMS(1),DIMS(4)
      do j=DIMS(2),DIMS(5)
         v(i,j,DIMS(6))=v(i,j,0)
         v(i,j,DIMS(3))=v(i,j,d6)

         v(i,j,DIMS(6))%z=v(i,j,0)%z + lattice(3)
         v(i,j,DIMS(3))%z=v(i,j,d6)%z - lattice(3)
      enddo; enddo

    end associate

  end subroutine

!-----------------------------------------------------------------------------------------
  subroutine compute_interatomic_distance_ON2(ac, mdf, myrank, num_ranks)
!-----------------------------------------------------------------------------------------
    type(analysis_context),intent(in out) :: ac
    type(mdframe),intent(in) :: mdf
    integer,intent(in) :: myrank, num_ranks

    integer :: i,j,k,n, it,jt, ia
    real(8) :: dr(3),r
    integer :: ir
    integer :: tid=0, num_threads=1, omp_get_thread_num, omp_get_num_threads

    call omp_set_num_threads(16)

!$omp parallel private(tid,i,j,dr,r,ir,it,jt)
    tid = omp_get_thread_num()
    num_threads = omp_get_num_threads()

    do i=1+tid, mdf%num_atoms, num_threads
       !if(myrank==0) print'(10i9)',myrank,i,mdf%num_atoms,tid,num_threads

       do j=1+myrank, mdf%num_atoms, num_ranks

          do ia=1,3
             dr(ia) = mdf%pos(j,ia) - mdf%pos(i,ia)

             if(dr(ia)>=0.5d0*mdf%lattice(ia)) dr(ia)=dr(ia)-mdf%lattice(ia)
             if(dr(ia)<-0.5d0*mdf%lattice(ia)) dr(ia)=dr(ia)+mdf%lattice(ia)
          enddo

          r = sum(dr*dr) 
          r = sqrt(r)
          ir = r*ac%box_dri + 1

          if(0.d0 < r .and. r < ac%box_rc) then
            it = mdf%itype(i)
            jt = mdf%itype(j)
!$omp atomic
            ac%gr(it,jt,ir)=ac%gr(it,jt,ir)+1.d0
          end if

       enddo

    enddo
!$omp end parallel

  end subroutine

!-----------------------------------------------------------------------------------------
  subroutine compute_interatomic_distance(nbrlist, v1, v2)
!-----------------------------------------------------------------------------------------
    type(nbrlist_type),allocatable :: nbrlist(:)
    type(voxel),intent(in) :: v1, v2
    integer :: i,j,k,n, it,jt, id
    real(8) :: rr(3)
    real(8),allocatable :: x(:),y(:),z(:),r(:)
    integer,allocatable :: ir(:)
    logical,allocatable :: l2b(:),l3b(:)

    type(base_atom_type) :: atom
    real(8) :: pos(3)

    n = v2%num_atoms
    allocate(l2b(n))

    do i=1, v1%num_atoms
       x = v2%x - v1%x(i)
       y = v2%y - v1%y(i)
       z = v2%z - v1%z(i)
       r = x*x + y*y + z*z
       r = sqrt(r)
       ir = r*DRI+1

       !print*,repeat('-',60)

       where(0.d0 < r .and. r < rcut) 
          l2b = .true.
       elsewhere
          l2b = .false.
       end where

       id = v1%id(i)
       do j=1, size(l2b)

          if(l2b(j)) then
             !atom = base_atom_type_ctor( pos=[v2%x(j),v2%y(j),v2%z(j)], &
             !                       itype=v2%itype(j), rr=r(j), ir=ir(j), id=v2%id(j) )
             atom = base_atom_type_ctor( [v2%x(j),v2%y(j),v2%z(j)], v2%itype(j), r(j), ir(j), v2%id(j) )
             nbrlist(id)%counter = nbrlist(id)%counter + 1
             nbrlist(id)%nbrs(nbrlist(id)%counter) = atom

          endif
       enddo

    enddo

  end subroutine

!-----------------------------------------------------------------------------------------
  subroutine print_voxels(voxels)
!-----------------------------------------------------------------------------------------
    type(voxel),intent(in) :: voxels(DIMS(1):DIMS(4),DIMS(2):DIMS(5),DIMS(3):DIMS(6))
    integer :: i,j,k, n

    print'(a)',repeat('-',60)
    n = 0
    do i=DIMS0(1),DIMS0(4)
    do j=DIMS0(2),DIMS0(5)
    do k=DIMS0(3),DIMS0(6)
       print'(3i3,3x,i9)',i,j,k,voxels(i,j,k)%num_atoms
       n = n + voxels(i,j,k)%num_atoms
    enddo; enddo; enddo

    print*,'total atoms: ', n
    print'(a)',repeat('-',60)

  end subroutine

!-----------------------------------------------------------------------------------------
  function get_voxels_from_mdcontext(frame, num_cells, lookup_elems) result(c)
!-----------------------------------------------------------------------------------------

    integer,intent(in) :: num_cells(:)
    type(mdframe),intent(in) :: frame
    type(string_array),allocatable,intent(in) :: lookup_elems(:)

    type(voxel),allocatable :: c(:,:,:)
    character(len=:),allocatable :: buf
 
    real(8) :: dr(3)
    integer :: i,j,k,ir(3)


    ! receiver needs to be indexed with the same range as passer
    allocate( c(DIMS(1):DIMS(4),DIMS(2):DIMS(5),DIMS(3):DIMS(6)) )

    do i=DIMS(1),DIMS(4)
    do j=DIMS(2),DIMS(5)
    do k=DIMS(3),DIMS(6)
       ! save original filename
       c(i,j,k)%filename = frame%filename

       c(i,j,k)%num_atoms=0 
       allocate(c(i,j,k)%x(0), c(i,j,k)%y(0), c(i,j,k)%z(0), c(i,j,k)%q(0), c(i,j,k)%elems(0))
       allocate(c(i,j,k)%id(0), c(i,j,k)%itype(0))
    enddo; enddo; enddo


    do i=1, frame%num_atoms
       dr(1:3) = frame%pos(i,1:3)/frame%lattice(1:3) 
       ir(1:3) = dr(1:3)*num_cells(1:3)

       !print'(i6,3i3,3f10.5)',i,ir(1:3), dr(1:3)
       associate( cc=>c(ir(1),ir(2),ir(3)) )
         cc%num_atoms = cc%num_atoms + 1
         cc%x = [cc%x,frame%pos(i,1)]
         cc%y = [cc%y,frame%pos(i,2)]
         cc%z = [cc%z,frame%pos(i,3)]
         cc%q = [cc%q,frame%q(i)]
         cc%elems = [cc%elems,frame%elems(i)]
         cc%itype = [cc%itype,get_index(lookup_elems, frame%elems(i)%str)]
         cc%id = [cc%id,i]
         !print'(i6,3i3,3x,50a)',i,ir(1:3),cc%elems
         !print'(i6,3i3,3x,50i3)',i,ir(1:3),cc%itype
       end associate

    enddo

  ! copy resident atoms to boundary voxels
    call cache_boundary_voxels(c, frame%lattice)

    !call print_voxels(c)

    return
  end function

!-----------------------------------------------------------------------------------------
  function get_mdframe_from_xyzfile(xyzpath, ac) result(c)
!-----------------------------------------------------------------------------------------
     type(mdframe) :: c
     type(analysis_context),optional,intent(in out) :: ac


     character(len=*),intent(in) :: xyzpath
     character(len=3) :: cbuf

     real(8) :: lattice(9)

     integer :: i,iunit

     ! save original filename
     c%filename = trim(adjustl(xyzpath))

     open(newunit=iunit,file=xyzpath,status='old',form='formatted')
     read(unit=iunit,fmt=*) c%num_atoms
     read(unit=iunit,fmt=*) c%lattice(1:6)

     !read(unit=iunit,fmt=*) lattice(1:9)
     !c%lattice(1)=lattice(1)
     !c%lattice(2)=lattice(5)
     !c%lattice(3)=lattice(9)
     !c%lattice(4:6)=90d0

     allocate(c%pos(c%num_atoms,3),c%v(c%num_atoms,3),c%f(c%num_atoms,3))
     allocate(c%elems(c%num_atoms),c%q(c%num_atoms), c%itype(c%num_atoms))
     
     !print'(a,i6,6f8.3)',xyzpath,c%num_atoms,c%lattice(1:6)
     do i=1, c%num_atoms
        read(unit=iunit,fmt=*) cbuf, c%pos(i,1), c%pos(i,2), c%pos(i,3)
        c%elems(i)%str=trim(adjustl(cbuf))
     enddo
     close(unit=iunit)

     if(present(ac)) then
       do i=1, c%num_atoms
          c%itype(i) = get_index(ac%elems, c%elems(i)%str)
       enddo
     endif
  end function

!-----------------------------------------------------------------------------------------
  function get_filelist_xyz(dirpath) result(c)
!-----------------------------------------------------------------------------------------
     type(string_array),allocatable :: c(:)

     character(*),intent(in) :: dirpath

     character(256) :: buf
     character(len=:),allocatable :: cmd, filename, filename_tmp

     integer :: iunit
    
     filename_tmp = dirpath//'/foo.txt'
     cmd = 'ls '//trim(adjustl(dirpath))//'/*.xyz | sort > '//filename_tmp
   
     call execute_command_line(cmd)

     ! allocate zero-element array
     allocate(c(0))

     open(newunit=iunit,file=filename_tmp,status='old',form='formatted')

     do while(.true.)
        read(unit=iunit,fmt='(a)',end=999) buf
        filename=trim(adjustl(buf))
        c = [c, string_array(str=filename)]
        !print*,filename
     enddo

     999 close(iunit, status='delete')
 
  end function

end module
