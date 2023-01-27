!------------------------------------------------------------------------------
module harmonic_potential_mod
!------------------------------------------------------------------------------
  use mpi_mod
  use base, only : myid, nprocs
  use utils, only : find_cmdline_argc, getstr, l2g, token

  implicit none

  type harmonic_potential_type

    logical :: apply = .false.
    character(len=:),allocatable :: filename 

    integer :: s_id, id_binding
    integer,allocatable :: id_except(:)

    real(8) :: s_const, s_pos(3)

  end type

  type(harmonic_potential_type) :: harmo_pot

  type feature_voxel_type
    real(8) :: boxsize, voxelsize, voxel_res(3), sigma
    integer :: nvoxels(3), nvoxel_total 
    real(8),allocatable :: feat(:,:,:), feature(:,:,:)

    contains  
       procedure :: init => feature_voxel_init_
       procedure :: get => get_feature_voxel_
  end type

  type(feature_voxel_type) :: fvox 

contains

!------------------------------------------------------------------------------
  subroutine show_message(myrank, num_ranks, message)
  !FIXME since F&P communicate through stdout, the output needs to be in order. let rank0 write. 
!------------------------------------------------------------------------------
   integer,intent(in) :: myrank, num_ranks
   character(len=:),allocatable,intent(in) :: message
   character(len=:),allocatable :: buf
   character(64) :: out_fmt
   integer :: rank, num_size, ierr, mystatus

   if (myrank==0) then

      write(out_fmt,'(i9)') len(message)
      write(*,fmt='(a'//trim(adjustl(out_fmt))//')') message

      do rank=1, num_ranks-1
         call MPI_Recv(num_size, 1, MPI_INTEGER, rank, 1112, MPI_COMM_WORLD, mystatus, ierr)
         !print*,'num_size', rank, num_size
         if (num_size > 0) then
            allocate(character(len=num_size) :: buf)
            call MPI_Recv(buf, num_size, MPI_CHARACTER, rank, 1113, MPI_COMM_WORLD, mystatus, ierr)
            write(out_fmt,'(i9)') num_size
            write(*,fmt='(a'//trim(adjustl(out_fmt))//')') buf
            deallocate(buf)
         endif
      enddo
   else
       !print*,myrank, message
       call MPI_Send(len(message), 1, MPI_INTEGER, 0, 1112, MPI_COMM_WORLD, mystatus, ierr)
       if (len(message) > 0) &
           call MPI_Send(message, len(message), MPI_CHARACTER, 0, 1113, MPI_COMM_WORLD, mystatus, ierr)
   endif

  end subroutine

!------------------------------------------------------------------------------
  subroutine feature_voxel_init_(self, boxsize_, voxelsize_, sigma_)
!------------------------------------------------------------------------------
     class(feature_voxel_type),intent(in out) :: self
     real(8),intent(in) :: boxsize_, voxelsize_, sigma_

     integer :: i1,i2,i3, idx
     real(8) :: dr(0:3)
     character(64) :: argv

     self%boxsize = boxsize_; self%voxelsize = voxelsize_; self%sigma = sigma_
     if(find_cmdline_argc('--fvox_configs',idx)) then
        call get_command_argument(idx+1,argv);  read(argv,*) self%boxsize 
        call get_command_argument(idx+2,argv);  read(argv,*) self%voxelsize 
        call get_command_argument(idx+3,argv);  read(argv,*) self%sigma 
     endif

     self%nvoxels(1:3)=self%boxsize/self%voxelsize
     self%voxel_res(1:3)=self%boxsize/self%nvoxels(1:3)
     self%nvoxel_total = (2*self%nvoxels(1)+1)*(2*self%nvoxels(2)+1)*(2*self%nvoxels(3)+1)
     allocate(self%feat(-self%nvoxels(1):self%nvoxels(1),-self%nvoxels(2):self%nvoxels(2),-self%nvoxels(3):self%nvoxels(3)))
     allocate(self%feature(-self%nvoxels(1):self%nvoxels(1),-self%nvoxels(2):self%nvoxels(2),-self%nvoxels(3):self%nvoxels(3)))

     self%feat = 0d0
     do i1=-self%nvoxels(1),self%nvoxels(1)
     do i2=-self%nvoxels(2),self%nvoxels(2)
     do i3=-self%nvoxels(3),self%nvoxels(3)
        dr(1:3) = ([i1,i2,i3]-0.5d0)*self%voxel_res(1:3)
        dr(0) = sqrt(sum(dr(1:3)*dr(1:3)))
        self%feat(i1,i2,i3) = exp(-(dr(0)/self%sigma)**2)
     enddo; enddo; enddo

  end subroutine

!------------------------------------------------------------------------------
  subroutine get_feature_voxel_(self, natoms, atype, pos, igd, nbplist) 
!------------------------------------------------------------------------------
     class(feature_voxel_type),intent(in out) :: self
     integer,intent(in) :: natoms, igd
     real(8),allocatable,intent(in) :: atype(:), pos(:,:)
     integer,allocatable,intent(in) :: nbplist(:,:) ! non-bonding list with 10A cutoff. note nbp(0,i) instead of nbp(i,0).

     character(10*self%nvoxel_total) :: strdata
     character(len=:),allocatable :: msg
     character(256) :: buf

     integer :: i,j,n,ii, ity,jty, i1,i2,i3, j1,j2,j3, k(3), amin(3),amax(3),bmin(3),bmax(3), ierr
     real(8) :: dr(0:3)

     msg = ""
     strdata = ""
     self%feature = 0d0
     do i=1, natoms
  
        if (igd == l2g(atype(i))) then
  
           self%feature = 0d0
  
           ity = nint(atype(i))
  
           do j1 = 1, nbplist(0,i)
              j = nbplist(j1,i)
              jty = nint(atype(j))
  
              if (jty == 1) cycle ! ignore H
  
              dr(1:3) = pos(j,1:3) - pos(i,1:3)
              do ii = 1, 3
                 k(ii) = int(dr(ii)/self%voxel_res(ii))
                 if (dr(ii)<0d0) k(ii)=k(ii)-1
              enddo
  
              do ii=1,3
                 amin(ii) = max(-self%nvoxels(ii)+k(ii),-self%nvoxels(ii))
                 amax(ii) = min( self%nvoxels(ii)+k(ii), self%nvoxels(ii))
                 bmin(ii) = max(-self%nvoxels(ii)-k(ii),-self%nvoxels(ii))
                 bmax(ii) = min( self%nvoxels(ii)-k(ii), self%nvoxels(ii))
              enddo
  
              self%feature(amin(1):amax(1),amin(2):amax(2),amin(3):amax(3)) = &
              self%feature(amin(1):amax(1),amin(2):amax(2),amin(3):amax(3)) + & 
                 self%feat(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3))
           enddo

           ii = 0
           do i1=-self%nvoxels(1),self%nvoxels(1)
           do i2=-self%nvoxels(2),self%nvoxels(2)
           do i3=-self%nvoxels(3),self%nvoxels(3)
             write(strdata(ii+1:ii+10),'(f9.5,a1)') self%feature(i1,i2,i3),' '
             ii = ii + 10
           enddo; enddo; enddo
           write(buf,'(a,i6,1x)') 'Feature: ', igd
           msg = msg//trim(buf)//strdata//new_line('A'); buf = ""
  
           write(buf,'(a,3i6)') 'Voxels: ', 2*self%nvoxels(1:3)+1;   
           msg = msg//trim(buf)//new_line('A'); buf = ""

           write(buf, '(a,i6,es15.5)') 'Vfree:', igd, self%feature(0,0,0)
           msg = msg//trim(buf)//new_line('A'); buf = ""

        endif
     enddo

     !FIXME since F&P communicate through stdout, the output needs to be in order. let rank0 write. 
     call show_message(myrank=myid, num_ranks=nprocs, message=msg)

    !call MPI_FINALIZE(i1)
    !stop
  end subroutine

!------------------------------------------------------------------------------
  subroutine show_harmonic_potential_params(hp)
!------------------------------------------------------------------------------
    type(harmonic_potential_type),intent(in) :: hp
    integer :: idx

    print'(a)',repeat('-',80)
    print'(a,l)','apply: ', hp%apply
    print'(a,a)','filename: ', hp%filename
    print'(a,i3,i6,3f10.5)','s_id,s_binding,s_pos: ', hp%s_id, hp%id_binding, hp%s_pos
    print'(a,10i6)','id_except: ', hp%id_except
    print'(a)',repeat('-',80)

  end subroutine

!------------------------------------------------------------------------------
  function harmonic_potential_ctor() result(harmo_pot)
!------------------------------------------------------------------------------
    integer :: idx, iunit
    character(256) :: argv, linein0
    character(len=:),allocatable :: linein 

    type(harmonic_potential_type) :: harmo_pot

    harmo_pot%filename="harmopot.in"

    if(find_cmdline_argc('--harmo_pot',idx).or.find_cmdline_argc('-harmo',idx)) then

      call get_command_argument(idx+1,argv)
      harmo_pot%filename = trim(adjustl(argv))

      open(newunit=iunit, file=harmo_pot%filename, form='formatted', status='old')

      do while (.true.)
        read(iunit,'(a)',end=10) linein0
        linein = trim(adjustl(linein0))

        if (getstr(linein, token) > 0 .and. token=='K') &
                call get_spring_params(linein, harmo_pot)
      enddo
      10 close(iunit)

    endif
  end function

!------------------------------------------------------------------------------
  subroutine get_spring_params(linein, hp)
!------------------------------------------------------------------------------
    character(len=:),allocatable,intent(in out) :: linein
    type(harmonic_potential_type),intent(in out) :: hp 

    real(8),allocatable :: K(:) ! spring constants
    integer :: i, ii

    if (getstr(linein, token) < 0) stop 'error while reading s_id'
    read(token, *) hp%s_id
    if (getstr(linein, token) < 0) stop 'error while reading s_const'
    read(token, *) hp%s_const
    do i=1, 3
       if (getstr(linein, token) < 0) stop 'error while reading s_pos'
       read(token, *) hp%s_pos(i)
    enddo
    if (getstr(linein, token) < 0) stop 'error while reading s_binding'
    read(token, *) hp%id_binding

    if (getstr(linein, token) < 0) stop 'error while reading num_except_ids'
    read(token, *) ii

    allocate(hp%id_except(ii))
    do i=1, ii
       if (getstr(linein, token) < 0) stop 'error while reading id_except'
       read(token, *) hp%id_except(i)
    enddo

    hp%apply = .true.

  end subroutine

!------------------------------------------------------------------------------
  subroutine apply_harmonic_potential(hp, natoms, atype, pos, f, nbrlist, bo, pe)
!------------------------------------------------------------------------------
    type(harmonic_potential_type) :: hp
    integer,intent(in) :: natoms
    real(8),allocatable,intent(in) :: atype(:), pos(:,:), bo(:,:,:)
    real(8),allocatable,intent(in out) :: f(:,:)
    integer,allocatable,intent(in) :: nbrlist(:,:)
    real(8),intent(in out) :: pe
    real(8) :: Eharmo

    character(len=:),allocatable :: nbr_str, msg
    character(256) :: buf

    integer :: i, j1, idx, igd, jgd, ierr
    real(8) :: dr(3), dr2, ff(3)
    real(8) :: sigma2 = 3d0**2

    msg = ""
    nbr_str = ""

    Eharmo = 0.d0
    do i=1, natoms
       igd = l2g(atype(i))

       dr(1:3) = pos(i,1:3) - hp%s_pos(1:3)
       dr2 = sum(dr(1:3)*dr(1:3))

       if (igd == hp%id_binding) then
          !print*,'igd: ', hp%s_id, dr2

          ff(1:3) = hp%s_const*dr(1:3)

          Eharmo = 0.5*hp%s_const*dr2
          f(i,1:3) = f(i,1:3) - ff(1:3)
          !print'(2i,3f)',i,igd,ff(1:3)

          nbr_str = itoa(igd)//" "
          do j1 = 1, nbrlist(i,0)
             if(bo(0,i,j1)>0.3d0) then
               jgd = l2g(atype(nbrlist(i,j1)))
               nbr_str = nbr_str//" "//itoa(jgd)
             endif
          enddo

          write(buf,'(2a)') 'Neighbors: ', nbr_str
          msg = msg//trim(buf)//new_line('A'); buf = ""

          write(buf,'(a,i6,f8.3,1x,6es15.5)') 'AtomSpringInfo: ', igd, hp%s_const, pos(i,1:3), hp%s_pos(1:3)
          msg = msg//trim(buf)//new_line('A'); buf = ""

          write(buf, '(a,es15.5)') 'Eharmo: ', Eharmo
          msg = msg//trim(buf)//new_line('A'); buf = ""

       else if(hp%s_id == 0) then ! apply bias potential
          if(.not. any(igd == hp%id_except)) then
                  !print*,'bias', any(igd == hp%id_except)
                  ff(1:3) = exp(-dr2*0.5d0/sigma2)*dr(1:3)
                  f(i,1:3) = f(i,1:3) + ff(1:3)
          endif
       endif

    enddo

     !FIXME since F&P communicate through stdout, the output needs to be in order. let rank0 write. 
    call show_message(myrank=myid, num_ranks=nprocs, message=msg)

    call MPI_ALLREDUCE(MPI_IN_PLACE, Eharmo, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)

    pe = pe + Eharmo

  contains
    function itoa(i) result(res)
      character(:),allocatable :: res
      integer,intent(in) :: i
      character(range(i)+2) :: tmp
      write(tmp,'(i0)') i
      res = trim(tmp)
    end function
  end subroutine
  
end module

!program main
!  use harmonic_potential 
!  harmo_pot = harmonic_potential_ctor() 
!  call show_harmonic_potential_params(harmo_pot)
!end program 
