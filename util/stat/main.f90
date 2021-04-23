!-----------------------------------------------------------------------------------------
program main
  use voxels_mod
!-----------------------------------------------------------------------------------------
  implicit none
  type(string_array),allocatable :: filenames(:)
  type(mdframe) :: oneframe
  type(mdframe),allocatable :: mdframes(:)
  type(analysis_context) :: ac
  type(voxel),allocatable :: voxels(:,:,:)
  real(8) :: cell_size(3)
  integer :: num_cells(3)

  integer :: n, i,j,k,l, i1,j1,k1, ity,jty,kty, idx, idx2, ii(3), jj(3)
  real(8) :: rr(0:3), rr1(0:3), rr2(0:3), cosine, theta
  character(len=256) :: argv
  character(len=:),allocatable :: dirpath
  logical :: is_there

  type(nbrlist_type),allocatable :: nbrlist(:)

  ! glob all files
  call getarg(1,argv)
  dirpath = trim(adjustl(argv))
  print'(a)', 'input data path: '//dirpath

  !inquire(directory=dirpath, exist=is_there)  ! for ifort
  !inquire(file=dirpath, exist=is_there)  ! for gfortran
  !if(.not. is_there) stop 'ERROR: input dir does not exist.'

  filenames = get_filelist_xyz(dirpath)

  ! pick one frame to setup parameters
  oneframe = get_mdframe_from_xyzfile(filenames(1)%str)
  ac = get_analysis_context_from_mdframe(oneframe)

  call get_voxel_params_from_mdframe(oneframe, num_cells, cell_size)

  do n=1, size(filenames)
     print*,'filename: ', filenames(n)%str

     ! increment the number of sampled mdframes
     ac%num_sample_frames = ac%num_sample_frames + 1

     oneframe = get_mdframe_from_xyzfile(filenames(n)%str, ac=ac)
     nbrlist = nbrlist_ctor_from_mdframe(oneframe)

     voxels = get_voxels_from_mdcontext(oneframe, num_cells, ac%elems)
    
     ! reset neighbor list
     do i=DIMS0(1), DIMS0(4)
     do j=DIMS0(2), DIMS0(5)
     do k=DIMS0(3), DIMS0(6)

        ii(1:3) = [i,j,k] + 1 + OFFSET

        do i1=-1, 1
        do j1=-1, 1
        do k1=-1, 1
           jj(1:3) = ii(1:3) + [i1,j1,k1]

           !print'(6i4,3x,3i4)',i,j,k, i1,j1,k1, i+i1,j+j1,k+k1
           call compute_interatomic_distance(nbrlist,  & 
                                             voxels(ii(1),ii(2),ii(3)), &
                                             voxels(jj(1),jj(2),jj(3)))
        enddo; enddo; enddo
     enddo; enddo; enddo

     do i = 1, size(nbrlist)

        associate(ia=>nbrlist(i)) 
        ity = ia%itype
        do j = 1, ia%counter
          jty = ia%nbrs(j)%itype
          idx = ia%nbrs(j)%ir
          ac%gr(ity,jty,idx)=ac%gr(ity,jty,idx)+1.d0
        enddo

        do j = 1, ia%counter-1
          jty = ia%nbrs(j)%itype

          rr1(1:3) = ia%nbrs(j)%pos(1:3)-ia%pos(1:3)
          rr1(0) = sqrt(sum(rr1(1:3)*rr1(1:3)))
          rr1(1:3) = rr1(1:3)/rr1(0)

          if(rr1(0)>ac%ba_rc(ity,jty)) cycle ! bond angle cutoff

          do k = j+1, ia%counter
             kty = ia%nbrs(k)%itype

             rr2(1:3) = ia%nbrs(k)%pos(1:3)-ia%pos(1:3)
             rr2(0) = sqrt(sum(rr2(1:3)*rr2(1:3)))
             rr2(1:3) = rr2(1:3)/rr2(0)

             if(rr2(0)>ac%ba_rc(ity,kty)) cycle ! bond angle cutoff
             !print'(3i3,2f8.3,2es10.2)',ity,jty,kty,rr1(0),rr2(0), ac%ba_rc(ity,jty), ac%ba_rc(ity,kty)
        
             cosine = sum(rr1(1:3)*rr2(1:3))
             theta = acos(cosine)*180d0/pi
             idx = int(theta) + 1

             ac%ba(jty,ity,kty,idx) = ac%ba(jty,ity,kty,idx)+0.5d0
             ac%ba(kty,ity,jty,idx) = ac%ba(kty,ity,jty,idx)+0.5d0

           enddo

        enddo
        end associate
     enddo

     deallocate(nbrlist)

  enddo

  call ac%save_stat()

end program
