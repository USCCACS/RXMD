!-----------------------------------------------------------------------------------------
program main
  use voxels_mod
!-----------------------------------------------------------------------------------------
  implicit none
  type(string_array),allocatable :: filenames(:)
  type(mdframe) :: mdf 
  type(mdframe),allocatable :: mdframes(:)
  type(analysis_context) :: ac
  type(voxel),allocatable :: voxels(:,:,:)
  real(8) :: cell_size(3)
  integer :: num_cells(3)

  integer :: n, i,j,k,l, ix,iy,iz, i1,j1,k1, ity,jty,kty, ir, idx, idx2, ii(3), jj(3)
  real(8) :: rr(0:3), rr1(0:3), rr2(0:3), cosine, theta
  character(len=256) :: argv
  character(len=:),allocatable :: dirpath
  logical :: is_there

  type(nbrlist_type),allocatable :: nbrlist(:)

  type(base_atom_type) :: atom

  ! glob all files
  call getarg(1,argv)
  dirpath = trim(adjustl(argv))
  print'(a)', 'input data path: '//dirpath

  !inquire(directory=dirpath, exist=is_there)  ! for ifort
  !inquire(file=dirpath, exist=is_there)  ! for gfortran
  !if(.not. is_there) stop 'ERROR: input dir does not exist.'

  filenames = get_filelist_xyz(dirpath)

  ! pick one frame to setup parameters
  mdf = get_mdframe_from_xyzfile(filenames(1)%str)
  ac = get_analysis_context_from_mdframe(mdf)

  do n=1, size(filenames)
     print*,'filename: ', filenames(n)%str

     ! increment the number of sampled mdframes
     ac%num_sample_frames = ac%num_sample_frames + 1

     mdf = get_mdframe_from_xyzfile(filenames(n)%str, ac=ac)
     nbrlist = nbrlist_ctor_from_mdframe(mdf)

     do i1=1, mdf%num_atoms
        ity = mdf%itype(i1)
        rr1(1:3) = mdf%pos(i1,1:3)
        do j1=1, mdf%num_atoms
           jty = mdf%itype(j1)
           do ix = -1, 1
           do iy = -1, 1
           do iz = -1, 1
              rr2(1:3) = mdf%pos(j1,1:3)+ix*mdf%hh(1:3,1)+iy*mdf%hh(1:3,2)+iz*mdf%hh(1:3,3)
              rr(1:3) = rr2(1:3)-rr1(1:3)
              rr(0)=sqrt(sum(rr(1:3)*rr(1:3)))

              if(abs(rr(0))<1d-3) cycle ! self-exclusion

              ir = rr(0)*DRI + 1

              if(rr(0) < RCUT) then
                 !ac%gr(ity,jty,ir)=ac%gr(ity,jty,ir)+1.d0

                 atom = base_atom_type(pos=[rr2(1),rr2(2),rr2(3)], &
                                       itype=jty, rr=rr(0), ir=ir, id=j1)

                 nbrlist(i1)%counter = nbrlist(i1)%counter + 1
                 nbrlist(i1)%nbrs(nbrlist(i1)%counter) = atom
              endif

           enddo; enddo; enddo
     enddo; enddo
    

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
