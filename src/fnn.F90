module fnn
  use iso_fortran_env, only: int32, int64, real32, real64, real128

  implicit none

  !integer,parameter :: rk = real64
  !integer,parameter :: rk = real128
  integer,parameter :: rk = real32

  !integer, parameter :: ik = int64
  integer, parameter :: ik = int32

  integer(ik),allocatable :: num_dims(:)

  real(rk),allocatable :: features(:,:) 

  integer(rk),parameter :: num_Mu = 9
  real(rk),parameter :: ml_Mu(num_Mu) = [1.0,2.0,2.86,4.06,4.96,5.74,6.42,7.02,7.58]

  integer(rk),parameter :: num_Eta = 3
  !real(rk),parameter :: ml_Eta(num_Eta) = [0.5,1.0,3.0,20.0]
  real(rk),parameter :: ml_Eta(num_Eta) = [0.5,1.0,3.0]

  !real(rk),parameter :: LJ_factor = 3.405
  !real(rk),parameter :: ml_Rc = 2.550*LJ_factor
  real(rk),parameter :: LJ_factor = 1.0
  real(rk),parameter :: ml_Rc = 8.0

  !integer(rk),parameter :: num_forcecomps = 3
  integer(rk),parameter :: num_forcecomps = 1
  integer(rk),parameter :: num_features = num_Mu*num_Eta*num_forcecomps

  type :: layer
    real(rk),allocatable :: b(:)
    real(rk),allocatable :: w(:,:)
  end type

  type :: network
    type(layer),allocatable :: layers(:)
  end type

  real(rk),allocatable :: infs(:,:)

contains

!------------------------------------------------------------------------------
subroutine get_feedforward_network(path) 
!------------------------------------------------------------------------------
   character(*),intent(in) :: path

   character(1),parameter :: cxyz(3) = ['x','y','z']
   type(network) :: netx,nety,netz

   netx = network_ctor(num_dims,path//'/'//cxyz(1))
   nety = network_ctor(num_dims,path//'/'//cxyz(2))
   netz = network_ctor(num_dims,path//'/'//cxyz(3))
    
end subroutine

!------------------------------------------------------------------------------
type(network) function network_ctor(dims, netdata_prefix) result(net)
!------------------------------------------------------------------------------
   implicit none

   integer(ik),intent(in) :: dims(:)
   character(*),intent(in) :: netdata_prefix
   character(3) :: arow, acol, alayer

   integer(ik) :: i, nrow, ncol, fileunit
   integer(ik) :: num_layers

   num_layers = size(dims)

   allocate(net%layers(num_layers))
   do i=1, num_layers-1
     nrow = dims(i)
     ncol = dims(i+1)

     allocate(net%layers(i)%b(ncol))
     allocate(net%layers(i)%w(nrow,ncol))
     print*,'i,nrow,ncol: ', i,nrow,ncol 
     write(alayer,'(i3)') i; alayer = adjustl(alayer)
     write(arow,'(i3)') nrow; arow = adjustl(arow)
     write(acol,'(i3)') ncol; acol = adjustl(acol)
     print*,'b: ', trim(netdata_prefix)//'_b_'//trim(alayer)//'_'//trim(acol)//'.bin', size(net%layers(i)%b)
     open(newunit=fileunit, access='stream', form='unformatted', status='old', &
          file=trim(netdata_prefix)//'_b_'//trim(alayer)//'_'//trim(acol)//'.bin')
     read(fileunit) net%layers(i)%b
     close(fileunit)

     print*,'w: ', trim(netdata_prefix)//'_w_'//trim(alayer)//'_'//trim(arow)//'_'//trim(acol)//'.bin'
     open(newunit=fileunit, access='stream', form='unformatted', status='old', &
          file=trim(netdata_prefix)//'_w_'//trim(alayer)//'_'//trim(arow)//'_'//trim(acol)//'.bin')
     read(fileunit) net%layers(i)%w
     close(fileunit)
   enddo

end function

!!------------------------------------------------------------------------------
!  function mynetwork_ctor(dims) result(netx)
!    implicit none
!!------------------------------------------------------------------------------
!    integer,intent(in) :: dims(:)
!    type(mynetwork) :: netx,nety,netz
!
!    integer(ik) :: i, j, k
!    integer(ik) :: ncol, nrow, num_layers=0, target_atom=1
!    integer(ik),parameter :: num_checks = 20, num_features = 27
!    integer(ik),parameter :: num_features3 = num_features*3
!    real(rk),allocatable :: x(:),y(:)
!    character(len=1) :: a1
!    real(rk) :: f(3)
!   
!    num_layers = size(dims)
!    print'(a,9i5)','number of layers from dims: ', dims
!
!!    !--- load features
!!    allocate(infs(num_checks, num_checks*num_features3))
!!    open(1,file='Al_data/input_feature.npy.txt')
!!    do i=1, num_checks
!!       read(1,*) infs(i, 1:num_features3)     
!!       !print*,infs(i, 1:num_features3)     
!!    enddo
!!    close(1)
!
!    !--- load network for x-force
!    netx = network_ctor(dims)
!
!    do i=1, num_layers-1
!      nrow = dims(i)
!      ncol = dims(i+1)
!
!      write(a1,'(i1)') i
!
!      open(1,file='Al_data/Al_b_x_'//a1//'.txt')
!      do j=1, ncol
!         read(1,*) netx%layers(i)%b(j)
!      enddo
!      close(1)
!
!      open(1,file='Al_data/Al_w_x_'//a1//'.txt')
!      do j=1, nrow 
!         read(1,*) netx%layers(i)%w(j,1:ncol)
!      enddo
!      close(1)
!    enddo
!
!    !--- load network for y-force
!    allocate(nety%layers(num_layers))
!    do i=1, num_layers-1
!      nrow = dims(i);  ncol = dims(i+1)
!
!      write(a1,'(i1)') i
!
!      allocate(nety%layers(i)%b(ncol))
!      open(1,file='Al_data/Al_b_y_'//a1//'.txt')
!      do j=1, ncol; read(1,*) nety%layers(i)%b(j);  enddo
!      close(1)
!
!      allocate(nety%layers(i)%w(nrow,ncol))
!      open(1,file='Al_data/Al_w_y_'//a1//'.txt')
!      do j=1, nrow;  read(1,*) nety%layers(i)%w(j,1:ncol);  enddo
!      close(1)
!    enddo
!
!    !--- load network for z-force
!    allocate(netz%layers(num_layers))
!    do i=1, num_layers-1
!      nrow = dims(i);  ncol = dims(i+1)
!
!      write(a1,'(i1)') i
!
!      allocate(netz%layers(i)%b(ncol))
!      open(1,file='Al_data/Al_b_z_'//a1//'.txt')
!      do j=1, ncol; read(1,*) netz%layers(i)%b(j);  enddo
!      close(1)
!
!      allocate(netz%layers(i)%w(nrow,ncol))
!      open(1,file='Al_data/Al_w_z_'//a1//'.txt')
!      do j=1, nrow;  read(1,*) netz%layers(i)%w(j,1:ncol);  enddo
!      close(1)
!    enddo
!
!    !--- initialize y with feature
!    do target_atom=1, 20
!
!       !--- z-direction force
!       if(allocated(y)) deallocate(y)
!       allocate(y(num_features))
!       y(1:num_features) = infs(target_atom,1:num_features)    
!   
!       !--- weight multiply plus bias 
!       do i=1, num_layers-1
!         nrow = dims(i);  ncol = dims(i+1)
!         !-- w*x + b
!         if(allocated(x)) deallocate(x);  allocate(x(ncol))
!         do k=1,ncol
!            x(k) = netx%layers(i)%b(k)+ & 
!                 sum(netx%layers(i)%w(1:nrow,k)*y(1:nrow)) 
!         enddo
!
!         !--- save prediction
!         f(1)=x(1)
!
!         !--- apply relu and update y
!         if(allocated(y)) deallocate(y); allocate(y(ncol))
!         y=max(x,0.0) 
!   
!       enddo
!   
!       !--- y-direction force
!       if(allocated(y)) deallocate(y); allocate(y(num_features))
!       y(1:num_features) = infs(target_atom,num_features+1:2*num_features)    
!   
!       !--- weight multiply plus bias 
!       do i=1, num_layers-1
!         nrow = dims(i);  ncol = dims(i+1)
!         !-- w*x + b
!         if(allocated(x)) deallocate(x);  allocate(x(ncol))
!         do k=1,ncol
!            x(k) = nety%layers(i)%b(k)+ & 
!                 sum(nety%layers(i)%w(1:nrow,k)*y(1:nrow)) 
!         enddo
!   
!         !--- save prediction
!         f(2)=x(1)
!
!         if(allocated(y)) deallocate(y); allocate(y(ncol))
!         y=max(x,0.0) 
!   
!       enddo
!   
!       !--- z-direction force
!       if(allocated(y)) deallocate(y); allocate(y(num_features))
!       y(1:num_features) = infs(target_atom,2*num_features+1:3*num_features)    
!   
!       do i=1, num_layers-1
!         nrow = dims(i);  ncol = dims(i+1)
!         !-- w*x + b
!         if(allocated(x)) deallocate(x);  allocate(x(ncol))
!         do k=1,ncol
!            x(k) = netz%layers(i)%b(k)+ & 
!                 sum(netz%layers(i)%w(1:nrow,k)*y(1:nrow)) 
!         enddo
!   
!         !--- save prediction
!         f(3)=x(1)
!
!         !--- apply relu and update y
!         if(allocated(y)) deallocate(y); allocate(y(ncol))
!         y=max(x,0.0)
!   
!       enddo
!   
!       print'(a,i6,3f15.8)' ,'f: ', target_atom, f(1:3)/50d0
!
!    enddo
!
!    stop 'foo'
!
!    return
!  end function
!
!  if(find_cmdline_argc('--network_dimensions',idx).or.find_cmdline_argc('-ndims',idx)) then
!      call get_command_argument(idx+1,argv)
!      read(argv, fmt=*) ii
!
!      num_dims = [num_features]
!      do ia=1, ii
!         call get_command_argument(idx+1+ia,argv)
!         read(argv, fmt=*) ib
!         num_dims = [num_dims, ib]
!      enddo
!      num_dims = [num_dims, num_forcecomps]
!  endif
!
!!  print'(/a)', repeat('=',60)
!!  print'(a30,2i6)', 'num_Mu, num_Eta : ', num_Mu, num_Eta
!!  print'(a30,i9)', 'num epochs : ', o%num_epochs
!!  print'(a30,a20,i5)', 'training datafile&size : ', trim(o%training_datafile), o%training_size
!!  print'(a30,a20,i5)', 'test datafile&size : ',  trim(o%test_datafile), o%test_size
!!  print'(a30,10i5)', 'network dimensions : ', o%num_dims
!!  print'(a30,f8.3)', 'learning rate : ', o%learning_rate
!!  print'(a30,a20)', 'network checkpoint file : ', trim(o%network_ckptfile)
!!  print'(a/)', repeat('=',60)
!
!  return
!end function

!!------------------------------------------------------------------------------
!function single_frame_ctor(num_atoms) result(o)
!  implicit none
!!------------------------------------------------------------------------------
!  integer(ik),intent(in) :: num_atoms
!  type(single_frame) :: o
!
!  allocate(o%aname(num_atoms), o%atype(num_atoms), o%q(num_atoms))
!  !allocate(o%pos(num_atoms,3), o%f(num_atoms,3),o%v(num_atoms,3))
!  allocate(o%pos(num_atoms,3), o%f(num_atoms,num_forcecomps),o%v(num_atoms,3))
!
!  o%atype = 0.0;  o%q = 0.0;  o%pos = 0.0;  o%f = 0.0;  o%v = 0.0
!
!  return
!end function

!!------------------------------------------------------------------------------
!function ml_load_data(dataset_path, num_frames) result(dataset)
!!------------------------------------------------------------------------------
!implicit none
!
!character(MAXSTRLENGTH), intent(in) :: dataset_path
!integer(ik),intent(in) :: num_frames 
!
!type(single_frame), allocatable :: dataset(:)
!integer(ik) :: fileunit, stat
!integer(ik) :: i, ia, num_atoms
!real(rk) :: lattice(6)
!
!allocate(dataset(num_frames))
!
!open(newunit=fileunit, file=dataset_path)
!
!do i=1, size(dataset) 
!   read(fileunit, fmt=*, iostat=stat) num_atoms
!   read(fileunit, fmt=*, iostat=stat) lattice(1:3)
!
!   lattice(1:3) = lattice(1:3)*LJ_factor
!   lattice(4:6) = 90.d0
!
!   dataset(i) = single_frame_ctor(num_atoms)
!
!   dataset(i)%num_atoms = num_atoms
!   dataset(i)%lattice = lattice
!
!   do ia=1, num_atoms
!      read(fileunit, fmt=*) dataset(i)%aname(ia), dataset(i)%pos(ia,1:3), dataset(i)%f(ia,1:num_forcecomps)
!      !write(6, fmt=*) dataset(i)%aname(ia), dataset(i)%pos(ia,1:3), dataset(i)%f(ia,1:3) 
!      dataset(i)%pos(ia,1:3) = dataset(i)%pos(ia,1:3)*LJ_factor
!
!      dataset(i)%f(ia,:) = dataset(i)%f(ia,:)*50
!   enddo
!
!enddo 
!   
!end function
!
!!------------------------------------------------------------------------------
!subroutine ml_initialize()
!!------------------------------------------------------------------------------
!  implicit none
!
!  integer(ik) :: iepoch, idata
!  real(rk) :: loss, cum_loss = 0.0, cum_loss0 = 0.0
!
!  type(network_type) :: network
!  type(cmdargs) :: cmd 
!
!  integer :: i, j
!
!  integer(ik) :: fileunit, n, num_layers
!  integer(ik), allocatable :: dims(:)
!
!  type(mynetwork) :: mynet
!
!  cmd = cmdargs_ctor()
!
!  mynet = mynetwork_ctor(cmd%num_dims)
!
!  dataset = ml_load_data(cmd%training_datafile, num_frames=cmd%training_size)
!
!  network = network_type(cmd%num_dims, 'relu')
!  !network = network_type(cmd%num_dims)
!
!  if (trim(cmd%network_ckptfile) /= 'n/a') & 
!     call network%load(cmd%network_ckptfile)
!
!  do iepoch = 1, cmd%num_epochs
!
!     cum_loss = 0.0
!     do idata = 1, size(dataset)
!    
!       associate(d => dataset(idata)) 
!
!         features = ml_get_features(single_frame_data = d) 
!  
!         loss = ml_train_model(network, num_atoms = d%num_atoms,  &
!                            x = features, y = d%f, &
!                            learning_rate = cmd%learning_rate) 
!       end associate
!
!       cum_loss = cum_loss + loss
!
!     enddo
!
!     if(mod(iepoch,10)==0) then
!        call network%save('network.ckpt')
!        write(*, fmt='(a,i6,es)') 'iepoch, loss: ', iepoch, cum_loss/size(dataset)
!     endif
!
!  enddo
!
!  associate(d => dataset(1))
!     features = ml_get_features(single_frame_data = d)
!     call ml_finalcheck(network, num_atoms = d%num_atoms, x = features, y = d%f)
!  end associate
!
!  return
!end subroutine
!
!!------------------------------------------------------------------------------
!function ml_get_features(single_frame_data) result(features)
!use atoms, only: pi
!!------------------------------------------------------------------------------
!implicit none
!type(single_frame),intent(in),optional :: single_frame_data
!real(rk),allocatable :: features(:,:)
!
!real(rk) :: rr(3), rr2, rij, dr_norm, dsum, eta, fr, rij_mu
!integer(ik) :: i, j, j1, l1, l2, ii, idx, ia
!
!!--- reset feature vector
!allocate(features(single_frame_data%num_atoms, num_features))
!features=0.d0
!
!!--- associate block to rename single_frame_data
!associate(sfd => single_frame_data)
!
!  do i=1, sfd%num_atoms
!  
!     do j=1, sfd%num_atoms
!  
!        if (i==j) cycle
!  
!        rr(1:3) = sfd%pos(i,1:3) - sfd%pos(j,1:3)
!        do ia=1,3
!           if(rr(ia)>=sfd%lattice(ia)*0.5d0) rr(ia)=rr(ia)-sfd%lattice(ia)
!           if(rr(ia)<-sfd%lattice(ia)*0.5d0) rr(ia)=rr(ia)+sfd%lattice(ia)
!        enddo
!  
!        rr2 = sum(rr(1:3)*rr(1:3))
!        rij = sqrt(rr2)
!  
!        if(rij < ml_Rc) then
!
!           rr(1:3) = rr(1:3)/rij
!
!           dr_norm = rij/ml_Rc
!           fr = 0.5*( 1.0 + cos(pi*dr_norm) ) 
!
!           do l1 = 1, num_Mu
!
!              rij_mu = rij - ml_Mu(l1)
!
!              do l2 = 1, num_Eta
!  
!                 eta = exp( -ml_Eta(l2) * rij_mu * rij_mu )
!  
!                 idx = (l1-1)*num_Eta*num_forcecomps + (l2-1)*num_forcecomps
!
!                 !features(i,idx+1:idx+3) = features(i,idx+1:idx+3) + rr(1:3)*eta*fr
!                 features(i,idx+1) = features(i,idx+1) + rr(1)*eta*fr
!              enddo
!           enddo
!
!        endif
!
!     enddo
!
!  enddo
!
!  !do i=1, sfd%num_atoms
!  !   print'(i6,30es13.5)',i,features(i,:)
!  !enddo
!  !stop 'foo'
!
!end associate
!
!return
!end function

end module
