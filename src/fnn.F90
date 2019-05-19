module fnnin_parser

  use utils, only : getstr, getfilenamebase
  use base, only : force_field_class
  use fileio, only : output

  use iso_fortran_env, only: int32, int64, real32, real64 

  implicit none

  !integer,parameter :: rk = real64
  integer,parameter :: rk = real32

  !integer, parameter :: ik = int64
  integer, parameter :: ik = int32

  type :: layer
    real(rk),allocatable :: b(:)
    real(rk),allocatable :: w(:,:)
  end type

  type :: network
    integer(ik),allocatable :: dims(:)
    type(layer),allocatable :: layers(:)
  end type

  type :: feature_stat
    real(rk),allocatable :: mean(:)
    real(rk),allocatable :: stddev(:)
  end type

  type :: model_params
    character(len=:),allocatable :: element
    real(rk) :: mass
    real(rk) :: scaling_factor

    integer(ik),allocatable :: hlayers(:)

    type(network),allocatable :: networks(:)
    type(feature_stat),allocatable :: fstat(:)
  end type

  type, extends(force_field_class) :: fnn_param
    real(rk),allocatable,dimension(:) :: rad_mu, rad_eta
    real(rk),allocatable,dimension(:) :: ang_mu, ang_eta, ang_zeta
    integer(ik),allocatable,dimension(:) :: ang_lambda

    real(rk) :: rad_rc, rad_damp, ang_rc, ang_damp

    type(model_params), allocatable :: models(:) 

    integer(ik) :: feature_size_rad, feature_size_ang, feature_size 

    contains 
       procedure :: print => fnn_param_print

  end type

  type(fnn_param),target :: fnn_param_obj

  interface get_tokens_and_append
      module procedure get_tokens_and_append_rv, get_tokens_and_append_iv, & 
                       get_tokens_and_append_rs, get_tokens_and_append_model
  end interface

  character(len=:),allocatable,private :: sbuf
  integer,private :: ibuf
  real(rk),private :: rbuf
  character(len=:),allocatable,private :: token

contains

  subroutine get_tokens_and_append_model(linein, models)
    character(len=:),allocatable,intent(in out) :: linein
    type(model_params),allocatable,intent(in out) :: models(:)
    type(model_params) :: mbuf

    if (getstr(linein, token) < 0) stop 'erro while reading element name'  
    mbuf%element = trim(adjustl(token))
    if (getstr(linein, token) < 0) stop 'erro while reading element mass'  
    read(token, *) mbuf%mass
    if (getstr(linein, token) < 0) stop 'erro while reading scaling factor'  
    read(token, *) mbuf%scaling_factor

    ! allocate zero-sized array
    if(.not.allocated(mbuf%hlayers)) allocate(mbuf%hlayers(0))
    do while( getstr(linein, token) > 0 )
       read(token, *) ibuf
       mbuf%hlayers = [mbuf%hlayers, ibuf]
    enddo

    ! allocate zero-sized array
    if(.not.allocated(models)) allocate(models(0)) 
    models = [models, mbuf]

    return
  end subroutine

  subroutine get_tokens_and_append_rv(linein, array)
    character(len=:),allocatable,intent(in out) :: linein
    real(4),allocatable,intent(in out) :: array(:)

    ! allocate zero-sized array
    if(.not.allocated(array)) allocate(array(0)) 

    do while (getstr(linein, token) > 0) 
       read(token,*) rbuf 
       array = [array, rbuf]
    enddo

    return
  end subroutine

  subroutine get_tokens_and_append_iv(linein, array)
    character(len=:),allocatable,intent(in out) :: linein
    integer(4),allocatable,intent(in out) :: array(:)

    ! allocate zero-sized array
    if(.not.allocated(array)) allocate(array(0)) 

    do while (getstr(linein, token) > 0) 
       read(token,*) ibuf 
       array = [array, ibuf]
    enddo

    return
  end subroutine

  subroutine get_tokens_and_append_rs(linein, scalar)
    character(len=:),allocatable,intent(in out) :: linein
    real(4),intent(in out) :: scalar

    if (getstr(linein, token) > 0) then
       read(token,*) rbuf 
       scalar = rbuf
    endif 

    return
  end subroutine

  function fnn_param_ctor(path) result(c)
    character(len=:),allocatable,intent(in) :: path
    character(256) :: linein0
    character(len=:),allocatable :: linein

    type(fnn_param) :: c 
    integer :: iunit

    open(newunit=iunit, file=path, status='old', form='formatted')

    ! allocate zero-sized array
    allocate(c%models(0), c%rad_mu(0),c%rad_eta(0))
    allocate(c%ang_mu(0),c%ang_eta(0),c%ang_zeta(0),c%ang_lambda(0))
    c%rad_rc = -1.0; c%ang_rc = -1.0; c%rad_damp = 1.0e9; c%ang_damp = 1.0e9
    c%feature_size_rad = 0; c%feature_size_ang = 0

    do while (.true.)
      read(iunit,'(a)',end=10) linein0
      linein = trim(adjustl(linein0))

      if(getstr(linein, token) > 0) then

        select case (token)
           case ('rad_eta')
             call get_tokens_and_append(linein, c%rad_eta)
           case ('rad_mu')
             call get_tokens_and_append(linein, c%rad_mu)
           case ('ang_eta')
             call get_tokens_and_append(linein, c%ang_eta)
           case ('ang_mu')
             call get_tokens_and_append(linein, c%ang_mu)
           case ('ang_lambda')
             call get_tokens_and_append(linein, c%ang_lambda)
           case ('ang_zeta')
             call get_tokens_and_append(linein, c%ang_zeta)

           case ('rad_rc')
             call get_tokens_and_append(linein, c%rad_rc)
           case ('rad_damp')
             call get_tokens_and_append(linein, c%rad_damp)
           case ('ang_rc')
             call get_tokens_and_append(linein, c%ang_rc)
           case ('ang_damp')
             call get_tokens_and_append(linein, c%ang_damp)
           case ('model')
             call get_tokens_and_append(linein, c%models)

           case default
        end select
     endif

    end do

    10 close(iunit)

    if (size(c%models)<0) stop 'ERROR: at least one model must be defined.'
    if (c%rad_rc<0) stop 'ERROR: radial feature cutoff must to be positive.'

    ! setup feature vector size. 
    c%feature_size_rad = 1
    if(size(c%rad_mu)>0) c%feature_size_rad = c%feature_size_rad * size(c%rad_mu)
    if(size(c%rad_eta)>0) c%feature_size_rad = c%feature_size_rad * size(c%rad_eta)

    ! upto here is radial feature
    if (c%ang_rc>0) then
       c%feature_size_ang = 1
       if(size(c%ang_mu)>0) c%feature_size_ang = c%feature_size_ang * size(c%ang_mu)
       if(size(c%ang_eta)>0) c%feature_size_ang = c%feature_size_ang * size(c%ang_eta)
       if(size(c%ang_lambda)>0) c%feature_size_ang = c%feature_size_ang * size(c%ang_lambda)
       if(size(c%ang_zeta)>0) c%feature_size_ang = c%feature_size_ang * size(c%ang_zeta)
    endif

    c%feature_size = c%feature_size_rad*size(c%models) + c%feature_size_ang

  end function

  subroutine fnn_param_print(this)
    class(fnn_param), intent(in) :: this
    integer :: i,j

    print'(a)',repeat('-',60)
       print'(a,2i9)', 'feature_size_rad(single & all pair):  ', &
               this%feature_size_rad, this%feature_size_rad*size(this%models)
       print'(a,i9)', 'feature_size_ang: ', this%feature_size_ang
       print'(a,i9)', 'feature_size: ', this%feature_size
    print'(a)',repeat('-',60)
       print'(a,20f6.2)', 'rad_mu: ', this%rad_mu
       print'(a,20f6.2)', 'rad_eta: ', this%rad_eta
       print'(a,20f6.2)', 'rad_rc: ', this%rad_rc
       print'(a,20f6.2)', 'rad_damp: ', this%rad_damp
    print'(a)',repeat('-',60)

    if(this%ang_rc>0.0) then
       print'(a)',repeat('-',60)
          print'(a,20f6.2)', 'ang_mu: ', this%ang_mu
          print'(a,20f6.2)', 'ang_eta: ', this%ang_eta
          print'(a,20i6)', 'ang_lambda: ', this%ang_lambda
          print'(a,20f6.2)', 'ang_zeta: ', this%ang_zeta
          print'(a,20f6.2)', 'ang_rc: ', this%ang_rc
          print'(a,20f6.2)', 'ang_damp: ', this%ang_damp
       print'(a)',repeat('-',60)
    endif

    do i = 1, size(this%models)
       associate(m => this%models(i))
          print'(a,2f6.2,10i6)', m%element, m%mass, m%scaling_factor, m%hlayers(:)
       end associate
    end do
    print'(a)',repeat('-',60)
    
  end subroutine


end module

!------------------------------------------------------------------------------
module fnn
!------------------------------------------------------------------------------

  use fnnin_parser

  use utils, only : pi, int_to_str
  use memory_allocator_mod
  use mpi_mod

  use base
  use lists_mod, only : getnonbondingmesh, linkedlist
  use communication_mod, only : copyatoms

  implicit none

  real(rk),allocatable :: features(:,:,:) 

  integer(ik),parameter :: num_forcecomps = 1
  integer(ik),parameter :: num_networks_per_atom = 3

  integer(ik) :: num_types = 0, num_pairs = 0

  integer,parameter :: NMINCELL_FNN=1

  integer,allocatable :: pair_types(:,:)

contains

!------------------------------------------------------------------------------
subroutine set_potentialtables_fnn()
!------------------------------------------------------------------------------
! TODO: implement 

end subroutine

!------------------------------------------------------------------------------
subroutine get_force_fnn(ff, natoms, atype, pos, f, q)
!------------------------------------------------------------------------------
class(force_field_class),pointer,intent(in out) :: ff
type(fnn_param),pointer :: fp => null()
integer,intent(in out) :: natoms
real(8),intent(in out),allocatable :: atype(:), pos(:,:), q(:), f(:,:)

type(network),allocatable :: networks(:)

integer(ik) :: i, j, k, n, ncol, nrow, num_layers
integer(ik) :: m1,k1,n1
integer(ik) :: na, nn, nl

real(rk),allocatable :: x(:,:),y(:,:)

! not sure if this is the best way, but binding force_field_class to fnn_parm
select type(ff); type is (fnn_param) 
   fp => ff
end select
num_layers = size(fp%models(1)%hlayers)

call COPYATOMS(imode=MODE_COPY_FNN, dr=lcsize(1:3), atype=atype, pos=pos)
call LINKEDLIST(atype, pos, lcsize, header, llist, nacell)

call get_features_fnn(natoms, atype, pos, features, fp) 


networks = fp%models(1)%networks ! for now for testing

network_loop : do nn=1, num_networks_per_atom

  y = features(nn,1:natoms,1:fp%feature_size)

  layer_loop : do nl=1, num_layers-1

#ifdef BLAS
    m1 = size(y,1); k1 = size(y,2)
    n1 = size(networks(nn)%layers(nl)%w,2)

    if(allocated(x)) deallocate(x); allocate(x,mold=y)
    do i=1,natoms
       x(i,:) = networks(nn)%layers(nl)%b(:)
    enddo

    call sgemm('n','n', m1,k1,n1, 1.0, y,m1, networks(nn)%layers(nl)%w,n1, 1.0, x,m1)
#else
    x = matmul(y,networks(nn)%layers(nl)%w)
    do i=1,natoms
       x(i,:) = networks(nn)%layers(nl)%b(:)
    enddo
#endif

    y = max(x,0.0) ! relu

  enddo layer_loop

  !--- update force
  f(1:natoms,nn) = x(1:natoms,1)

enddo network_loop

end subroutine

!------------------------------------------------------------------------------
subroutine mddriver_fnn(mdbase, num_mdsteps) 
!------------------------------------------------------------------------------
type(mdbase_class),intent(in out) :: mdbase
integer,intent(in) :: num_mdsteps
real(4) :: ke

integer :: i,ity,nstep


!--- set force model
do nstep=0, num_mdsteps-1

   if(mod(nstep,fstep)==0) &
        call OUTPUT(atype, pos, v, q, GetFileNameBase(DataDir,current_step+nstep))

   ke=0.d0
   do i=1, NATOMS
      ity=nint(atype(i))
      ke = ke + hmas(ity)*sum(v(i,1:3)*v(i,1:3))
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, ke, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
   ke = ke/GNATOMS

   call get_force_fnn(mdbase%ff, natoms, atype, pos, f, q)
   if(mod(nstep,pstep)==0) print'(a,i6,es15.5)','step: ', nstep, ke

enddo

return
end subroutine

!------------------------------------------------------------------------------
subroutine get_features_fnn(num_atoms, atype, pos, features, fp)
!------------------------------------------------------------------------------
integer,intent(in) :: num_atoms
real(8),intent(in),allocatable :: atype(:), pos(:,:)
real(rk),allocatable,intent(in out) :: features(:,:,:)
type(fnn_param),intent(in) :: fp

real(rk) :: rr(3), rr2, rij, dsum 
integer(ik) :: i, j, k, i1, j1, k1, l1, l2, l3, l4, ii
integer(ik) :: c1,c2,c3,ic(3),c4,c5,c6,ity,jty,inxn

integer(ik) :: nnbr, lnbr(MAXNEIGHBS)
integer(ik) :: idx, idx_stride, l1_stride, l2_stride, l3_stride
real(rk) :: r_ij(0:3), r_kj(0:3), eta_ij, eta_kj, fc_ij, fc_kj, rij_mu, rkj_mu
real(rk) :: cos_ijk, theta_ijk, lambda_ijk, zeta_G3a, zeta_G3b, zeta_const

real(rk) :: G3_mu_eta, G3a_xyz(3), G3a_c1, G3a_c2, G3a, G3b_xyz(3), G3b

features = 0.0

nbrlist(:,0) = 0

idx_stride = fp%feature_size_rad
l1_stride = size(fp%rad_eta)

!$omp parallel do default(shared) collapse(3) & 
!$omp private(c1,c2,c3,ic,c4,c5,c6,n,n1,m,m1,nty,mty,inxn,rr,rr2,rij,fr_ij,rij_mu,eta_ij,idx) 
do c1=0, cc(1)-1
do c2=0, cc(2)-1
do c3=0, cc(3)-1

  i = header(c1, c2, c3)
  do i1=1, nacell(c1, c2, c3)
     ity = nint(atype(i))

     !print'(3i6,i6,3f10.5)',c1,c2,c3,m,pos(m,1:3)

     do c4 = -1, 1
     do c5 = -1, 1
     do c6 = -1, 1
        ic(1:3) = [c1+c4, c2+c5, c3+c6]

        j = header(ic(1),ic(2),ic(3))
        do j1=1, nacell(ic(1), ic(2), ic(3))

           if(i/=j) then
             jty = nint(atype(j))
             inxn = pair_types(ity, jty)

             rr(1:3) = pos(j,1:3) - pos(i,1:3)
             rr2 = sum(rr(1:3)*rr(1:3))
             rij = sqrt(rr2)

             if(rij<fp%rad_rc) then
     !print'( f10.5, 2(3i6,i6,3f10.5) )',rij, c1,c2,c3,m,pos(m,1:3), c4,c5,c6,n,pos(n,1:3)

                if(rij<fp%ang_rc) then
                   nbrlist(i, 0) = nbrlist(i, 0) + 1
                   nbrlist(i, nbrlist(i, 0)) = j
                endif

                rr(1:3) = rr(1:3)/rij
                fc_ij = 0.5*( 1.0 + cos(pi*rij/fp%rad_damp) ) 

                do l1 = 1, size(fp%rad_mu)

                   rij_mu = rij - fp%rad_mu(l1)

                   do l2 = 1, size(fp%rad_eta)

                      eta_ij = exp( -fp%rad_eta(l2) * rij_mu * rij_mu )

                      idx = (inxn-1)*idx_stride + (l1-1)*l1_stride + l2

                      features(1:3,i,idx) = features(1:3,i,idx) + eta_ij*fc_ij*rr(1:3)
                   enddo

                enddo
             endif

           endif

           j=llist(j)
        enddo
     enddo; enddo; enddo

     i = llist(i)
  enddo
enddo; enddo; enddo
!$omp end parallel do 


! return if the angular cutoff is not given.
if(fp%ang_rc < 0.0) return

idx_stride = fp%feature_size_rad*size(fp%models)
l1_stride = size(fp%ang_zeta)*size(fp%ang_lambda)*size(fp%ang_eta)
l2_stride = size(fp%ang_zeta)*size(fp%ang_lambda)
l3_stride = size(fp%ang_zeta)

do j=1, num_atoms

   do i1=1, nbrlist(j,0)-1
      i = nbrlist(j,i1)

      r_ij(1:3) = pos(i,1:3) - pos(j,1:3)
      r_ij(0) = sqrt( sum(r_ij(1:3)*r_ij(1:3)) )

      fc_ij = 0.5*( 1.0 + cos(pi*r_ij(0)/fp%ang_damp) ) 

      do k1=i1+1, nbrlist(i,0)
         k = nbrlist(i,k1)

         r_kj(1:3) = pos(k,1:3) - pos(j,1:3)
         r_kj(0) = sqrt( sum(r_kj(1:3)*r_kj(1:3)) )

         fc_kj = 0.5*( 1.0 + cos(pi*r_kj(0)/fp%ang_damp) ) 

         cos_ijk = sum( r_ij(1:3)*r_kj(1:3) ) / ( r_ij(0) * r_kj(0) )
         theta_ijk = acos(cos_ijk)

         G3a_c1 = (r_ij(0)-r_kj(0)*cos_ijk)/(r_ij(0)*r_kj(0))/r_ij(0)
         G3a_c2 = (r_kj(0)-r_ij(0)*cos_ijk)/(r_ij(0)*r_kj(0))/r_kj(0)

         G3a_xyz(1:3) = r_ij(1:3)*G3a_c1 + r_kj(1:3)*G3a_c2
         G3b_xyz(1:3) = r_ij(1:3)/r_kj(0) + r_kj(1:3)/r_ij(0)

         do l1=1, size(fp%ang_mu)

            rij_mu = r_ij(0) - fp%ang_mu(l1)
            rkj_mu = r_kj(0) - fp%ang_mu(l1)

            do l2=1, size(fp%ang_eta)

               eta_ij = exp( -fp%ang_eta(l2) * rij_mu * rij_mu )
               eta_kj = exp( -fp%ang_eta(l2) * rkj_mu * rkj_mu )

               G3_mu_eta = eta_ij*eta_kj*fc_ij*fc_kj

               do l3=1, size(fp%ang_lambda)

                  lambda_ijk = 1.0 + fp%ang_lambda(l3)*cos_ijk 

                  do l4=1, size(fp%ang_zeta)

                     zeta_const = 2**(1-fp%ang_zeta(l4)) 
                     zeta_G3a = zeta_const * (lambda_ijk**fp%ang_zeta(l4))
                     zeta_G3b = zeta_const * (lambda_ijk**(fp%ang_zeta(l4)-1))

                     idx = fp%feature_size_rad + &
                         (l1-1)*l1_stride + (l2-1)*l2_stride + (l3-1)*l3_stride + l4

                     features(1:3,i,idx) = features(1:3,i,idx) + &
                         (G3a_xyz(1:3)*zeta_G3a + G3b_xyz(1:3)*zeta_G3b)*G3_mu_eta

         enddo; enddo; enddo; enddo
          
      enddo
   enddo

enddo

do i=1,num_atoms
   ity = int(atype(i))
   do j = 1, 3 ! xyz-loop
      features(j,i,:)=features(j,i,:) - fp%models(ity)%fstat(j)%mean(:)
      features(j,i,:)=features(j,i,:) / fp%models(ity)%fstat(j)%stddev(:)
   enddo 
enddo


!do i=1, num_atoms
!   print'(i6,i6,a1,3x,30i6)',i,nbrlist(i,0),',', nbrlist(i,1:10)
!enddo
!stop 'foo'

return
end subroutine

!------------------------------------------------------------------------------
subroutine mean_stddev_loader(mean, stddev, feature_size, path, suffix) 
!------------------------------------------------------------------------------
real(rk),allocatable,intent(in out) :: mean(:), stddev(:)
integer(ik),intent(in) :: feature_size
character(len=:),allocatable,intent(in) :: path
character(len=1),intent(in) :: suffix !<- ['x','y','z']

integer(ik) :: funit_m, funit_s
logical :: exist_m, exist_s
character(len=:),allocatable :: filename_m, filename_s

filename_m = path//'feature_mean_'//int_to_str(feature_size)//'.'//suffix
filename_s = path//'feature_stddev_'//int_to_str(feature_size)//'.'//suffix

if(.not.allocated(mean)) allocate(mean(feature_size))
if(.not.allocated(stddev)) allocate(stddev(feature_size))

inquire(file=filename_m, exist=exist_m)
inquire(file=filename_s, exist=exist_s)

! mean and stddev must exist, otherwise no feature vector standardization.
if(exist_m .and. exist_s) then
  open(newunit=funit_m, file=filename_m, access='stream', form='unformatted', status='old')
  read(funit_m) mean
  close(funit_m)
  open(newunit=funit_s, file=filename_s, access='stream', form='unformatted', status='old')
  read(funit_s) stddev
  close(funit_s)
  write(*,fmt='(a15,a30,a15,a30)') 'mean: ', filename_m, ', stddev: ', filename_s
else
  print'(a)', repeat('-',60)
  print'(a)', 'INFO: '//filename_m//' & '//filename_s//' not be used.'
  print'(a)', repeat('-',60)

  mean=0.0
  stddev=1.0
endif

end subroutine

!------------------------------------------------------------------------------
type(network) function network_ctor(dims, path, suffix) result(net)
!------------------------------------------------------------------------------

integer(ik),intent(in) :: dims(:)
character(len=:),allocatable,intent(in) :: path
character(len=1),intent(in) :: suffix !<- ['x','y','z']

character(len=:),allocatable :: filename_b, filename_w
character(len=:),allocatable :: arow, acol, alayer

integer(ik) :: i, nrow, ncol, fileunit
integer(ik) :: num_layers

net%dims = dims
num_layers = size(net%dims)

allocate(net%layers(num_layers))

do i=1, num_layers-1
  nrow = dims(i)
  ncol = dims(i+1)

  allocate(net%layers(i)%b(ncol))
  allocate(net%layers(i)%w(nrow,ncol))
  !print*,'i,nrow,ncol: ', i,nrow,ncol 

  alayer = int_to_str(i)
  arow = int_to_str(nrow)
  acol = int_to_str(ncol)

  filename_b = trim(path)//'b_'//alayer//'_'//acol//'.'//suffix
  open(newunit=fileunit, file=filename_b, access='stream', form='unformatted', status='old')
  read(fileunit) net%layers(i)%b
  close(fileunit)

  filename_w = trim(path)//'w_'//alayer//'_'//arow//'_'//acol//'.'//suffix
  open(newunit=fileunit, file=filename_w, access='stream', form='unformatted', status='old')
  read(fileunit) net%layers(i)%w

  write(*, fmt='(a30,2i6,a30,i6)') &
       'w: '//filename_w, shape(net%layers(i)%w), ' b: '//filename_b, size(net%layers(i)%b)

  close(fileunit)
enddo

end function

!------------------------------------------------------------------------------------------
subroutine get_cutoff_fnn(rcut, rcut2, maxrcut, radial_cutoff)
!------------------------------------------------------------------------------------------
real(8),allocatable,intent(in out) :: rcut(:), rcut2(:)
real(8),intent(in out) :: maxrcut
real(rk),intent(in) :: radial_cutoff

integer :: ity,jty,inxn

!--- get the cutoff length 
call allocator(rcut, 1, num_pairs)
call allocator(rcut2, 1, num_pairs)
call allocator(pair_types, 1, num_types, 1, num_types)

inxn=0
do ity=1, num_types
do jty=ity, num_types
   inxn = inxn + 1
   pair_types(ity,jty) = inxn

   rcut(inxn)  = radial_cutoff
   rcut2(inxn) = radial_cutoff*radial_cutoff
   !print'(a,3i6,2f10.5)','ity, jty, inxn: ', ity, jty, inxn, rcut(inxn), rcut2(inxn)

   pair_types(jty,ity) = pair_types(ity,jty) 
enddo
enddo

maxrcut = maxval(rcut)

end subroutine

end module
