module fnnin_parser

  use utils, only : getstr
  use base, only : force_field_class

  use iso_fortran_env, only: int32, int64, real32, real64, real128
  !integer,parameter :: rk = real64
  !integer,parameter :: rk = real128
  integer,parameter :: rk = real32

  !integer, parameter :: ik = int64
  integer, parameter :: ik = int32

  type :: layer
    real(rk),allocatable :: b(:)
    real(rk),allocatable :: w(:,:)
  end type

  type :: network
    type(layer),allocatable :: layers(:)
  end type

  type :: model_params
    character(len=:),allocatable :: element
    real(rk) :: mass
    integer(ik),allocatable :: hlayers(:)

    type(network),allocatable :: networks(:)
  end type

  type, extends(force_field_class) :: fnn_param
    real(rk),allocatable,dimension(:) :: rad_mu, rad_eta
    real(rk),allocatable,dimension(:) :: ang_mu, ang_eta, ang_zeta
    integer(ik),allocatable,dimension(:) :: ang_lambda

    real(rk) :: rad_rc, rad_damp, ang_rc, ang_damp

    type(model_params), allocatable :: models(:) 

    contains 
       procedure :: print => fnn_param_print

  end type

  interface get_tokens_and_append
      module procedure :: get_tokens_and_append_rv, get_tokens_and_append_iv, & 
                          get_tokens_and_append_rs, get_tokens_and_append_model
  end interface

  character(len=:),allocatable,private :: sbuf
  integer,private :: ibuf
  real(rk),private :: rbuf
  character(len=:),allocatable,private :: token

contains

  subroutine get_tokens_and_append_model(linein, models)
    implicit none
    character(len=:),allocatable,intent(in out) :: linein
    type(model_params),allocatable,intent(in out) :: models(:)
    type(model_params) :: mbuf

    if (getstr(linein, token) < 0) stop 'erro while reading element name'  
    mbuf%element = trim(adjustl(token))
    if (getstr(linein, token) < 0) stop 'erro while reading element mass'  
    read(token, *) mbuf%mass

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
    implicit none
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
    implicit none
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
    implicit none
    character(len=:),allocatable,intent(in out) :: linein
    real(4),intent(in out) :: scalar

    if (getstr(linein, token) > 0) then
       read(token,*) rbuf 
       scalar = rbuf
    endif 

    return
  end subroutine

  function fnn_param_ctor(path) result(c)
    implicit none
    character(len=:),allocatable,intent(in) :: path
    character(256) :: linein0
    character(len=:),allocatable :: linein

    type(fnn_param) :: c 
    integer :: iunit

    open(newunit=iunit, file=path, status='old', form='formatted')

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

  end function

  subroutine fnn_param_print(this)
    implicit none
    class(fnn_param), intent(in) :: this
    integer :: i,j

    print'(a)',repeat('-',60)
       print'(a,20f6.2)', 'rad_eta: ', this%rad_eta
       print'(a,20f6.2)', 'rad_mu: ', this%rad_mu
    print'(a)',repeat('-',60)
       print'(a,20f6.2)', 'ang_eta: ', this%ang_eta
       print'(a,20f6.2)', 'ang_mu: ', this%ang_mu
       print'(a,20i6)', 'ang_lambda: ', this%ang_lambda
       print'(a,20f6.2)', 'ang_zeta: ', this%ang_zeta
    print'(a)',repeat('-',60)
       print'(a,20f6.2)', 'rad_rc: ', this%rad_rc
       print'(a,20f6.2)', 'rad_damp: ', this%rad_damp
       print'(a,20f6.2)', 'ang_rc: ', this%ang_rc
       print'(a,20f6.2)', 'ang_damp: ', this%ang_damp
    print'(a)',repeat('-',60)

    do i = 1, size(this%models)
       associate(m => this%models(i))
          print*, m%element, m%mass, m%hlayers(:)
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

  integer(ik),allocatable :: num_dims(:)

  real(rk),allocatable :: features(:,:) 

  !integer(ik) :: num_Mu = 0, num_Eta = 0
  !real(rk),allocatable :: ml_Mu(:), ml_Eta(:)

  real(rk),parameter :: LJ_factor = 1.0
  !real(rk) :: ml_Rc = 0.0

  integer(ik),parameter :: num_forcecomps = 1
  integer(ik) :: num_features = 0

  integer(ik),parameter :: num_networks = 3

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
implicit none
class(force_field_class),allocatable,target,intent(in out) :: ff
type(fnn_param),pointer :: fp => null()
integer,intent(in out) :: natoms
real(8),intent(in out),allocatable :: atype(:), pos(:,:), q(:), f(:,:)

type(network),allocatable :: networks(:)

integer(ik) :: i, j, k, n, ncol, nrow, num_layers
integer(ik) :: m1,k1,n1
integer(ik) :: na, nn, nl

real(rk),allocatable :: x(:,:),y(:,:)

call COPYATOMS(imode=MODE_COPY_FNN, dr=lcsize(1:3), atype=atype, pos=pos)

call LINKEDLIST(atype, pos, lcsize, header, llist, nacell)

!call neighborlist(NMINCELL_FNN, atype, pos, pair_types, skip_check=.true.)

! not sure if this is the best way, but binding force_field_class to fnn_parm
select type(ff); type is (fnn_param) 
   fp => ff
end select
num_layers = size(fp%models(1)%hlayers)

call get_features_fnn(natoms, atype, pos, features, fp) 

networks = fp%models(1)%networks ! for now for testing

network_loop : do nn=1, num_networks

  y = features(1:natoms,1:num_features)
  !print*,'ncol,nrow:', num_dims(1), num_dims(2)

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

    !print'(a,4i6)','ncol,nrow:', num_dims(nl), num_dims(nl+1),shape(y)
    y = max(x,0.0) ! relu

  enddo layer_loop

  !--- update force
  f(1:natoms,nn) = x(1:natoms,1)

enddo network_loop

end subroutine

!------------------------------------------------------------------------------
subroutine mddriver_fnn(mdbase, num_mdsteps) 
!------------------------------------------------------------------------------
implicit none
type(mdbase_class),intent(in out) :: mdbase
integer,intent(in) :: num_mdsteps

integer :: i

!--- set force model
do i=1, num_mdsteps
   call get_force_fnn(mdbase%ff, natoms, atype, pos, f, q)
enddo

return
end subroutine

!------------------------------------------------------------------------------
subroutine get_features_fnn(num_atoms, atype, pos, features, fp)
!------------------------------------------------------------------------------
implicit none
integer,intent(in) :: num_atoms
real(8),intent(in),allocatable :: atype(:), pos(:,:)
real(rk),allocatable,intent(in out) :: features(:,:)
type(fnn_param),intent(in) :: fp

real(rk) :: rr(3), rr2, rij, dr_norm, dsum, eta, fr, rij_mu
integer(ik) :: i, j, j1, l1, l2, ii, idx, ia
integer(ik) :: c1,c2,c3,ic(3),c4,c5,c6,n,n1,m,m1,nty,mty,inxn

features = 0.0

!$omp parallel do default(shared) collapse(3) & 
!$omp private(c1,c2,c3,ic,c4,c5,c6,n,n1,m,m1,nty,mty,inxn,rr,rr2,rij,dr_norm,fr,rij_mu,eta,idx) 
do c1=0, cc(1)-1
do c2=0, cc(2)-1
do c3=0, cc(3)-1

  m = header(c1, c2, c3)
  do m1=1, nacell(c1, c2, c3)
     mty = nint(atype(m))

     do c4 = -1, 1
     do c5 = -1, 1
     do c6 = -1, 1
        ic(1:3) = [c1+c4, c2+c5, c3+c6]

        n = header(ic(1),ic(2),ic(3))
        do n1=1, nacell(ic(1), ic(2), ic(3))

           if(n/=m) then
             nty = nint(atype(n))
             inxn = pair_types(mty, nty)

             rr(1:3) = pos(n,1:3) - pos(m,1:3)
             rr2 = sum(rr(1:3)*rr(1:3))
             rij = sqrt(rr2)

             if(rr2<rc2(inxn)) then

                rr(1:3) = rr(1:3)/rij
                dr_norm = rij/fp%rad_rc
                fr = 0.5*( 1.0 + cos(pi*dr_norm) ) 

                do l1 = 1, size(fp%rad_mu)

                   rij_mu = rij - fp%rad_mu(l1)

                   do l2 = 1, size(fp%rad_eta)
                      eta = exp( -fp%rad_eta(l2) * rij_mu * rij_mu )
                      idx = (l1-1)*size(fp%rad_eta)*num_forcecomps + (l2-1)*num_forcecomps
                      features(n,idx+1) = features(n,idx+1) + eta*fr 
                   enddo

                enddo
             endif

           endif

           n=llist(n)
        enddo
     enddo; enddo; enddo

     m = llist(m)
  enddo
enddo; enddo; enddo
!$omp end parallel do 

!do i=1, num_atoms
!   print'(i6,30es13.5)',i,features(i,:)
!enddo
!stop 'foo'

return
end subroutine

!------------------------------------------------------------------------------
subroutine load_weight_and_bais_fnn(networks, num_dims, path) 
!------------------------------------------------------------------------------
   implicit none
   type(network),allocatable,intent(in out) :: networks(:)
   integer(ik),intent(in) :: num_dims(:)
   character(len=:),allocatable,intent(in) :: path

   character(1),parameter :: cxyz(3) = ['x','y','z']

   print*,'In get_feedforward_network, path: ', path
   networks(1) = network_ctor(num_dims,path,cxyz(1))
   networks(2) = network_ctor(num_dims,path,cxyz(2))
   networks(3) = network_ctor(num_dims,path,cxyz(3))
    
end subroutine

!------------------------------------------------------------------------------
type(network) function network_ctor(dims, path, suffix) result(net)
!------------------------------------------------------------------------------
implicit none

integer(ik),intent(in) :: dims(:)
character(len=:),allocatable,intent(in) :: path
character(len=1),intent(in) :: suffix

character(len=:),allocatable :: filename_b, filename_w
character(len=:),allocatable :: arow, acol, alayer

integer(ik) :: i, nrow, ncol, fileunit
integer(ik) :: num_layers

num_layers = size(dims)

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
implicit none
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
   print'(a,3i6,2f10.5)','ity, jty, inxn: ', ity, jty, inxn, rcut(inxn), rcut2(inxn)

   pair_types(jty,ity) = pair_types(ity,jty) 
enddo
enddo

maxrcut = maxval(rcut)

end subroutine

end module
