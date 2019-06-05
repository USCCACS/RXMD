module fnnin_parser

  use utils, only : getstr, getfilenamebase, l2g
  use base, only : force_field_class
  use fileio, only : output, writebin
  use velocity_modifiers_mod, only : vkick

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
function sort_pos_by_type(num_atoms, num_types, atype, pos, v, q) result(anum)
!------------------------------------------------------------------------------
integer(ik),intent(in)  :: num_atoms, num_types
real(8),intent(in out),allocatable:: atype(:),pos(:,:),v(:,:),q(:)
real(8),allocatable :: atype0(:),pos0(:,:),v0(:,:),q0(:)

integer(ik) :: anum(num_types+1), anum1(num_types)

integer :: i, ii, ity, nsum

!FIXME should be a better way to sort these arrays here.

allocate(atype0(num_atoms),pos0(num_atoms,3),v0(num_atoms,3),q0(num_atoms)) 

anum=0
do i=1, num_atoms
   ity = nint(atype(i))+1
   anum(ity) = anum(ity)+1
enddo

anum1=0; nsum=0
do i=1, num_atoms

   ity = nint(atype(i))
   anum1(ity) = anum1(ity)+1

   ii = anum(ity)+anum1(ity)

   atype0(ii) = atype(i)
   pos0(ii,1:3) = pos(i,1:3)
   v0(ii,1:3) = v(i,1:3)
   q0(ii) = q(i)
enddo

atype(1:num_atoms)=atype0(1:num_atoms)
pos(1:num_atoms,1:3)=pos0(1:num_atoms,1:3)
v(1:num_atoms,1:3)=v0(1:num_atoms,1:3)
q(1:num_atoms)=q0(1:num_atoms)

deallocate(atype0,pos0,v0,q0)

end function

!------------------------------------------------------------------------------
subroutine get_force_fnn1(ff, num_atoms, atype, pos, f, q)
!------------------------------------------------------------------------------
class(force_field_class),pointer,intent(in out) :: ff
type(fnn_param),pointer :: fp => null()
integer,intent(in out) :: num_atoms 
real(8),intent(in out),allocatable :: atype(:), pos(:,:), q(:), f(:,:)

integer(ik) :: i, ii, ity, nn, nl, ncol, nrow, i1,i2, nsum
integer(ik),allocatable :: atom_per_type(:)

real(rk),allocatable :: x(:,:),y(:,:)

! not sure if this is the best way, but binding force_field_class to fnn_parm
select type(ff); type is (fnn_param) 
   fp => ff
end select

if(.not.allocated(atom_per_type)) allocate(atom_per_type(0:size(fp%models)))
atom_per_type = sort_pos_by_type(num_atoms, size(fp%models), atype, pos, v, q)

call COPYATOMS(imode=MODE_COPY_FNN, dr=lcsize(1:3), atype=atype, pos=pos)
call LINKEDLIST(atype, pos, lcsize, header, llist, nacell)

call get_features_fnn(num_atoms, atype, pos, features, fp) 

nsum = 0
do ity = 1, size(fp%models)

   nsum = nsum + atom_per_type(ity-1)
   i1 = nsum + 1
   i2 = nsum + atom_per_type(ity)
   !print*,'atom_per_type,nsum,i1,i2: ', atom_per_type, nsum, i1, i2

   do nn=1, num_networks_per_atom  ! fx,fy,fz loop
 
     y = features(nn,1:fp%feature_size,i1:i2)
 
     associate(n=>fp%models(ity)%networks(nn), m=>fp%models(ity)) 
        do nl=1, size(n%dims)-1

           x = matmul(n%layers(nl)%w,y)
           do i=1, size(x,dim=2)
              x(:,i)=x(:,i)+n%layers(nl)%b(:)
           enddo

           y = max(x,0.0) ! relu
           !print*,'ity,nl,shape(x),shape(%w),shape(y): ', ity,nl,shape(x),shape(n%layers(nl)%w),shape(y)
        enddo 
     end associate

     f(i1:i2,nn) = x(1,:) ! update force. 
   enddo

! TODO: obtained force is scaled by the scaling_factor assuming that the trained
! weight matrix & biased are also scaled.  better to have a check mechanism on their consistency. 
   f(i1:i2,1:3)=f(i1:i2,1:3)/fp%models(ity)%scaling_factor

enddo

end subroutine

!------------------------------------------------------------------------------
subroutine get_force_fnn(ff, num_atoms, atype, pos, f, q)
!------------------------------------------------------------------------------
class(force_field_class),pointer,intent(in out) :: ff
type(fnn_param),pointer :: fp => null()
integer,intent(in out) :: num_atoms 
real(8),intent(in out),allocatable :: atype(:), pos(:,:), q(:), f(:,:)

integer(ik) :: i, ity, nn, nl, ncol, nrow

real(rk),allocatable :: x(:),y(:)

! not sure if this is the best way, but binding force_field_class to fnn_parm
select type(ff); type is (fnn_param) 
   fp => ff
end select

call COPYATOMS(imode=MODE_COPY_FNN, dr=lcsize(1:3), atype=atype, pos=pos)
call LINKEDLIST(atype, pos, lcsize, header, llist, nacell)

call get_features_fnn(num_atoms, atype, pos, features, fp) 

do i=1, num_atoms 
  
   ity = atype(i)

   do nn=1, num_networks_per_atom  ! fx,fy,fz loop
 
     y = features(nn,1:fp%feature_size,i)
 
     associate(m=>fp%models(ity), n=>fp%models(ity)%networks(nn)) 
        do nl=1, size(n%dims)-1
           x = matmul(n%layers(nl)%w,y) + n%layers(nl)%b

           !print'(a,3i4,4i6)','i,ity,nn,shape(w),shape(b),shape(x): ', &
           !        i,ity,nn,shape(n%layers(nl)%w),shape(n%layers(nl)%b),shape(x)

           y = max(x,0.0) ! relu
        enddo 
     end associate
 
     f(i,nn) = x(1) ! update force
   enddo

! TODO: obtained force is scaled by the scaling_factor assuming that the trained
! weight matrix & biased are also scaled.  better to have a check mechanism on their consistency. 
   f(i,1:3)=f(i,1:3)/fp%models(ity)%scaling_factor

enddo

end subroutine

!------------------------------------------------------------------------------
subroutine mddriver_fnn(mdbase, num_mdsteps) 
use utils, only : UTEMP, UTEMP0
use velocity_modifiers_mod, only : gaussian_dist_velocity, adjust_temperature, scale_temperature
!------------------------------------------------------------------------------
type(mdbase_class),intent(in out) :: mdbase
integer,intent(in) :: num_mdsteps
real(8) :: ctmp,cpu0,cpu1,cpu2,comp=0.d0

integer :: i,ity

call get_force_fnn(mdbase%ff, natoms, atype, pos, f, q)

call cpu_time(cpu0)

!--- set force model
do nstep=0, num_mdsteps-1

  if(mod(nstep,pstep)==0) call print_e_fnn(atype, v, q)

  if(mod(nstep,fstep)==0) &
        call OUTPUT(atype, pos, v, q, GetFileNameBase(DataDir,current_step+nstep))

  if(mod(nstep,sstep)==0.and.mdmode==4) &
      v(1:NATOMS,1:3)=vsfact*v(1:NATOMS,1:3)

   if(mod(nstep,sstep)==0.and.mdmode==5) then
      ctmp = (treq*UTEMP0)/( GKE*UTEMP )
      v(1:NATOMS,1:3)=sqrt(ctmp)*v(1:NATOMS,1:3)
   endif

   if(mod(nstep,sstep)==0.and.(mdmode==0.or.mdmode==6)) &
      call gaussian_dist_velocity(atype, v)

!--- element-wise velocity scaling
   if(mod(nstep,sstep)==0.and.mdmode==7) &
      call scale_temperature(atype, v)

   if(mod(nstep,sstep)==0.and.mdmode==8) &
      call adjust_temperature(atype, v)

!--- update velocity & position
   call vkick(1.d0, atype, v, f)
   pos(1:natoms,1:3)=pos(1:natoms,1:3)+dt*v(1:natoms,1:3)

!--- migrate atoms after positions are updated
   call COPYATOMS(imode=MODE_MOVE_FNN,dr=[0.d0, 0.d0, 0.d0],atype=atype,pos=pos,v=v,f=f,q=q)

   call cpu_time(cpu1)
   call get_force_fnn(mdbase%ff, natoms, atype, pos, f, q)
   call cpu_time(cpu2)
   comp = comp + (cpu2-cpu1)

!--- update velocity
   call vkick(1.d0, atype, v, f)

enddo

!--- save the final configurations
call OUTPUT(atype, pos, v, q,  GetFileNameBase(DataDir,current_step+nstep))

!--- update rxff.bin in working directory for continuation run
call WriteBIN(atype, pos, v, q, GetFileNameBase(DataDir,-1))


call cpu_time(cpu2)
if(myid==0) print'(a,2f12.5)','comp, total (sec): ', comp, cpu2-cpu0

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
integer(ik) :: idx, idx_stride, l1_stride, l2_stride, l3_stride, l4_stride
real(rk) :: r_ij(0:3), r_kj(0:3), r_ij_norm(3), r_kj_norm(3), eta_ij, eta_kj, fc_ij, fc_kj, rij_mu, rkj_mu
real(rk) :: cos_ijk, lambda_ijk, rijk_inv, zeta_G3a, zeta_G3b, zeta_G3b_0, zeta_const

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

             rr(1:3) = pos(i,1:3) - pos(j,1:3)
             rr2 = sum(rr(1:3)*rr(1:3))
             rij = sqrt(rr2)

             if(rij<fp%rad_rc) then

                if(rij<fp%ang_rc .and. ity /= jty) then
                   nbrlist(i, 0) = nbrlist(i, 0) + 1
                   nbrlist(i, nbrlist(i, 0)) = j
                endif

                rr(1:3) = rr(1:3)/rij
                !fc_ij = 0.5*( 1.0 + cos(pi*rij/fp%rad_damp) ) 
                fc_ij = 0.5*( 1.0 + cos(3.14*rij/fp%rad_damp) )  ! for debugging. to be removed. 

                do l1 = 1, size(fp%rad_mu)

                   rij_mu = rij - fp%rad_mu(l1)

                   do l2 = 1, size(fp%rad_eta)

                      eta_ij = exp( -fp%rad_eta(l2) * rij_mu * rij_mu )

                      idx = (jty-1)*idx_stride + (l1-1)*l1_stride + l2
                      features(1:3,idx,i) = features(1:3,idx,i) + eta_ij*fc_ij*rr(1:3)
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

l1_stride = size(fp%ang_lambda)*size(fp%ang_zeta)*size(fp%ang_eta)
l2_stride = size(fp%ang_zeta)*size(fp%ang_eta)
l3_stride = size(fp%ang_eta)
l4_stride = 1

do j=1, num_atoms

   do i1=1, nbrlist(j,0)-1
      i = nbrlist(j,i1)

      r_ij(1:3) = pos(j,1:3) - pos(i,1:3)
      r_ij(0) = sqrt( sum(r_ij(1:3)*r_ij(1:3)) )
      r_ij_norm(1:3) = r_ij(1:3)/r_ij(0)

      !fc_ij = 0.5*( 1.0 + cos(pi*r_ij(0)/fp%ang_damp) ) 
      fc_ij = 0.5*( 1.0 + cos(3.14*r_ij(0)/fp%ang_damp) )  ! for debugging. to be removed. 

      do k1=i1+1, nbrlist(j,0)
         k = nbrlist(j,k1)

         r_kj(1:3) = pos(j,1:3) - pos(k,1:3)
         r_kj(0) = sqrt( sum(r_kj(1:3)*r_kj(1:3)) )
         r_kj_norm(1:3) = r_kj(1:3)/r_kj(0)

         rijk_inv = 1.0/(r_ij(0) * r_kj(0))

         !fc_kj = 0.5*( 1.0 + cos(pi*r_kj(0)/fp%ang_damp) ) 
         fc_kj = 0.5*( 1.0 + cos(3.14*r_kj(0)/fp%ang_damp) )  ! for debugging. to be removed. 

         cos_ijk = sum( r_ij(1:3)*r_kj(1:3) ) * rijk_inv

         G3a_c1 = r_ij(0)-r_kj(0)*cos_ijk
         G3a_c2 = r_kj(0)-r_ij(0)*cos_ijk

         G3a_xyz(1:3) = (r_ij_norm(1:3)*G3a_c1 + r_kj_norm(1:3)*G3a_c2)*rijk_inv*fc_ij*fc_kj
         G3b_xyz(1:3) = (r_ij_norm(1:3) + r_kj_norm(1:3))*fc_ij*fc_kj

! l1: mu, l2:lambda, l3:zeta, l4:eta
         do l1=1, size(fp%ang_mu)

            rij_mu = r_ij(0) - fp%ang_mu(l1)
            rkj_mu = r_kj(0) - fp%ang_mu(l1)

            do l2=1, size(fp%ang_lambda)

               lambda_ijk = 1.0 + fp%ang_lambda(l2)*cos_ijk 

               do l3=1, size(fp%ang_zeta)

                  zeta_const = 2.0**(1.0-fp%ang_zeta(l3))
                  zeta_G3a = fp%ang_zeta(l3) * fp%ang_lambda(l2) * zeta_const * (lambda_ijk**(fp%ang_zeta(l3)-1))
                  zeta_G3b_0 = zeta_const * (lambda_ijk**fp%ang_zeta(l3))

                  do l4=1, size(fp%ang_eta)

                     zeta_G3b = - 2.0*fp%ang_eta(l4) * zeta_G3b_0

                     eta_ij = exp( -fp%ang_eta(l4) * rij_mu * rij_mu )
                     eta_kj = exp( -fp%ang_eta(l4) * rkj_mu * rkj_mu )

                     G3_mu_eta = eta_ij*eta_kj

                     idx = idx_stride + 1 + &
                         (l1-1)*l1_stride + (l2-1)*l2_stride + (l3-1)*l3_stride + (l4-1)*l4_stride 

                     features(1:3,idx,j) = features(1:3,idx,j) + &
                         G3a_xyz(1:3)*zeta_G3a*G3_mu_eta + G3b_xyz(1:3)*zeta_G3b*G3_mu_eta

         enddo; enddo; enddo; enddo
          
      enddo
   enddo

enddo

do i=1, num_atoms
   ity = int(atype(i))
   do j = 1, 3 ! xyz-loop
      features(j,:,i)=(features(j,:,i)-fp%models(ity)%fstat(j)%mean(:))/fp%models(ity)%fstat(j)%stddev(:)
   enddo 
enddo


return
end subroutine

!------------------------------------------------------------------------------
subroutine mean_stddev_loader(mean, stddev, feature_size, path, suffix, verbose) 
!------------------------------------------------------------------------------
real(rk),allocatable,intent(in out) :: mean(:), stddev(:)
integer(ik),intent(in) :: feature_size
character(len=:),allocatable,intent(in) :: path
character(len=1),intent(in) :: suffix !<- ['x','y','z']
logical,optional :: verbose

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

  open(newunit=funit_m, file=filename_m, access='stream', form='formatted', status='old')
  if(myid==0) read(funit_m,*) mean
  call MPI_BCAST(mean, size(mean), MPI_FLOAT, 0, MPI_COMM_WORLD, ierr) ! TODO: support only MPI_FLOAT for now
  close(funit_m)

  open(newunit=funit_s, file=filename_s, access='stream', form='formatted', status='old')
  if(myid==0) read(funit_s,*) stddev
  call MPI_BCAST(stddev, size(stddev), MPI_FLOAT, 0, MPI_COMM_WORLD, ierr) ! TODO: support only MPI_FLOAT for now
  close(funit_s)

  if(present(verbose) .and. verbose) &
     write(*,fmt='(a15,a30,a15,a30)') 'mean: ', filename_m, ', stddev: ', filename_s
else
  if(present(verbose) .and. verbose) then
     print'(a)', repeat('-',60)
     print'(a)', 'INFO: '//filename_m//' & '//filename_s//' not be used.'
     print'(a)', repeat('-',60)
  endif

  mean=0.0
  stddev=1.0
endif

end subroutine

!------------------------------------------------------------------------------
type(network) function network_ctor(dims, path, suffix, verbose) result(net)
!------------------------------------------------------------------------------

integer(ik),intent(in) :: dims(:)
character(len=:),allocatable,intent(in) :: path
character(len=1),intent(in) :: suffix !<- ['x','y','z']
logical,optional :: verbose

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
  allocate(net%layers(i)%w(ncol,nrow))
  !print*,'i,nrow,ncol: ', i,nrow,ncol 

  alayer = int_to_str(i)
  arow = int_to_str(nrow)
  acol = int_to_str(ncol)

  filename_b = trim(path)//'b_'//alayer//'_'//acol//'.'//suffix
  open(newunit=fileunit, file=filename_b, access='stream', form='formatted', status='old')
  if(myid==0) read(fileunit,*) net%layers(i)%b
  call MPI_BCAST(net%layers(i)%b, size(net%layers(i)%b), MPI_FLOAT, 0, MPI_COMM_WORLD, ierr) ! TODO: support only MPI_FLOAT for now
  close(fileunit)

  filename_w = trim(path)//'w_'//alayer//'_'//arow//'_'//acol//'.'//suffix
  open(newunit=fileunit, file=filename_w, access='stream', form='formatted', status='old')
  if(myid==0) read(fileunit,*) net%layers(i)%w
  call MPI_BCAST(net%layers(i)%w, size(net%layers(i)%w), MPI_FLOAT, 0, MPI_COMM_WORLD, ierr) ! TODO: support only MPI_FLOAT for now
  close(fileunit)

  if(present(verbose) .and. verbose) &
     write(*, fmt='(a30,2i6,a30,i6)') &
        'w: '//filename_w, shape(net%layers(i)%w), ' b: '//filename_b, size(net%layers(i)%b)

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

!-------------------------------------------------------------------------------------------
subroutine print_e_fnn(atype, v, q)
use mpi_mod
use base, only : hh, hhi, natoms, gnatoms, mdbox, myid, ierr, hmas, wt0
use atoms
use memory_allocator_mod
! calculate the kinetic energy and sum up all of potential energies, then print them.
!-------------------------------------------------------------------------------------------
implicit none

real(8),allocatable,intent(in) :: atype(:), q(:)
real(8),allocatable,intent(in) :: v(:,:)

integer :: i,ity,cstep
real(8) :: tt=0.d0

ke=0.d0
do i=1, NATOMS
   ity=nint(atype(i))
   ke = ke + hmas(ity)*sum(v(i,1:3)*v(i,1:3))
enddo

call MPI_ALLREDUCE (MPI_IN_PLACE, ke, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
ke = ke/GNATOMS
tt = ke*UTEMP
GKE = ke ! FIXME for ctmp = (treq*UTEMP0)/( GKE*UTEMP )

if(myid==0) then
   
   cstep = nstep + current_step 

   write(6,'(a,i9,es13.5,2f8.2)') 'MDstep: ', cstep,ke, tt,GetTotalMemory()*1e-9

endif

end subroutine

end module
