module fnn
  use iso_fortran_env, only: int32, int64, real32, real64, real128

  use utils
  use base
  use lists_mod
  use communication_mod
  use memory_allocator_mod

  implicit none

  !integer,parameter :: rk = real64
  !integer,parameter :: rk = real128
  integer,parameter :: rk = real32

  !integer, parameter :: ik = int64
  integer, parameter :: ik = int32

  integer(ik),allocatable :: num_dims(:)

  real(rk),allocatable :: features(:,:) 

  integer(ik),parameter :: num_Mu = 9
  real(rk),parameter :: ml_Mu(num_Mu) = [1.0,2.0,2.86,4.06,4.96,5.74,6.42,7.02,7.58]

  integer(ik),parameter :: num_Eta = 3
  !real(rk),parameter :: ml_Eta(num_Eta) = [0.5,1.0,3.0,20.0]
  real(rk),parameter :: ml_Eta(num_Eta) = [0.5,1.0,3.0]

  !real(rk),parameter :: LJ_factor = 3.405
  !real(rk),parameter :: ml_Rc = 2.550*LJ_factor
  real(rk),parameter :: LJ_factor = 1.0
  real(rk),parameter :: ml_Rc = 8.0

  !integer(rk),parameter :: num_forcecomps = 3
  integer(ik),parameter :: num_forcecomps = 1
  integer(ik),parameter :: num_features = num_Mu*num_Eta*num_forcecomps

  integer(ik),parameter :: num_networks = 3

  type :: layer
    real(rk),allocatable :: b(:)
    real(rk),allocatable :: w(:,:)
  end type

  type :: network
    type(layer),allocatable :: layers(:)
  end type

  type(network),allocatable :: networks(:)

  real(rk),allocatable :: infs(:,:)

  integer,parameter :: num_types=1, num_pairs=num_types*(num_types+1)/2
  integer,parameter :: NMINCELL_FNN=1

  integer,allocatable :: pair_types(:,:)

contains

!------------------------------------------------------------------------------
subroutine get_force_fnn(networks, natoms, atype, pos, f, q)
!------------------------------------------------------------------------------
implicit none
type(network),intent(in),allocatable :: networks(:)
integer,intent(in out) :: natoms
real(8),intent(in out),allocatable :: atype(:), pos(:,:), q(:), f(:,:)

integer(ik) :: i, j, k, n, ncol, nrow, num_layers
integer(ik) :: na, nn, nl

real(rk),allocatable :: x(:),y(:)

call COPYATOMS(imode=MODE_COPY_FNN, dr=lcsize(1:3), atype=atype, pos=pos)

call LINKEDLIST(atype, pos, lcsize, header, llist, nacell)

!call neighborlist(NMINCELL_FNN, atype, pos, pair_types, skip_check=.true.)

call get_features_fnn(natoms, atype, pos, features) 

num_layers = size(num_dims)

do na=1, natoms

   network_loop : do nn=1, num_networks

      if(allocated(y)) deallocate(y)
      allocate(y(num_features))

      y(1:num_features) = features(na,1:num_features)    
   
      layer_loop : do nl=1, num_layers-1

        nrow = num_dims(nl)
        ncol = num_dims(nl+1)
     
        !-- w*x + b
        if(allocated(x)) deallocate(x);  allocate(x(ncol))
     
        do k=1,ncol
          x(k) = networks(nn)%layers(nl)%b(k) + & 
                sum(networks(nn)%layers(nl)%w(1:nrow,k)*y(1:nrow)) 
        enddo

     
        !--- apply relu and update y
        if(allocated(y)) deallocate(y); allocate(y(ncol))
        y = max(x,0.0) ! relu
     
      enddo layer_loop

      if(size(x)==1) then
         !if(na==1) print'(a,3i6,f8.5)','f(na,nn): ', na,nn,nl,x(1)
         f(na,nn) = x(1)
      else
         print*,'ERROR: the last layer size is not 1'
         stop 
      endif

   enddo network_loop

   !print'(a,i6,6f10.5)','na,pos(na,1:3),f(na,1:3): ', na,pos(na,1:3),f(na,1:3)

enddo

end subroutine

!------------------------------------------------------------------------------
subroutine mddriver_fnn(num_mdsteps) 
!------------------------------------------------------------------------------
implicit none
integer,intent(in) :: num_mdsteps

integer :: i

!--- set force model
do i=1, num_mdsteps
   call get_force_fnn(networks, natoms, atype, pos, f, q)
enddo

return
end subroutine

!------------------------------------------------------------------------------
subroutine set_name_and_mass_fnn(atom_mass, atom_name)
!------------------------------------------------------------------------------
implicit none

real(8),allocatable,intent(in out) :: atom_mass(:)
character(2),allocatable,intent(in out) :: atom_name(:)

if(.not.allocated(atom_mass)) allocate(atom_mass(num_types))
if(.not.allocated(atom_name)) allocate(atom_name(num_types))

! aluminum name & mass 
atom_mass(1) = 27d0
atom_name(1) = 'Al'

print'(a,a3,f8.3,2i6)','atmname, mass: ', atom_name, atom_mass, &
       size(atom_name), size(atom_mass)

end subroutine

!------------------------------------------------------------------------------
subroutine get_features_fnn(num_atoms, atype, pos, features) 
!------------------------------------------------------------------------------
implicit none
integer,intent(in) :: num_atoms
real(8),intent(in),allocatable :: atype(:), pos(:,:)
real(rk),allocatable,intent(in out) :: features(:,:)

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
        ic(1:3) = [c1, c2, c3] + [c4, c5, c6]

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
                dr_norm = rij/ml_Rc
                fr = 0.5*( 1.0 + cos(pi*dr_norm) ) 

                do l1 = 1, num_Mu

                   rij_mu = rij - ml_Mu(l1)

                   do l2 = 1, num_Eta
                      eta = exp( -ml_Eta(l2) * rij_mu * rij_mu )
                      idx = (l1-1)*num_Eta*num_forcecomps + (l2-1)*num_forcecomps
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
subroutine load_weight_and_bais_fnn(networks, path) 
!------------------------------------------------------------------------------
   implicit none
   type(network),allocatable,intent(in out) :: networks(:)
   character(len=:),allocatable,intent(in) :: path

   character(1),parameter :: cxyz(3) = ['x','y','z']

   print*,'In get_feedforward_network, path: ', path
   networks(1) = network_ctor(num_dims,path//'/'//cxyz(1))
   networks(2) = network_ctor(num_dims,path//'/'//cxyz(2))
   networks(3) = network_ctor(num_dims,path//'/'//cxyz(3))
    
end subroutine

!------------------------------------------------------------------------------
type(network) function network_ctor(dims, netdata_prefix) result(net)
!------------------------------------------------------------------------------
implicit none

integer(ik),intent(in) :: dims(:)
character(len=*),intent(in) :: netdata_prefix

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

  filename_b = trim(netdata_prefix)//'_b_'//alayer//'_'//acol//'.net'
  open(newunit=fileunit, file=filename_b, access='stream', form='unformatted', status='old')
  read(fileunit) net%layers(i)%b
  close(fileunit)

  filename_w = trim(netdata_prefix)//'_w_'//alayer//'_'//arow//'_'//acol//'.net'
  open(newunit=fileunit, file=filename_w, access='stream', form='unformatted', status='old')
  read(fileunit) net%layers(i)%w

  write(*, fmt='(a30,2i6,a30,i6)') &
       'w: '//filename_w, shape(net%layers(i)%w), ' b: '//filename_b, size(net%layers(i)%b)

  close(fileunit)
enddo

end function

!------------------------------------------------------------------------------------------
subroutine get_cutoff_fnn(rcut, rcut2, maxrcut, natoms_per_type)
!------------------------------------------------------------------------------------------
implicit none
real(8),allocatable,intent(in out) :: rcut(:), rcut2(:)
real(8),intent(in out) :: maxrcut
integer(8),allocatable,intent(in out),optional :: natoms_per_type(:)

integer :: ity,jty,inxn

!--- get the cutoff length 
call allocator(rcut, 1, num_pairs)
call allocator(rcut2, 1, num_pairs)
call allocator(pair_types, 1, num_types, 1, num_types)

do ity=1, num_types
do jty=ity, num_types
   pair_types(ity,jty) = ity + (jty-1)*num_types
   inxn = pair_types(ity,jty) 

   rcut(inxn)  = ml_Rc
   rcut2(inxn) = ml_Rc*ml_Rc
   print'(a,3i6,2f10.5)','ity, jty, inxn: ', ity, jty, inxn, rcut(inxn), rcut2(inxn)

   pair_types(jty,ity) = pair_types(ity,jty) 
enddo
enddo

maxrcut = maxval(rcut)

end subroutine

end module
