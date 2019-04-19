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

  integer,parameter :: num_types=1, num_pairs=num_types*(num_types+1)/2
  integer,parameter :: NMINCELL_FNN=1

  integer,allocatable :: pair_types(:,:)


contains

!------------------------------------------------------------------------------------------
subroutine mdcontext_fnn(atype, pos, v, f, q)
!------------------------------------------------------------------------------------------
implicit none

real(8),intent(in out),allocatable,dimension(:) :: atype, q
real(8),intent(in out),allocatable,dimension(:,:) :: pos, v, f

integer :: i,ity
real(8) :: dns, mm

call set_name_and_mass_fnn(mass, atmname)

!--- FIXME set all atomtype 1 for now
do i=1, size(atype)
   ity=nint(atype(i))
   if(ity>0) then 
      atype(i)=1.d0+l2g(atype(i))*1d-13
   endif
enddo

!--- set md dirver function 
mddriver_func => mddriver_fnn

!--- set force field parameters
call get_feedforward_network(str_gen('DAT'))

!--- set cutoff distance
call get_cutoff_fnn(rc, rc2, maxrc)

!--- features(natoms, num_features)
call allocator(features,1, NBUFFER, 1, num_features) 

!=============================================================================
! TODO this part should be merged into the basic context. 
!=============================================================================
!--- get number of atoms per each type. 
call allocator(natoms_per_type, 1, size(atmname))
do i=1, NATOMS
   ity=nint(atype(i))
   natoms_per_type(ity)=natoms_per_type(ity)+1
enddo
call MPI_ALLREDUCE(MPI_IN_PLACE, natoms_per_type, size(natoms_per_type), &
                   MPI_INTEGER8, MPI_SUM,  MPI_COMM_WORLD, ierr)

!--- dt/2*mass, mass/2
call allocator(dthm, 1, size(atmname))
call allocator(hmas, 1, size(atmname))
do ity=1, size(mass)
   if(mass(ity) > 0.d0) then
      dthm(ity) = dt*0.5d0/mass(ity)
      hmas(ity) = 0.5d0*mass(ity)
   endif
enddo

!--- get density 
mm = 0.d0
do i=1, NATOMS
   ity = nint(atype(i))
   mm = mm + mass(ity)
enddo
call MPI_ALLREDUCE(MPI_IN_PLACE, mm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
dns=mm/MDBOX*UDENS

!--- update box-related variables based on the cutoff distance
call update_box_params(vprocs, vid, hh, lata, latb, latc, maxrc, &
                       cc, lcsize, hhi, mdbox, lbox, obox)

!--- Linked List & Near Neighb Parameters
call allocator(nbrlist,1,NBUFFER,0,MAXNEIGHBS)
call allocator(llist,1,NBUFFER)
call allocator(header,-MAXLAYERS,cc(1)-1+MAXLAYERS,-MAXLAYERS,cc(2)-1+MAXLAYERS,-MAXLAYERS,cc(3)-1+MAXLAYERS)
call allocator(nacell,-MAXLAYERS,cc(1)-1+MAXLAYERS,-MAXLAYERS,cc(2)-1+MAXLAYERS,-MAXLAYERS,cc(3)-1+MAXLAYERS)

!=============================================================================
! TODO this part should be merged into the basic context. 
!=============================================================================

!--- FNN specific output 
if(myid==0) then

   write(6,'(a)') repeat('-',60)
   write(6,'(a30,3f10.4)') "density [g/cc]:",dns
   write(6,'(a30,3i6)')  '# of linkedlist cell:', cc(1:3)
   write(6,'(a30,f10.3,2x,3f10.2)') "maxrc, lcsize [A]:", &
        maxrc,lata/cc(1)/vprocs(1),latb/cc(2)/vprocs(2),latc/cc(3)/vprocs(3)
   write(6,'(a30,2i6)') "NMINCELL, MAXNEIGHBS10:", NMINCELL, MAXNEIGHBS10

   print'(a30 $)','# of atoms per type:'
   do ity=1, num_pairs
      if(natoms_per_type(ity)>0) print'(i12,a2,i2 $)',natoms_per_type(ity),' -',ity
   enddo
   print*

   write(6,'(a)') repeat('-',60)
   write(6,'(a30,a12)') "DataDir :", trim(DataDir)
   write(6,'(a30,2(a12,1x))') &
         "FFPath, ParmPath:", trim(FFPath),trim(ParmPath)
   write(6,'(a)') repeat('-',60)

endif


end subroutine

!------------------------------------------------------------------------------
subroutine get_force_fnn(natoms, atype, pos, f, q)
!------------------------------------------------------------------------------
implicit none
integer,intent(in out) :: natoms
real(8),intent(in out),allocatable :: atype(:), pos(:,:), q(:), f(:,:)

integer :: i, j, k
real(8) :: dr(3)

!do i=1, natoms
!   print'(a,i6,4es25.15)','i,atype(i),pos(i,1:3): ', i,atype(i),pos(i,1:3)
!enddo

call COPYATOMS(imode=MODE_COPY_FNN, dr=lcsize(1:3), atype=atype, pos=pos)

call LINKEDLIST(atype, pos, lcsize, header, llist, nacell)

!call neighborlist(NMINCELL_FNN, atype, pos, pair_types, skip_check=.true.)

call get_features(natoms, atype, pos, features) 

end subroutine

!------------------------------------------------------------------------------
subroutine mddriver_fnn(num_mdsteps) 
!------------------------------------------------------------------------------
implicit none
integer,intent(in) :: num_mdsteps

integer :: i

!do i=1, NATOMS
!   print*,i,atype(i),pos(i,1:3)
!enddo

!--- set force model
call get_force_fnn(natoms, atype, pos, f, q)

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
subroutine get_features(num_atoms, atype, pos, features) 
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
subroutine get_feedforward_network(path) 
!------------------------------------------------------------------------------
   character(len=:),allocatable,intent(in) :: path

   character(1),parameter :: cxyz(3) = ['x','y','z']
   type(network) :: netx,nety,netz

   print*,'In get_feedforward_network, path: ', path
   netx = network_ctor(num_dims,path//'/'//cxyz(1))
   nety = network_ctor(num_dims,path//'/'//cxyz(2))
   netz = network_ctor(num_dims,path//'/'//cxyz(3))
    
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
  print*,'i,nrow,ncol: ', i,nrow,ncol 

  alayer = int_to_str(i)
  arow = int_to_str(nrow)
  acol = int_to_str(ncol)

  filename_b = trim(netdata_prefix)//'_b_'//alayer//'_'//acol//'.net'
  print*,'b: ', filename_b, size(net%layers(i)%b)
  open(newunit=fileunit, file=filename_b, access='stream', form='unformatted', status='old')
  read(fileunit) net%layers(i)%b
  close(fileunit)

  filename_w = trim(netdata_prefix)//'_w_'//alayer//'_'//arow//'_'//acol//'.net'
  print*,'w: ', filename_w, shape(net%layers(i)%w)
  open(newunit=fileunit, file=filename_w, access='stream', form='unformatted', status='old')
  read(fileunit) net%layers(i)%w
  close(fileunit)
enddo

end function

!------------------------------------------------------------------------------------------
subroutine get_cutoff_fnn(rcut, rcut2, maxrcut)
!------------------------------------------------------------------------------------------
implicit none
real(8),allocatable,intent(in out) :: rcut(:), rcut2(:)
real(8),intent(in out) :: maxrcut

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

end module
