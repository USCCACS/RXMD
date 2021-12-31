!------------------------------------------------------------------------------
module rxmdnn
!------------------------------------------------------------------------------
  use iso_c_binding

  use utils, only : pi, int_to_str
  use memory_allocator_mod
  use mpi_mod

  use base
  use lists_mod, only : getnonbondingmesh, linkedlist
  use communication_mod, only : copyatoms
  use msd_mod, only : msd_data, msd_add_initial_pos, msd_measure, msd_save
  use velocity_modifiers_mod
  use fileio
 
  implicit none

  type rxmdnn_model_params 
    character(len=:),allocatable :: filename
    character(len=:),allocatable :: element
    real(8) :: mass
  end type

  type, extends(force_field_class) :: rxmdnn_param
    type(rxmdnn_model_params), allocatable :: models(:) 
    contains 
       procedure :: print => rxmdnn_param_print
  end type

  type(rxmdnn_param),target :: rxmdnn_param_obj

  character(len=:),allocatable,private :: token

  ! timing for 1-feature calc & 2-force inference
  real(8),save,private :: tstart(0:3)=0.0, tfinish(0:3)=0.0

  real(8) :: maxrc_rxmdnn=0.d0

  real(c_double),allocatable,target :: nbrdist(:)
  type(c_ptr) :: nbrdist_ptr

  interface 

    ! create NN model, pass w&b to GPU, allocate nbrdist on CPU
    subroutine init_rxmdnn() bind(c,name="init_rxmdnn")
    end subroutine

    ! create NN model, pass w&b to GPU, allocate nbrdist on CPU
    subroutine init_rxmdnn_hybrid(natoms) bind(c,name="init_rxmdnn_hybrid")
       import :: c_int
       integer(c_int),value :: natoms
    end subroutine

    subroutine predict_rxmdnn_hybrid(natoms, maxnbrs, nbrdist_ptr) bind(c,name="predict_rxmdnn_hybrid")
       import :: c_int, c_ptr
       integer(c_int),value :: natoms, maxnbrs
       type(c_ptr),value :: nbrdist_ptr
    end subroutine
  
    subroutine predict_rxmdnn() bind(c,name="predict_rxmdnn")
    ! compute nbrlist (on CPU?), feature and predict E&F (return them to CPU?)

    end subroutine

    subroutine get_maxrc_rxmdnn(maxrc) bind(c,name="get_maxrc_rxmdnn")
        import :: c_double
        real(c_double) :: maxrc
    end subroutine
  
  end interface

contains

!------------------------------------------------------------------------------
subroutine allocate_nbrdist_rxmdnn()
!------------------------------------------------------------------------------

  allocate(nbrdist(NBUFFER*MAXNEIGHBS))

  return 

end subroutine


!------------------------------------------------------------------------------
subroutine rxmdnn_param_print(this)
!------------------------------------------------------------------------------
  class(rxmdnn_param), intent(in) :: this
  integer :: ia

  write(*,fmt='(a)') repeat('=',80)
  do ia=1, size(this%models)

     associate(m=>this%models(ia))
        print'(a,i3,2a,1x,f10.3,3x,a,1x,f8.3)', 'element,mass,filename: ', & 
           get_index_of_model(m%element,this%models),'-',m%element, m%mass, m%filename
     end associate
    
  enddo
  write(*,fmt='(a)') repeat('=',80)

end subroutine

!------------------------------------------------------------------------------
function get_index_of_model(element, models) result(idx)
!------------------------------------------------------------------------------
  type(rxmdnn_model_params),allocatable,intent(in) :: models(:)
  character(len=:),allocatable,intent(in) :: element
  integer :: idx
  do idx=1, size(models)
     if(models(idx)%element == element) return
  enddo
  idx = -1
  return
end function

!------------------------------------------------------------------------------
function rxmdnn_param_ctor(path) result(c)
!------------------------------------------------------------------------------
  character(len=:),allocatable,intent(in) :: path
  character(256) :: linein0
  character(len=:),allocatable :: linein 

  integer :: iunit
  type(rxmdnn_param) :: c 

  open(newunit=iunit, file=path, status='old', form='formatted')

  ! find how many models exist. 
  do while (.true.)
    read(iunit,'(a)',end=10) linein0
    linein = trim(adjustl(linein0))

    if(getstr(linein, token) > 0) then
       if(token=='model') call get_tokens_and_append_model(linein, c%models)
    endif
  end do
  10 rewind(iunit)

  if (size(c%models)<=0) stop 'ERROR: at least one model must be defined.'

  return

contains

  !------------------------------------------------------------------------------
  subroutine get_tokens_and_append_model(linein, models)
  !------------------------------------------------------------------------------
    character(len=:),allocatable,intent(in out) :: linein
    type(rxmdnn_model_params),allocatable,intent(in out) :: models(:)
    type(rxmdnn_model_params) :: mbuf
  
    if (getstr(linein, token) < 0) stop 'error while reading element name'
    mbuf%element = trim(adjustl(token))
    if (getstr(linein, token) < 0) stop 'error while reading element mass'
    read(token, *) mbuf%mass
    if (getstr(linein, token) < 0) stop 'error while reading model filename'
    mbuf%filename = trim(adjustl(token))
    read(token, *) mbuf%filename
  
    ! allocate zero-sized array
    if(.not.allocated(models)) allocate(models(0)) 
    models = [models, mbuf]
  
    return
  end subroutine
  

end function

!------------------------------------------------------------------------------
subroutine get_force_rxmdnn(ff, num_atoms, atype, pos, f, q)
!------------------------------------------------------------------------------
class(force_field_class),pointer,intent(in out) :: ff

integer,intent(in out) :: num_atoms 
real(8),intent(in out),allocatable :: atype(:), pos(:,:), q(:), f(:,:)

real(8) :: coo_i(3),  coo_j(3, MAXNEIGHBS), E_i, f3r(3, NBUFFER), F_i(3, num_atoms)
integer :: type_i, index_i, type_j(MAXNEIGHBS), index_j(MAXNEIGHBS)
integer :: i,j,j1,ity,n_j,stat,idx

call COPYATOMS(imode = MODE_COPY_FNN, dr=lcsize(1:3), atype=atype, pos=pos, ipos=ipos)
call LINKEDLIST(atype, pos, lcsize, header, llist, nacell)

call get_nbrlist_rxmdnn(num_atoms, atype, pos, maxrc) 

nbrdist_ptr = c_loc(nbrdist(1))
call predict_rxmdnn_hybrid(num_atoms, MAXNEIGHBS, nbrdist_ptr) 


CALL COPYATOMS(imode=MODE_CPBK, dr=dr_zero, atype=atype, pos=pos, f=f, q=q)

end subroutine

!------------------------------------------------------------------------------
subroutine mddriver_rxmdnn(mdbase, num_mdsteps) 
use utils, only : UTEMP, UTEMP0
use velocity_modifiers_mod, only : gaussian_dist_velocity, adjust_temperature, scale_temperature
!------------------------------------------------------------------------------
type(mdbase_class),intent(in out) :: mdbase
integer,intent(in) :: num_mdsteps
real(8) :: ctmp,cpu0,cpu1,cpu2,comp=0.d0

character(len=:),allocatable :: filebase

integer :: i,ity

if(reset_velocity_random) call gaussian_dist_velocity(atype, v)

call get_force_rxmdnn(mdbase%ff, natoms, atype, pos, f, q)

call cpu_time(cpu0)

!--- set force model
do nstep=0, num_mdsteps-1

  if(mod(nstep,pstep)==0) call print_e_rxmdnn(atype, v, q)

  call cpu_time(tstart(0))

  if(mod(nstep,fstep)==0) call OUTPUT(filebase, atype, pos, v, q, v)

  if(mod(nstep,sstep)==0.and.mdmode==4) &
      v(1:NATOMS,1:3)=vsfact*v(1:NATOMS,1:3)

   if(mod(nstep,sstep)==0.and.mdmode==5) &
      call scale_to_target_temperature(atype, v, treq)

   if(mod(nstep,sstep)==0.and.(mdmode==0.or.mdmode==6)) &
      call gaussian_dist_velocity(atype, v)

!--- element-wise velocity scaling
   if(mod(nstep,sstep)==0.and.mdmode==7) &
      call scale_temperature(atype, v)

   if(mod(nstep,sstep)==0.and.mdmode==8) &
      call adjust_temperature(atype, v)

   if(mod(nstep,sstep)==0.and.mdmode==9) &
      call maximally_preserving_BD(atype, v, vsfact) 

!--- MSD measurements
   call msd_add_initial_pos(msd_data, nstep, NATOMS, pos, ipos)
   call msd_measure(msd_data, nstep, NATOMS, atype, pos, ipos)

!--- total force may not be zero with FNN. fix linear momentum every pstep.
   if(mod(nstep,pstep)==0) call linear_momentum(atype, v)

!--- update velocity & position
   call vkick(1.d0, atype, v, f)

   pos(1:natoms,1:3)=pos(1:natoms,1:3)+dt*v(1:natoms,1:3)

!--- migrate atoms after positions are updated
   call COPYATOMS(imode=MODE_MOVE_FNN,dr=dr_zero,atype=atype,pos=pos, &
                  v=v,f=f,q=q,ipos=ipos)

   call cpu_time(cpu1)
   call get_force_rxmdnn(mdbase%ff, natoms, atype, pos, f, q)
   call cpu_time(cpu2)
   comp = comp + (cpu2-cpu1)

!--- update velocity
   call vkick(1.d0, atype, v, f)

   call cpu_time(tfinish(0))

enddo

!--- save the final configurations
call OUTPUT(filebase, atype, pos, v, q, v)

!--- update rxff.bin in working directory for continuation run
call WriteBIN(atype, pos, v, q, filebase)

!--- save result if msd_data%is_msd == true
call msd_save(msd_data)

call cpu_time(cpu2)
if(myid==0) print'(a,2f12.5)','comp, total (sec): ', comp, cpu2-cpu0

return
end subroutine

!------------------------------------------------------------------------------
subroutine get_nbrlist_rxmdnn(num_atoms, atype, pos, rcmax)
!------------------------------------------------------------------------------
integer,intent(in) :: num_atoms
real(8),intent(in),allocatable :: atype(:), pos(:,:)
real(8),intent(in) :: rcmax

real(8) :: rr(3), rr2, rij, dsum 
integer :: i, j, k, i1, j1, k1, l1, l2, l3, l4, ii, idx
integer :: c1,c2,c3,ic(3),c4,c5,c6,ity,jty,kty,inxn

nbrlist(:,0) = 0
nbrdist(:) = 0.d0

call cpu_time(tstart(2))

!!$omp parallel do default(shared) collapse(3) & 
!!$omp private(c1,c2,c3,ic,c4,c5,c6,rr,rr2,rij) 
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

             rr(1:3) = pos(i,1:3) - pos(j,1:3)
             rr2 = sum(rr(1:3)*rr(1:3))
             rij = sqrt(rr2)

             if(rij<rcmax) then 
               nbrlist(i, 0) = nbrlist(i, 0) + 1
               nbrlist(i, nbrlist(i, 0)) = j

               ii = nbrlist(i,0)
               idx = 4*((i-1)*MAXNEIGHBS+ii-1)
               !print'(2i6,f8.3,4f10.5,a)',i,ii,rij,nbrdist(idx+1:idx+4),' before'
               nbrdist(idx+1) = rij
               nbrdist(idx+2:idx+4) = rr(1:3)
               !print'(2i6,f8.3,4f10.5,i6,1x,a)',i,ii,rij,nbrdist(idx+1:idx+4),idx,' after'
             endif

           endif

           j=llist(j)
        enddo
     enddo; enddo; enddo

     i = llist(i)
  enddo
enddo; enddo; enddo
!!$omp end parallel do 
call cpu_time(tfinish(2))

!do i=1, natoms
!   print'(i6,i6,f8.3,a2)',i,nbrlist(i,0),rcmax,': '
!   do j1=1, MAXNEIGHBS
!      idx = 4*((i-1)*MAXNEIGHBS+j1-1)
!      print'(4f6.2 $)',nbrdist(idx+1:idx+4)
!      if(mod(j1-1,8)==7) print'(a1,i6)', ',',j1
!   enddo
!   print*
!enddo
!stop 'foo'

return

end subroutine

!-------------------------------------------------------------------------------------------
subroutine print_e_rxmdnn(atype, v, q)
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
real(8) :: tt=0.d0, Etotal

!ke=0.d0
!do i=1, NATOMS
!   ity=nint(atype(i))
!   ke = ke + hmas(ity)*sum(v(i,1:3)*v(i,1:3))
!enddo
!
!call MPI_ALLREDUCE (MPI_IN_PLACE, ke, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!ke = ke/GNATOMS
!tt = ke*UTEMP
!GKE = ke ! FIXME for ctmp = (treq*UTEMP0)/( GKE*UTEMP )
!
!call MPI_ALLREDUCE (MPI_IN_PLACE, Epot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!Epot = Epot/GNATOMS
!Etotal = Epot + ke
!
!if(myid==0) then
!   
!   cstep = nstep + current_step 
!
!   write(6,'(a,i9,es13.5,f10.3,3x,4f10.5)') &
!        'MDstep,Etotal,T(K),onestep(sec),force_calc,nbr_calc: ', cstep, Epot + ke, tt, tfinish(0:3)-tstart(0:3)
!endif

end subroutine

end module
