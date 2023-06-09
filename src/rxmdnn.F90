!------------------------------------------------------------------------------
module rxmdnn
!------------------------------------------------------------------------------
  use iso_c_binding

  use utils, only : pi, int_to_str, token
  use memory_allocator_mod
  use mpi_mod

  use base
  use lists_mod, only : getnonbondingmesh, linkedlist
  use communication_mod, only : copyatoms
  use msd_mod, only : msd_data, msd_add_initial_pos, msd_measure, msd_save
  use velocity_modifiers_mod
  use fileio

  use qeq_mod
 
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

  !character(len=:),allocatable,private :: token

  ! timing for 1-feature calc & 2-force inference
  real(8),save,private :: tstart(0:3)=0.0, tfinish(0:3)=0.0

  real(8) :: maxrc_rxmdnn=0.d0

  type nbrdist_type
     !real(c_float),allocatable:: rij(:)
     !integer(c_int),allocatable :: jtype(:)

     ! FIXME dynamics allocation got released in torch side?
     real(c_float) :: rij(4*NBUFFER*MAXNEIGHBS) 
     integer(c_int) :: jtype(NBUFFER*MAXNEIGHBS) 
  end type

  !type(nbrdist_type),allocatable :: nbrdists(:)
  type(nbrdist_type) :: nbrdists

  integer,allocatable :: ndst_counts(:)

  type(c_ptr) :: nbrdist_ptr, nbrtype_ptr
  type(c_ptr) :: pos_ptr, type_ptr, force_ptr

#ifdef RXMDNN

  interface 
    ! initialize rxmdtorch 
    subroutine init_rxmdtorch(myrank) bind(c,name="init_rxmdtorch")
       import :: c_int 
       integer(c_int),value :: myrank
    end subroutine

    subroutine predict_rxmdtorch(natoms, nglobal, nbuffer, pos_ptr, type_ptr, force_ptr, energy) bind(c,name="get_nn_force_torch")
       import :: c_int, c_double, c_ptr
       integer(c_int),value :: natoms, nbuffer, nglobal
       type(c_ptr),value :: pos_ptr, type_ptr, force_ptr
       real(c_double) :: energy
    end subroutine

    subroutine get_maxrc_rxmdnn(maxrc) bind(c,name="get_maxrc_rxmdnn")
        import :: c_double
        real(c_double) :: maxrc
    end subroutine
  end interface

contains

#else

contains

    subroutine init_rxmdtorch(myrank) 
       integer(c_int),value :: myrank
    end subroutine

    subroutine predict_rxmdtorch(natoms, nbuffer, maxnbrs, pos_ptr, type_ptr, force_ptr, energy) 
       integer(c_int),value :: natoms, maxnbrs, nbuffer
       type(c_ptr),value :: pos_ptr, type_ptr, force_ptr
       real(c_double) :: energy
    end subroutine

    subroutine get_maxrc_rxmdnn(maxrc) 
        real(c_double) :: maxrc
    end subroutine

#endif


!------------------------------------------------------------------------------
subroutine allocate_nbrdist_rxmdnn(num_atoms)
!------------------------------------------------------------------------------
  integer,intent(in) :: num_atoms

  ! FIXME dynamics allocation got released in torch side?
  !allocate(nbrdists%rij(4*MAXNEIGHBS*num_atoms)) 
  !allocate(nbrdists%jtype(MAXNEIGHBS*num_atoms))
  nbrdists%rij = 0.d0
  nbrdists%jtype = 0

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
       if(token=='aenet'.or.token=='torch') call get_tokens_and_append_model(linein, c%models)
       if(token=='allegro') call get_tokens_and_append_model_allegro(linein, c%models)
    endif
  end do
  10 rewind(iunit)

  if (size(c%models)<=0) stop 'ERROR: at least one model must be defined.'

  return

contains
  !------------------------------------------------------------------------------
  subroutine get_tokens_and_append_model_allegro(linein, models)
  !------------------------------------------------------------------------------
    character(len=:),allocatable,intent(in out) :: linein
    type(rxmdnn_model_params),allocatable,intent(in out) :: models(:)
    type(rxmdnn_model_params),allocatable :: mbuf(:)
    character(len=:),allocatable :: modelpath
    integer :: i, n_elements
    real(8) :: mass
  
    if (getstr(linein, token) < 0) stop 'error while reading modelpath'
    modelpath = trim(adjustl(token))

    if (getstr(linein, token) < 0) stop 'error while reading n_elements'
    read(token, *) n_elements

    ! allocate zero-sized array
    if(.not.allocated(models)) allocate(models(0)) 

    allocate(mbuf(n_elements))
    do i = 1, n_elements
       if (getstr(linein, token) < 0) stop 'error while reading element'
       mbuf(i)%element = trim(adjustl(token))
       if (getstr(linein, token) < 0) stop 'error while reading mass'
       read(token, *) mbuf(i)%mass
       mbuf(i)%filename = modelpath
       models = [models, mbuf(i)]
    enddo

    return
  end subroutine

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
subroutine get_force_rxmdnn(ff, num_atoms, atype, pos, f, q, energy)
!------------------------------------------------------------------------------
class(force_field_class),pointer,intent(in out) :: ff

integer,intent(in out) :: num_atoms 
real(8),intent(in out),allocatable, target :: atype(:), pos(:,:), q(:), f(:,:)
real(8),intent(in out) :: energy

integer :: i,j,j1,ity,n_j,stat,idx

call COPYATOMS(imode = MODE_COPY_FNN, dr=lcsize(1:3), atype=atype, pos=pos, q=q, ipos=ipos)
call LINKEDLIST(atype, pos, lcsize, header, llist, nacell)

f=0.d0

! for Fortran/C interface
pos_ptr = c_loc(pos(1,1))
type_ptr = c_loc(atype(1))
force_ptr = c_loc(f(1,1))
call predict_rxmdtorch(NATOMS, copyptr(6) , NBUFFER, pos_ptr, type_ptr, force_ptr, energy) 

CALL COPYATOMS(imode=MODE_CPBK, dr=dr_zero, atype=atype, pos=pos, f=f, q=q)

! don't forget to convert energy to kcal/mol
f(1:num_atoms,:) = f(1:num_atoms,:)*Eev_kcal
energy = energy*Eev_kcal

end subroutine

!------------------------------------------------------------------------------
subroutine mddriver_rxmdnn(mdbase, num_mdsteps) 
use utils, only : UTEMP, UTEMP0
use velocity_modifiers_mod, only : gaussian_dist_velocity, adjust_temperature, scale_temperature
!------------------------------------------------------------------------------
type(mdbase_class),intent(in out) :: mdbase
integer,intent(in) :: num_mdsteps
real(8) :: ctmp,cpu0,cpu1,cpu2,comp=0.d0
real(8) :: nn_energy

character(len=:),allocatable :: filebase

integer :: i,ity


if(reset_velocity_random.or.current_step==0) call gaussian_dist_velocity(atype, v)

call get_force_rxmdnn(mdbase%ff, natoms, atype, pos, f, q, nn_energy)

filebase = GetFileNameBase(DataDir,-1)

call cpu_time(cpu0)


!--- set force model
do nstep=0, num_mdsteps-1

  if(mod(nstep,pstep)==0) call print_e_rxmdnn(atype, v, q, nn_energy)

  call cpu_time(tstart(0))

  if(mod(nstep,fstep)==0) then
     filebase = GetFileNameBase(DataDir,current_step+nstep)
     call OUTPUT(filebase, atype, pos, v, q, v)
  endif

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

   !print*,'ff_type_flag ', ff_type_flag, ff_type_flag == TYPE_NNQEQ

   call cpu_time(cpu1)
   call get_force_rxmdnn(mdbase%ff, natoms, atype, pos, f, q, nn_energy)
   if (ff_type_flag == TYPE_NNQEQ) call QEq(atype, pos, q)
   call cpu_time(cpu2)
   comp = comp + (cpu2-cpu1)

!--- update velocity
   call vkick(1.d0, atype, v, f)

   call cpu_time(tfinish(0))

enddo

!--- save the final configurations
call OUTPUT(filebase, atype, pos, v, q, v)

!--- update rxff.bin in working directory for continuation run
filebase = GetFileNameBase(DataDir,-1)
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
nbrdists%rij = 0.d0
nbrdists%jtype = 0

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
               idx = (i-1)*4*MAXNEIGHBS+(ii-1)*4

               !print'(3i6,f8.3,4f10.5,a)',i,ity,ii,rij,nbrdists(ity)%rij(idx+1:idx+4),' before'
               nbrdists%jtype(idx+1) = jty
               !print'(a,3i4,i9,i3)','i,j,ii,idx,jtype: ',i,j,ii,idx, nbrdists%jtype(idx+1)
               nbrdists%rij(idx+1) = rij
               nbrdists%rij(idx+2:idx+4) = rr(1:3)
               !print'(3i6,f8.3,4f10.5,i6,1x,a)',i,ity,ii,rij,nbrdists(ity)%rij(idx+1:idx+4),idx,' after'

               !print'(a,3i4,i9,i3,4f)','i,j,ii,idx,jtype,rij: ',i,j,ii,idx, nbrdists%jtype(idx+1), nbrdists%rij(idx+1:idx+4)
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

return

end subroutine

!-------------------------------------------------------------------------------------------
subroutine print_e_rxmdnn(atype, v, q, nn_energy)
use mpi_mod
use base, only : hh, hhi, natoms, gnatoms, mdbox, myid, ierr, hmas, wt0
use atoms
use memory_allocator_mod
! calculate the kinetic energy and sum up all of potential energies, then print them.
!-------------------------------------------------------------------------------------------
implicit none

real(8),allocatable,intent(in) :: atype(:), q(:)
real(8),allocatable,intent(in) :: v(:,:)
real(8),intent(in out) :: nn_energy

integer :: i,ity,cstep
real(8) :: tt=0.d0, Etotal

call MPI_ALLREDUCE (MPI_IN_PLACE, nn_energy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

ke=0.d0
do i=1, NATOMS
   ity=nint(atype(i))
   ke = ke + hmas(ity)*sum(v(i,1:3)*v(i,1:3))
enddo

call MPI_ALLREDUCE (MPI_IN_PLACE, ke, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
tt = ke/GNATOMS*UTEMP
Etotal = ke + nn_energy 

if(myid==0) then
   cstep = nstep + current_step 
   write(6,'(a,i9,3es13.5,f10.3)') 'MDstep,Etotal,KE,PE,T(K): ', cstep, Etotal, ke, nn_energy, tt
endif

end subroutine

end module

