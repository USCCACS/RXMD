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

  type(c_ptr) :: pos_ptr, type_ptr, force_ptr, nbrlist_ptr, fvar_ptr

  type nn_stat_type
    integer :: num_models
    real(8) :: emean, evar, fvar_max
    real(8),allocatable :: fvar(:,:)
  end type

  type(nn_stat_type) :: nn_stat

#ifdef RXMDNN

  interface 
    ! initialize rxmdtorch 
    subroutine init_rxmdtorch(myrank) bind(c,name="init_rxmdtorch")
       import :: c_int 
       integer(c_int),value :: myrank
    end subroutine

    subroutine predict_force_rxmdnn(natoms, nglobal, nbuffer, pos_ptr, type_ptr, force_ptr, nbrlist_ptr, energy) &
                                 bind(c,name="get_nn_force_torch")
       import :: c_int, c_double, c_ptr
       integer(c_int),value :: natoms, nbuffer, nglobal
       type(c_ptr),value :: pos_ptr, type_ptr, force_ptr, nbrlist_ptr
       real(c_double) :: energy 
    end subroutine

    subroutine update_current_model_rxmdnn(id) bind(c,name="update_current_model_rxmdnn")
        import :: c_int
        integer(c_int),value :: id 
    end subroutine

    subroutine get_num_models_rxmdnn(num_models) bind(c,name="get_num_models_rxmdnn")
        import :: c_int
        integer(c_int) :: num_models
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

    subroutine predict_force_rxmdnn(natoms, nbuffer, maxnbrs, pos_ptr, type_ptr, force_ptr, nbrlist_ptr, energy) 
       integer(c_int),value :: natoms, maxnbrs, nbuffer
       type(c_ptr),value :: pos_ptr, type_ptr, force_ptr, nbrlist_ptr
       real(c_double) :: energy 
    end subroutine

    subroutine update_current_model_rxmdnn(id)
        import :: c_int
        integer(c_int),value :: id 
    end subroutine

    subroutine get_num_models_rxmdnn(num_models)
        import :: c_int
        integer(c_int) :: num_models
    end subroutine

    subroutine get_maxrc_rxmdnn(maxrc) 
        real(c_double) :: maxrc
    end subroutine

#endif

!------------------------------------------------------------------------------
subroutine update_nn_stat(ns, num_models, num_atoms, esum, e2sum, fsum, f2sum)
!------------------------------------------------------------------------------
  type(nn_stat_type),intent(in out) :: ns 
  integer,intent(in) :: num_models, num_atoms
  real(8),intent(in) :: esum, e2sum, fsum(num_atoms,3), f2sum(num_atoms,3)

  real(8) :: fvar_max

  ns%num_models = num_models

  ! energy mean & variance
  ns%emean = esum/num_models
  ns%evar = e2sum/num_models - ns%emean**2

  ! force mean & variance
  if(.not.allocated(ns%fvar)) allocate(ns%fvar(NBUFFER,3))

  ns%fvar(1:num_atoms,1:3) = f2sum(1:num_atoms,1:3)/num_models - &
                             (fsum(1:num_atoms,1:3)/num_models)**2

  fvar_max = maxval(ns%fvar)
  call MPI_ALLREDUCE(MPI_IN_PLACE, fvar_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  ns%fvar_max = fvar_max

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
subroutine get_force_rxmdnn(ff, num_atoms, atype, pos, f, q, nn_stat)
use nnip_modifier, only : shortrep, shparams, springpot, spparams
!------------------------------------------------------------------------------
!use param_dftd

class(force_field_class),pointer,intent(in out) :: ff

integer,intent(in out) :: num_atoms 
real(8),intent(in out),allocatable, target :: atype(:), pos(:,:), q(:), f(:,:)
type(nn_stat_type),intent(in out) :: nn_stat 

integer :: i,j,j1,ity,jty, n_j,stat,idx

integer,allocatable,dimension(:),target :: nbrlist_nn

integer :: num_models, id
real(8) :: Enn, esum, e2sum, fsum(num_atoms, 3), f2sum(num_atoms, 3)
real(8) :: Erep

if(.not.allocated(nbrlist_nn)) allocate(nbrlist_nn(NBUFFER*MAXNEIGHBS))

call COPYATOMS(imode = MODE_COPY_FNN, dr=lcsize(1:3), atype=atype, pos=pos, q=q, ipos=ipos)
call LINKEDLIST(atype, pos, lcsize, header, llist, nacell)
call get_nbrlist_rxmdnn(pos, maxrc)

!--- packing
idx=0
do i=1, num_atoms
   idx=idx+1
   nbrlist_nn(idx)=nbrlist(i,0)
   !print'(a4,i6 $)','ni: ', nbrlist_nn(idx)
   do j1=1, nbrlist(i,0)
      idx=idx+1
      nbrlist_nn(idx)=nbrlist(i,j1)
      !print'(i6 $)', nbrlist_nn(idx)
   enddo
   !print*
enddo

call get_num_models_rxmdnn(num_models) 

esum=0.d0; e2sum=0.d0
fsum=0.d0; f2sum=0.d0 
do id = 1, num_models

  !print'(a,i3)','=== INFO === current_model id: ', id-1
  call update_current_model_rxmdnn(id-1)

  f=0.d0; Enn=0.d0

! for Fortran/C interface
  pos_ptr = c_loc(pos(1,1))
  type_ptr = c_loc(atype(1))
  force_ptr = c_loc(f(1,1))
  nbrlist_ptr = c_loc(nbrlist_nn(1))
  if(NATOMS>0) &
    call predict_force_rxmdnn(NATOMS, copyptr(6), NBUFFER, pos_ptr, type_ptr, force_ptr, nbrlist_ptr, Enn) 
  !print'(a,3es15.5)','myid,energy:',Enn

  call MPI_ALLREDUCE(MPI_IN_PLACE, Enn, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  esum = esum + Enn; e2sum = e2sum + Enn**2

  CALL COPYATOMS(imode=MODE_CPBK, dr=dr_zero, atype=atype, pos=pos, f=f, q=q)

  fsum(1:num_atoms,1:3) = fsum(1:num_atoms,1:3) + f(1:num_atoms, 1:3)
  f2sum(1:num_atoms,1:3) = f2sum(1:num_atoms,1:3) + f(1:num_atoms, 1:3)**2

enddo

! don't forget to convert energy to kcal/mol
esum = esum*Eev_kcal; e2sum = e2sum*Eev_kcal**2
fsum = fsum*Eev_kcal; f2sum = f2sum*Eev_kcal**2

! predicted force 
f(1:num_atoms,1:3) = fsum(1:num_atoms,1:3)/num_models

call update_nn_stat(nn_stat, num_models, num_atoms, esum, e2sum, fsum, f2sum)

if (shparams%flag) &
   call shortrep(myid, num_atoms, atype, pos, nbrlist, f, shparams, NBUFFER, MAXNEIGHBS)

if (spparams%flag) &
   call springpot(myid, num_atoms, atype, pos, nbrlist, f, spparams, NBUFFER, MAXNEIGHBS)


end subroutine

!------------------------------------------------------------------------------
subroutine mddriver_rxmdnn(mdbase, num_mdsteps) 
use utils, only : UTEMP, UTEMP0
use velocity_modifiers_mod, only : gaussian_dist_velocity, adjust_temperature, elementwise_scaling_temperature
use nnip_modifier, only : ceiling_array, ceparams, shortrep, shparams
!------------------------------------------------------------------------------
type(mdbase_class),intent(in out) :: mdbase
integer,intent(in) :: num_mdsteps
real(8) :: ctmp,cpu0,cpu1,cpu2,comp=0.d0
type(nn_stat_type) :: nn_stat

character(len=:),allocatable :: filebase

integer :: i,ity

if(reset_velocity_random.or.current_step==0) call gaussian_dist_velocity(atype, v)

call get_force_rxmdnn(mdbase%ff, natoms, atype, pos, f, q, nn_stat)

filebase = GetFileNameBase(DataDir,-1)

call cpu_time(cpu0)


!--- set force model
do nstep=0, num_mdsteps-1

  call cpu_time(tstart(0))

  if(mod(nstep,pstep)==0) call print_e_rxmdnn(atype, v, q, nn_stat)

  if(mod(nstep,fstep)==0) then
     filebase = GetFileNameBase(DataDir,current_step+nstep)
     call OUTPUT(filebase, atype, pos, v, q, v)
     !call OUTPUT(filebase, atype, pos, v, q, nn_stat%fvar)
  endif

  if(mod(nstep,sstep)==0 .and. is_tramp) &
     treq = ramp_to_target_temperature(myid, nstep, num_mdsteps) 

  if(mod(nstep,sstep)==0.and.mdmode==4) &
      v(1:NATOMS,1:3)=vsfact*v(1:NATOMS,1:3)

  if(mod(nstep,sstep)==0.and.mdmode==5) &
     call scale_to_target_temperature(atype, v, treq)

  if(mod(nstep,sstep)==0.and.(mdmode==0.or.mdmode==6)) &
     call gaussian_dist_velocity(atype, v)

!--- element-wise velocity scaling
  if(mod(nstep,sstep)==0.and.mdmode==7) &
     call elementwise_scaling_temperature(atype, v, treq)

  if(mod(nstep,sstep)==0.and.mdmode==8) &
     call adjust_temperature(atype, v)

  if(mod(nstep,sstep)==0.and.mdmode==9) &
     call maximally_preserving_BD(atype, v, vsfact) 

  if(is_vfceiling) &
     call velocity_and_force_ceiling(myid, atype, v, f, treq)

!--- MSD measurements
  call msd_add_initial_pos(msd_data, nstep, NATOMS, pos, ipos)
  call msd_measure(msd_data, nstep, NATOMS, atype, pos, ipos)

!--- total force may not be zero with FNN. fix linear momentum every pstep.
  if(mod(nstep,pstep)==0) call linear_momentum(atype, v)

!--- update velocity & position
  if (mod(nstep,sstep)==0 .and. ceparams%flag) &
      call ceiling_array(myid, NATOMS, f, 0.5d0*dt, ceparams%max_velocity, atype, dthm)

  call vkick(1.d0, atype, v, f)

  if (mod(nstep,sstep)==0.and.ceparams%flag) &
      call ceiling_array(myid, NATOMS, v, dt, ceparams%max_disp)

  pos(1:natoms,1:3)=pos(1:natoms,1:3)+dt*v(1:natoms,1:3)

!--- migrate atoms after positions are updated
  call COPYATOMS(imode=MODE_MOVE_FNN,dr=dr_zero,atype=atype,pos=pos, &
                  v=v,f=f,q=q,ipos=ipos)

   !print*,'ff_type_flag ', ff_type_flag, ff_type_flag == TYPE_NNQEQ

  call cpu_time(cpu1)
  call get_force_rxmdnn(mdbase%ff, natoms, atype, pos, f, q, nn_stat)
  if (ff_type_flag == TYPE_NNQEQ) call QEq(atype, pos, q)
  call cpu_time(cpu2)
  comp = comp + (cpu2-cpu1)

!--- update velocity
  if (mod(nstep,sstep)==0.and.ceparams%flag) &
      call ceiling_array(myid, NATOMS, f, 0.5d0*dt, ceparams%max_velocity, atype, dthm)

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
subroutine get_nbrlist_rxmdnn(pos, rcmax)
!------------------------------------------------------------------------------
real(8),intent(in),allocatable :: pos(:,:)
real(8),intent(in) :: rcmax

real(8) :: rr(3), rr2, rij, dsum 
integer :: i, j, k, i1, j1, k1, l1, l2, l3, l4, ii, idx
integer :: c1,c2,c3,ic(3),c4,c5,c6,ity,jty,kty,inxn

call cpu_time(tstart(2))

nbrlist(:,0) = 0

!!$omp parallel do default(shared) collapse(3) & 
!!$omp private(c1,c2,c3,ic,c4,c5,c6,rr,rr2,rij) 
do c1=0, cc(1)-1
do c2=0, cc(2)-1
do c3=0, cc(3)-1

  i = header(c1, c2, c3)
  do i1=1, nacell(c1, c2, c3)

     !print'(3i6,i6,3f10.5)',c1,c2,c3,m,pos(m,1:3)

     do c4 = -1, 1
     do c5 = -1, 1
     do c6 = -1, 1
        ic(1:3) = [c1+c4, c2+c5, c3+c6]

        j = header(ic(1),ic(2),ic(3))
        j = header(ic(1),ic(2),ic(3))
        do j1=1, nacell(ic(1), ic(2), ic(3))

           if(i/=j) then

             rr(1:3) = pos(i,1:3) - pos(j,1:3)
             rr2 = sum(rr(1:3)*rr(1:3))
             rij = sqrt(rr2)

             if(rij<rcmax) then 
               nbrlist(i, 0) = nbrlist(i, 0) + 1
               nbrlist(i, nbrlist(i, 0)) = j
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
subroutine print_e_rxmdnn(atype, v, q, ns)
use mpi_mod
use base, only : hh, hhi, natoms, gnatoms, mdbox, myid, ierr, hmas, wt0
use atoms
use memory_allocator_mod
! calculate the kinetic energy and sum up all of potential energies, then print them.
!-------------------------------------------------------------------------------------------
implicit none

real(8),allocatable,intent(in) :: atype(:), q(:)
real(8),allocatable,intent(in) :: v(:,:)
type(nn_stat_type),intent(in out) :: ns 

integer :: i,ity,cstep
real(8) :: tt=0.d0, Enn, Etotal
real(8),save :: time1 = 0d0

if (time1 == 0.d0) time1 = MPI_WTIME()

Enn = ns%emean

ke=0.d0
do i=1, NATOMS
   ity=nint(atype(i))
   ke = ke + hmas(ity)*sum(v(i,1:3)*v(i,1:3))
enddo

call MPI_ALLREDUCE (MPI_IN_PLACE, ke, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
tt = ke/GNATOMS*UTEMP
Etotal = ke + Enn 

if(myid==0) then
   cstep = nstep + current_step 
   write(6,'(a35,i9,3es13.5,2f10.3)') 'MDstep,Etotal,KE,PE,T(K),time(s): ', cstep, Etotal, ke, Enn, tt, MPI_WTIME()-time1
   if(ns%num_models > 1) &
           write(6,'(a35,i9,3es15.5,f10.3)') 'MDstep,Estdev(frac),Emean,Fvar_max,T(K): ', &
           cstep,sqrt(ns%evar)/abs(ns%emean), ns%emean, ns%fvar_max, tt
endif

time1 = MPI_WTIME()

end subroutine

end module
