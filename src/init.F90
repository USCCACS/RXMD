!------------------------------------------------------------------------------------------
module init
!------------------------------------------------------------------------------------------

  use base
  use atoms
  use mpi_mod
  use qeq_mod, only : QEq
  use pqeq_mod, only : PQEq, initialize_eField, initialize_pqeq
  use force_mod, only : force_reaxff
  use memory_allocator_mod
  use fileio, only : ReadBIN

  use fnn, only : features, get_cutoff_fnn, ml_eta, ml_mu, ml_rc, networks, num_networks, &
                  num_dims, num_mu, num_eta, num_features, num_forcecomps, num_pairs, & 
                  mddriver_fnn, getnonbondingmesh, load_weight_and_bais_fnn, set_name_and_mass_fnn

  use reaxff_param_mod, only : chi, eta, mddriver_reaxff, &
      get_cutoff_bondorder, set_potentialtables_reaxff, get_forcefield_params_reaxff

contains

!------------------------------------------------------------------------------------------
subroutine mdcontext_base(atype, pos, v, f, q)
!------------------------------------------------------------------------------------------
implicit none

real(8),intent(in out),allocatable,dimension(:) :: atype, q
real(8),intent(in out),allocatable,dimension(:,:) :: pos, v, f

integer :: i,j,k, ity, l(3)

!--- an error trap
if(vprocs(1)*vprocs(2)*vprocs(3) /= nprocs ) then
  if(myid==0) write(6,'(a60,3i3,i5)')  &
  "ERROR: requested/allocated # of procs are not consistent: ", vprocs(1:3), nprocs
  call MPI_FINALIZE(ierr)
  stop
endif

!--- TODO make this stat a class
!--- allocate & initialize Array size Stat variables
call allocator(maxas, 1, (ntime_step/pstep)+1, 1, nmaxas)

!--- setup the vector ID and parity for processes, in x, y and z order.
vID(1)=mod(myid,vprocs(1))
vID(2)=mod(myid/vprocs(1),vprocs(2))
vID(3)=myid/(vprocs(1)*vprocs(2))
myparity(1)=mod(vID(1),2)
myparity(2)=mod(vID(2),2)
myparity(3)=mod(vID(3),2)

k=0
do i=1,3
!--- Add (1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1) to vector ID applying PBC. 
   do j=1,-1,-2

!--- save vector ID.
      l(1:3) = vID(1:3)

!--- apply PBC.
      l(i)=mod( (vID(i) + j + vprocs(i)), vprocs(i) )
     
!--- get direction index.
      k=k+1 ![123456] 

!--- convert vector ID to sequential ID.
      target_node(k) = l(1)+l(2)*vprocs(1)+l(3)*vprocs(1)*vprocs(2) 
   enddo         
                   
enddo    

!--- read atom and MD box data
call allocator(atype,1,NBUFFER)
call allocator(q,1,NBUFFER)
call allocator(pos,1,NBUFFER,1,3)
call allocator(v,1,NBUFFER,1,3)
call allocator(f,1,NBUFFER,1,3)
!--- index array for returning reaction force
call allocator(frcindx,1,NBUFFER)

call ReadBIN(atype, pos, v, q, f, trim(DataDir)//"/rxff.bin")

!--- get global number of atoms
GNATOMS = NATOMS ! Convert 4 byte to 8 byte
call MPI_ALLREDUCE(MPI_IN_PLACE, GNATOMS, 1, MPI_INTEGER8, MPI_SUM,  MPI_COMM_WORLD, ierr)

!--- TODO where to put the spring force extention? 
!--- for spring force
if (isSpring) then
  call allocator(ipos,1,NBUFFER,1,3)
  ipos(1:NATOMS,1:3)=pos(1:NATOMS,1:3)

  write(6,'(a)') repeat('-',60)
  write(6,'(a,f8.2 $)') 'springConst [kcal/mol] ', springConst 
  do ity=1, size(hasSpringForce) 
     if(hasSpringForce(ity)) print'(i3 $)',ity
  enddo
  print*
  write(6,'(a)') repeat('-',60)
endif

!--- print out parameters and open data file
if(myid==0) then
   write(6,'(a)') repeat('-',60)
   write(6,'(a30,i9,a3,i9)') "req/alloc # of procs:", vprocs(1)*vprocs(2)*vprocs(3), "  /",nprocs
   write(6,'(a30,3i9)')      "req proc arrengement:", vprocs(1),vprocs(2),vprocs(3)
   write(6,'(a30,es12.2)')   "time step[fs]:",dt*UTIME
   write(6,'(a30,i3, i10, i10)') "MDMODE CURRENTSTEP NTIMESTPE:", &
                                  mdmode, current_step, ntime_step
   write(6,'(a30,f12.3,f8.3,i9)') 'treq,vsfact,sstep:',treq*UTEMP0, vsfact, sstep
   write(6,'(a30,2i6)') 'fstep,pstep:', fstep,pstep
   write(6,'(a30,i24,i24)') "NATOMS GNATOMS:", NATOMS, GNATOMS
   write(6,'(a30,3f12.3)') "LBOX:",LBOX(1:3)
   write(6,'(a30,3f15.3)') "Hmatrix [A]:",HH(1:3,1,0)
   write(6,'(a30,3f15.3)') "Hmatrix [A]:",HH(1:3,2,0)
   write(6,'(a30,3f15.3)') "Hmatrix [A]:",HH(1:3,3,0)
   print'(a30,3f12.3)', 'lata,latb,latc:', lata,latb,latc
   print'(a30,3f12.3)', 'lalpha,lbeta,lgamma:', lalpha,lbeta,lgamma
   write(6,'(a30,2i9)') "NBUFFER, MAXNEIGHBS:", NBUFFER, MAXNEIGHBS
   write(6,'(a)') repeat('-',60)
   write(6,'(a30,a12)') "DataDir :", trim(DataDir)
   write(6,'(a30,2(a12,1x))') &
         "FFPath, ParmPath:", trim(FFPath),trim(ParmPath)
   write(6,'(a)') repeat('-',60)

endif

if(is_reaxff) then
  call mdcontext_reaxff()
  if(myid==0) print*,'get_mdcontext_func : mdcontext_reaxff'
else if(is_fnn) then
  call mdcontext_fnn()
  if(myid==0) print*,'get_mdcontext_func : mdcontext_fnn'
else
  call mdcontext_reaxff()
  if(myid==0) print*,'get_mdcontext_func : mdcontext_reaxff'
endif


end subroutine

!------------------------------------------------------------------------------------------
subroutine get_drived_properties(get_cutoff_func, set_potentialtables_func, &
                                 atmname, mass, natoms_per_type, &
                                 dthm, hmas, dns)
!------------------------------------------------------------------------------------------
implicit none

interface cutoff_func_interface
  subroutine get_cutoff(rcut, rcut2, maxrcut, natoms_per_type)
    real(8),allocatable,intent(in out) :: rcut(:), rcut2(:)
    real(8),intent(in out) :: maxrcut
    integer(8),allocatable,intent(in out),optional :: natoms_per_type(:)
  end subroutine
end interface

interface potentialtable_interface
  subroutine set_potentialtable()
  end subroutine
end interface

procedure(get_cutoff) :: get_cutoff_func
procedure(set_potentialtable),optional :: set_potentialtables_func

character(2),allocatable,intent(in) :: atmname(:)
real(8),allocatable,intent(in) :: mass(:)

integer(8),allocatable,intent(in out) :: natoms_per_type(:)
real(8),allocatable,intent(in out) :: dthm(:),hmas(:)
real(8),intent(in out) :: dns

integer :: i, ity

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
dns = 0.d0
do i=1, NATOMS
   ity = nint(atype(i))
   dns = dns + mass(ity)
enddo
call MPI_ALLREDUCE(MPI_IN_PLACE, dns, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, ierr)
dns=dns/mdbox*UDENS

!--- setup cutoff distance, rc, rc2, and maxrc. 
call get_cutoff_func(rc, rc2, maxrc, natoms_per_type)

!--- update box-related variables based on the cutoff distance
call update_box_params(vprocs, vid, hh, lata, latb, latc, maxrc, &
                       cc, lcsize, hhi, mdbox, lbox, obox)

!--- Linked List & Near Neighb Parameters
call allocator(nbrlist,1,NBUFFER,0,MAXNEIGHBS)
call allocator(llist,1,NBUFFER)
call allocator(header,-MAXLAYERS,cc(1)-1+MAXLAYERS,-MAXLAYERS,cc(2)-1+MAXLAYERS,-MAXLAYERS,cc(3)-1+MAXLAYERS)
call allocator(nacell,-MAXLAYERS,cc(1)-1+MAXLAYERS,-MAXLAYERS,cc(2)-1+MAXLAYERS,-MAXLAYERS,cc(3)-1+MAXLAYERS)

!--- setup potential tables
if(present(set_potentialtables_func)) call set_potentialtables_func()

if(myid==0) then
   write(6,'(a)') repeat('-',60)
   write(6,'(a30,3f10.4)') "density [g/cc]:",dns
   write(6,'(a30,3i6)')  '# of linkedlist cell:', cc(1:3)
   write(6,'(a30,f10.3,2x,3f10.2)') "maxrc, lcsize [A]:", &
        maxrc,lata/cc(1)/vprocs(1),latb/cc(2)/vprocs(2),latc/cc(3)/vprocs(3)

   print'(a30 $)','# of atoms per type:'
   do ity=1, num_pairs
      if(natoms_per_type(ity)>0) print'(i12,a2,i2 $)',natoms_per_type(ity),' -',ity
   enddo
   print*
   write(6,'(a)') repeat('-',60)
endif


return
end subroutine

!------------------------------------------------------------------------------------------
subroutine mdcontext_fnn()
!------------------------------------------------------------------------------------------
implicit none

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

!--- features(natoms, num_features)
call allocator(features,1, NBUFFER, 1, num_features) 

!--- use three types of networks, [x,y,z]
allocate(networks(num_networks))

!--- set force field parameters
call load_weight_and_bais_fnn(networks, str_gen('DAT'))

call get_drived_properties(get_cutoff_bondorder, set_potentialtables_reaxff, &
                           atmname, mass, natoms_per_type, &
                           dthm, hmas, dns)

!--- FNN specific output 
if(myid==0) then

   write(6,'(a)') repeat('-',60)
   write(6,*) 'num_dims: ', num_dims
   write(6,*) 'num_Mu, ml_Mu: ', num_Mu, ml_Mu
   write(6,*) 'num_Eta, ml_Eta', num_Eta, ml_Eta
   write(6,*) 'ml_Rc: ', ml_Rc
   write(6,*) 'num_forcecomps: ', num_forcecomps
   write(6,*) 'num_features: ', num_features
   write(6,*) 'num_networks: ', num_networks
   write(6,'(a)') repeat('-',60)

endif


end subroutine

!------------------------------------------------------------------------------------------
subroutine mdcontext_reaxff()
! This subroutine takes care of setting up initial system configuration.
! Unit conversion of parameters (energy, length & mass) are also done here.
!------------------------------------------------------------------------------------------
implicit none

integer :: i,j,k, ity
real(8) :: dns

wt0 = MPI_WTIME()

!--- allocate variales used during reaxff energy & force calculations
call mdvariables_allocator_reaxff()

!--- set charge model
if(isPQEq) then
  charge_model_func => PQEq
  rctap = rctap0_pqeq

  call allocator(spos,1,NBUFFER,1,3)

  call initialize_pqeq(chi,eta)
  if(isEfield) call initialize_eField(myid)

else
  charge_model_func => QEq
  rctap = rctap0
endif

!--- set taper function for the vdw and coulomb terms
rctap2 = rctap**2
CTap(0:7)=(/1.d0, 0.d0, 0.d0, 0.d0,   -35.d0/(rctap)**4, &
          84.d0/(rctap)**5, -70.d0/(rctap)**6, &
          20.d0/(rctap)**7 /)

!--- set force model
force_model_func => force_reaxff

!--- set md dirver function 
mddriver_func => mddriver_reaxff

!--- set force field parameters
call get_forcefield_params_reaxff(FFPath)

!--- get force field dependent parameters
call get_drived_properties(get_cutoff_func=get_cutoff_bondorder, &
                           set_potentialtables_func=set_potentialtables_reaxff, &
                           atmname=atmname, mass=mass, natoms_per_type=natoms_per_type, &
                           dthm=dthm, hmas=hmas, dns=dns)

!--- setup 10[A] radius mesh to avoid visiting unecessary cells 
call GetNonbondingMesh()

!--- ReaxFF specific output 
if(myid==0) then

   write(6,'(a)') repeat('-',60)
   write(6,'(a30,2i6)') "NMINCELL, MAXNEIGHBS10:", NMINCELL, MAXNEIGHBS10
   write(6,'(a30,i6,es10.1,i6,i6)') "isQEq,QEq_tol,NMAXQEq,qstep:", &
                                     isQEq,QEq_tol,NMAXQEq,qstep
   write(6,'(a30,f8.3,f8.3)') 'Lex_fqs,Lex_k:',Lex_fqs,Lex_k
   write(6,'(a)') repeat('-',60)

   write(6,'(a)')  &
   "nstep  TE  PE  KE: 1-Ebond 2-(Elnpr,Eover,Eunder) 3-(Eval,Epen,Ecoa) 4-(Etors,Econj) 5-Ehbond 6-(Evdw,EClmb,Echarge)"

endif

end subroutine

end module
