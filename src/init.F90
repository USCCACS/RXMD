!------------------------------------------------------------------------------------------
module init
!------------------------------------------------------------------------------------------

  use base, only : hh,hhi,lbox,obox, mdbox
  use qeq_mod
  use pqeq_mod
  use force_mod
  use reaxff_param_mod
  use velocity_modifiers_mod
  use memory_allocator_mod
  use fileio
  use fnn
  use lists_mod

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

if(is_reaxff) then
  call mdcontext_reaxff(atype, pos, v, f, q)
  print*,'get_mdcontext_func : mdcontext_reaxff'
else if(is_fnn) then
  call mdcontext_fnn(atype, pos, v, f, q)
  print*,'get_mdcontext_func : mdcontext_fnn'
else
  call mdcontext_reaxff(atype, pos, v, f, q)
  print*,'get_mdcontext_func : mdcontext_reaxff'
endif


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

end subroutine

!------------------------------------------------------------------------------------------
subroutine mdcontext_reaxff(atype, pos, v, f, q)
! This subroutine takes care of setting up initial system configuration.
! Unit conversion of parameters (energy, length & mass) are also done here.
!------------------------------------------------------------------------------------------
implicit none

real(8),intent(in out),allocatable,dimension(:) :: atype, q
real(8),intent(in out),allocatable,dimension(:,:) :: pos, v, f

integer :: i,j,k, ity
real(8) :: mm, dns, mat(3,3)

wt0 = MPI_WTIME()

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

!--- set force model
force_model_func => force_reaxff

!--- set md dirver function 
mddriver_func => mddriver_reaxff

!--- set force field parameters
call get_forcefield_params_reaxff(FFPath)

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



CTap(0:7)=(/1.d0, 0.d0, 0.d0, 0.d0,   -35.d0/(rctap)**4, &
          84.d0/(rctap)**5, -70.d0/(rctap)**6, &
          20.d0/(rctap)**7 /)


!--- determine cutoff distances only for exsiting atom pairs. cutoff cleanup using natoms_per_type()
call get_bondorder_cutoff(rc, rc2, maxrc, natoms_per_type)

!--- setup potential table. need the cutoff distance. 
call set_potentialtables_reaxff()

!--- update box-related variables based on the cutoff distance
call update_box_params(vprocs, vid, hh, lata, latb, latc, maxrc, &
                       cc, lcsize, hhi, mdbox, lbox, obox)

!--- Linked List & Near Neighb Parameters
call allocator(nbrlist,1,NBUFFER,0,MAXNEIGHBS)
call allocator(llist,1,NBUFFER)
call allocator(header,-MAXLAYERS,cc(1)-1+MAXLAYERS,-MAXLAYERS,cc(2)-1+MAXLAYERS,-MAXLAYERS,cc(3)-1+MAXLAYERS)
call allocator(nacell,-MAXLAYERS,cc(1)-1+MAXLAYERS,-MAXLAYERS,cc(2)-1+MAXLAYERS,-MAXLAYERS,cc(3)-1+MAXLAYERS)

call mdvariables_allocator_reaxff()

!--- setup 10[A] radius mesh to avoid visiting unecessary cells 
call GetNonbondingMesh()

!--- allocate & initialize Array size Stat variables
call allocator(maxas, 1, (ntime_step/pstep)+1, 1, nmaxas)

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
   write(6,'(a30,3f10.4)') "density [g/cc]:",dns
   write(6,'(a30,3i6)')  '# of linkedlist cell:', cc(1:3)
   write(6,'(a30,f10.3,2x,3f10.2)') "maxrc, lcsize [A]:", &
        maxrc,lata/cc(1)/vprocs(1),latb/cc(2)/vprocs(2),latc/cc(3)/vprocs(3)
   write(6,'(a30,3i6)')  '# of linkedlist cell (NB):', nbcc(1:3)
   write(6,'(a30,3f10.2)') "lcsize [A] (NB):", &
        lata/nbcc(1)/vprocs(1),latb/nbcc(2)/vprocs(2),latc/nbcc(3)/vprocs(3)
   write(6,'(a30,2i6)') "MAXNEIGHBS, MAXNEIGHBS10:", MAXNEIGHBS,MAXNEIGHBS10
   write(6,'(a30,i6,i9)') "NMINCELL, NBUFFER:", NMINCELL, NBUFFER
   write(6,'(a30,a12)') "DataDir :", trim(DataDir)

   print'(a30 $)','# of atoms per type:'
   do ity=1, nso
      if(natoms_per_type(ity)>0) print'(i12,a2,i2 $)',natoms_per_type(ity),' -',ity
   enddo
   print*

   write(6,'(a)') repeat('-',60)
   write(6,'(a30,2(a12,1x))') &
         "FFPath, ParmPath:", trim(FFPath),trim(ParmPath)
   write(6,'(a30,i6,es10.1,i6,i6)') "isQEq,QEq_tol,NMAXQEq,qstep:", &
                                     isQEq,QEq_tol,NMAXQEq,qstep
   write(6,'(a30,f8.3,f8.3)') 'Lex_fqs,Lex_k:',Lex_fqs,Lex_k

   print'(a)', "----------------------------------------------------------------"
   write(6,'(a)')  &
   "nstep  TE  PE  KE: 1-Ebond 2-(Elnpr,Eover,Eunder) 3-(Eval,Epen,Ecoa) 4-(Etors,Econj) 5-Ehbond 6-(Evdw,EClmb,Echarge)"

endif

end subroutine

end module
