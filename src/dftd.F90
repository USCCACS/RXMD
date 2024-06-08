!-----------------------------------------------------------------------
module param_dftd
use element_name_manager
!-----------------------------------------------------------------------
! type declaration and initialization of variables for DFT-D
!
!  DFT-D second version : S.Grimme, J.Comput.Chem. 27, (2006), 1787-1799
!-----------------------------------------------------------------------
implicit none

integer :: init
!c VDW radii are in ang. determined by
!c atomic ROHF/TZV calculations and then taken as the radius of the 0.01 density contour
!c (scaling done below)
!c data obatained using the 0.005 contour and scale=1.0 gave too large differences
!c between the atoms i.e 1.1 1.5 1.42 1.36 1.3 for H,C,N,O,F
real(8), dimension(103) :: vander =                               &
   (/ 0.91d0,0.92d0, & !c H, He
      0.75d0,1.28d0,1.35d0,1.32d0,1.27d0,1.22d0,1.17d0,1.13d0, & !c Li-Ne
      1.04d0,1.24d0,1.49d0,1.56d0,1.55d0,1.53d0,1.49d0,1.45d0, & !c Na-Ar
      1.35d0,1.34d0, & !c K, Ca
      1.42d0,1.42d0,1.42d0,1.42d0,1.42d0, & !c Sc-Zn
      1.42d0,1.42d0,1.42d0,1.42d0,1.42d0, &
      1.50d0,1.57d0,1.60d0,1.61d0,1.59d0,1.57d0, & !c Ga-Kr
      1.48d0,1.46d0, & !c Rb, Sr
      1.49d0,1.49d0,1.49d0,1.49d0,1.49d0, & !c Y-Cd
      1.49d0,1.49d0,1.49d0,1.49d0,1.49d0, &
      1.52d0,1.64d0,1.71d0,1.72d0,1.72d0,1.71d0, & !c In, Sn, Sb, Te, I, Xe
      (1.d+99,init=1,49) /)

!c C6 parameters in J mol^-1 nm^6
real(8), dimension(103) :: c6 = &
   (/ 0.14d0, 0.08d0, &
      1.61d0, 1.61d0, 3.13d0, 1.75d0, 1.23d0, 0.70d0, 0.75d0, 0.63d0, &
      5.71d0, 5.71d0,10.79d0, 9.23d0, 7.84d0, 5.57d0, 5.07d0, 4.61d0, &
     10.80d0,10.80d0, &
     10.80d0,10.80d0,10.80d0,10.80d0,10.80d0, &
     10.80d0,10.80d0,10.80d0,10.80d0,10.80d0, &
     16.99d0,17.10d0,16.37d0,12.64d0,12.47d0,12.01d0, &
     24.67d0,24.67d0,&
     24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,&
     24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,& 
     37.32d0,38.71d0,38.44d0,31.74d0,31.50d0,29.99d0, &
     (1.d-14,init=1,49) /)

real(8) :: alp = 20.0d0
real(8) :: scalerad = 1.10d0   ! vander <- vander * scalerad
real(8) :: scalec6             ! DFT-functional depenent scale factor
real(8),parameter :: Joule2Cal = 4.184d0

real(8),  allocatable, dimension(:,:) :: c6ij, rrij

!-----energy, force, and stress
real(8) :: edisp                                 ! dispersion energy
real(8), allocatable, dimension(:,:) :: fdftd    ! force
real(8), dimension(3,3) :: strdftd = 0.d0        ! stress
real(8), allocatable, dimension(:) :: atom_dftd  ! for EDA

!-----cutoff distance
real(8)  :: cutoff

save

end module

!-----------------------------------------------------------------------
!subroutine setdftd(myid, jgga, ntype, nhk1, nhk2, zatom, evdj, hrdev, audang, avogad, ldcmd)
subroutine setdftd(myid, jgga)
!-----------------------------------------------------------------------
! setup DFT-D : an empirical correction for the vdW interaction
!-----------------------------------------------------------------------
use param_dftd

implicit none
integer :: myid, jgga 
integer :: ntype, it, jt, iz, jz, status

real(8) :: const

ntype = size(elements%e)

if(myid==0) print*, 'INFO: Parameters for DFT-D'
do it = 1, ntype
   iz = elements%e(it)%atomic
   if(myid==0) print'(a,i3,2es10.2)','  it,c6,R0: ', it, c6(iz), vander(iz)
end do

!------allocate memory
allocate( c6ij(ntype,ntype), rrij(ntype,ntype), stat=status )

if( jgga == 2 ) then      !c PBE
    scalec6 = 0.75d0
else if( jgga == 3 ) then !c RPBE
    scalec6 = 0.75d0      ! Caution!!  This is not optimized value!
else if( jgga == 4 ) then !c revPBE
    scalec6 = 0.75d0      ! Caution!!  This is not optimized value!
else
    print*, 'not yet supported in DFT-D'
end if

const  = scalec6*1d6/Joule2Cal   ! given in (J/mol)*nm^6 -> (kcal/mol)*A^6
c6     = c6     * const
vander = vander * scalerad 

cutoff = 0d0
do it =  1, ntype
   iz = elements%e(it)%atomic
   do jt = it, ntype
   jz = elements%e(jt)%atomic

   c6ij(it,jt) = sqrt(c6(iz)*c6(jz))
   rrij(it,jt) = vander(iz) + vander(jz)

   c6ij(jt,it) = c6ij(it,jt)
   rrij(jt,it) = rrij(it,jt)
   print'(a,2a3,2es10.2)','  ity,jty,c6ij,rrij: ', elements%e(it)%name,elements%e(jt)%name,c6ij(it,jt),rrij(it,jt)

   if (cutoff<rrij(it,jt) .and. rrij(it,jt)<1d9) cutoff = rrij(it,jt)
enddo;  enddo

print'(a,es10.2)','  cutoff: ', cutoff
cutoff = cutoff * 4.d0
cutoff = cutoff * cutoff

return
end

!!-----------------------------------------------------------------------
!subroutine get_force_dftd(num_atoms, atype, pos, f, energy)
!!-----------------------------------------------------------------------
!use param_dftd
!use base, only : nbrlist
!
!integer,intent(in out) :: num_atoms 
!real(8),intent(in out),allocatable, target :: atype(:), pos(:,:)
!real(8),intent(in out),allocatable, target :: f(:,:)
!real(8),intent(in out) :: energy 
!
!integer :: i,j,ity,jty
!real(8) :: rr(3), rr1, rr2, rr3, rr6, fdmp, edmp
!
!print*,'in get_force_dftd'
!
!energy=0
!do i=1, num_atoms
!   ity = nint(atype(i))
!   do j1=1, nbrlist(i,0)
!      j=nbrlist(i,j1)
!      jty = nint(atype(j))
!      print*,i,j,ity,jty
!
!      rr(1:3)=pos(j,1:3)-pos(i,1:3)
!      rr2 = sum(rr(1:3)*rr(1:3))
!      rr1 = sqrt(rr2)
!      rr3 = rr1*rr2
!      rr6 = rr3*rr3
!
!      edmp = exp(-alp*(rr1/rrij(ity,jty)-1d0))
!      fdmp = 1d0/(1d0+edmp)
!
!      print'(4i3,3es10.2)',i,j,ity,jty,rr1,edmp,fdmp
!
!      energy = energy + c6ij(ity,jty)/rr6*fmp
!
!   enddo
!enddo
!end subroutine

!------------------------------------------------------------------------------
subroutine get_force_rxmdnn(ff, num_atoms, atype, pos, f, q, energy, evar)
!------------------------------------------------------------------------------
use param_dftd

class(force_field_class),pointer,intent(in out) :: ff

integer,intent(in out) :: num_atoms 
real(8),intent(in out),allocatable, target :: atype(:), pos(:,:), q(:), f(:,:)
real(8),intent(in out) :: energy, evar

integer :: i,j,j1,ity,jty, n_j,stat,idx

integer,allocatable,dimension(:),target :: nbrlist_nn

real(8) :: rr(3), frc(3), fcoef, eij, rr1, rr2, rr3, rr6, rbr, fdmp, edmp, Edftd, sdftd
character(16) :: argv

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

!print'(a,2i6)','myid,maxnbr: ',myid, maxval(nbrlist(:,0))

f=0.d0; q=0.d0

! for Fortran/C interface
pos_ptr = c_loc(pos(1,1))
type_ptr = c_loc(atype(1))
force_ptr = c_loc(f(1,1))
nbrlist_ptr = c_loc(nbrlist_nn(1))
fvar_ptr = c_loc(q(1)) ! use q to store force uncertainty
if(NATOMS>0) &
   call predict_rxmdtorch(NATOMS, copyptr(6) , NBUFFER, pos_ptr, type_ptr, force_ptr, nbrlist_ptr, energy, evar, fvar_ptr) 
!print'(a,3es15.5)','myid,energy,evar,fvar:',energy,evar,fvar

CALL COPYATOMS(imode=MODE_CPBK, dr=dr_zero, atype=atype, pos=pos, f=f, q=q)

! don't forget to convert energy to kcal/mol
f(1:num_atoms,:) = f(1:num_atoms,:)*Eev_kcal
energy = energy*Eev_kcal

! DFT-D
Edftd=0d0; sdftd=1d0
if(find_cmdline_argc('--dftd',idx)) then

  call get_command_argument(idx+1,argv)
  read(argv,*) sdftd

  do i = 1, num_atoms
     ity = nint(atype(i))
     do j1 = 1, nbrlist(i,0)
        j = nbrlist(i,j1)
        jty = nint(atype(j))
  
        rr(1:3) = pos(j,1:3)-pos(i,1:3)
        rr2 = sum(rr(1:3)*rr(1:3))
        if (cutoff < rr2) continue
        rr1 = sqrt(rr2)
        rr3 = rr1*rr2
        rr6 = rr3*rr3
        rbr = rr1/rrij(ity,jty)
  
        edmp = exp(-alp*(rbr-1d0))
        fdmp = 1d0/(1d0+edmp)
        eij = -sdftd*c6ij(ity,jty)/rr6*fdmp
        Edftd = Edftd + eij
  
        fcoef = (eij/rr2)*(-6d0 + alp*rbr*edmp*fdmp)
        !print'(a,2i6,2a3,6es10.2,2f12.1)',&
        !        'i,j,ity,jty,c6ij,rrij,rr1,rc,edmp,fdmp,fcoef,eij: ', &
        !         i,j,elements%e(ity)%name,elements%e(jty)%name,&
        !         c6ij(ity,jty),rrij(ity,jty),&
        !         rr1,sqrt(cutoff),edmp,fdmp, fcoef,eij
  
        frc(1:3) = fcoef*rr(1:3)
        f(i,1:3) = f(i,1:3) - frc(1:3)
        !f(j,1:3) = f(j,1:3) + frc(1:3)
     enddo
  enddo

endif


end subroutine

