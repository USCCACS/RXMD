module nnmm_mod

use base 
use atoms, only : qsfp, qsfv
use utils, only : get_boxparameters, update_box_params, getstr
use memory_allocator_mod

implicit none

integer,parameter :: NN_LAYER_TYPE=0, MM_LAYER_TYPE=1

type unitcell_type
  real(8) :: lattice(6)
  integer :: num_atoms
  real(8),allocatable :: pos(:), q(:)
  character(2),allocatable :: atype(:)
  integer,allocatable :: itype(:)
end type

! carries unitcell data for NN & MM domains
type(unitcell_type) :: NN, MM

character(len=:),allocatable,private :: token

type LJ_params
  real(8) :: eps, s11, s12, s13, s22, s23, s33
  contains 
    procedure :: show => show_MM_params
    procedure :: set => set_MM_params
end type

type shift
  character(len=:),allocatable :: atype 
  real(8) :: zshift
  contains 
    procedure :: show => show_shift_params
end type

type remove
  real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
  contains 
    procedure :: show => show_remove_params
end type

type stripe
  character(len=:),allocatable :: dtype 
  integer :: dir, num_stripes
  real(8) :: shift
  contains 
    procedure :: show => show_stripe_params
end type

type circular
  character(len=:),allocatable :: dtype 
  real(8) :: center(2), radius
  contains 
    procedure :: show => show_circ_params
end type

type rectangular
  character(len=:),allocatable :: dtype 
  real(8) :: xmin=0d0, ymin=0d0
  real(8) :: xmax=0d0, ymax=0d0
  contains 
    procedure :: show => show_rectangular_params
end type

type layer
  integer :: ltype, num_ucells
  real(8) :: zmin, zmax
  type(remove),allocatable :: removes(:)
  type(stripe),allocatable :: stripes(:)
  type(circular),allocatable :: circs(:)
  type(rectangular),allocatable :: rects(:)
  type(unitcell_type) :: ucell
  contains 
    procedure :: show => show_layer_params
end type

type nnmm_params
  type(LJ_params) :: MM_params
  type(shift),allocatable :: shifts(:)

  integer :: lx, ly
  type(layer),allocatable :: layers(:)

  contains 
    procedure :: show => show_nnmm_params
end type

type(nnmm_params) :: nnmmp

contains

!------------------------------------------------------------------------------
subroutine get_mm_force_and_mix(dp, PEmm, num_atoms, atype, pos, f, q, nbrlist)
!------------------------------------------------------------------------------
type(nnmm_params),intent(in) :: dp
integer,intent(in) :: num_atoms 
real(8),allocatable,intent(in) :: atype(:), pos(:,:), q(:)
integer,allocatable,intent(in)  :: nbrlist(:,:)

real(8),intent(in out) :: PEmm 
real(8),allocatable,intent(in out) :: f(:,:)

integer :: i, j, j1, ity, jty, num_frac_mm
real(8) :: f_mm(3), f_nn(3), ff(3), frac_mm, PElj 
real(8) :: eps, sigma(3,3), rr(3), rr1, rr2, rsig, rsig2, rsig6, rsig12

real(8),parameter :: rcut_mm = 3d0
integer :: qq

call dp%MM_params%set(eps,sigma)

PEmm = 0.d0
do i = 1, num_atoms
   ity = nint(atype(i))

   f_mm = 0.d0
   f_nn = 0.d0
   frac_mm = 0.d0
   num_frac_mm = 0
   do j1 = 1, nbrlist(i,0)
      j = nbrlist(i,j1)
      jty = nint(atype(j))

      rr(1:3) = pos(i,1:3) - pos(j,1:3)
      rr2 = sum(rr(1:3)*rr(1:3))
      rr1 = sqrt(rr2)

      rsig = sigma(ity,jty)/rr1
      rsig2 = rsig**2
      rsig6 = rsig**6
      rsig12 = rsig**12

      PElj = 4d0*eps*(rsig12-rsig6)
      !ff(1:3) = (-24d0)*eps*(2d0*rsig12-rsig6)/rr2*rr(1:3)
      ff(1:3) = 24d0*eps*rsig6/rr2*rr(1:3) ! only attractive term

      PEmm = PEmm + PElj
      f_mm(1:3) = f_mm(1:3) - ff(1:3)

      ff(1:3) = -48d0*eps*rsig12/rr2*rr(1:3)
      f_nn(1:3) = f_nn(1:3) - ff(1:3)

      !print'(a,2i3,3f10.5,1x,3f10.5,1x,f10.5)','ity,jty,sig(ity,jty),eps,rr1,ff(1:3),PElj: ',&
      !        ity,jty,sigma(ity,jty),eps,rr1,ff(1:3),PElj

      if(rr1 < rcut_mm) then
        num_frac_mm = num_frac_mm + 1
        qq = mod(nint(q(j)*1d9),10)
        if (qq == MM_LAYER_TYPE) frac_mm = frac_mm + 1d0
      endif

   enddo
   if(num_frac_mm>0) frac_mm = frac_mm/num_frac_mm
   !print'(a,i3,f10.5,i6)','i,frac,num_frac_mm: ', i, frac_mm, num_frac_mm

   f(i,1:3) = f_nn(1:3) + (1d0-frac_mm)*f(i,1:3) + frac_mm*f_mm(1:3)

enddo

end subroutine

!-------------------------------------------------------------------------------------------
function to_be_removed(lattice, rr0, la) result(flag)
!-------------------------------------------------------------------------------------------
   real(8),intent(in) :: lattice(6), rr0(3)
   type(layer),intent(in) :: la

   real(8) :: rr(3), zlength
   integer :: ia
   logical :: flag

   flag = .false.
   rr = rr0
   zlength = la%zmax - la%zmin

   if( (rr(3)<la%zmin) .or. (la%zmax<rr(3)) ) return 

   rr(3) = rr(3) - la%zmin

   do ia = 1, size(la%removes)

     if( (la%removes(ia)%xmin*lattice(1)<rr(1)) .and. (rr(1)<la%removes(ia)%xmax*lattice(1)) .and. &
         (la%removes(ia)%ymin*lattice(2)<rr(2)) .and. (rr(2)<la%removes(ia)%ymax*lattice(2)) .and. &
         (la%removes(ia)%zmin*zlength<rr(3))    .and. (rr(3)<la%removes(ia)%zmax*zlength)) flag = .true.
   enddo

end function

!-------------------------------------------------------------------------------------------
function apply_domain_operators(lattice, atype, rr0, la, shifts) result(rr)
!-------------------------------------------------------------------------------------------
   character(2),intent(in) :: atype
   real(8),intent(in) :: lattice(6), rr0(3)

   type(layer),intent(in) :: la
   type(shift),allocatable,intent(in) :: shifts(:)

   real(8) :: rr(3), dr(0:2), lstripe
   integer :: ia, ib, shift_sign

   rr = rr0

   ! must be within this layer
   if( (rr0(3)<la%zmin) .or. (la%zmax<rr0(3)) ) return 

!print'(a,3f10.5)','before, rr: ', rr
   do ia = 1, size(la%circs)

      shift_sign = 1
      if(la%circs(ia)%dtype == 'down') shift_sign = -1

      dr(1:2) = la%circs(ia)%center(1:2)*lattice(1:2) - rr(1:2)
      dr(0) = sqrt(sum(dr(1:2)*dr(1:2)))
      if(dr(0) < la%circs(ia)%radius) then

         do ib = 1, size(shifts)
            if(trim(adjustl(atype)) == shifts(ib)%atype) then
!print*,trim(adjustl(atype)), shifts(ib)%atype, trim(adjustl(atype))==shifts(ib)%atype, shift_sign*shifts(ib)%zshift
                    rr(3) = rr(3) + shift_sign*shifts(ib)%zshift
            endif
         enddo

      endif

   enddo

   do ia = 1, size(la%rects)

      shift_sign = 1
      if(la%rects(ia)%dtype == 'down') shift_sign = -1

      if((rr(1)<la%rects(ia)%xmin*lattice(1)) .or. (la%rects(ia)%xmax*lattice(1)<rr(1))) cycle
      if((rr(2)<la%rects(ia)%ymin*lattice(2)) .or. (la%rects(ia)%ymax*lattice(2)<rr(2))) cycle
      do ib = 1, size(shifts)
         if(trim(adjustl(atype)) == shifts(ib)%atype) then
!print*,trim(adjustl(atype)), shifts(ib)%atype, trim(adjustl(atype))==shifts(ib)%atype, shift_sign*shifts(ib)%zshift
                 rr(3) = rr(3) + shift_sign*shifts(ib)%zshift
         endif
      enddo
   enddo

   do ia = 1, size(la%stripes)

      shift_sign = 1
      if(la%stripes(ia)%dtype == 'down') shift_sign = -1

      ib = la%stripes(ia)%dir
      lstripe = lattice(ib)/la%stripes(ia)%num_stripes
      if (mod(int( (rr(ib) + la%stripes(ia)%shift*lattice(ib) )/lstripe),2) /= 0) cycle

      do ib = 1, size(shifts)
         if(trim(adjustl(atype)) == shifts(ib)%atype) then
!print*,trim(adjustl(atype)), shifts(ib)%atype, trim(adjustl(atype))==shifts(ib)%atype, shift_sign*shifts(ib)%zshift
                 rr(3) = rr(3) + shift_sign*shifts(ib)%zshift
         endif
      enddo
   enddo
!print'(a,3f10.5)',' after, rr: ', rr
   
end function

!-------------------------------------------------------------------------------------------
subroutine nnmd_setup_system(atype, pos, v, q, f, atmname, dp, NN, MM)
!-------------------------------------------------------------------------------------------
real(8),allocatable,dimension(:) :: atype,q
real(8),allocatable,dimension(:,:) :: pos,v,f

character(2),allocatable :: atmname(:)
type(nnmm_params),intent(in) :: dp
type(unitcell_type),intent(in)  :: NN, MM

type(unitcell_type)  :: NNMM

integer :: iunit, ix, iy, iz, il, ia, ib, ity, num_atoms, num_local, nc, layer_type
real(8) :: lattice(6), rr(3), mat(3,3), rmin(3), rmax(3), runit(3)

!--- allocate arrays
if(.not.allocated(atype)) call allocator(atype,1,NBUFFER)
if(.not.allocated(q)) call allocator(q,1,NBUFFER)
if(.not.allocated(pos)) call allocator(pos,1,NBUFFER,1,3)
if(.not.allocated(v)) call allocator(v,1,NBUFFER,1,3)
if(.not.allocated(f)) call allocator(f,1,NBUFFER,1,3)
if(.not.allocated(qsfp)) call allocator(qsfp,1,NBUFFER)
if(.not.allocated(qsfv)) call allocator(qsfv,1,NBUFFER)
f(:,:)=0.0

! get global number of atoms and system dimensions
num_atoms = 0
do il = 1, size(dp%layers)

   NNMM = dp%layers(il)%ucell

   ! REMARK: assuming all layers have the same cross section lx*NNMM%lattice(1) & ly*NNMM%lattice(2)
   lattice = [dp%lx*NNMM%lattice(1), dp%ly*NNMM%lattice(2), dp%layers(il)%zmax, 90d0, 90d0, 90d0]
   num_atoms = num_atoms + NNMM%num_atoms*dp%layers(il)%num_ucells
enddo

num_atoms = num_atoms * dp%lx*dp%ly

!print*,'num_atoms : ', num_atoms

!if(myid==0) then
!open(newunit=iunit,file="spto.xyz",form="formatted")
!write(iunit,'(i9)') num_atoms
!write(iunit,'(6f10.3)')  lattice
!endif

runit(1:3) = lattice(1:3)/vprocs(1:3)
rmin = runit*vid
rmax = runit*(vid+1)

!print'(a,i3,3i4,9f8.3)','myid,vprocs,runit,rmin,rmax:', myid,vid,runit,rmin,rmax

num_atoms = 0
num_local = 0

pos=0.d0

do il = 1, size(dp%layers)

   NNMM = dp%layers(il)%ucell

   do ix = 0, dp%lx-1
   do iy = 0, dp%ly-1
   do iz = 0, dp%layers(il)%num_ucells - 1

      do ia = 0, NNMM%num_atoms-1

         rr(1) = (NNMM%pos(3*ia+1) + ix)*NNMM%lattice(1)
         rr(2) = (NNMM%pos(3*ia+2) + iy)*NNMM%lattice(2)
         rr(3) = (NNMM%pos(3*ia+3) + iz)*NNMM%lattice(3) + dp%layers(il)%zmin
         rr = rr + 1d-9  ! tiny shift
         !if(myid==0) print'(2i6,3i3,3f10.5)',dp%layers(il)%ltype,ia,ix,iy,iz,NNMM%pos(3*ia+1:3*ia+3)

         if(to_be_removed(lattice, rr, dp%layers(il))) cycle

         num_atoms = num_atoms + 1

         rr = apply_domain_operators(lattice, NNMM%atype(ia+1), rr, dp%layers(il), dp%shifts)

         if(rr(1) <= rmin(1) .or. rmax(1) < rr(1)) cycle
         if(rr(2) <= rmin(2) .or. rmax(2) < rr(2)) cycle
         if(rr(3) <= rmin(3) .or. rmax(3) < rr(3)) cycle

         num_local = num_local + 1
         pos(num_local,1:3) = rr(1:3)

         do ity = 1, size(atmname)
            if(atmname(ity) == trim(adjustl(NNMM%atype(ia+1)))) atype(num_local) = ity + num_atoms*1d-13
         enddo

         ! store atomic charge and NN/MM region info in q array
         q(num_local) = NNMM%q(ia+1) + 1d-9*dp%layers(il)%ltype 

         !if(myid==0) write(iunit,'(a3,3f10.3,i3)') NNMM%atype(ia+1), rr(1:3), layer_type
      enddo
   enddo; enddo; enddo
enddo

NATOMS = num_local

!if(myid==0) close(iunit)

lata=lattice(1);   latb=lattice(2);  latc=lattice(3)
lalpha=lattice(4); lbeta=lattice(5); lgamma=lattice(6)

call get_boxparameters(mat,lata,latb,latc,lalpha,lbeta,lgamma)
do ia=1, 3
do ib=1, 3
   HH(ia,ib,0)=mat(ia,ib)
enddo; enddo
call update_box_params(vprocs, vid, hh, lata, latb, latc, maxrc, &
                       cc, lcsize, hhi, mdbox, lbox, obox)

!print*,'myid,natoms,vid,obox,lbox,mdbox: ',myid,natoms,vid, obox, lbox, mdbox

end subroutine

!------------------------------------------------------------------------------
subroutine nnmm_set_unitcells(NN, MM)
!------------------------------------------------------------------------------
  type(unitcell_type),intent(in out) :: NN,MM

  real(8) :: avex, avey 
  
  NN%lattice = [3.902d0,3.902d0,4.156d0,90d0,90d0,90d0]
  NN%pos = [ 0.000d0,   0.000d0,   0.000d0, & 
             0.500d0,   0.500d0,   0.5377d0, & 
             0.500d0,   0.500d0,   0.1118d0, & 
             0.000d0,   0.500d0,   0.6174d0, & 
             0.500d0,   0.000d0,   0.6174d0 ]
  NN%q = [0.0d0,0.0d0,0.0d0,0.0d0,0.0d0]
  NN%atype = ["Pb", "Ti", "O ", "O ", "O "]

  NN%num_atoms = size(NN%atype)

  MM%lattice = [3.945130d0, 3.945130d0, 3.945130d0, 90d0,  90d0,  90d0]
  MM%pos = [0.0d0,0.0d0,0.0d0, 0.5d0,0.5d0,0.5d0, 0.5d0,0.0d0,0.5d0, &
            0.5d0,0.5d0,0.0d0, 0.0d0,0.5d0,0.5d0]
  MM%q = [0.0d0,0.0d0,0.0d0,0.0d0,0.0d0]
! since NN model is for Pb/Ti/O, relabel Sr to Pb
  !MM%atype = ["Sr", "Ti", "O ", "O ", "O "]
  MM%atype = ["Pb", "Ti", "O ", "O ", "O "] 
  MM%num_atoms = size(MM%atype)

  !===================================================================
  ! REMARK: use MM lattice to apply mismatch strain to NN
  !===================================================================
  !NN%lattice(1:2) = MM%lattice(1:2)
  avex = (NN%lattice(1) + MM%lattice(1))*0.5d0
  avey = (NN%lattice(2) + MM%lattice(2))*0.5d0
  NN%lattice(1) = avex;  MM%lattice(1) = avex
  NN%lattice(2) = avey;  MM%lattice(2) = avey
  
  return
end subroutine

!------------------------------------------------------------------------------
subroutine show_nnmm_params(this)
!------------------------------------------------------------------------------
  class(nnmm_params),intent(in) :: this
  integer :: ia

  print'(a)',repeat('-',60)
  call this%MM_params%show()
  print*

  do ia=1, size(this%shifts)
     call this%shifts(ia)%show()
  enddo 
  print*

  do ia = 1, size(this%layers)
     print'(a10,i3,a5)','=== layer ', ia, ' === '
     call this%layers(ia)%show()
     print*
  enddo
  print'(a)',repeat('-',60)
end subroutine

!------------------------------------------------------------------------------
subroutine show_layer_params(this)
!------------------------------------------------------------------------------
  class(layer),intent(in) :: this
  integer :: ia

  print'(a40,2i6,2f10.3)', 'ltype,num_ucells,zmin,zmax: ', &
          this%ltype,this%num_ucells,this%zmin,this%zmax

  do ia = 1, size(this%removes)
     call this%removes(ia)%show()
  enddo 
  do ia = 1, size(this%stripes)
     call this%stripes(ia)%show()
  enddo 
  do ia = 1, size(this%circs)
     call this%circs(ia)%show()
  enddo 
  do ia = 1, size(this%rects)
     call this%rects(ia)%show()
  enddo 
end subroutine

!------------------------------------------------------------------------------
subroutine set_MM_params(this, eps, sigma)
!------------------------------------------------------------------------------
  class(LJ_params),intent(in) :: this
  real(8) :: eps, sigma(3,3)
  eps = this%eps
  sigma(1,1) = this%s11; sigma(2,2) = this%s22; sigma(3,3) = this%s33
  sigma(1,2) = this%s12; sigma(2,1) = sigma(1,2)
  sigma(1,3) = this%s13; sigma(3,1) = sigma(1,3)
  sigma(2,3) = this%s23; sigma(3,2) = sigma(2,3)
end subroutine

!------------------------------------------------------------------------------
subroutine show_MM_params(this)
!------------------------------------------------------------------------------
  class(LJ_params),intent(in) :: this
  print'(a40,es12.4,2x,6es12.4)','MM(eps,sigma[11,12,13,22,23,33]): ', &
          this%eps, this%s11, this%s12, this%s13, this%s22, this%s23, this%s33
end subroutine

!------------------------------------------------------------------------------
subroutine show_shift_params(this)
!------------------------------------------------------------------------------
  class(shift),intent(in) :: this
  print'(a40,a8,3es12.2,2x,es12.2)','atom shift(type,zshift): ', this%atype, this%zshift
end subroutine

!------------------------------------------------------------------------------
subroutine show_remove_params(this)
!------------------------------------------------------------------------------
  class(remove),intent(in) :: this
  print'(a40,6es12.2)','remove(xmin,xmax,ymin,ymax,zmin,zmax): ', &
          this%xmin, this%xmax, this%ymin, this%ymax, this%zmin, this%zmax
end subroutine

!------------------------------------------------------------------------------
subroutine show_stripe_params(this)
!------------------------------------------------------------------------------
  class(stripe),intent(in) :: this
  print'(a40,a8,i6,es12.2,i6)','stripe_domain(dir,shift,num_stripes): ', &
  this%dtype, this%dir, this%shift, this%num_stripes
end subroutine

!------------------------------------------------------------------------------
subroutine show_circ_params(this)
!------------------------------------------------------------------------------
  class(circular),intent(in) :: this
  print'(a40,a8,2es12.2,2x,es12.2)','circular_domain(type,center,radius): ', &
  this%dtype, this%center(1:2), this%radius
end subroutine

!------------------------------------------------------------------------------
subroutine show_rectangular_params(this)
!------------------------------------------------------------------------------
  class(rectangular),intent(in) :: this
  print'(a40,a8,4es12.2)','rectangular_domain(type,min,max): ', &
  this%dtype, this%xmin, this%ymin, this%xmax, this%ymax
end subroutine

!------------------------------------------------------------------------------
subroutine get_tokens_and_append_MM_params(linein, ljp)
!------------------------------------------------------------------------------
  character(len=:),allocatable,intent(in out) :: linein
  type(LJ_params),intent(in out) :: ljp

  if (getstr(linein, token) < 0) stop 'error while reading MMparam eps'
  read(token, *) ljp%eps
  if (getstr(linein, token) < 0) stop 'error while reading MMparam s11'
  read(token, *) ljp%s11
  if (getstr(linein, token) < 0) stop 'error while reading MMparam s12'
  read(token, *) ljp%s12
  if (getstr(linein, token) < 0) stop 'error while reading MMparam s13'
  read(token, *) ljp%s13
  if (getstr(linein, token) < 0) stop 'error while reading MMparam s22'
  read(token, *) ljp%s22
  if (getstr(linein, token) < 0) stop 'error while reading MMparam s23'
  read(token, *) ljp%s23
  if (getstr(linein, token) < 0) stop 'error while reading MMparam s33'
  read(token, *) ljp%s33

  return
end subroutine

!------------------------------------------------------------------------------
subroutine get_tokens_and_append_shift_per_atomtype(linein, shifts)
!------------------------------------------------------------------------------
  character(len=:),allocatable,intent(in out) :: linein
  type(shift),allocatable,intent(in out) :: shifts(:)
  type(shift) :: s

  if (getstr(linein, token) < 0) stop 'error while reading atom type for shift op'
  s%atype = token
  if (getstr(linein, token) < 0) stop 'error while reading shift value'
  read(token, *) s%zshift

  ! allocate zero-sized array
  if(.not.allocated(shifts)) allocate(shifts(0)) 

  shifts = [shifts, s]

  return
end subroutine

!------------------------------------------------------------------------------
subroutine get_tokens_and_append_remove_domain(linein, removes)
!------------------------------------------------------------------------------
  character(len=:),allocatable,intent(in out) :: linein
  type(remove),allocatable,intent(in out) :: removes(:)
  type(remove) :: remove

  if (getstr(linein, token) < 0) stop 'error while reading remove xmin'
  read(token, *) remove%xmin
  if (getstr(linein, token) < 0) stop 'error while reading remove xmax'
  read(token, *) remove%xmax
  if (getstr(linein, token) < 0) stop 'error while reading remove ymin'
  read(token, *) remove%ymin
  if (getstr(linein, token) < 0) stop 'error while reading remove ymax'
  read(token, *) remove%ymax
  if (getstr(linein, token) < 0) stop 'error while reading remove zmin'
  read(token, *) remove%zmin
  if (getstr(linein, token) < 0) stop 'error while reading remove zmax'
  read(token, *) remove%zmax

  ! allocate zero-sized array
  if(.not.allocated(removes)) allocate(removes(0)) 

  removes = [removes, remove]

  return
end subroutine

!------------------------------------------------------------------------------
subroutine get_tokens_and_append_stripe_domain(linein, stripes)
!------------------------------------------------------------------------------
  character(len=:),allocatable,intent(in out) :: linein
  type(stripe),allocatable,intent(in out) :: stripes(:)
  type(stripe) :: stripe 

  if (getstr(linein, token) < 0) stop 'error while reading stripe domain type'
  stripe%dtype = token
  if (getstr(linein, token) < 0) stop 'error while reading stripe domain dir'
  read(token, *) stripe%dir
  if (getstr(linein, token) < 0) stop 'error while reading stripe domain shift'
  read(token, *) stripe%shift
  if (getstr(linein, token) < 0) stop 'error while reading num_stripe domains'
  read(token, *) stripe%num_stripes

  ! allocate zero-sized array
  if(.not.allocated(stripes)) allocate(stripes(0)) 

  stripes = [stripes, stripe]

  return
end subroutine

!------------------------------------------------------------------------------
subroutine get_tokens_and_append_rectangular_domain(linein, rects)
!------------------------------------------------------------------------------
  character(len=:),allocatable,intent(in out) :: linein
  type(rectangular),allocatable,intent(in out) :: rects(:)
  type(rectangular) :: rect

  if (getstr(linein, token) < 0) stop 'error while reading rec domain type'
  rect%dtype = token
  if (getstr(linein, token) < 0) stop 'error while reading rectangular domain xmin'
  read(token, *) rect%xmin
  if (getstr(linein, token) < 0) stop 'error while reading rectangular domain ymin'
  read(token, *) rect%ymin

  if (getstr(linein, token) < 0) stop 'error while reading rectangular domain xmax'
  read(token, *) rect%xmax
  if (getstr(linein, token) < 0) stop 'error while reading rectangular domain ymax'
  read(token, *) rect%ymax

  ! allocate zero-sized array
  if(.not.allocated(rects)) allocate(rects(0)) 

  rects = [rects, rect]

  return
end subroutine

!------------------------------------------------------------------------------
subroutine get_tokens_and_append_circular_domain(linein, circs)
!------------------------------------------------------------------------------
  character(len=:),allocatable,intent(in out) :: linein
  type(circular),allocatable,intent(in out) :: circs(:)
  type(circular) :: circ

  if (getstr(linein, token) < 0) stop 'error while reading domain type'
  circ%dtype = token
  if (getstr(linein, token) < 0) stop 'error while reading domain center(1)'
  read(token, *) circ%center(1)
  if (getstr(linein, token) < 0) stop 'error while reading domain center(2)'
  read(token, *) circ%center(2)
  if (getstr(linein, token) < 0) stop 'error while reading radius'
  read(token, *) circ%radius

  ! allocate zero-sized array
  if(.not.allocated(circs)) allocate(circs(0)) 

  circs = [circs, circ]

  return
end subroutine

!------------------------------------------------------------------------------
function layer_ctor(ltype, num_ucells, ucell, zmin, zmax) result(c)
!------------------------------------------------------------------------------
  integer,intent(in) :: ltype, num_ucells
  real(8),intent(in) :: zmin, zmax
  type(unitcell_type) :: ucell

  type(layer) :: c

  allocate(c%circs(0), c%rects(0))

  c%ltype = ltype
  c%num_ucells = num_ucells
  c%ucell = ucell
  c%zmin = zmin
  c%zmax = zmax

  return
end function

!------------------------------------------------------------------------------
subroutine get_tokens_and_append_nnmm_layers(linein, lx, ly, layers)
!------------------------------------------------------------------------------
  character(len=:),allocatable,intent(in out) :: linein
  integer,intent(in out) :: lx, ly
  type(layer),allocatable,intent(in out) :: layers(:)

  integer :: i, ltype, num_ucells
  type(layer) :: c

  real(8) :: zlength = 0.d0

  if(.not.allocated(layers)) allocate(layers(0)) 

  do while (getstr(linein, token) > 0) 
     select case(token)
       case("XY")
         if (getstr(linein, token) < 0) stop 'error while reading X unitcells in layer'
         read(token, *) lx
         if (getstr(linein, token) < 0) stop 'error while reading Y unitcells in layer'
         read(token, *) ly

       case("MM")
         if (getstr(linein, token) < 0) stop 'error while reading MM layer'
         read(token, *) num_ucells
         
         c = layer_ctor(ltype=MM_LAYER_TYPE, num_ucells=num_ucells, ucell=MM, &
                 zmin=zlength, zmax=zlength+num_ucells*MM%lattice(3))
         layers = [layers, c]

         zlength = c%zmax

       case("NN")
         if (getstr(linein, token) < 0) stop 'error while reading NN layer'
         read(token, *) num_ucells

         c = layer_ctor(ltype=NN_LAYER_TYPE, num_ucells=num_ucells, ucell=NN, &
                 zmin=zlength, zmax=zlength+num_ucells*NN%lattice(3))
         layers = [layers, c]

         zlength = c%zmax
     end select
  enddo

  !print*,'XY: ', lx, ly
  !do i=1, size(layers)
  !   print'(a,2i6,2f10.3)','layers: ', layers(i)%ltype, layers(i)%num_ucells, layers(i)%zmin, layers(i)%zmax
  !enddo

end subroutine

!-------------------------------------------------------------------------------------------
function nnmm_params_ctor(filename) result(dp)
!-------------------------------------------------------------------------------------------
!character(len=:),allocatable,intent(in) :: filename
character(*) :: filename
character(256) :: linein0
character(len=:),allocatable :: linein

integer :: iunit, lid

type(nnmm_params) :: dp

! box dimensions are used in the layers construction
call nnmm_set_unitcells(NN,MM)

allocate(dp%shifts(0), dp%layers(0))

open(newunit=iunit, file=filename, form='formatted', status='old')

do while (.true.)
  read(iunit,'(a)',end=10) linein0
  linein = trim(adjustl(linein0))

  if (getstr(linein, token) > 0) then

     if(token=='MM') call get_tokens_and_append_MM_params(linein, dp%MM_params)
     if(token=='shift') call get_tokens_and_append_shift_per_atomtype(linein, dp%shifts)

     if(token=='layers') call get_tokens_and_append_nnmm_layers(linein, dp%lx, dp%ly, dp%layers)

     if(token=='remove') then
       if (getstr(linein, token) < 0) stop 'error while reading layer index for remove'
       read(token, *) lid
       call get_tokens_and_append_remove_domain(linein, dp%layers(lid)%removes)
     endif

     if(token=='stripe') then
       if (getstr(linein, token) < 0) stop 'error while reading layer index for stripe domain'
       read(token, *) lid
       call get_tokens_and_append_stripe_domain(linein, dp%layers(lid)%stripes)
     endif

     if(token=='circular'.or.token=='circ') then
       if (getstr(linein, token) < 0) stop 'error while reading layer index for circular domain'
       read(token, *) lid
       call get_tokens_and_append_circular_domain(linein, dp%layers(lid)%circs)
     endif

     if(token=='rectangular'.or.token=='rec') then
       if (getstr(linein, token) < 0) stop 'error while reading layer index for rectangular domain'
       read(token, *) lid
       call get_tokens_and_append_rectangular_domain(linein, dp%layers(lid)%rects)
     endif
  endif

end do
10 close(iunit)

if(myid==0) call dp.show()

end function

end module
