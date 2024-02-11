module fileio

  use mpi_mod
  use utils, only : l2g, get_boxparameters
  use base
  use atoms
  use memory_allocator_mod

type :: xyz_aggregator
  logical :: use_aggregator = .false.
  integer :: num_stack, num_total_atoms
  integer :: stack_counter = 0
  real(8),allocatable :: buf_pos(:), buf_type(:)
  contains
    procedure :: init => xyz_aggregator_init 
    procedure :: gather => xyz_aggregator_gather_atom_data
    procedure :: save => xyz_aggregator_save_into_file
    procedure :: reset => xyz_aggregator_reset
end type

type(xyz_aggregator) :: xyz_agg

contains

!----------------------------------------------------------------------------------------
subroutine xyz_aggregator_init(this, num_total_atoms, num_stack)
!----------------------------------------------------------------------------------------
  class(xyz_aggregator),intent(in out) :: this
  integer,intent(in) :: num_stack
  integer(8),intent(in) :: num_total_atoms

  this%stack_counter = 0

  this%use_aggregator = .true.
  this%num_stack = num_stack
  this%num_total_atoms = num_total_atoms

  allocate(this%buf_pos(3*this%num_total_atoms*this%num_stack))
  allocate(this%buf_type(this%num_total_atoms*this%num_stack))
  this%buf_pos=0.d0; this%buf_type=0.d0

  if(myid==0) then
    print*,'stack_counter,use_aggregator,num_stack,num_total_atoms: ', &
          this%stack_counter,this%use_aggregator,this%num_stack,this%num_total_atoms
    print*,'shape(buf_pos), shape(buf_type): ',shape(this%buf_pos), shape(this%buf_type)
  endif

  return
end subroutine

!----------------------------------------------------------------------------------------
subroutine xyz_aggregator_gather_atom_data(this,num_atoms,atype,pos)
!----------------------------------------------------------------------------------------
  class(xyz_aggregator),intent(in out) :: this
  integer,intent(in) :: num_atoms
  real(8),allocatable,intent(in) :: atype(:)
  real(8),allocatable,intent(in) :: pos(:,:)

  integer :: i,ii,ii3,igid,offset 

  offset = this%stack_counter*this%num_total_atoms

  print*,'entering gather(): ', this%stack_counter,offset, size(this%buf_type)

  do i=1, num_atoms
     igid = l2g(atype(i))
     ii = offset+igid  ! offset is zero-indexed, but igid is 1-indexed
     this%buf_type(ii) = atype(i)

     ii3 = ii*3 ! 1-indexed multiplied by 3, so the 1st element is ii3-2
     this%buf_pos(ii3-2:ii3) = pos(i,1:3)
  enddo

  call mpi_allreduce(mpi_in_place, this%buf_pos(3*offset+1), 3*this%num_total_atoms, &
                     mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
  call mpi_allreduce(mpi_in_place, this%buf_type(offset+1), this%num_total_atoms, &
                     mpi_double_precision, mpi_sum, mpi_comm_world, ierr)

  this%stack_counter = this%stack_counter + 1

  return
end subroutine

!----------------------------------------------------------------------------------------
subroutine xyz_aggregator_save_into_file(this,myrank,filename)
!----------------------------------------------------------------------------------------
  class(xyz_aggregator),intent(in out) :: this
  character(len=:),allocatable,intent(in) :: filename
  integer,intent(in) :: myrank

  integer :: i,ii,ii3, istack, iunit, offset

  if(myrank/=0) return  ! only master rank save data

  open(newunit=iunit,file=filename,form='formatted')

  do istack=0, this%num_stack-1

    write(iunit,'(i9)') this%num_total_atoms
    write(iunit,'(3f12.5,3f8.3,a,i5)')  lata,latb,latc,lalpha,lbeta,lgamma, ' stack ', istack

    offset = istack*this%num_total_atoms

    do i=1, this%num_total_atoms
      ii = offset + i
      ii3 = ii*3
       
      ity = nint(this%buf_type(ii))
      !print*,'i,ity,atmname(ity),buf_pos(ii3-2:ii)',istack,i,ity,atmname(ity),this%buf_type(ii),this%buf_pos(ii3-2:ii3)

      write(iunit,'(a,3f10.5,i6)') atmname(ity), this%buf_pos(ii3-2:ii3), l2g(this%buf_type(ii))
    enddo

  enddo

  close(iunit)

end subroutine

!----------------------------------------------------------------------------------------
subroutine xyz_aggregator_reset(this)
!----------------------------------------------------------------------------------------
  class(xyz_aggregator),intent(in out) :: this

  ! reset everything
  this%stack_counter = 0
  this%buf_pos = 0.d0
  this%buf_type = 0.d0

end subroutine

!----------------------------------------------------------------------------------------
subroutine OUTPUT(fileNameBase, atype, pos, v, q, f, energy)
!----------------------------------------------------------------------------------------
implicit none

character(len=:),allocatable,intent(in) :: fileNameBase
real(8),allocatable,intent(in) :: atype(:), q(:)
real(8),allocatable,intent(in) :: pos(:,:),v(:,:)
real(8),allocatable,intent(in),optional :: f(:,:)
real(8),optional :: energy

character(len=:),allocatable :: filename

integer :: idx

if(xyz_num_stack>1) then

  call xyz_agg%gather(NATOMS,atype,pos)

  if(xyz_agg%stack_counter == xyz_agg%num_stack) then

    if(isBinary) call WriteBIN(atype,pos,v,q,fileNameBase)

    filename = fileNameBase//"-agg.xyz" 
    call xyz_agg%save(myid, filename)
    call xyz_agg%reset()

  endif

else

  if(isBinary) call WriteBIN(atype,pos,v,q,fileNameBase)
  if(isBondFile) call WriteBND(fileNameBase)
  if(isPDB) call WritePDB(fileNameBase)
  if(isXYZ) then
     if( find_cmdline_argc('--xyz_pto',idx)) then
       call WriteXYZ_PTO(fileNameBase, atype, pos, v, q, f)
     else
       call WriteXYZ(fileNameBase, atype, pos, v, q, f, energy)
     endif
  endif

endif

return

Contains 

!--------------------------------------------------------------------------
subroutine WriteBND(fileNameBase)
!--------------------------------------------------------------------------
implicit none

character(*),intent(in) :: fileNameBase

integer :: i, ity, j, j1, jty, m
real(8) :: bndordr(MAXNEIGHBS)
integer :: igd,jgd,bndlist(0:MAXNEIGHBS)

integer (kind=MPI_OFFSET_KIND) :: offset
integer (kind=MPI_OFFSET_KIND) :: fileSize
integer :: localDataSize
integer :: fh ! file handler

integer :: BNDLineSize, baseCharPerAtom
integer,parameter :: MaxBNDLineSize=4096
character(MaxBNDLineSize) :: BNDOneLine
real(8),parameter :: BNDcutoff=0.3d0

character(len=:),allocatable :: BNDAllLines

integer :: scanbuf

integer :: ti,tj,tk
call system_clock(ti,tk)

! precompute the total # of neighbors
m=0
do i=1, NATOMS
   do j1 = 1, nbrlist(i,0)
!--- don't count if BO is less than BNDcutoff.
       if(BO(0,i,j1) > BNDcutoff) then 
           m=m+1
       endif
   enddo
enddo

200 format(i12.12,1x,3f12.3,1x,2i3,20(1x,i12.12,f6.3)) 

! get local datasize based on above format and the total # of neighbors
baseCharPerAtom=12+1+3*12+1+2*3 +1 ! last 1 for newline
localDataSize=NATOMS*(baseCharPerAtom)+m*(1+12+6)

if( (baseCharPerAtom+MAXNEIGHBS*(1+12+6)) > MaxBNDLineSize) then
    print'(a,i6,2i12)', 'ERROR: MaxBNDLineSize is too small @ WriteBND', &
                    myid, baseCharPerAtom+MAXNEIGHBS*(1+12+6), MaxBNDLineSize
endif

call MPI_File_Open(MPI_COMM_WORLD,trim(fileNameBase)//".bnd", &
    MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

! offset will point the end of local write after the scan
call MPI_Scan(localDataSize,scanbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

! since offset is MPI_OFFSET_KIND and localDataSize is integer, use an integer as buffer
offset=scanbuf

! nprocs-1 rank has the total data size
call MPI_Bcast(scanbuf,1,MPI_INTEGER,nprocs-1,MPI_COMM_WORLD,ierr)
fileSize=scanbuf

call MPI_File_set_size(fh, fileSize, ierr)

! set offset at the beginning of the local write
offset=offset-localDataSize

call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)

allocate(character(len=localDataSize) :: BNDAllLines)
BNDALLLines=""

BNDLineSize=0
do i=1, NATOMS
   ity = nint(atype(i))
!--- get global ID for i-atom
   igd = l2g(atype(i))

!--- count the number bonds to be shown.
   bndlist(0)=0
   do j1 = 1, nbrlist(i,0)
      j = nbrlist(i,j1)
      jty = nint(atype(j))

!--- get global ID for j-atom
      jgd = l2g(atype(j))

!--- if bond order is less than 0.3, ignore the bond.
      if( BO(0,i,j1) < 0.3d0) cycle

      bndlist(0) = bndlist(0) + 1
      bndlist(bndlist(0)) = jgd
      bndordr(bndlist(0)) = BO(0,i,j1)
   enddo

   BNDOneLine=""
   write(BNDOneLine,200) igd, pos(i,1:3),nint(atype(i)),bndlist(0), &
         (bndlist(j1),bndordr(j1),j1=1,bndlist(0))

   ! remove space and add new_line
   BNDOneLine=trim(adjustl(BNDOneLine))//NEW_LINE('A')
   BNDLineSize=BNDLineSize+len(trim(BNDOneLine))

   BNDAllLines=trim(BNDAllLines)//trim(BNDOneLine)
enddo

if(localDataSize>0) then
   call MPI_File_Write(fh,BNDAllLines,localDataSize, &
        MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
endif

deallocate(BNDAllLines) 

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
call MPI_File_Close(fh,ierr)

call system_clock(tj,tk)
it_timer(20)=it_timer(20)+(tj-ti)

return
end subroutine

!--------------------------------------------------------------------------
subroutine WritePDB(fileNameBase)
!--------------------------------------------------------------------------
implicit none

character(*),intent(in) :: fileNameBase

integer :: i, ity, igd 
real(8) :: tt=0.d0, ss=0.d0

integer (kind=MPI_OFFSET_KIND) :: offset
integer (kind=MPI_OFFSET_KIND) :: fileSize
integer :: localDataSize
integer :: fh ! file handler

integer,parameter :: PDBLineSize=67
character(PDBLineSize) :: PDBOneLine

character(len=:),allocatable :: PDBAllLines

integer :: scanbuf

integer :: ti,tj,tk
call system_clock(ti,tk)

! get local datasize
localDataSize=NATOMS*PDBLineSize

call MPI_File_Open(MPI_COMM_WORLD,trim(fileNameBase)//".pdb", &
     MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

! offset will point the end of local write after the scan
call MPI_Scan(localDataSize,scanbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

! since offset is MPI_OFFSET_KIND and localDataSize is integer, use an integer as buffer
offset=scanbuf

! nprocs-1 rank has the total data size
call MPI_Bcast(scanbuf,1,MPI_INTEGER,nprocs-1,MPI_COMM_WORLD,ierr)
fileSize=scanbuf

call MPI_File_set_size(fh, fileSize, ierr)

! set offset at the beginning of the local write
offset=offset-localDataSize

allocate(character(len=localDataSize) :: PDBAllLines)
PDBAllLines=""

call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)

do i=1, NATOMS

  ity = nint(atype(i))
!--- calculate atomic temperature 
  tt = hmas(ity)*sum(v(i,1:3)*v(i,1:3))
  tt = tt*UTEMP*1d-2 !scale down to use two decimals in PDB format 

!--- sum up diagonal atomic stress components 
  ss = sum(astr(1:3))/3.d0*USTRS

  tt = q(i)

  igd = l2g(atype(i))
  write(PDBOneLine,100)'ATOM  ',0, atmname(ity), igd, pos(i,1:3), tt, ss

  PDBOneLine(PDBLineSize:PDBLineSize)=NEW_LINE('A')
  PDBAllLines=trim(PDBAllLines)//trim(PDBOneLine)

enddo

if(localDataSize>0) then
    call MPI_File_Write(fh,PDBAllLines,localDataSize, &
         MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
endif

deallocate(PDBAllLines)

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
call MPI_File_Close(fh,ierr)

100 format(A6,I5,1x,A2,i12,4x,3f8.3,f6.2,f6.2)

call system_clock(tj,tk)
it_timer(21)=it_timer(21)+(tj-ti)


end subroutine

!--------------------------------------------------------------------------
subroutine WriteXYZ_PTO(fileNameBase, atype, pos, v, q, f)
!--------------------------------------------------------------------------
implicit none
real(8),allocatable,intent(in) :: atype(:), q(:)
real(8),allocatable,intent(in) :: pos(:,:),v(:,:)
real(8),allocatable,intent(in),optional :: f(:,:)

character(*),intent(in) :: fileNameBase

integer :: i, ity, idx1, idx0, igd

integer (kind=MPI_OFFSET_KIND) :: offset
integer (kind=MPI_OFFSET_KIND) :: fileSize
integer :: localDataSize
integer :: fh ! file handler

integer :: OneLineSize, MetaDataSize
character(60) :: a60
character(len=:),allocatable :: OneLine,AllLines

integer :: scanbuf

integer :: ti,tj,tk

integer :: j1,j,jty,nTi, ibond
real(8),allocatable :: dpol(:,:)
real(8) :: centroid(3), dr(3), drsq

integer,parameter :: mxbond = 6
real(8),parameter :: rcut = 3.d0, rcutsq = rcut*rcut

integer :: octa_parity(2)
integer,allocatable :: rotation_flag(:)
real(8),parameter :: aaxis = 3.902d0

integer(8) :: GNATOMS_Ti

call system_clock(ti,tk)

allocate(dpol(NATOMS,3))
dpol = 0.d0
allocate(rotation_flag(NATOMS))
rotation_flag = 0

nTi=0
do i=1, NATOMS !i atom loop

   ity = nint(atype(i))

   if(atmname(ity) == "Ti") then ! if i atom is Ti compute polarization

      nTi=nTi+1
      ibond=0
      centroid = 0.0d0 !reset centroid positions

      octa_parity(1:2)=nint(pos(i,1:2)/aaxis)

      do j1 = 1, nbrlist(i,0) !j atom loop

         j = nbrlist(i,j1)
         jty = nint(atype(j))

         if(atmname(jty) == "O") then !is atom j O ?
           dr(1:3) = pos(i,1:3) - pos(j,1:3)
           drsq = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)

           if(drsq < rcutsq) then ! Am I bonded ?
             centroid(1:3) = centroid(1:3) + pos(j,1:3)
             ibond=ibond+1

             if(dr(2)>1d0.and.abs(dr(1))>0.3d0) rotation_flag(i) = nint(sign(1d0,dr(1)))
           endif 

         endif 

      enddo !end j loop


      if(ibond<mxbond) print'(a,3i9)', "WARNING Broken Symmetry LT 6 Ti-O Bonds: myid,ibond,mxbond ", myid, ibond, mxbond
      if(ibond>mxbond) print'(a,3i9)', "WARNING Broken Symmetry GT 6 Ti-O Bonds: myid,ibond,mxbond ", myid, ibond, mxbond

      centroid(1:3) = centroid(1:3)/dble(ibond)
      dpol(i,1:3) = pos(i,1:3) - centroid(1:3)

   endif !end Ti check

enddo !end i atom loop

! GNATOMS, 6 lattice parameters, two newlines.
MetaDataSize = 9 + 60 + 2
write(a60,'(3f12.5,3f8.3)')  lata,latb,latc,lalpha,lbeta,lgamma

! name + pos(i,1:3) + gid + dpol(i,1:3) + rot_flag + rot_angle + newline
OneLineSize = 3 + 36 + 9 + 36 + 2 + 12 + 1 

! get local datasize
if(find_cmdline_argc('--xyz_pto_tionly',idx)) then
  localDataSize=nTi*OneLineSize
else
  localDataSize=NATOMS*OneLineSize
endif
if(myid==0) localDataSize = localDataSize + MetaDataSize

call MPI_File_Open(MPI_COMM_WORLD,trim(fileNameBase)//".xyz", &
     MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

! offset will point the end of local write after the scan
call MPI_Scan(localDataSize,scanbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

! since offset is MPI_OFFSET_KIND and localDataSize is integer, use an integer as buffer
offset=scanbuf

! nprocs-1 rank has the total data size
call MPI_Bcast(scanbuf,1,MPI_INTEGER,nprocs-1,MPI_COMM_WORLD,ierr)
fileSize=scanbuf

call MPI_File_set_size(fh, fileSize, ierr)

! set offset at the beginning of the local write
offset=offset-localDataSize

call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)

allocate(character(OneLineSize) :: OneLine)
allocate(character(localDataSize) :: AllLines)

GNATOMS_Ti = nTi
call MPI_ALLREDUCE(MPI_IN_PLACE, GNATOMS_Ti, 1, MPI_INTEGER8, MPI_SUM,  MPI_COMM_WORLD, ierr)

idx0=1 ! allline index
if(myid==0) then

  if(find_cmdline_argc('--xyz_pto_tionly',idx)) then
    write(AllLines(idx0:idx0+8),'(i9)') GNATOMS_Ti; idx0=idx0+9
  else
    write(AllLines(idx0:idx0+8),'(i9)') GNATOMS; idx0=idx0+9
  endif

  write(AllLines(idx0:idx0),'(a1)') new_line('A'); idx0=idx0+1
  write(AllLines(idx0:idx0+59),'(a60)') a60; idx0=idx0+60
  write(AllLines(idx0:idx0),'(a1)') new_line('A'); idx0=idx0+1
endif

do i=1, NATOMS
  idx1 = 1 ! oneline index

  ity = nint(atype(i))
  igd = l2g(atype(i))

  if(find_cmdline_argc('--xyz_pto_tionly',idx) .and. atmname(ity) /= "Ti") cycle

! element name
  write(OneLine(idx1:idx1+2),'(a3)') atmname(ity); idx1=idx1+3

! atom positions. with high precision for dielectric calculation
  write(OneLine(idx1:idx1+35),'(3f12.5)') pos(i,1:3); idx1=idx1+36

! global Id
  write(OneLine(idx1:idx1+8),'(i9)') igd; idx1=idx1+9

  write(OneLine(idx1:idx1+35),'(3es12.4)') dpol(i,1:3); idx1=idx1+36

! rotation flag
  octa_parity(1:2) = nint(pos(i,1:2)/aaxis)
  if( rotation_flag(i) == 0 ) then
    write(OneLine(idx1:idx1+1),'(i2)') 0; idx1=idx1+2
  else
    if( mod((octa_parity(1) + octa_parity(2)),2)*2 - 1 == rotation_flag(i) ) then
      write(OneLine(idx1:idx1+1),'(i2)') 1; idx1=idx1+2
    else
      write(OneLine(idx1:idx1+1),'(i2)') 2; idx1=idx1+2
    endif 
  endif

! rotation angle in x&y plane. use dr(3) to store the norm of dpol in xy 
  dr(1:3) = dpol(i,1:3)
  dr(3) = sqrt(sum(dr(1:2)*dr(1:2)))
  if(dr(3) < 1d-6) then 
     write(OneLine(idx1:idx1+11),'(es12.4)') 0.d0; idx1=idx1+12
  else
     dr(1:2) = dr(1:2)/dr(3) 
     if(dr(2) > 0.d0) then
        write(OneLine(idx1:idx1+11),'(es12.4)') acos(dr(1)); idx1=idx1+12
     else
        write(OneLine(idx1:idx1+11),'(es12.4)') 3.1415926d0+acos(-1d0*dr(1)); idx1=idx1+12
     endif
  endif



  write(OneLine(idx1:idx1),'(a1)') new_line('A'); idx1=idx1+1

  write(AllLines(idx0:idx0+OneLineSize-1),'(a)') OneLine; idx0=idx0+OneLineSize

enddo

if(localDataSize>0) then
    call MPI_File_Write(fh,AllLines,localDataSize, &
         MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
endif

deallocate(AllLines, OneLine)
deallocate(dpol, rotation_flag)

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
call MPI_File_Close(fh,ierr)

call system_clock(tj,tk)
it_timer(21)=it_timer(21)+(tj-ti)

end subroutine


!--------------------------------------------------------------------------
subroutine WriteXYZ(fileNameBase, atype, pos, v, q, f, energy)
!--------------------------------------------------------------------------
implicit none
real(8),allocatable,intent(in) :: atype(:), q(:)
real(8),allocatable,intent(in) :: pos(:,:),v(:,:)
real(8),allocatable,intent(in),optional :: f(:,:)
real(8),intent(in),optional :: energy

character(*),intent(in) :: fileNameBase

integer :: i, ity, idx1, idx0 , igd

integer (kind=MPI_OFFSET_KIND) :: offset
integer (kind=MPI_OFFSET_KIND) :: fileSize
integer :: localDataSize
integer :: fh ! file handler

integer :: OneLineSize, MetaDataSize
character(256) :: a256
character(len=:),allocatable :: OneLine,AllLines

integer :: scanbuf

integer :: ti,tj,tk

character(len=:),allocatable :: header
integer :: header_size

call system_clock(ti,tk)


write(a256, '(a10,f0.2,8f10.2,a2 $)') 'Lattice="', lata,0d0,0d0, 0d0,latb,0d0, 0d0,0d0,latc, '" '
header=adjustl(trim(a256))

if(present(energy)) then
  write(a256, '(a10,f0.5 $)') 'energy=', energy
  header=header//adjustl(trim(a256))
endif
header=header//'Properties=species:S:1:pos:R:3:charge:R:1:forces:R:3:id:I:1 pbc="T T T" config_type=md'
header_size=len(adjustl(trim(header)))
!print*,header_size, header

MetaDataSize = 9 + header_size + 2

if(isPQEq) then
  ! name + pos(i,1:3) + q(i) + gid + spos(i,1:3) + newline
  OneLineSize = 3 + 60 + 20 + 60 + 9 +  1 
else if(present(f)) then
  ! name + pos(i,1:3) + q(i) + gid + f(i,1:3) + newline
  OneLineSize = 3 + 60 + 20 + 60 + 9 + 1 
else
  ! name + pos(i,1:3) + q(i) + gid + newline
  OneLineSize = 3 + 36 + 8 + 9 + 1 
endif

! get local datasize
localDataSize=NATOMS*OneLineSize
if(myid==0) localDataSize = localDataSize + MetaDataSize

call MPI_File_Open(MPI_COMM_WORLD,trim(fileNameBase)//".xyz", &
     MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

! offset will point the end of local write after the scan
call MPI_Scan(localDataSize,scanbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

! since offset is MPI_OFFSET_KIND and localDataSize is integer, use an integer as buffer
offset=scanbuf

! nprocs-1 rank has the total data size
call MPI_Bcast(scanbuf,1,MPI_INTEGER,nprocs-1,MPI_COMM_WORLD,ierr)
fileSize=scanbuf

call MPI_File_set_size(fh, fileSize, ierr)

! set offset at the beginning of the local write
offset=offset-localDataSize

call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)

allocate(character(OneLineSize) :: OneLine)
allocate(character(localDataSize) :: AllLines)

idx0=1 ! allline index
if(myid==0) then
  write(AllLines(idx0:idx0+8),'(i9)') GNATOMS; idx0=idx0+9
  write(AllLines(idx0:idx0),'(a1)') new_line('A'); idx0=idx0+1
  write(AllLines(idx0:idx0+header_size-1),'(a)') header; idx0=idx0+header_size
  write(AllLines(idx0:idx0),'(a1)') new_line('A'); idx0=idx0+1
endif

do i=1, NATOMS
  idx1 = 1 ! oneline index

  ity = nint(atype(i))
  igd = l2g(atype(i))

! element name
  write(OneLine(idx1:idx1+2),'(a3)') atmname(ity); idx1=idx1+3

! atom positions. with high precision for dielectric calculation
  if(isPQEq.or.present(f)) then
     write(OneLine(idx1:idx1+59),'(3es20.12)') pos(i,1),pos(i,2),pos(i,3); idx1=idx1+60
     write(OneLine(idx1:idx1+19),'(es20.12)') q(i); idx1=idx1+20
  else
     write(OneLine(idx1:idx1+35),'(3f12.5)') pos(i,1:3); idx1=idx1+36
     write(OneLine(idx1:idx1+7),'(3f8.3)') q(i); idx1=idx1+8
  endif

! shell charge positions with high precision for dielectric calculation
  if(isPQEq) then
     write(OneLine(idx1:idx1+59),'(3es20.12)') spos(i,1:3); idx1=idx1+60
  else if(present(f)) then
     write(OneLine(idx1:idx1+59),'(3es20.12)') f(i,1),f(i,2),f(i,3); idx1=idx1+60
  endif

! global Id
  write(OneLine(idx1:idx1+8),'(i9)') igd; idx1=idx1+9

  write(OneLine(idx1:idx1),'(a1)') new_line('A'); idx1=idx1+1

  write(AllLines(idx0:idx0+OneLineSize-1),'(a)') OneLine; idx0=idx0+OneLineSize

enddo

if(localDataSize>0) then
    call MPI_File_Write(fh,AllLines,localDataSize, &
         MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
endif

deallocate(AllLines, OneLine)

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
call MPI_File_Close(fh,ierr)

call system_clock(tj,tk)
it_timer(21)=it_timer(21)+(tj-ti)


end subroutine

end subroutine OUTPUT

!--------------------------------------------------------------------------
subroutine ReadXYZ(atype, rreal, v, q, f, fileName)
!--------------------------------------------------------------------------
implicit none

character(*),intent(in) :: fileName
real(8),allocatable,dimension(:),intent(inout) :: atype,q
real(8),allocatable,dimension(:,:),intent(inout) :: rreal,v,f

real(8) :: dbuf6(6), mat(3,3), mati(3,3)

integer :: ti,tj,tk

integer :: i,j,i1, ntot, num_unit, funit
integer(8) :: iigd
real(8),allocatable :: pos0(:)
integer,allocatable :: atype0(:)
character(2) :: c2

call system_clock(ti,tk)

if(myid==0) then
  open(newunit=funit, file=trim(filename), status='old', form='formatted')

  read(funit,fmt=*) num_unit
  allocate(atype0(num_unit), pos0(3*num_unit))

  read(funit,fmt=*) dbuf6(1:6)

  do i = 1, num_unit
     read(funit,fmt=*) c2,pos0(i*3-2:i*3) 

     do j=1, 3
        call assert( pos0((i-1)*3+j)>0.d0, 'ERROR: negative coordinate found in ReadXYZ(): '//trim(filename))
     enddo

!TODO: use atmname to get the element-to-integer mapping.
!      atmname has to be set based on the input files, e.g. ffield or fnn.in.
     do i1=1, size(atmname)
       if(c2(1:2) == atmname(i1)) exit
     enddo
     atype0(i)=i1
  enddo

  close(funit)
endif

call MPI_Bcast(num_unit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if(myid>0) allocate(pos0(3*num_unit), atype0(num_unit))

call MPI_Bcast(dbuf6,size(dbuf6),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
lata=dbuf6(1); latb=dbuf6(2); latc=dbuf6(3)
lalpha=dbuf6(4); lbeta=dbuf6(5); lgamma=dbuf6(6)

call MPI_Bcast(atype0,size(atype0),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast(pos0,size(pos0),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!--- allocate arrays
if(.not.allocated(atype)) call allocator(atype,1,NBUFFER)
if(.not.allocated(q)) call allocator(q,1,NBUFFER)
if(.not.allocated(rreal)) call allocator(rreal,1,NBUFFER,1,3)
if(.not.allocated(v)) call allocator(v,1,NBUFFER,1,3)
if(.not.allocated(f)) call allocator(f,1,NBUFFER,1,3)
if(.not.allocated(qsfp)) call allocator(qsfp,1,NBUFFER)
if(.not.allocated(qsfv)) call allocator(qsfv,1,NBUFFER)
f(:,:)=0.0

iigd = num_unit*myid ! for global ID
do i=1,num_unit
   pos(i,1:3) = pos0(i*3-2:i*3)
   atype(i) = dble(atype0(i)) + (iigd+i)*1d-13
enddo 

!--- assuming XYZ format is in the real coordinates given for one domain. 
!    To run multiple domain job, normalize the coords, update H-matrix, 
!    update the box origin, then scale back to the real coords. 
call get_boxparameters(mat,lata,latb,latc,lalpha,lbeta,lgamma)
call matinv(mat,mati)
call xu2xs_inplace(mati,[0.d0,0.d0,0.d0],num_unit,pos)

NATOMS=num_unit

!--- update to glocal cell parameters
lata=lata*vprocs(1)
latb=latb*vprocs(2)
latc=latc*vprocs(3)

call get_boxparameters(mat,lata,latb,latc,lalpha,lbeta,lgamma)
do i=1, 3
do j=1, 3
   HH(i,j,0)=mat(i,j)
enddo; enddo
call update_box_params(vprocs, vid, hh, lata, latb, latc, maxrc, &
                       cc, lcsize, hhi, mdbox, lbox, obox)

call xs2xu_inplace(hh,obox,num_unit,pos)

call system_clock(tj,tk)
it_timer(22)=it_timer(22)+(tj-ti)

return
end

!--------------------------------------------------------------------------
subroutine ReadPSTO(atype, rreal, v, q, f, fileName)
!--------------------------------------------------------------------------
implicit none

character(*),intent(in) :: fileName
real(8),allocatable,dimension(:),intent(inout) :: atype,q
real(8),allocatable,dimension(:,:),intent(inout) :: rreal,v,f

integer :: i,i1

integer (kind=MPI_OFFSET_KIND) :: offset, offsettmp
integer (kind=MPI_OFFSET_KIND) :: fileSize
integer :: localDataSize, metaDataSize, scanbuf
integer :: fh ! file handler

integer :: nmeta
integer,allocatable :: idata(:)
real(8),allocatable :: dbuf(:)
real(8) :: ddata(6), d10(10)

real(8) :: rnorm(NBUFFER,3), mat(3,3)
integer :: j

!=== # of unit cells ===
!integer :: mx=4,my=4,mz=4
!integer :: mx=20,my=20,mz=20
!integer :: mx=50,my=50,mz=50
integer :: mx=24,my=24,mz=24

character(MAXSTRLENGTH) :: argv
integer :: idx 

integer :: ix,iy,iz,nn, ii
integer(8) :: iigd

integer :: ti,tj,tk

integer,parameter :: ntot=5
real(8) :: pos0(ntot*3)
integer :: atype0(ntot)

real(8) :: drnd(3)
integer :: num_seed
integer,allocatable :: seeds(:)

!--- initialize random seed with MPI rank
call random_seed(size=num_seed)
allocate(seeds(num_seed))
seeds=myid
call random_seed(put=seeds)

call system_clock(ti,tk)

if(find_cmdline_argc('--mxyz',idx)) then
  call get_command_argument(idx+1,argv);  read(argv,*) mx
  call get_command_argument(idx+2,argv);  read(argv,*) my
  call get_command_argument(idx+3,argv);  read(argv,*) mz
endif

!--- allocate arrays
if(.not.allocated(atype)) call allocator(atype,1,NBUFFER)
if(.not.allocated(q)) call allocator(q,1,NBUFFER)
if(.not.allocated(rreal)) call allocator(rreal,1,NBUFFER,1,3)
if(.not.allocated(v)) call allocator(v,1,NBUFFER,1,3)
if(.not.allocated(f)) call allocator(f,1,NBUFFER,1,3)
if(.not.allocated(qsfp)) call allocator(qsfp,1,NBUFFER)
if(.not.allocated(qsfv)) call allocator(qsfv,1,NBUFFER)
f(:,:)=0.0d0

pos0=(/ &
1.000000000d0,   1.000000000d0,   1.000000000d0, &
0.500000000d0,   0.500000000d0,   0.537699952d0, &
0.500000000d0,   0.500000000d0,   0.111800048d0, &
1.000000000d0,   0.500000000d0,   0.617399904d0, &
0.500000000d0,   1.000000000d0,   0.617399904d0  &
/)

atype0=(/1, 3, 4, 4, 4/)  ! Pb-1, Sr-2, Ti-3, O-4

!--- local unit cell parameters
lata=3.90200d0
latb=3.90200d0
latc=4.15600d0
lalpha=90.0000d0
lbeta=90.0000d0
lgamma=90.0000d0

iigd = mx*my*mz*ntot*myid ! for global ID
nn=0
do ix=0,mx-1
do iy=0,my-1
do iz=0,mz-1
   do ii=1,ntot
      nn=nn+1
      rreal(nn,1:3) = pos0(3*ii-2:3*ii)+(/ix,iy,iz/)  ! repeat unit cell
      rreal(nn,1:3) = rreal(nn,1:3)+vID(1:3)*(/mx,my,mz/) ! adding the box origin
      rreal(nn,1:3) = rreal(nn,1:3)*(/lata,latb,latc/) ! real coords
      call random_number(drnd)
      drnd(1:3) = (2.d0*drnd(1:3)-1.d0)*0.001d0
      !rreal(nn,1:3) = rreal(nn,1:3) + drnd(1:3)
      rreal(nn,1:3) = rreal(nn,1:3) - 1d-3
      atype(nn) = dble(atype0(ii)) + (iigd+nn)*1d-13
   enddo
enddo; enddo; enddo
NATOMS=nn

!--- update to glocal cell parameters
lata=lata*mx*vprocs(1)
latb=latb*my*vprocs(2)
latc=latc*mz*vprocs(3)

call get_boxparameters(mat,lata,latb,latc,lalpha,lbeta,lgamma)
do i=1, 3
do j=1, 3
   HH(i,j,0)=mat(i,j)
enddo; enddo
call update_box_params(vprocs, vid, hh, lata, latb, latc, maxrc, &
                       cc, lcsize, hhi, mdbox, lbox, obox)

call system_clock(tj,tk)
it_timer(22)=it_timer(22)+(tj-ti)

return
end


!--------------------------------------------------------------------------
subroutine ReadH2O(atype, rreal, v, q, f, fileName)
!--------------------------------------------------------------------------
implicit none

character(*),intent(in) :: fileName
real(8),allocatable,dimension(:),intent(inout) :: atype,q
real(8),allocatable,dimension(:,:),intent(inout) :: rreal,v,f

integer :: i,i1

integer (kind=MPI_OFFSET_KIND) :: offset, offsettmp
integer (kind=MPI_OFFSET_KIND) :: fileSize
integer :: localDataSize, metaDataSize, scanbuf
integer :: fh ! file handler

integer :: nmeta
integer,allocatable :: idata(:)
real(8),allocatable :: dbuf(:)
real(8) :: ddata(6), d10(10)

real(8) :: rnorm(NBUFFER,3), mat(3,3)
integer :: j

!=== # of unit cells ===
integer :: mx=2,my=1,mz=1

integer :: ix,iy,iz,ntot, ii
integer(8) :: iigd

integer :: ti,tj,tk

integer,parameter :: nH2O=24
real(8) :: pos0(nH2O*3)
integer :: atype0(nH2O)

call system_clock(ti,tk)

!--- allocate arrays
if(.not.allocated(atype)) call allocator(atype,1,NBUFFER)
if(.not.allocated(q)) call allocator(q,1,NBUFFER)
if(.not.allocated(rreal)) call allocator(rreal,1,NBUFFER,1,3)
if(.not.allocated(v)) call allocator(v,1,NBUFFER,1,3)
if(.not.allocated(f)) call allocator(f,1,NBUFFER,1,3)
if(.not.allocated(qsfp)) call allocator(qsfp,1,NBUFFER)
if(.not.allocated(qsfv)) call allocator(qsfv,1,NBUFFER)
f(:,:)=0.0

pos0=(/ &
5.000000D-01, 1.666262D-01, 7.617751D-01,&
5.000000D-01, 1.666262D-01, 8.729921D-01,&
6.550273D-01, 2.179985D-01, 7.240675D-01,&
0.000001D+00, 3.333738D-01, 6.381704D-01,&
0.000001D+00, 4.364829D-01, 6.760142D-01,&
1.543963D-01, 2.815157D-01, 6.760142D-01,&
0.000001D+00, 6.666262D-01, 7.617751D-01,&
0.000001D+00, 6.666262D-01, 8.729921D-01,&
1.550273D-01, 7.179985D-01, 7.240675D-01,&
0.999999D+00, 3.333738D-01, 2.617751D-01,&
8.456037D-01, 2.815157D-01, 2.240675D-01,&
0.999999D+00, 3.333738D-01, 3.729921D-01,&
5.000000D-01, 8.333738D-01, 2.617751D-01,&
5.000000D-01, 8.333738D-01, 3.729921D-01,&
3.456037D-01, 7.815157D-01, 2.240675D-01,&
5.000000D-01, 8.333738D-01, 6.381704D-01,&
4.995793D-01, 9.364829D-01, 6.760142D-01,&
6.543963D-01, 7.815157D-01, 6.760142D-01,&
0.999999D+00, 6.666262D-01, 1.381704D-01,&
0.999999D+00, 5.635171D-01, 1.760142D-01,&
8.449727D-01, 7.179985D-01, 1.760142D-01,&
5.000000D-01, 1.666262D-01, 1.381704D-01,&
4.995793D-01, 6.351712D-02, 1.760142D-01,&
3.449727D-01, 2.179985D-01, 1.760142D-01 &
/)

atype0=(/2,1,1, 2,1,1, 2,1,1, 2,1,1, 2,1,1, 2,1,1, 2,1,1, 2,1,1/)

!--- local unit cell parameters
lata=4.5181d0
latb=lata*sqrt(3.d0)
latc=7.346d0
lalpha=90.0000d0
lbeta=90.0000d0
lgamma=90.0000d0

iigd = mx*my*mz*nH2O*myid ! for global ID
ntot=0
do ix=0,mx-1
do iy=0,my-1
do iz=0,mz-1
   do ii=1,nH2O
      ntot=ntot+1
      rreal(ntot,1:3) = pos0(3*ii-2:3*ii)+(/ix,iy,iz/)  ! repeat unit cell
      rreal(ntot,1:3) = rreal(ntot,1:3)+vID(1:3)*(/mx,my,mz/) ! adding the box origin
      rreal(ntot,1:3) = rreal(ntot,1:3)*(/lata,latb,latc/) ! real coords
      atype(ntot) = dble(atype0(ii)) + (iigd+ntot)*1d-13
   enddo 
enddo; enddo; enddo
NATOMS=ntot

!--- update to glocal cell parameters
lata=lata*mx*vprocs(1)
latb=latb*my*vprocs(2)
latc=latc*mz*vprocs(3)

call get_boxparameters(mat,lata,latb,latc,lalpha,lbeta,lgamma)
do i=1, 3
do j=1, 3
   HH(i,j,0)=mat(i,j)
enddo; enddo
call update_box_params(vprocs, vid, hh, lata, latb, latc, maxrc, &
                       cc, lcsize, hhi, mdbox, lbox, obox)

call system_clock(tj,tk)
it_timer(22)=it_timer(22)+(tj-ti)

return
end

!--------------------------------------------------------------------------
subroutine ReadSiOH(atype, rreal, v, q, f, fileName)
!--------------------------------------------------------------------------
implicit none

character(*),intent(in) :: fileName
real(8),allocatable,dimension(:),intent(inout) :: atype,q
real(8),allocatable,dimension(:,:),intent(inout) :: rreal,v,f

integer :: i,i1

integer (kind=MPI_OFFSET_KIND) :: offset, offsettmp
integer (kind=MPI_OFFSET_KIND) :: fileSize
integer :: localDataSize, metaDataSize, scanbuf
integer :: fh ! file handler

integer :: nmeta
integer,allocatable :: idata(:)
real(8),allocatable :: dbuf(:)
real(8) :: ddata(6), d10(10)

real(8) :: rnorm(NBUFFER,3), mat(3,3)
integer :: j

!=== # of unit cells ===
!integer :: mx=256,my=256,mz=4
!integer :: mx=1,my=1,mz=1
!integer :: mx=32,my=32,mz=32
!integer :: mx=4,my=4,mz=4
!integer :: mx=16,my=16,mz=32
!integer :: mx=64,my=64,mz=64
!integer :: mx=48,my=48,mz=48
!integer :: mx=36,my=36,mz=36
!integer :: mx=35,my=35,mz=35
!integer :: mx=28,my=28,mz=28
!integer :: mx=30,my=30,mz=30
!integer :: mx=24,my=24,mz=24
integer :: mx=36,my=36,mz=36

integer :: ix,iy,iz,nn, ii
integer(8) :: iigd

integer :: ti,tj,tk

integer,parameter :: ntot=195
real(8) :: pos0(ntot*3)
integer :: atype0(ntot)

real(8) :: drnd(3)
integer :: num_seed
integer,allocatable :: seeds(:)

!--- initialize random seed with MPI rank
call random_seed(size=num_seed)
allocate(seeds(num_seed))
seeds=myid
call random_seed(put=seeds)

call system_clock(ti,tk)

!--- allocate arrays
if(.not.allocated(atype)) call allocator(atype,1,NBUFFER)
if(.not.allocated(q)) call allocator(q,1,NBUFFER)
if(.not.allocated(rreal)) call allocator(rreal,1,NBUFFER,1,3)
if(.not.allocated(v)) call allocator(v,1,NBUFFER,1,3)
if(.not.allocated(f)) call allocator(f,1,NBUFFER,1,3)
if(.not.allocated(qsfp)) call allocator(qsfp,1,NBUFFER)
if(.not.allocated(qsfv)) call allocator(qsfv,1,NBUFFER)
f(:,:)=0.0d0

pos0=(/ &
2.896559E-01,9.010904E-01,1.853931E-01, 3.284581E-01,9.100838E-01,2.208674E-01,&
3.633743E-01,8.542179E-01,2.208674E-01, 4.730954E-01,4.851874E-01,4.816755E-01,&
7.831386E-01,3.815818E-01,8.181271E-01, 7.841347E-01,8.345101E-01,8.294751E-03,&
4.210080E-01,5.890756E-01,6.276347E-02, 8.186435E-01,1.274824E-01,2.424698E-01,&
9.604279E-01,5.989141E-02,7.948055E-01, 7.745654E-01,5.371070E-01,1.020711E-01,&
5.086421E-01,1.000226E-02,1.866714E-01, 6.300482E-01,9.137201E-01,7.859194E-01,&
3.381796E-01,8.144284E-01,9.454980E-01, 4.007597E-01,8.910814E-01,8.871297E-01,&
6.898484E-01,9.067020E-01,6.949365E-01, 7.150693E-01,9.430006E-01,1.835806E-01,&
2.766885E-01,7.390428E-01,7.370425E-03, 3.825250E-01,1.422157E-01,1.838318E-02,&
6.049462E-01,6.134800E-01,7.167335E-01, 5.098488E-01,1.095387E-01,2.343418E-01,&
9.890532E-01,9.093238E-01,9.473161E-01, 8.635606E-01,5.985550E-01,9.495005E-01,&
3.553503E-01,3.621865E-02,3.704168E-02, 5.373235E-01,6.853075E-01,6.666658E-01,&
9.232759E-01,6.328259E-01,8.615342E-01, 6.391770E-01,4.062857E-01,6.617748E-01,&
5.176866E-01,6.125224E-01,9.488865E-03, 8.805516E-01,1.893870E-01,1.726099E-01,&
5.984725E-01,5.087002E-01,6.813437E-01, 7.097055E-01,8.688037E-01,2.649550E-01,&
5.852917E-01,5.747511E-01,9.268552E-01, 3.466217E-01,6.587114E-01,2.060747E-02,&
6.102737E-01,9.733487E-01,1.600488E-01, 7.769045E-01,7.238509E-01,3.218695E-02,&
4.698903E-01,1.856655E-01,6.569082E-01, 7.090914E-01,6.508966E-01,7.095283E-01,&
3.713505E-01,6.573180E-01,6.297700E-01, 7.732309E-01,3.484686E-02,2.031148E-01,&
6.389623E-01,8.948461E-01,5.945991E-01, 8.594670E-01,6.984814E-01,7.988840E-01,&
7.944160E-01,7.197825E-01,7.115469E-01, 8.900807E-01,8.670972E-01,9.840691E-01,&
9.618669E-01,9.757115E-01,8.643601E-01, 5.746862E-01,1.717205E-01,8.110785E-01,&
7.768524E-01,6.118316E-01,1.674961E-02, 9.897773E-01,9.659317E-01,5.918342E-01,&
9.663126E-01,2.183416E-02,6.871722E-01, 8.880904E-01,2.212976E-01,5.438771E-01,&
5.666917E-01,2.074291E-01,7.063006E-01, 4.545161E-01,7.257991E-01,6.053247E-01,&
3.788093E-01,9.808763E-01,9.457685E-01, 7.590583E-01,8.233711E-01,7.122782E-01,&
7.645411E-01,8.905862E-01,1.008070E-01, 6.990060E-01,8.658996E-01,9.430050E-01,&
5.107072E-01,8.783041E-01,9.026108E-01, 4.587294E-01,2.218041E-01,2.991219E-02,&
6.210815E-02,9.610189E-01,1.501334E-02, 6.451592E-01,1.588904E-01,6.472993E-01,&
4.008399E-01,4.893338E-01,1.386982E-02, 4.930919E-01,7.109935E-01,5.029526E-01,&
5.531571E-01,1.800233E-01,9.911525E-01, 7.416187E-01,2.010803E-01,2.633929E-01,&
8.710486E-01,9.017191E-02,3.335449E-01, 8.394341E-01,4.648740E-01,2.529003E-01,&
9.497607E-01,7.539658E-01,5.686438E-01, 4.806395E-01,9.434077E-01,2.720277E-01,&
2.973899E-01,5.753444E-01,6.301227E-01, 9.544116E-01,5.499997E-01,7.942898E-01,&
9.675689E-01,2.757018E-01,5.965484E-01, 3.160196E-01,6.737609E-01,2.655250E-01,&
5.780100E-01,6.133282E-01,8.241823E-01, 4.257382E-01,4.660718E-03,1.150892E-01,&
1.986427E-02,6.755438E-01,8.916163E-01, 6.144657E-01,9.190183E-01,8.939318E-01,&
4.256239E-01,8.322072E-01,6.281890E-01, 6.057205E-01,2.568390E-02,9.136587E-01,&
6.122239E-01,1.346789E-01,9.091428E-01, 3.773772E-01,9.256900E-01,6.672892E-01,&
3.864387E-01,1.380905E-01,6.025954E-01, 4.050902E-01,2.787312E-02,6.249908E-01,&
9.538102E-01,8.625610E-01,5.898874E-01, 3.121571E-01,5.186360E-01,5.347804E-01,&
4.278099E-02,1.259515E-01,8.231410E-01, 2.984451E-01,5.113200E-01,7.207668E-01,&
4.565571E-01,9.140708E-01,3.758934E-01, 4.208843E-01,5.752568E-01,1.733197E-01,&
4.041795E-02,6.461991E-02,3.346638E-01, 3.653323E-02,7.073202E-01,6.193164E-01,&
2.943749E-01,1.724462E-01,6.508252E-01, 7.348493E-01,7.361405E-03,6.767735E-01,&
9.827535E-01,4.447695E-01,7.816432E-01, 3.937478E-01,5.990458E-01,2.776857E-01,&
4.786979E-01,6.448668E-01,3.294016E-01, 7.251099E-01,9.553030E-02,6.081733E-01,&
8.756580E-01,1.225915E-01,8.318812E-01, 2.400876E-01,9.076028E-01,3.538787E-01,&
9.560741E-01,7.601678E-02,4.038507E-01, 7.855981E-01,4.401962E-01,1.608259E-01,&
4.546469E-01,2.416283E-01,1.385770E-01, 6.901990E-01,5.790918E-01,9.582157E-01,&
4.850583E-01,2.162330E-01,2.413446E-01, 4.190563E-01,2.350195E-01,3.280114E-01,&
7.275494E-02,5.792757E-01,2.638648E-01, 5.782224E-01,2.682066E-01,2.579288E-01,&
3.460393E-01,9.241968E-01,3.774716E-01, 2.270407E-01,2.400341E-01,7.044641E-01,&
7.207597E-01,8.001341E-01,3.502985E-01, 9.885522E-01,6.500119E-01,2.931208E-01,&
8.220983E-01,1.460230E-01,5.913700E-01, 8.210295E-01,1.896173E-01,9.041591E-01,&
7.205995E-01,1.493645E-01,9.301202E-01, 1.369419E-01,7.129716E-01,6.623844E-01,&
2.478784E-01,1.179310E-02,5.805950E-02, 3.808442E-01,4.242854E-01,8.493069E-01,&
3.758382E-01,9.147423E-01,7.788336E-01, 2.886746E-01,4.263777E-01,7.913482E-01,&
8.614164E-01,7.156365E-01,6.228500E-01, 8.604704E-01,4.806478E-01,3.599740E-01,&
5.418759E-01,4.703822E-01,9.209730E-01, 8.774398E-01,2.267100E-01,9.946778E-01,&
2.301355E-01,7.278814E-01,2.256577E-01, 3.621561E-01,5.209354E-01,3.500545E-01,&
5.939486E-02,4.972579E-01,1.852525E-01, 1.051101E-01,7.434405E-01,8.756937E-01,&
8.220401E-01,7.523909E-01,3.457173E-01, 9.108433E-01,6.872574E-01,3.673464E-01,&
1.257596E-01,3.804269E-01,9.278830E-01, 9.839576E-02,9.514671E-01,5.873459E-01,&
6.822895E-01,2.916038E-01,2.765987E-01, 2.753930E-01,9.293491E-01,6.261838E-01,&
1.260493E-01,1.700427E-01,8.798536E-01, 3.728272E-01,1.831091E-01,4.986538E-01,&
3.562509E-02,4.098399E-01,6.899745E-01, 1.384538E-01,3.462313E-02,5.156590E-02,&
5.570552E-02,8.373294E-01,8.991911E-01, 6.844452E-01,3.776806E-02,5.198634E-01,&
9.011644E-01,2.700973E-01,9.602314E-02, 6.402284E-01,9.389418E-01,4.919515E-01,&
4.411661E-01,3.165469E-01,9.786233E-01, 2.508962E-01,7.592453E-01,1.174561E-01,&
3.963012E-02,4.198085E-01,8.729506E-01, 6.351798E-01,7.318370E-01,3.735427E-01,&
1.099922E-01,5.582236E-02,1.575476E-01, 7.247754E-01,3.981602E-01,7.289840E-01,&
1.388353E-01,1.132990E-01,9.734019E-01, 4.432136E-01,4.180779E-01,9.385761E-01,&
5.275180E-01,7.138659E-01,3.983762E-01, 1.256044E-01,7.144048E-02,2.666783E-01,&
4.773155E-01,8.067879E-01,3.728665E-01, 7.828880E-01,2.772001E-01,8.478717E-01,&
9.340268E-02,5.159661E-01,3.545894E-01, 8.673314E-01,5.879349E-01,3.844986E-01,&
2.091274E-01,4.240206E-01,8.676743E-01, 9.592589E-01,6.625397E-03,4.917686E-01,&
2.230493E-01,8.188328E-01,2.895099E-01, 7.405778E-02,4.000256E-01,5.175708E-01,&
2.578153E-01,2.711233E-01,3.685286E-01, 3.563688E-01,2.593042E-01,4.159720E-01,&
8.127136E-02,3.942987E-01,1.398042E-01, 1.132795E-01,2.728742E-01,9.127418E-01,&
1.564464E-01,1.750041E-01,2.757450E-01, 6.715995E-01,3.892918E-01,5.576335E-01,&
1.208811E-01,4.076306E-01,3.632859E-02, 7.329602E-01,3.560607E-01,4.737680E-01,&
5.804217E-01,3.177917E-01,6.941749E-01, 7.076547E-01,3.781770E-01,2.080472E-01,&
1.413328E-01,2.731903E-01,6.422246E-01, 1.383902E-01,7.331808E-01,7.710320E-01,&
8.142595E-01,2.894255E-01,5.029664E-01, 2.772015E-01,3.316053E-01,7.368999E-01,&
9.537991E-01,7.266697E-01,4.610852E-01, 5.999115E-02,4.230914E-01,4.094235E-01,&
9.985171E-01,3.193229E-01,1.226266E-01, 7.054952E-01,8.749797E-01,4.305852E-01,&
1.735232E-01,2.800678E-01,2.986322E-01, 1.983702E-01,9.970196E-01,3.040427E-01,&
9.435719E-01,1.694680E-01,4.616827E-01, 5.385193E-01,9.488189E-01,4.448314E-01,&
1.981691E-01,1.874020E-01,7.981469E-01, 1.832686E-01,7.219048E-01,9.495778E-01,&
8.858276E-01,3.956711E-01,7.879858E-01, 1.778784E-01,7.958976E-01,6.061625E-01,&
6.852379E-01,3.152792E-01,3.842456E-01, 1.993977E-01,6.264346E-01,6.342948E-01,&
8.357823E-01,3.603480E-01,1.032717E-01, 1.501349E-01,6.510423E-01,2.308655E-01,&
9.830321E-02,3.321805E-01,3.593327E-01, 3.709520E-01,4.737921E-01,4.505308E-01,&
2.015521E-01,8.943016E-01,4.555886E-01, 5.789244E-02,3.384479E-01,6.087263E-01,&
3.693252E-01,3.631289E-01,4.515184E-01, 7.829931E-01,4.410501E-01,4.268626E-01,&
9.529140E-01,4.254272E-01,3.855712E-01, 1.907147E-01,8.957682E-01,5.653515E-01,&
1.573753E-01,3.423661E-01,2.054531E-01 &
/)

atype0=(/ &
1,2,1,2,3,3,3,3,3,2,3,2,2,3,3,3,3,2,3,2,&
3,2,3,2,3,3,2,2,2,2,3,2,2,2,2,2,2,2,2,2,&
3,2,2,2,3,3,2,3,3,3,2,2,2,2,2,3,2,2,2,2,&
2,2,2,2,3,2,3,2,2,2,2,2,2,3,2,2,3,3,3,2,&
2,2,2,2,3,2,2,2,2,2,3,3,2,3,2,3,3,3,2,2,&
3,2,3,2,2,3,3,2,2,3,2,3,2,2,2,3,2,3,2,2,&
3,2,2,3,2,3,3,2,3,2,3,2,2,3,2,2,3,3,2,2,&
2,2,2,2,2,3,3,3,2,2,2,2,2,2,2,2,2,3,3,2,&
2,2,2,3,2,2,2,2,2,2,2,3,2,2,3,2,2,2,2,2,&
2,2,2,2,2,2,2,3,2,3,2,2,2,3,2 &
/)

!--- local unit cell parameters
lata=14.32d0
latb=14.32d0
latc=14.32d0
lalpha=90.0000d0
lbeta=90.0000d0
lgamma=90.0000d0

iigd = mx*my*mz*ntot*myid ! for global ID
nn=0
do ix=0,mx-1
do iy=0,my-1
do iz=0,mz-1
   do ii=1,ntot
      nn=nn+1
      rreal(nn,1:3) = pos0(3*ii-2:3*ii)+(/ix,iy,iz/)  ! repeat unit cell
      rreal(nn,1:3) = rreal(nn,1:3)+vID(1:3)*(/mx,my,mz/) ! adding the box origin
      rreal(nn,1:3) = rreal(nn,1:3)*(/lata,latb,latc/) ! real coords
      call random_number(drnd)
      drnd(1:3) = (2.d0*drnd(1:3)-1.d0)*0.001d0
      rreal(nn,1:3) = rreal(nn,1:3) + drnd(1:3)
      atype(nn) = dble(atype0(ii)) + (iigd+nn)*1d-13
   enddo
enddo; enddo; enddo
NATOMS=nn

!--- update to glocal cell parameters
lata=lata*mx*vprocs(1)
latb=latb*my*vprocs(2)
latc=latc*mz*vprocs(3)

call get_boxparameters(mat,lata,latb,latc,lalpha,lbeta,lgamma)
do i=1, 3
do j=1, 3
   HH(i,j,0)=mat(i,j)
enddo; enddo
call update_box_params(vprocs, vid, hh, lata, latb, latc, maxrc, &
                       cc, lcsize, hhi, mdbox, lbox, obox)

call system_clock(tj,tk)
it_timer(22)=it_timer(22)+(tj-ti)

return
end



!--------------------------------------------------------------------------
subroutine ReadBIN(atype, rreal, v, q, f, fileName)
!--------------------------------------------------------------------------
implicit none

character(*),intent(in) :: fileName
real(8),allocatable,dimension(:) :: atype,q
real(8),allocatable,dimension(:,:) :: rreal,v,f

integer :: i,i1

integer (kind=MPI_OFFSET_KIND) :: offset, offsettmp
integer (kind=MPI_OFFSET_KIND) :: fileSize
integer :: localDataSize, metaDataSize, scanbuf
integer :: fh ! file handler

integer :: nmeta
integer,allocatable :: idata(:)
real(8),allocatable :: dbuf(:)
real(8) :: ddata(6), d10(10)

real(8) :: rnorm(NBUFFER,3), mat(3,3)
integer :: j

integer :: ti,tj,tk
call system_clock(ti,tk)


! Meta Data: 
!  Total Number of MPI ranks and MPI ranks in xyz (4 integers)
!  Number of resident atoms per each MPI rank (nprocs integers) 
!  current step (1 integer) + lattice parameters (6 doubles)

nmeta=4+nprocs+1
allocate(idata(nmeta))
metaDataSize = 4*nmeta + 8*6

if(myid==0) write(6,'(2a)') 'INFO: Opening file in ReadBIN(): ', trim(fileName)
call MPI_File_Open(MPI_COMM_WORLD,trim(fileName),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)

if(ierr > 0) then
   if(myid==0) write(6,'(2a)') 'MPI_File_Open Error in ReadBIN(): ', trim(fileName)
   call MPI_Finalize(ierr)
   stop
endif

! read metadata at the beginning of file
offsettmp=0
call MPI_File_Seek(fh,offsettmp,MPI_SEEK_SET,ierr)
call MPI_File_Read(fh,idata,nmeta,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)

offsettmp=4*nmeta
call MPI_File_Seek(fh,offsettmp,MPI_SEEK_SET,ierr)
call MPI_File_Read(fh,ddata,6,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

NATOMS = idata(4+myid+1)
current_step = idata(nmeta)
deallocate(idata)
lata=ddata(1); latb=ddata(2); latc=ddata(3)
lalpha=ddata(4); lbeta=ddata(5); lgamma=ddata(6)

! Get local datasize: 10 doubles for each atoms
localDataSize = 8*NATOMS*10

! offset will point the end of local write after the scan
call MPI_Scan(localDataSize,scanbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

! Since offset is MPI_OFFSET_KIND and localDataSize is integer, use an integer as buffer
offset = scanbuf + metaDataSize

! nprocs-1 rank has the total data size
fileSize = offset
!call MPI_Bcast(fileSize,1,MPI_INTEGER,nprocs-1,MPI_COMM_WORLD,ierr)
!call MPI_File_set_size(fh, fileSize, ierr)

! set offset at the beginning of the local write
offset=offset-localDataSize
call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)

allocate(dbuf(10*NATOMS))
call MPI_File_Read(fh,dbuf,10*NATOMS,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

if(.not.allocated(atype)) call allocator(atype,1,NBUFFER)
if(.not.allocated(q)) call allocator(q,1,NBUFFER)
if(.not.allocated(rreal)) call allocator(rreal,1,NBUFFER,1,3)
if(.not.allocated(v)) call allocator(v,1,NBUFFER,1,3)
if(.not.allocated(f)) call allocator(f,1,NBUFFER,1,3)
if(.not.allocated(qsfp)) call allocator(qsfp,1,NBUFFER)
if(.not.allocated(qsfv)) call allocator(qsfv,1,NBUFFER)
f(:,:)=0.0

do i=1, NATOMS
    i1=10*(i-1)
    rnorm(i,1:3)=dbuf(i1+1:i1+3)
    v(i,1:3)=dbuf(i1+4:i1+6)
    q(i)=dbuf(i1+7)
    atype(i)=dbuf(i1+8)
    qsfp(i)=dbuf(i1+9)
    qsfv(i)=dbuf(i1+10)
enddo
deallocate(dbuf)

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
call MPI_File_Close(fh,ierr)

call get_boxparameters(mat,lata,latb,latc,lalpha,lbeta,lgamma)
do i=1, 3
do j=1, 3
   HH(i,j,0)=mat(i,j)
enddo; enddo

! Check: box params are updated first time here but cutoff distance is determined 
! based on ffield params. give a flag? 
call update_box_params(vprocs, vid, hh, lata, latb, latc, maxrc, &
                       cc, lcsize, hhi, mdbox, lbox, obox)

call xs2xu(hh,obox,NATOMS,rnorm,rreal)

call system_clock(tj,tk)
it_timer(22)=it_timer(22)+(tj-ti)

return
end

!--------------------------------------------------------------------------
subroutine WriteBIN(atype, rreal, v, q, fileNameBase)
!--------------------------------------------------------------------------
implicit none

real(8),allocatable,intent(in) :: atype(:), q(:)
real(8),allocatable,intent(in) :: rreal(:,:),v(:,:)
character(len=:),allocatable,intent(in) :: fileNameBase

integer :: i,j

integer (kind=MPI_OFFSET_KIND) :: offset, offsettmp
integer :: localDataSize, metaDataSize, scanbuf
integer :: fh ! file handler

integer :: nmeta
integer,allocatable :: ldata(:),gdata(:)
real(8) :: ddata(6)
real(8),allocatable :: dbuf(:)

real(8) :: rnorm(NBUFFER,3)

integer :: ti,tj,tk
call system_clock(ti,tk)

call xu2xs(hhi,obox,NATOMS,rreal,rnorm)

if(.not. isBinary) return

! Meta Data: 
!  Total Number of MPI ranks and MPI ranks in xyz (4 integers)
!  Number of resident atoms per each MPI rank (nprocs integers) 
!  current step (1 integer) + lattice parameters (6 doubles)
nmeta=4+nprocs+1
metaDataSize = 4*nmeta + 8*6

! Get local datasize: 10 doubles for each atoms
localDataSize = 8*NATOMS*10

call MPI_File_Open(MPI_COMM_WORLD,trim(fileNameBase)//".bin",MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

! offset will point the end of local write after the scan
call MPI_Scan(localDataSize,scanbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

! Since offset is MPI_OFFSET_KIND and localDataSize is integer, use an integer as buffer
offset = scanbuf + metaDataSize

! save metadata at the beginning of file
allocate(ldata(nmeta),gdata(nmeta))
ldata(:)=0
ldata(4+myid+1)=NATOMS
call MPI_ALLREDUCE(MPI_IN_PLACE,ldata,nmeta,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
gdata = ldata
gdata(1)=nprocs
gdata(2:4)=vprocs
gdata(nmeta)=nstep+current_step

ddata(1)=lata; ddata(2)=latb; ddata(3)=latc
ddata(4)=lalpha; ddata(5)=lbeta; ddata(6)=lgamma

if(myid==0) then
   offsettmp=0
   call MPI_File_Seek(fh,offsettmp,MPI_SEEK_SET,ierr)
   call MPI_File_Write(fh,gdata,nmeta,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)

   offsettmp=4*nmeta
   call MPI_File_Seek(fh,offsettmp,MPI_SEEK_SET,ierr)
   call MPI_File_Write(fh,ddata,6,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
endif
deallocate(ldata,gdata)

! set offset at the beginning of the local write
offset=offset-localDataSize
call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)

allocate(dbuf(10*NATOMS))
do i=1, NATOMS
   j = (i - 1)*10
   dbuf(j+1:j+3)=rnorm(i,1:3)
   dbuf(j+4:j+6)=v(i,1:3)
   dbuf(j+7)=q(i)
   dbuf(j+8)=atype(i)
   dbuf(j+9)=qsfp(i)
   dbuf(j+10)=qsfv(i)
enddo
call MPI_File_Write(fh,dbuf,10*NATOMS,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
deallocate(dbuf)

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
call MPI_File_Close(fh,ierr)

call system_clock(tj,tk)
it_timer(23)=it_timer(23)+(tj-ti)

return
end

end module
