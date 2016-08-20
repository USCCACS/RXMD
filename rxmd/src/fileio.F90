!----------------------------------------------------------------------------------------
subroutine OUTPUT()
use atoms; use parameters
!----------------------------------------------------------------------------------------
implicit none
integer :: i, ity, j, j1, jty, m, n,n3, cstep, iscan
integer :: l2g
real(8) :: ri(3), rj(3), bndordr(MAXNEIGHBS), tt=0.d0, ss=0.d0
integer :: igd,jgd,bndlist(0:MAXNEIGHBS)
character(8) :: fname0
character(6) :: a6
character(9) :: a9
character(128) :: FileName

integer (kind=MPI_OFFSET_KIND) :: offset
integer (kind=MPI_OFFSET_KIND) :: fileSize
integer :: localDataSize
integer :: fh ! file handler

integer,parameter :: PDBLineSize=67
character(PDBLineSize) :: PDBOneLine

integer :: BNDLineSize
integer,parameter :: MaxBNDLineSize=512
character(MaxBNDLineSize) :: BNDOneLine
real(8),parameter :: BNDcutoff=0.3d0

integer :: scanbuf


!--- binary ------------------------------------------------------------------
if(isBinary) then
  call xu2xs()
  call coio_write(1)
  call xs2xu()
endif
!------------------------------------------------------------------ binary ---

write(a6(1:6),'(i6.6)') myid
write(a9(1:9),'(i9.9)') nstep + current_step

!--- BondFile -------------------------------------------------------------
if(isBondFile) then

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

200 format(i9,1x,3f9.3,1x,2i3,20(i9,f6.3)) 

    ! get local datasize based on above format and the total # of neighbors
    localDataSize=NATOMS*(9+1+3*9+1+2*3 +1)+m*(9+6)

    call MPI_File_Open(MPI_COMM_WORLD,trim(DataDir)//"/"//a9//".bnd", &
        MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

    ! offset will point the end of local write after the scan
    call MPI_Scan(localDataSize,scanbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

    ! since offset is MPI_OFFSET_KIND and localDataSize is integer, use an integer as buffer
    offset=scanbuf

    ! nprocs-1 rank has the total data size
    fileSize=offset
    call MPI_Bcast(fileSize,1,MPI_INTEGER,nprocs-1,MPI_COMM_WORLD,ierr)
    call MPI_File_set_size(fh, fileSize, ierr)

    ! set offset at the beginning of the local write
    offset=offset-localDataSize

    !open(10,file=trim(FileName)//".bnd")

   do i=1, NATOMS
      ity = atype(i)
!--- get global ID for i-atom
      igd = l2g(atype(i))

!--- count the number bonds to be shown.
      bndlist(0)=0
      do j1 = 1, nbrlist(i,0)
         j = nbrlist(i,j1)
         jty = atype(j)

!--- get global ID for j-atom
         jgd = l2g(atype(j))

!--- if bond order is less than 0.3, ignore the bond.
         if( BO(0,i,j1) < 0.3d0) cycle

         bndlist(0) = bndlist(0) + 1
         bndlist(bndlist(0)) = jgd
         bndordr(bndlist(0)) = BO(0,i,j1)
      enddo

        BNDOneLine=""
        write(BNDOneLine,200) igd, pos(1:3,i),int(atype(i)),bndlist(0), &
            (bndlist(j1),bndordr(j1),j1=1,bndlist(0))
        ! remove space and add new_line
        BNDOneLine=trim(BNDOneLine)//NEW_LINE('A')
        BNDLineSize=len(trim(BNDOneLine))

        call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)
        call MPI_File_Write(fh,BNDOneLine,BNDLineSize, &
            MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)

        offset=offset+BNDLineSize

   enddo

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   call MPI_File_Close(fh,ierr)

endif
!------------------------------------------------------------ BondFile ----

!--- PDB ---------------------------------------------------------------------
if(isPDB) then

   ! get local datasize
   localDataSize=NATOMS*PDBLineSize

   call MPI_File_Open(MPI_COMM_WORLD,trim(DataDir)//"/"//a9//".pdb", &
       MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

   ! offset will point the end of local write after the scan
   call MPI_Scan(localDataSize,scanbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

   ! since offset is MPI_OFFSET_KIND and localDataSize is integer, use an integer as buffer
   offset=scanbuf

   ! nprocs-1 rank has the total data size
   fileSize=offset
   call MPI_Bcast(fileSize,1,MPI_INTEGER,nprocs-1,MPI_COMM_WORLD,ierr)
   call MPI_File_set_size(fh, fileSize, ierr)

   ! set offset at the beginning of the local write
   offset=offset-localDataSize

   do i=1, NATOMS

     ity = atype(i)
!--- calculate atomic temperature 
     tt = hmas(ity)*sum(v(1:3,i)*v(1:3,i))
     tt = tt*UTEMP*1d-2 !scale down to use two decimals in PDB format 

!--- sum up diagonal atomic stress components 
#ifdef STRESS
     ss = sum(astr(1:3,i))/3.d0
#endif
     ss = ss*USTRS

     ss = q(i)*10 ! 10x atomic charge

     igd = l2g(atype(i))
     select case(ity)
       case(1) 
         write(PDBOneLine,100)'ATOM  ',0, 'C', igd, pos(1:3,i), tt, ss
       case(2) 
         write(PDBOneLine,100)'ATOM  ',0, 'H', igd, pos(1:3,i), tt, ss
       case(3) 
         write(PDBOneLine,100)'ATOM  ',0, 'O', igd, pos(1:3,i), tt, ss
       case(4) 
         write(PDBOneLine,100)'ATOM  ',0, 'N', igd, pos(1:3,i), tt, ss
       case(5) 
         write(PDBOneLine,100)'ATOM  ',0, 'S', igd, pos(1:3,i), tt, ss
       case(6) 
         write(PDBOneLine,100)'ATOM  ',0,'Si', igd, pos(1:3,i), tt, ss
       case(7) 
         write(PDBOneLine,100)'ATOM  ',0,'Al', igd, pos(1:3,i), tt, ss
     end select
       PDBOneLine(PDBLineSize:PDBLineSize)=NEW_LINE('A')

       call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)
       call MPI_File_Write(fh,PDBOneLine,PDBLineSize, &
           MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)

       offset=offset+PDBLineSize
   enddo

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   call MPI_File_Close(fh,ierr)
endif
100 format(A6,I5,1x,A2,i12,4x,3f8.3,f6.2,f6.2)
!-------------------------------------------------------------------- PDB ----

return
end subroutine

!--------------------------------------------------------------------------
subroutine coio_read()
use atoms; use parameters
!--------------------------------------------------------------------------
implicit none
integer :: i,j,k,k1,c1,n1
character(6) :: a6
character(9) :: a9

integer (kind=MPI_OFFSET_KIND) :: offset, offsettmp
integer (kind=MPI_OFFSET_KIND) :: fileSize
integer :: localDataSize, metaDataSize, scanbuf
integer :: fh ! file handler

integer :: nmeta
integer,allocatable :: idata(:)
real(8) :: ddata(6), d10(10)

character(128) :: FileName

write(a6(1:6),'(i6.6)') myid
write(a9(1:9),'(i9.9)') current_step

! Meta Data: 
!  Total Number of MPI ranks and MPI ranks in xyz (4 integers)
!  Number of resident atoms per each MPI rank (nprocs integers) 
!  current step (1 integer) + lattice parameters (6 doubles)
nmeta=4+nprocs+1
allocate(idata(nmeta))
metaDataSize = 4*nmeta + 8*6

call MPI_File_Open(MPI_COMM_WORLD,trim(DataDir)//"/rxff.bin",MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)

! read metadata at the beginning of file
offsettmp=0
call MPI_File_Seek(fh,offsettmp,MPI_SEEK_SET,ierr)
call MPI_File_Read(fh,idata,nmeta,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)

offsettmp=4*nmeta
call MPI_File_Seek(fh,offsettmp,MPI_SEEK_SET,ierr)
call MPI_File_Read(fh,ddata,6,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

NATOMS = idata(4+myid+1)
current_step = idata(nmeta)
lata=ddata(1); latb=ddata(2); latc=ddata(3)
lalpha=ddata(4); lbeta=ddata(5); lgamma=ddata(6)
!print*,'idata: ', idata(1:nmeta)
!print*,'ddata: ', ddata(1:6)

! Get local datasize: 10 doubles for each atoms
localDataSize = 8*NATOMS*10

! offset will point the end of local write after the scan
call MPI_Scan(localDataSize,scanbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

! Since offset is MPI_OFFSET_KIND and localDataSize is integer, use an integer as buffer
offset = scanbuf + metaDataSize

! nprocs-1 rank has the total data size
fileSize = offset
call MPI_Bcast(fileSize,1,MPI_INTEGER,nprocs-1,MPI_COMM_WORLD,ierr)
call MPI_File_set_size(fh, fileSize, ierr)

! set offset at the beginning of the local write
offset=offset-localDataSize
call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)

do i=1, NATOMS
   call MPI_File_Read(fh,d10,10,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

   pos(1:3,i)=d10(1:3)
   v(1:3,i)=d10(4:6)
   q(i)=d10(7)
   atype(i)=d10(8)
   qsfp(i)=d10(9)
   qsfv(i)=d10(10)

   offset=offset+10*8 ! 10 x 8bytes
   call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)
enddo

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
call MPI_File_Close(fh,ierr)

deallocate(idata)

return
end

!--------------------------------------------------------------------------
subroutine coio_write(imode)
use atoms; use parameters
!--------------------------------------------------------------------------
implicit none
integer,intent(in) :: imode 
integer :: i,j,k,k1,n1
character(6) :: a6
character(9) :: a9
character(128) :: FileName

integer (kind=MPI_OFFSET_KIND) :: offset, offsettmp
integer (kind=MPI_OFFSET_KIND) :: fileSize
integer :: localDataSize, metaDataSize, scanbuf
integer :: fh ! file handler

integer :: nmeta
integer,allocatable :: ldata(:),gdata(:)
real(8) :: ddata(6), d10(10)
real(8),allocatable :: dbuf(:)

! Meta Data: 
!  Total Number of MPI ranks and MPI ranks in xyz (4 integers)
!  Number of resident atoms per each MPI rank (nprocs integers) 
!  current step (1 integer) + lattice parameters (6 doubles)
nmeta=4+nprocs+1
metaDataSize = 4*nmeta + 8*6

! Get local datasize: 10 doubles for each atoms
localDataSize = 8*NATOMS*10

if(imode==-1) then
  FileName=trim(DataDir)//"/"//"rxff.bin"
else
  write(a9(1:9),'(i9.9)') nstep+current_step
  FileName=trim(DataDir)//"/"//a9//".bin"
endif

call MPI_File_Open(MPI_COMM_WORLD,trim(FileName),MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

! offset will point the end of local write after the scan
call MPI_Scan(localDataSize,scanbuf,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

! Since offset is MPI_OFFSET_KIND and localDataSize is integer, use an integer as buffer
offset = scanbuf + metaDataSize

! nprocs-1 rank has the total data size
fileSize = offset
call MPI_Bcast(fileSize,1,MPI_INTEGER,nprocs-1,MPI_COMM_WORLD,ierr)
call MPI_File_set_size(fh, fileSize, ierr)

! save metadata at the beginning of file
allocate(ldata(nmeta),gdata(nmeta))
ldata(:)=0
ldata(4+myid+1)=NATOMS
call MPI_ALLREDUCE(ldata,gdata,nmeta,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
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
   dbuf(j+1:j+3)=pos(1:3,i)
   dbuf(j+4:j+6)=v(1:3,i)
   dbuf(j+7)=q(i)
   dbuf(j+8)=atype(i)
   dbuf(j+9)=qsfp(i)
   dbuf(j+10)=qsfv(i)

   !offset=offset+10*8 ! 10 x 8bytes
   !call MPI_File_Seek(fh,offset,MPI_SEEK_SET,ierr)
enddo
call MPI_File_Write(fh,dbuf,10*NATOMS,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
deallocate(dbuf)

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
call MPI_File_Close(fh,ierr)

return
end
