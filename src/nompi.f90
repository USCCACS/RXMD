module nompi

   use iso_fortran_env, only: int32, int64, real32, real64, real128

   integer,parameter :: mpi_status_size = 1
   integer :: mpi_sum=1, mpi_max=2, mpi_min=3

   integer :: mpi_comm_world, mpi_info, mpi_file, mpi_datatype, mpi_status, mpi_in_place
   integer :: mpi_mode_rdonly, mpi_mode_wronly, mpi_mode_create, mpi_seek_set
   integer :: mpi_info_null, mpi_status_ignore

   integer,parameter :: mpi_double = real64
   integer,parameter :: mpi_comm = int32
   integer,parameter :: mpi_character = int32
   integer,parameter :: mpi_integer = int32
   integer,parameter :: mpi_integer8 = int64
   integer,parameter :: mpi_double_precision = real64
   integer,parameter :: mpi_offset_kind = real64

   interface mpi_file_write
      module procedure :: mpi_file_write_r8, mpi_file_write_i4, mpi_file_write_c
   end interface

   interface mpi_send
      module procedure :: mpi_send_vector, mpi_send_scalar
   end interface

   interface mpi_allreduce
      module procedure :: mpi_allreduce_inplace, mpi_allreduce_outplace, mpi_allreduce_scalar
      module procedure :: mpi_allreduce_inplace_2d, mpi_allreduce_outplace_2d
   end interface

contains

!all mpi objects (e.g., mpi_datatype, mpi_comm) are of type integer in fortran.
!https://www.mpich.org/static/docs/v3.2/www3/mpi_comm_rank.html

   subroutine mpi_init( ierr )
     integer :: ierr
   end subroutine

   subroutine mpi_comm_rank(comm, rank, ierror)
     !type(mpi_comm), intent(in) :: comm
     integer, intent(in) :: comm
     integer, intent(out) :: rank
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_comm_size(comm, size, ierror)
     !type(mpi_comm), intent(in) :: comm
     integer, intent(in) :: comm
     integer, intent(out) :: size
     integer, optional, intent(out) :: ierror
   end subroutine


   subroutine mpi_finalize(ierror)
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_barrier(comm, ierror)
     !type(mpi_comm), intent(in) :: comm
     integer, intent(in) :: comm
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_send_vector(buf, count, datatype, dest, tag, comm, ierror)
     type(*), dimension(:), intent(in) :: buf
     integer, intent(in) :: count, dest, tag
     !type(mpi_datatype), intent(in) :: datatype
     !type(mpi_comm), intent(in) :: comm
     integer, intent(in) :: comm, datatype
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_allreduce_outplace(sbuf, rbuf, count, datatype, mpiop, comm, ierror)
     type(*), dimension(:) :: sbuf, rbuf
     integer, intent(in) :: count, mpiop
     integer, intent(in) :: comm, datatype
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_allreduce_inplace(inplace, buf, count, datatype, mpiop, comm, ierror)
     type(*), dimension(:) :: buf
     integer, intent(in) :: inplace, count, mpiop
     integer, intent(in) :: comm, datatype
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_allreduce_outplace_2d(sbuf, rbuf, count, datatype, mpiop, comm, ierror)
     type(*), dimension(:,:) :: sbuf, rbuf
     integer, intent(in) :: count, mpiop
     integer, intent(in) :: comm, datatype
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_allreduce_inplace_2d(inplace, buf, count, datatype, mpiop, comm, ierror)
     type(*), dimension(:,:) :: buf
     integer, intent(in) :: inplace, count, mpiop
     integer, intent(in) :: comm, datatype
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_allreduce_scalar(sbuf, rbuf, count, datatype, mpiop, comm, ierror)
     type(*) :: sbuf, rbuf
     integer, intent(in) :: count, mpiop
     integer, intent(in) :: comm, datatype
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_send_scalar(buf, count, datatype, dest, tag, comm, ierror)
     type(*), intent(in) :: buf
     integer, intent(in) :: count, dest, tag
     integer, intent(in) :: comm, datatype
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_recv(buf, count, datatype, source, tag, comm, status, ierror)
     type(*), dimension(:) :: buf
     integer, intent(in) :: count, source, tag
     integer, intent(in) :: comm, datatype
     integer :: status(mpi_status_size)
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_scan(sendbuf, recvbuf, count, datatype, op, comm, ierror)
     integer :: sendbuf, recvbuf
     integer :: count, datatype, op, comm, ierror
   end subroutine


   subroutine mpi_probe(source, tag, comm, status, ierror)
     integer :: source, tag, comm, status(mpi_status_size), ierror
   end subroutine


   subroutine mpi_get_count(status, datatype, count, ierror)
     integer :: status(mpi_status_size), datatype, count, ierror
   end subroutine

   subroutine mpi_bcast(buffer, count, datatype, root, comm, ierror)
     type(*), dimension(..) :: buffer
     integer, intent(in) :: count, root
     integer, intent(in) :: comm, datatype
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_file_set_size(fh, size, ierror)
     !type(mpi_file), intent(in) :: fh
     integer, intent(in) :: fh
     integer(kind=mpi_offset_kind), intent(in) :: size
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_file_open(comm, filename, amode, info, fh, ierror)
     integer, intent(in) :: comm
     character(len=*), intent(in) :: filename
     integer, intent(in) :: amode
     integer, intent(in) :: info, fh
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_file_write_r8(fh, buf, count, datatype, status, ierror)
     real(8), intent(in), dimension(:) :: buf
     integer, intent(in) :: fh, count, datatype, status
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_file_write_i4(fh, buf, count, datatype, status, ierror)
     integer, intent(in), dimension(:) :: buf
     integer, intent(in) :: fh, count, datatype, status
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_file_write_c(fh, buf, count, datatype, status, ierror)
     character(len=:),allocatable, intent(in) :: buf
     integer, intent(in) :: fh, count, datatype, status
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_file_read(fh, buf, count, datatype, status, ierror)
     type(*), dimension(:) :: buf
     integer, intent(in) :: fh, count, datatype, status 
     integer, optional, intent(out) :: ierror
   end subroutine

   subroutine mpi_file_seek(fh, offset, whence, ierror)
     integer :: fh, whence, ierror
     integer(kind=mpi_offset_kind) :: offset
   end subroutine

   subroutine mpi_file_close(fh, ierror)
     integer,  intent(inout) :: fh
     integer, optional, intent(out) :: ierror
   end subroutine

   real(8) function mpi_wtime() result(t)
     !t = cpu_time() 
   end function

end module
