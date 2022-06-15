module rxmdnn
  use iso_c_binding
  implicit none

  interface 

    subroutine init(natoms) bind(c,name="init_rxmdtorch")
    ! create NN model, pass w&b to GPU, allocate nbrdist on CPU
       import :: c_int
       integer(c_int),value :: natoms
    end subroutine
  
  end interface

contains
  
end module

program main
  use rxmdnn
  real(8) :: maxrc=0.d0

  integer(c_int) :: natoms=10, nnbrs=20

  real(c_float),allocatable,target :: nbrdist(:)
  type(c_ptr) :: nbrdist_ptr

  integer :: ity

  allocate(nbrdist(natoms*nnbrs))

  call random_number(nbrdist) 
  nbrdist_ptr = c_loc(nbrdist(1))

  call init(natoms)

end program
