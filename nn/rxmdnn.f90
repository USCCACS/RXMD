module rxmdnn
  use iso_c_binding
  implicit none

  interface 

    subroutine init() bind(c,name="init_rxmdnn")
    ! create NN model, pass w&b to GPU, allocate nbrdist on CPU
    end subroutine

    subroutine init_hybrid(natoms) bind(c,name="init_rxmdnn_hybrid")
    ! create NN model, pass w&b to GPU, allocate nbrdist on CPU
       import :: c_int
       integer(c_int),value :: natoms
    end subroutine
  
    subroutine predict() bind(c,name="predict_rxmdnn")
    ! compute nbrlist (on CPU?), feature and predict E&F (return them to CPU?)
    end subroutine

    subroutine predict_hybrid(natoms, maxnbrs, nbrdist_ptr) bind(c,name="predict_rxmdnn_hybrid")
    ! compute nbrlist (on CPU?), feature and predict E&F (return them to CPU?)
       import :: c_ptr, c_int
       integer(c_int),value :: natoms, maxnbrs
       type(c_ptr),value :: nbrdist_ptr
    end subroutine

    subroutine get_maxrc(maxrc) bind(c,name="get_maxrc_rxmdnn")
        import :: c_double
        real(c_double) :: maxrc
    end subroutine
  
  end interface

contains
  
end module

program main
  use rxmdnn
  real(8) :: maxrc=0.d0

  integer(c_int) :: natoms=10, nnbrs=20

  real(c_double),allocatable,target :: nbrdist(:)
  type(c_ptr) :: nbrdist_ptr

  allocate(nbrdist(natoms*nnbrs))

  call random_number(nbrdist) 
  nbrdist_ptr = c_loc(nbrdist(1))

  call init_hybrid(natoms)
  call predict_hybrid(natoms,natoms,nbrdist_ptr)
  call get_maxrc(maxrc)
  print*,'maxrc: ', maxrc

end program
