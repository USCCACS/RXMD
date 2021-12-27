module rxmdnn
  use iso_c_binding
  implicit none

  interface 

    subroutine init() bind(c,name="init_rxmdnn")
    ! create NN model, pass w&b to GPU, allocate nbrdist on CPU
    end subroutine
  
    subroutine predict() bind(c,name="predict_rxmdnn")

    ! compute nbrlist (on CPU?), feature and predict E&F (return them to CPU?)
    end subroutine

    subroutine get_maxrc(maxrc) bind(c,name="get_maxrc_rxmdnn")
        import :: c_float
        real(c_float) :: maxrc
    end subroutine
  
  end interface

contains
  
end module

program main
  use rxmdnn
  real(4) :: maxrc=0.d0

  call init()
  call predict()
  call get_maxrc(maxrc)
  print*,'maxrc: ', maxrc

end program
