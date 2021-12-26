module rxmdnn
  use iso_c_binding
  implicit none

  interface 

    subroutine init_rxmdnn() bind(c,name="init_rxmdnn")
    end subroutine
  
    subroutine predict_rxmdnn() bind(c,name="predict_rxmdnn")
    end subroutine
  
  end interface

contains
  
end module

program main
  use rxmdnn

  call init_rxmdnn()
  call predict_rxmdnn()

end program
