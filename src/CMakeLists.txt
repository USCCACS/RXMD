set(src_rxmd
    bo.F90
    cg.F90
    cmdline.F90
    comm.F90
    constant.F90
    ensembles.F90
    fileio.F90
    fnn.F90
    harmopot.F90
    init.F90
    lists.F90
    memory.F90
    module.F90
    msd.F90
    nnmm.F90
    pot.F90
    pqeq.F90
    qeq.F90
    reaxff.F90
    remd_pot_ANN_original.F90
    rxmdnn.F90
    shortrep.F90
    stepmod.F90
    symmfunc.F90
    utils.F90
    velocity.F90)

if(ENABLE_MPI)
  add_library(mpi_wrapper INTERFACE)
  target_link_libraries(mpi_wrapper INTERFACE MPI::MPI_Fortran)
else()
  add_library(mpi_wrapper nompi.F90)
  target_compile_definitions(mpi_wrapper INTERFACE NOMPI)
endif()

add_library(lib_rxmd ${src_rxmd})
target_link_libraries(lib_rxmd PUBLIC mpi_wrapper)

add_executable(rxmd main.F90)
set_target_properties(rxmd PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(rxmd PUBLIC lib_rxmd)
if(RXMD_ENABLE_TORCH)
  target_link_libraries(rxmd PUBLIC torch_rxmd)
endif()
