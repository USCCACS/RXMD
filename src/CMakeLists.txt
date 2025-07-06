# Source files for rxmd executable (complete list from Makefile)
set(RXMD_SOURCES
  module.F90
  utils.F90
  memory.F90
  velocity.F90
  pot.F90
  fileio.F90
  comm.F90
  init.F90
  bo.F90
  qeq.F90
  reaxff.F90
  main.F90
  pqeq.F90
  cmdline.F90
  lists.F90
  msd.F90
  ensembles.F90
  constant.F90
  rxmdnn.F90
  nnmm.F90
  stepmod.F90
  harmopot.F90
  nnip_modifier.F90
)



## Build main rxmd executable
#add_executable(rxmd ${RXMD_SOURCES})
#
## Set target properties
#set_target_properties(rxmd PROPERTIES
#  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
#)
#
## Link MPI libraries
#target_link_libraries(rxmd PRIVATE MPI::MPI_Fortran)


# Build rxmdnn executable (with C++ support)
add_executable(rxmdnn ${RXMD_SOURCES})

set_target_properties(rxmdnn PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
  OUTPUT_NAME rxmdnn
  LINKER_LANGUAGE Fortran
)

target_link_libraries(rxmdnn PRIVATE
  MPI::MPI_Fortran
  stdc++
  rxmdtorch
  ${TORCH_LIBRARIES}
)

target_compile_definitions(rxmdnn PRIVATE RXMDNN)

target_include_directories(rxmdnn PRIVATE ${TORCH_INCLUDE_DIRS})
