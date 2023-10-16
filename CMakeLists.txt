######################################################################
# CMake version and policies
######################################################################
cmake_minimum_required(VERSION 3.17.0)


######################################################################
# RXMD project
######################################################################
project(
  RXMD
  VERSION 0.1.1
  Description "RXMD: Linear Scalable Parallel ReaxFF Molecular Dynamics Simulator"
  LANGUAGES Fortran CXX)

#####################################################################
# Build and install options
######################################################################

#--------------------------------------------------------------------
# PyTorch, MPI related build options
#--------------------------------------------------------------------
include(CMakeDependentOption)
option(RXMD_ENABLE_MPI "enable MPI" ON)
option(RXMD_ENABLE_TORCH "enable PyTorch" OFF)
cmake_dependent_option(RXMD_ENABLE_IPEX "use Intel extensions for PyTorch" RXMD_ENABLE_TORCH OFF)

#####################################################################
# Set compiler specific options and flags
######################################################################

#-------------------------------------------------------------------
# Check MPI installation
#-------------------------------------------------------------------
set(ENABLE_MPI OFF)
if(RXMD_ENABLE_MPI)
  find_package(MPI Required Fortran)

  if(NOT MPI_FOUND)
    message(
      FATAL_ERROR
        "MPI support not found! Provide MPI compiler wrappers or build without MPI by passing '-DRXMD_ENABLEMPI=OFF' to cmake."
    )
  endif()
  set(ENABLE_MPI ON)
endif()

#-------------------------------------------------------------------
# Check Torch installation
#-------------------------------------------------------------------
add_library(rxmd_torch)
set(ENABLE_TORCH OFF)
if(RXMD_ENABLE_TORCH)
  find_package(Torch REQUIRED)
  if(NOT TORCH_FOUND)
    message(
      FATAL_ERROR
        "Torch support not found. Provide path to TorchConfig.cmake")
    )
  endif()  
  set(ENABLE_TORCH ON)
endif()

if(RXMD_ENABLE_IPEX)
  find_package(IPEX REQUIRED)
  if(NOT IPEX_FOUND)
    FATAL_ERROR
      "Intel Extensios for PyTorch are not found. Provide path to IPEXConfig.cmake"
  endif()
  set(ENABLE_IPEX ON)
endif() 

######################################################################
message(STATUS "Ready to parse RXMD source tree")
######################################################################
add_subdirectory(src)
add_subdirectory(rxmdtorch)
