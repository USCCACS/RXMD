######################################################################
# CMake version and policies
######################################################################
cmake_minimum_required(VERSION 3.17.0)


######################################################################
# RXMD project
######################################################################
project(
  rxmd
  VERSION 0.1.1
  DESCRIPTION "RXMD: Linear Scalable Parallel ReaxFF Molecular Dynamics Simulator"
  LANGUAGES Fortran CXX)

set(PROJECT_CMAKE ${rxmd_SOURCE_DIR}/CMake)
list(APPEND CMAKE_MODULE_PATH ${rxmd_SOURCE_DIR}/CMake)

#####################################################################
# Build and install options
######################################################################

#--------------------------------------------------------------------
# PyTorch, MPI related build options
#--------------------------------------------------------------------
include(CMakeDependentOption)
option(RXMD_ENABLE_MPI "enable MPI" ON)
option(RXMD_ENABLE_TORCH "enable PyTorch" OFF)
cmake_dependent_option(RXMD_ENABLE_IPEX OFF "use Intel extensions for PyTorch" RXMD_ENABLE_TORCH OFF)

#####################################################################
# Set compiler specific options and flags
######################################################################

#-------------------------------------------------------------------
# Check MPI installation
#-------------------------------------------------------------------
set(ENABLE_MPI OFF)
if(RXMD_ENABLE_MPI)
  find_package(MPI REQUIRED Fortran)

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
add_library(rxmd_torch INTERFACE)

if(RXMD_ENABLE_IPEX)
  find_package(IPEX REQUIRED)
  if(NOT IPEX_FOUND)
    message(
      FATAL_ERROR
      "Intel Extensios for PyTorch are not found. Provide path to IPEXConfig.cmake"
    )
  endif()
elseif(RXMD_ENABLE_TORCH)
  find_package(Torch REQUIRED)
  if(NOT TORCH_FOUND)
    message(
      FATAL_ERROR
        "Torch support not found. Provide path to TorchConfig.cmake"
    )
  endif()  
endif()

if(IPEX_FOUND)
# if(NOT "x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xIntelLLVM")
#   message(FATAL_ERROR "Support for SYCL requires LLVM-based Intel oneAPI compiler (icpx).")
#  endif()
#  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 2023.0)
#    message(FATAL_ERROR "RXMD with IPEX Requires Intel oneAPI 2023.0 or higher!")
  endif()
  add_library(SYCL::host INTERFACE IMPORTED)
  add_library(SYCL::device INTERFACE IMPORTED)
  find_package(
    IntelDPCPP
    REQUIRED
    CONFIGS
    IntelDPCPPConfig-modified.cmake
    PATHS
    ${PROJECT_CMAKE}
    NO_DEFAULT_PATH)
  target_link_libraries(SYCL::host INTERFACE OneAPI::DPCPP-host)
  target_link_libraries(SYCL::device INTERFACE OneAPI::DPCPP-device)
endif()
######################################################################
message(STATUS "Ready to parse RXMD source tree")
######################################################################
add_subdirectory(src)
add_subdirectory(rxmdtorch)
