cmake_minimum_required(VERSION 3.16)
project(rxmd VERSION 1.0.0 LANGUAGES Fortran CXX)

# Set Fortran standard
set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler-specific configurations based on make.inc
if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # Intel Fortran compiler
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -C")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -C -traceback")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # GNU Fortran compiler (Sophia configuration)
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fallow-argument-mismatch -mcmodel=medium")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -fcheck=all -fallow-argument-mismatch")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
  # IBM XL Fortran compiler
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -traceback")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -C -traceback")
endif()

message(STATUS "CMAKE_Fortran_COMPILER_ID: ${CMAKE_Fortran_COMPILER_ID}")

# Find MPI
find_package(MPI REQUIRED)

# Set MPI Fortran compiler
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})


# Find PyTorch for rxmdnn target
find_package(Torch REQUIRED)

# Find rxmdnn
add_subdirectory(src)

# Add rxmdtorch subdirectory if it exists
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/rxmdtorch/CMakeLists.txt")
  add_subdirectory(rxmdtorch)
endif()


# Set install paths relative to project root - use absolute paths
set(INSTALL_BIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/bin")

# Create install directories
file(MAKE_DIRECTORY ${INSTALL_BIN_DIR})

# Install targets
install(TARGETS rxmdnn RUNTIME DESTINATION ${INSTALL_BIN_DIR})

