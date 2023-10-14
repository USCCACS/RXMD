cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(rxmdtorch LANGUAGES CXX)

option(RXMD_ENABLE_IPEX "Enable use of IPEX on Intel Accelerators when using Torch" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(${PROJECT_NAME} src/rxmdtorch.cpp)

#####################################
# Torch
#####################################
find_package(Torch REQUIRED)
target_link_libraries(${PROJECT_NAME} "${TORCH_LIBRARIES}")

if(RXMD_ENABLE_IPEX)
  find_package(IPEX REQUIRED)
  target_link_libraries(${PROJECT_NAME} "${TORCH_IPEX_LIBRARIES}")
endif()

#####################################

#set(CMAKE_CXX_FLAGS "-O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fs\
anitize=address")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-point\
er -fsanitize=address")
