cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(rxmdtorch LANGUAGES CXX)

find_package(Python COMPONENTS Interpreter Development)
find_package(Torch REQUIRED) # Polaris
find_package(IPEX REQUIRED)  # Sunspot

#set(CMAKE_CXX_FLAGS "-O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
#set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(${PROJECT_NAME} src/rxmdtorch.cpp)

#target_link_libraries(${PROJECT_NAME} "${TORCH_LIBRARIES}") # Polaris
target_link_libraries(${PROJECT_NAME} "${TORCH_IPEX_LIBRARIES}") # Sunspot
