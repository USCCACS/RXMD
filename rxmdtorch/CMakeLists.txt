####################################
# Set C++ for building torch
####################################

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(TORCH_SRCS src/rxmdtorch.cpp)
add_library(torch_rxmd ${TORCH_SRCS})

if(IPEX_FOUND)
  target_link_libraries(torch_rxmd "${TORCH_IPEX_LIBRARIES}")
  target_link_libraries(torch_rxmd PUBLIC SYCL::host)
  target_link_libraries(torch_rxmd PUBLIC SYCL::device)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_IPEX_CXX_FLAGS}")
elseif(TORCH_FOUND)
  target_link_libraries(torch_rxmd "${TORCH_LIBRARIES}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fs\
anitize=address")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-point\
er -fsanitize=address")
