set(TORCH_SRCS src/rxmdtorch.cpp)
add_library(torch_rxmd ${TORCH_SRCS})
target_link_libraries(torch_rxmd PRIVATE rxmd_torchlibs)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fs\
anitize=address")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-point\
er -fsanitize=address")
