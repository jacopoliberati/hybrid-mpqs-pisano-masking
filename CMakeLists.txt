cmake_minimum_required(VERSION 3.16)
project(HybridMPQS CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)

add_executable(hybrid_mpqs_core src/main.cpp)

target_compile_options(hybrid_mpqs_core PRIVATE -O3 -mavx2 -march=native)
target_link_libraries(hybrid_mpqs_core PRIVATE OpenMP::OpenMP_CXX)
