cmake_minimum_required(VERSION 3.1)
project(lnrr)
set (${PROJECT_NAME}_VERSION_MAJOR 1)
set (${PROJECT_NAME}_VERSION_MINOR 0)
set(CMAKE_CXX_STANDARD 17) # avoid using EIGEN_MAKE_ALIGNED_OPERATOR_NEW
# add_compile_options(-march=native)
add_compile_options(-O3)
# add_compile_options(-g)
add_compile_options(-std=c++17)
# add_compile_options(-fopenmp)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

find_package(Eigen3 3.1 REQUIRED)
message(STATUS "[fgt] Eigen3 version: ${EIGEN3_VERSION}")

add_library(Library-C++
    src/scan_to_model.cpp
)
target_include_directories(Library-C++
    PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${EIGEN3_INCLUDE_DIR}
    )
