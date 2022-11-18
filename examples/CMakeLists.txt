project(lnrr-examples CXX)
cmake_minimum_required(VERSION 3.16)

find_package(lnrr REQUIRED)
find_package(Fgt REQUIRED)
find_package(Ceres REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_executable(lnrr_example
    lnrr.cpp
    )
target_link_libraries(lnrr_example
    Lnrr::Library-C++
    )
