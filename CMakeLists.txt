cmake_minimum_required(VERSION 3.16)
project(lnrr CXX)
set(LNRR_VERSION 0.1.0)
set(PROJECT_SOVERSION 0)

set(CMAKE_CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED ON)
# add_compile_options(-O3)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

option(BUILD_TESTS "Build test suite" ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(OpenMP REQUIRED)
find_package(Eigen3 3.3.7 REQUIRED)
find_package(Ceres 2.1 REQUIRED)
find_package(Fgt REQUIRED)

add_subdirectory(src)

if(BUILD_TESTS)
  include(CTest)
  add_subdirectory(tests)
endif()

# if(NOT CMAKE_BUILD_TYPE)
#   set(default_cmake_build_type "RelWithDebInfo")
#   message(STATUS "Setting build type to '${default_cmake_build_type}' as none was specified.")
#   set(CMAKE_BUILD_TYPE "${default_cmake_build_type}" CACHE
#       STRING "Choose the type of build." FORCE)
#   # Set the possible values of build type for cmake-gui
#   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
#     "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
# endif()
# message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
# if(CMAKE_BUILD_TYPE MATCHES Debug)
#   message(STATUS "Compiling Debug mode")
#   add_compile_options(-DMODE_DEBUG)
# endif(CMAKE_BUILD_TYPE MATCHES Debug)



# include(CMakePackageConfigHelpers)
# write_basic_package_version_file(${PROJECT_BINARY_DIR}/lnrr-config-version.cmake VERSION ${LNRR_VERSION} COMPATIBILITY AnyNewerVersion)
# configure_file(cmake/lnrr-config.cmake ${PROJECT_BINARY_DIR}/lnrr-config.cmake @ONLY)
# install(FILES ${PROJECT_BINARY_DIR}/lnrr-config.cmake ${PROJECT_BINARY_DIR}/lnrr-config-version.cmake DESTINATION lib/cmake/lnrr)


# if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
#     target_compile_options(Library-C++
#         PUBLIC
#         -std=c++11
#         PRIVATE
#         -Wall
#         -pedantic
#         -Wno-nested-anon-types
#         )
# elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
#     target_compile_options(Library-C++
#         PUBLIC
#         -std=c++11
#         -Wno-deprecated-declarations
#         PRIVATE
#         -Wall
#         -pedantic
#         -Wno-unknown-pragmas
#         )
# endif()
# install(TARGETS Library-C++ DESTINATION lib EXPORT lnrr-targets)
# install(DIRECTORY include/lnrr DESTINATION include)
# install(EXPORT lnrr-targets DESTINATION lib/cmake/lnrr NAMESPACE Lnrr::)

set(EXECUTABLE_OUTPUT_PATH ..)
add_executable(test_lnrr
    tests/test_lnrr.cpp
    )
target_link_libraries(test_lnrr
    Library-C++
    )

# # Create uninstall command
# configure_file("${CMAKE_SOURCE_DIR}/cmake/uninstall_target.cmake.in"
#                "${CMAKE_BINARY_DIR}/uninstall_target.cmake" IMMEDIATE @ONLY)
# add_custom_target(uninstall "${CMAKE_COMMAND}" -P
#                   "${CMAKE_BINARY_DIR}/uninstall_target.cmake")
