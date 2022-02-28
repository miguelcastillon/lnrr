cmake_minimum_required(VERSION 3.1)
set(LNRR_VERSION 0.1.0)
set(PROJECT_SOVERSION 0)

project(lnrr CXX)

message(STATUS "[lnrr] Compiler: ${CMAKE_CXX_COMPILER_ID}")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(Eigen3 REQUIRED)
message(STATUS "[lnrr] Eigen3 version: ${EIGEN3_VERSION}")

add_library(Library-C++
    src/scan_to_model.cpp
    )
target_include_directories(Library-C++ INTERFACE
    "$<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>"
    "$<INSTALL_INTERFACE:include>"
    PRIVATE include ${EIGEN3_INCLUDE_DIR}
    )
set_target_properties(Library-C++ PROPERTIES
    OUTPUT_NAME lnrr
    VERSION ${LNRR_VERSION}
    SOVERSION ${PROJECT_SOVERSION}
    )
target_compile_definitions(Library-C++
    PUBLIC
    FGT_WITH_OPENMP
    )
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    target_compile_options(Library-C++
        PUBLIC
        -std=c++11
        PRIVATE
        -Wall
        -pedantic
        # The no-nested-anon-types is for nanoflann
        -Wno-nested-anon-types
        )
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    target_compile_options(Library-C++
        PUBLIC
        -std=c++11
        # For some deprecations inside Eigen
        -Wno-deprecated-declarations
        PRIVATE
        -Wall
        -pedantic
        # For pragma omp's
        -Wno-unknown-pragmas
        )
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
    target_compile_options(Library-C++
        PUBLIC
        # Eigen does relative includes
        /wd4464
        # This is the famous "unreferenced inline function has been removed"
        /wd4514
        # Yo, I don't care if that function wasn't inlined
        /wd4710
        # Padding
        /wd4820
        # Automatic inline expansion
        /wd4711
        )
endif()
install(DIRECTORY include/lnrr DESTINATION include)
install(TARGETS Library-C++
    DESTINATION lib
    EXPORT lnrr-targets
    )

include(GenerateExportHeader)
generate_export_header(Library-C++)

include(CMakePackageConfigHelpers)
configure_file(cmake/lnrr-config.cmake
    ${PROJECT_BINARY_DIR}/lnrr-config.cmake
    @ONLY
    )
write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/lnrr-config-version.cmake
    VERSION ${LNRR_VERSION}
    COMPATIBILITY AnyNewerVersion
    )
install(FILES
    ${PROJECT_BINARY_DIR}/lnrr-config.cmake
    ${PROJECT_BINARY_DIR}/lnrr-config-version.cmake
    DESTINATION lib/cmake/lnrr
    )
install(EXPORT lnrr-targets
    DESTINATION lib/cmake/lnrr
    NAMESPACE Lnrr::
    )