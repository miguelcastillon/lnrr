cmake_minimum_required(VERSION 3.1)
set(LNRR_VERSION 0.1.0)
set(PROJECT_SOVERSION 0)

project(lnrr CXX)

message(STATUS "[lnrr] Compiler: ${CMAKE_CXX_COMPILER_ID}")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(Eigen3 REQUIRED)
find_package(Ceres REQUIRED)
find_package(Fgt REQUIRED)
message(STATUS "[lnrr] Eigen3 version: ${EIGEN3_VERSION}")

include(CMakePackageConfigHelpers)
write_basic_package_version_file(${PROJECT_BINARY_DIR}/lnrr-config-version.cmake VERSION ${LNRR_VERSION} COMPATIBILITY AnyNewerVersion)
configure_file(cmake/lnrr-config.cmake ${PROJECT_BINARY_DIR}/lnrr-config.cmake @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/lnrr-config.cmake ${PROJECT_BINARY_DIR}/lnrr-config-version.cmake DESTINATION lib/cmake/lnrr)

add_library(Library-C++
    STATIC
    src/scan_to_model.cpp
    )
set_target_properties(Library-C++ PROPERTIES
    OUTPUT_NAME lnrr
    VERSION ${LNRR_VERSION}
    SOVERSION ${PROJECT_SOVERSION}
    )
target_include_directories(Library-C++ PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
    ${EIGEN3_INCLUDE_DIR} 
    ${CERES_INCLUDE_DIR}
    )
target_link_libraries(Library-C++
    PUBLIC
    Fgt::Library-C++
    ${CERES_LIBRARIES}
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
install(TARGETS Library-C++ DESTINATION lib EXPORT lnrr-targets)
install(DIRECTORY include/lnrr DESTINATION include)
install(EXPORT lnrr-targets DESTINATION lib/cmake/lnrr NAMESPACE Lnrr::)

# include(GenerateExportHeader)
# generate_export_header(Library-C++)