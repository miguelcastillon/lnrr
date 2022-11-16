if(CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG)
    message(WARNING "Compiling Debug mode. Information will be printed but performance will be slow. Compile on Release for best performance.")
endif(CMAKE_BUILD_TYPE MATCHES Debug)

add_library(Library-C++
    STATIC
    scan_to_model.cpp
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