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
    include
    ${EIGEN3_INCLUDE_DIR} 
    ${CERES_INCLUDE_DIR}
    )
target_link_libraries(Library-C++
    PUBLIC
    Fgt::Library-C++
    ${CERES_LIBRARIES}
    )