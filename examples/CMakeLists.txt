if(CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

add_executable(lnrr
    lnrr.cpp
    )
target_link_libraries(lnrr
    Library-C++
    )
