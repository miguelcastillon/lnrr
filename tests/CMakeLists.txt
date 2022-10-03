enable_testing()
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.1
)
FetchContent_MakeAvailable(googletest)
add_library(GTest::GTest INTERFACE IMPORTED)
target_link_libraries(GTest::GTest INTERFACE gtest_main)

include(GoogleTest)

################## 

add_executable(conversions_test conversions_test.cpp)
target_include_directories(conversions_test PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
    ${EIGEN3_INCLUDE_DIR} 
    ${CERES_INCLUDE_DIR}
    )
target_link_libraries(conversions_test 
  GTest::gtest_main
  )
gtest_discover_tests(conversions_test)

################## 

add_executable(E_step_test E_step_test.cpp)
target_include_directories(E_step_test PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
    ${EIGEN3_INCLUDE_DIR} 
    ${CERES_INCLUDE_DIR}
    )
target_link_libraries(E_step_test 
  GTest::gtest_main
  Library-C++
  )
gtest_discover_tests(E_step_test)

################## 



