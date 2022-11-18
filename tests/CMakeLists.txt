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

add_executable(jacobian_G_test jacobian_G_test.cpp)
target_include_directories(jacobian_G_test PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
    ${EIGEN3_INCLUDE_DIR} 
    ${CERES_INCLUDE_DIR}
    )
target_link_libraries(jacobian_G_test 
  GTest::gtest_main
  Library-C++
  )
gtest_discover_tests(jacobian_G_test)

################## 

add_executable(jacobian_YR_T_test jacobian_YR_T_test.cpp)
target_include_directories(jacobian_YR_T_test PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
    ${EIGEN3_INCLUDE_DIR} 
    ${CERES_INCLUDE_DIR}
    )
target_link_libraries(jacobian_YR_T_test 
  GTest::gtest_main
  Library-C++
  )
gtest_discover_tests(jacobian_YR_T_test)

################## 

add_executable(M_step_test M_step_test.cpp)
target_include_directories(M_step_test PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
    ${EIGEN3_INCLUDE_DIR} 
    ${CERES_INCLUDE_DIR}
    )
target_link_libraries(M_step_test 
  GTest::gtest_main
  Library-C++
  )
gtest_discover_tests(M_step_test)

################## 
