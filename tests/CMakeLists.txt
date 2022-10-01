include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.1
)
FetchContent_MakeAvailable(googletest)
add_library(GTest::GTest INTERFACE IMPORTED)
target_link_libraries(GTest::GTest INTERFACE gtest_main)

add_executable(lnrr_test lnrr_test.cpp)

target_link_libraries(lnrr_test
 PRIVATE
  GTest::GTest
  )

add_test(lnrr_gtests lnrr_test)