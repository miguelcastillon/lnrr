#include <gtest/gtest.h>

TEST(MultiplyTests, TestIntegerOne_One) { EXPECT_EQ(5 * 5, 20); }

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}