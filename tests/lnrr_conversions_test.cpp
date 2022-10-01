#include <gtest/gtest.h>
#include <lnrr/utils/conversions.h>

struct RotationEquivalence {
    double roll;
    double pitch;
    double yaw;
    Eigen::Quaterniond q;
    Eigen::Matrix3d matrix;
};

RotationEquivalence data1 = {0.0, 0.0, 0.0,
                             Eigen::Quaterniond(1.0, 0.0, 0.0, 0.0),
                             Eigen::Matrix3d::Identity()};
RotationEquivalence data2 = {
    0.0, 0.0, M_PI_2, Eigen::Quaterniond(0.7071068, 0.0, 0.0, 0.7071068),
    (Eigen::Matrix3d() << 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
        .finished()};
// RotationEquivalence data3 = {0.0, M_PI_2, 0.0,
//                              Eigen::Quaterniond(1.0, 0.0, 0.0, 0.0),
//                              Eigen::Matrix3d::Identity()};
// RotationEquivalence data4 = {M_PI_2, 0.0, 0.0,
//                              Eigen::Quaterniond(1.0, 0.0, 0.0, 0.0),
//                              Eigen::Matrix3d{1, 0, 0, 0, 0, -1, 0, 1, 0}};

std::vector<RotationEquivalence> data = {data1, data2};

TEST(ConversionsTests, TestRpy2Quaternion) {
    for (size_t i = 0; i < data.size(); i++) {
        Eigen::Quaterniond q_function = lnrr::conversions::rpy2Quaternion(
            data[i].roll, data[i].pitch, data[i].yaw);
        EXPECT_NEAR(
            0.0,
            (data[i].q.toRotationMatrix() - q_function.toRotationMatrix())
                .cwiseAbs()
                .sum(),
            1e-6)
            << "i: " << i << std::endl
            << "roll: " << data[i].roll << std::endl
            << "pitch: " << data[i].pitch << std::endl
            << "yaw: " << data[i].yaw << std::endl
            << "q: " << data[i].q << std::endl
            << "q_function: " << q_function << std::endl;
    }
}

TEST(ConversionsTests, TestQuaternion2rpy) {
    for (size_t i = 0; i < data.size(); i++) {
        Eigen::Vector3d rpy = lnrr::conversions::quaternion2rpy(data[i].q);
        EXPECT_NEAR(data[i].roll, rpy[0], 1e-6)
            << "i: " << i << std::endl
            << "q: " << data[i].q << std::endl
            << "roll: " << data[i].roll << std::endl
            << "roll_function: " << rpy[0] << std::endl;
        EXPECT_NEAR(data[i].pitch, rpy[1], 1e-6)
            << "i: " << i << std::endl
            << "q: " << data[i].q << std::endl
            << "pitch: " << data[i].pitch << std::endl
            << "pitch_function: " << rpy[1] << std::endl;
        EXPECT_NEAR(data[i].yaw, rpy[2], 1e-6)
            << "i: " << i << std::endl
            << "q: " << data[i].q << std::endl
            << "yaw: " << data[i].yaw << std::endl
            << "yaw_function: " << rpy[2] << std::endl;
    }
}

TEST(ConversionsTests, TestRpy2RotationMatrix) {
    for (size_t i = 0; i < data.size(); i++) {
        Eigen::Matrix3d matrix_function = lnrr::conversions::rpy2RotationMatrix(
            data[i].roll, data[i].pitch, data[i].yaw);
        EXPECT_NEAR(0.0, (data[i].matrix - matrix_function).cwiseAbs().sum(),
                    1e-6)
            << "i: " << i << std::endl
            << "roll: " << data[i].roll << std::endl
            << "pitch: " << data[i].pitch << std::endl
            << "yaw: " << data[i].yaw << std::endl
            << "matrix: " << data[i].matrix << std::endl
            << "matrix_function: " << matrix_function << std::endl;
    }
}

TEST(ConversionsTests, TestRotationMatrix2rpy) {
    for (size_t i = 0; i < data.size(); i++) {
        Eigen::Vector3d rpy =
            lnrr::conversions::rotationMatrix2rpy(data[i].matrix);
        EXPECT_NEAR(data[i].roll, rpy[0], 1e-6)
            << "i: " << i << std::endl
            << "matrix: " << data[i].matrix << std::endl
            << "roll: " << data[i].roll << std::endl
            << "roll_function: " << rpy[0] << std::endl;
        EXPECT_NEAR(data[i].pitch, rpy[1], 1e-6)
            << "i: " << i << std::endl
            << "matrix: " << data[i].matrix << std::endl
            << "pitch: " << data[i].pitch << std::endl
            << "pitch_function: " << rpy[1] << std::endl;
        EXPECT_NEAR(data[i].yaw, rpy[2], 1e-6)
            << "i: " << i << std::endl
            << "matrix: " << data[i].matrix << std::endl
            << "yaw: " << data[i].yaw << std::endl
            << "yaw_function: " << rpy[2] << std::endl;
    }
}