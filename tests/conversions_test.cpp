#include <gtest/gtest.h>
#include <lnrr_se3/utils/conversions.h>

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
    0.0, 0.0, M_PI_4, Eigen::Quaterniond(0.9238795, 0.0, 0.0, 0.3826834),
    (Eigen::Matrix3d() << 0.7071068, -0.7071068, 0.0000000, 0.7071068,
     0.7071068, 0.0000000, 0.0000000, 0.0000000, 1.0000000)
        .finished()};
RotationEquivalence data3 = {
    0.0, M_PI_4, 0.0, Eigen::Quaterniond(0.9238795, 0.0, 0.3826834, 0.0),
    (Eigen::Matrix3d() << 0.7071068, 0.0000000, 0.7071068, 0.0000000, 1.0000000,
     0.0000000, -0.7071068, 0.0000000, 0.7071068)
        .finished()};
RotationEquivalence data4 = {
    M_PI_4, 0.0, 0.0, Eigen::Quaterniond(0.9238795, 0.3826834, 0.0, 0.0),
    (Eigen::Matrix3d() << 1.0000000, 0.0000000, 0.0000000, 0.0000000, 0.7071068,
     -0.7071068, 0.0000000, 0.7071068, 0.7071068)
        .finished()};
RotationEquivalence data5 = {
    M_PI_4, M_PI_4, M_PI_4,
    Eigen::Quaterniond(0.8446232, 0.1913417, 0.4619398, 0.1913417),
    (Eigen::Matrix3d() << 0.5, -0.1464466, 0.8535534, 0.5, 0.8535534,
     -0.1464466, -0.7071068, 0.5, 0.5)
        .finished()};
RotationEquivalence data6 = {
    0.122173, -0.5410521, 3.0194196,
    Eigen::Quaterniond(0.0424344, 0.2698338, 0.0424344, 0.9610351),
    (Eigen::Matrix3d() << -0.8507781, -0.0586615, 0.5222408, 0.1044624,
     -0.9927973, 0.0586615, 0.5150381, 0.1044624, 0.8507781)
        .finished()};

std::vector<RotationEquivalence> data = {data1, data2, data3,
                                         data4, data5, data6};

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
            << "q.w(): " << data[i].q.w() << ", "
            << "q.vec(): " << data[i].q.vec() << std::endl
            << "q_function.w(): " << q_function.w() << ", "
            << "q_function.vec(): " << q_function.vec() << std::endl;
    }
}

TEST(ConversionsTests, TestQuaternion2rpy) {
    for (size_t i = 0; i < data.size(); i++) {
        Eigen::Vector3d rpy = lnrr::conversions::quaternion2rpy(data[i].q);
        EXPECT_NEAR(data[i].roll, rpy[0], 1e-6)
            << "i: " << i << std::endl
            << "q.w(): " << data[i].q.w() << ", "
            << "q.vec(): " << data[i].q.vec() << std::endl
            << "roll: " << data[i].roll << std::endl
            << "roll_function: " << rpy[0] << std::endl;
        EXPECT_NEAR(data[i].pitch, rpy[1], 1e-6)
            << "i: " << i << std::endl
            << "q.w(): " << data[i].q.w() << ", "
            << "q.vec(): " << data[i].q.vec() << std::endl
            << "pitch: " << data[i].pitch << std::endl
            << "pitch_function: " << rpy[1] << std::endl;
        EXPECT_NEAR(data[i].yaw, rpy[2], 1e-6)
            << "i: " << i << std::endl
            << "q.w(): " << data[i].q.w() << ", "
            << "q.vec(): " << data[i].q.vec() << std::endl
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