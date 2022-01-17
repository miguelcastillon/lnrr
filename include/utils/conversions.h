#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

namespace cld {
namespace conversions {

template <typename T>
Eigen::Quaternion<T> rpy2Quaternion(const T& roll, const T& pitch,
                                    const T& yaw) {
    Eigen::AngleAxis<T> rollAngle(roll, Eigen::Matrix<T, 3, 1>::UnitX());
    Eigen::AngleAxis<T> pitchAngle(pitch, Eigen::Matrix<T, 3, 1>::UnitY());
    Eigen::AngleAxis<T> yawAngle(yaw, Eigen::Matrix<T, 3, 1>::UnitZ());

    Eigen::Quaternion<T> q = yawAngle * pitchAngle * rollAngle;
    return q;
}

template <typename T>
Eigen::Matrix<T, 3, 1> quaternion2rpy(const Eigen::Quaternion<T>& q) {
    return q.toRotationMatrix().eulerAngles(2, 1, 0).colwise().reverse();
}

template <typename T>
Eigen::Matrix<T, 3, 3> rpy2RotationMatrix(const T& roll, const T& pitch,
                                          const T& yaw) {
    return rpy2Quaternion(roll, pitch, yaw).toRotationMatrix();
}

template <typename T>
Eigen::Matrix<T, 3, 1>
rotationMatrix2rpy(const Eigen::Matrix<T, 3, 3>& rotation_matrix) {
    return rotation_matrix.eulerAngles(2, 1, 0).colwise().reverse();
}

} // namespace conversions
} // namespace cld
