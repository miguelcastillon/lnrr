#pragma once

#include <lnrr/utils/types.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

namespace lnrr {

inline Matrix computeG(const double& beta, const int& number_lines) {
    double k = -2.0 * beta * beta;
    Matrix G = Matrix(number_lines, number_lines);
    for (int i = 0; i < number_lines; i++)
        for (int j = 0; j < number_lines; j++)
            G(i, j) = std::abs(i - j);
    G = (G.array().pow(2) / k).exp();
    return G;
}

template <typename T>
MatrixT<T> computeJacobianG(const MatrixT<T>& G) {
    assert(G.rows() == G.cols());
    int number_lines = G.rows();
    MatrixT<T> JG = MatrixT<T>::Zero(number_lines * 6, number_lines * 6);
    for (Eigen::Index block = 0; block < 6; block++)
        JG.block(block * number_lines, block * number_lines, number_lines,
                 number_lines) = G;

    return JG;
}

template <typename T>
MatrixT<T> computeJacobianPointComposition(
    std::vector<RigidTransform<T>> T_lines, MatrixX3T<T>& moving,
    const Eigen::Matrix<int, Eigen::Dynamic, 1>& line_sizes) {
    int number_lines = T_lines.size();
    int number_points_M = moving.rows();
    assert(number_lines > 0);
    assert(number_points_M > 0);

    MatrixT<T> J = MatrixT<T>::Zero(3 * number_points_M, 6 * number_lines);

    int m = 0;
    for (int l = 0; l < number_lines; l++) {
        Matrix3T<T> rotation_matrix_l = T_lines[l].rot_mat;
        MatrixT<T> jacobian_block(3, 6);
        jacobian_block.block(0, 0, 3, 3) = rotation_matrix_l;
        for (size_t k = 0; k < line_sizes[l]; k++) {
            // Skew-symmetric matrix
            Matrix3T<T> moving_m_x = Matrix3T<T>::Zero();
            moving_m_x(0, 1) = -moving(m, 2);
            moving_m_x(1, 0) = moving(m, 2);
            moving_m_x(0, 2) = moving(m, 1);
            moving_m_x(2, 0) = -moving(m, 1);
            moving_m_x(1, 2) = -moving(m, 0);
            moving_m_x(2, 1) = moving(m, 0);
            jacobian_block.block(0, 3, 3, 3) = -rotation_matrix_l * moving_m_x;

            for (size_t i = 0; i < jacobian_block.rows(); i++)
                for (size_t j = 0; j < jacobian_block.cols(); j++)
                    J(i * number_points_M + m, j * number_lines + l) =
                        jacobian_block(i, j);

            m++;
        }
    }
    return J;
}

template <typename T>
std::vector<RigidTransform<T>> computeTransformations(const MatrixT<T>& G,
                                                      const MatrixX6T<T>& W) {
    assert(G.rows() == G.cols());
    assert(G.rows() == W.rows());
    assert(W.cols() == 6);

    int number_lines = G.rows();

    std::vector<RigidTransform<T>> T_lines;
    T_lines.resize(number_lines);

    MatrixT<T> GW = G * W;
    for (size_t l = 0; l < number_lines; l++) {
        T_lines[l].xyz(0) = GW(l, 0);
        T_lines[l].xyz(1) = GW(l, 1);
        T_lines[l].xyz(2) = GW(l, 2);
        Vector3T<T> rotation;
        rotation(0) = GW(l, 3);
        rotation(1) = GW(l, 4);
        rotation(2) = GW(l, 5);
        T angle = rotation.norm();
        if (abs(angle) > T(0.0)) {
            Vector3T<T> axis = rotation / angle;
            T_lines[l].rot_mat = Matrix3T<T>(Eigen::AngleAxis<T>(angle, axis));
        } else
            T_lines[l].rot_mat.setIdentity();
    }
    return T_lines;
}

template <typename T>
MatrixX3T<T> computeTransformedMoving(
    const MatrixX3T<T>& moving, const MatrixT<T>& G, const MatrixX6T<T>& W,
    const Eigen::Matrix<int, Eigen::Dynamic, 1>& line_sizes) {
    int number_lines = G.rows();

    std::vector<RigidTransform<T>> T_lines = computeTransformations(G, W);
    MatrixX3T<T> transformed_moving(moving.rows(), 3);
    int m = 0;
    // #pragma omp parallel private(m)
    for (size_t l = 0; l < number_lines; l++) {
        for (size_t j = 0; j < line_sizes[l]; j++) {
            Vector3T<T> y_m;
            y_m(0) = moving(m, 0);
            y_m(1) = moving(m, 1);
            y_m(2) = moving(m, 2);
            Vector3T<T> result = T_lines[l].rot_mat * y_m + T_lines[l].xyz;
            transformed_moving(m, 0) = result(0);
            transformed_moving(m, 1) = result(1);
            transformed_moving(m, 2) = result(2);
            m++;
        }
    }
    assert(m == moving.rows());
    return transformed_moving;
}

} // namespace lnrr
