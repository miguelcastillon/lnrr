#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <chrono>

namespace lnrr {

class lnrr_error : public std::runtime_error {
public:
    explicit lnrr_error(const std::string what_arg)
        : std::runtime_error(what_arg) {}
};

typedef Eigen::MatrixXd Matrix;
typedef Eigen::SparseMatrix<double> SparseMatrix;
typedef Eigen::Matrix3d Matrix3;
typedef Eigen::MatrixX3d MatrixX3;
typedef Eigen::VectorXd Vector;
typedef Eigen::Vector3d Vector3;

struct RigidTransform {
    Vector3 xyz;
    Matrix3 rot_mat;
};

struct RigidTransformRPY {
    Vector3 xyz;
    Vector3 rpy;
};

struct Result {
    Matrix points;
    std::vector<RigidTransformRPY> line_transforms;
    double sigma2;
    std::chrono::microseconds runtime;
    size_t iterations;
};

struct Probabilities {
    Vector p1;
    Vector pt1;
    Matrix px;
    double l;
};
} // namespace lnrr