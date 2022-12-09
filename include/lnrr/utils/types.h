#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <chrono>

#include <pcl/io/pcd_io.h>

namespace lnrr {

class lnrr_error : public std::runtime_error {
public:
    explicit lnrr_error(const std::string what_arg)
        : std::runtime_error(what_arg) {}
};

typedef Eigen::MatrixXd Matrix;
typedef Eigen::SparseMatrix<double> Sparse;
typedef Eigen::Matrix3d Matrix3;
typedef Eigen::MatrixX3d MatrixX3;
typedef Eigen::Matrix<double, Eigen::Dynamic, 6> MatrixX6;
typedef Eigen::VectorXd Vector;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VectorInt;
typedef Eigen::Vector3d Vector3;

template <typename T>
using MatrixT = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template <typename T>
using SparseT = Eigen::SparseMatrix<T>;
template <typename T>
using Matrix3T = Eigen::Matrix<T, 3, 3>;
template <typename T>
using MatrixX3T = Eigen::Matrix<T, Eigen::Dynamic, 3>;
template <typename T>
using MatrixX6T = Eigen::Matrix<T, Eigen::Dynamic, 6>;
template <typename T>
using VectorT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <typename T>
using Vector3T = Eigen::Matrix<T, 3, 1>;

typedef pcl::PointCloud<pcl::PointXYZ> PointCloud;
typedef pcl::PointCloud<pcl::PointXYZ>::Ptr PointCloudPtr;

template <typename T>
struct RigidTransform {
    Vector3T<T> xyz;
    Matrix3T<T> rot_mat;
};

struct Result {
    MatrixX3 points;
    std::vector<RigidTransform<double>> line_transforms;
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