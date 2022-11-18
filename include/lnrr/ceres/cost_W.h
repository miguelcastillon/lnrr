#pragma once

#include <ceres/ceres.h>
#include <lnrr/utils/conversions.h>
#include <lnrr/utils/operations.h>
#include <lnrr/utils/types.h>

namespace lnrr {
class CostFunctionW {
public:
    CostFunctionW(const Matrix& moving, const Matrix& G,
                  const Matrix& jacobianG, const Vector& PX_vec,
                  const Vector& P1, const VectorInt& line_sizes,
                  const double& lambda)
        : moving_(moving),
          G_(G),
          jacobianG_(jacobianG),
          PX_vec_(PX_vec),
          P1_(P1),
          line_sizes_(line_sizes),
          lambda_(lambda) {}

    template <typename T>
    bool operator()(T const* const* parameters, T* residual) const {
        int number_lines_scan = G_.rows();
        const Eigen::Map<const MatrixX6T<T>> W_matrix(parameters[0],
                                                      number_lines_scan, 6);
        MatrixT<T> G_T = G_.cast<T>();

        std::vector<RigidTransform<T>> T_lines =
            computeTransformations(G_T, MatrixX6T<T>(W_matrix));

        MatrixX3T<T> moving_T = moving_.cast<T>();
        MatrixT<T> F =
            computeJacobianPointComposition(T_lines, moving_T, line_sizes_);
        MatrixX3T<T> moving_transformed = computeTransformedMoving(
            moving_T, G_T, MatrixX6T<T>(W_matrix), line_sizes_);
        MatrixX6T<T> lambdaGW = T(lambda_) * G_T * W_matrix;
        const Eigen::Map<const VectorT<T>> lambdaGW_vec(lambdaGW.data(),
                                                        lambdaGW.size());
        const MatrixT<T> dP1Yr =
            P1_.cast<T>().asDiagonal() * moving_transformed;
        const Eigen::Map<const VectorT<T>> dP1Yr_vec(dP1Yr.data(),
                                                     dP1Yr.size());
        const VectorT<T> res = (F * jacobianG_.cast<T>()).transpose() *
                                   (-PX_vec_.cast<T>() + dP1Yr_vec) +
                               lambdaGW_vec;
        for (Eigen::Index i = 0; i < res.rows(); i++)
            residual[i] = res(i);
        return true;
    }

private:
    MatrixX3 moving_;
    Matrix G_;
    Matrix jacobianG_;
    Vector PX_vec_;
    Vector P1_;
    VectorInt line_sizes_;
    double lambda_;
};
} // namespace lnrr