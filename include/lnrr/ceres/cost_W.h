#pragma once

#include <ceres/ceres.h>
#include <lnrr/utils/conversions.h>
#include <lnrr/utils/operations.h>
#include <lnrr/utils/types.h>

namespace lnrr {
class CostFunctionW {
public:
    CostFunctionW(const Matrix& moving, const Matrix& G,
                  const Sparse& jacobianG, const Vector& PX_vec,
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
#ifdef DEBUG
        auto tic = std::chrono::high_resolution_clock::now();
#endif
        int number_lines_scan = G_.rows();
        const Eigen::Map<const MatrixX6T<T>> W_matrix(parameters[0],
                                                      number_lines_scan, 6);
        MatrixT<T> G_T = G_.cast<T>();
#ifdef DEBUG
        auto toc = std::chrono::high_resolution_clock::now();
        double time_A =
            std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                .count();

        tic = std::chrono::high_resolution_clock::now();
#endif
        // This line is the second bottleneck (10% - 40%)
        std::vector<RigidTransform<T>> T_lines =
            computeTransformations(G_T, MatrixX6T<T>(W_matrix));
#ifdef DEBUG
        toc = std::chrono::high_resolution_clock::now();
        double time_B =
            std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                .count();

        tic = std::chrono::high_resolution_clock::now();
#endif
        MatrixX3T<T> moving_T = moving_.cast<T>();
#ifdef DEBUG
        toc = std::chrono::high_resolution_clock::now();
        double time_C =
            std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                .count();

        tic = std::chrono::high_resolution_clock::now();
#endif
        SparseT<T> F =
            computeJacobianPointComposition(T_lines, moving_T, line_sizes_);
#ifdef DEBUG
        toc = std::chrono::high_resolution_clock::now();
        double time_D =
            std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                .count();

        tic = std::chrono::high_resolution_clock::now();
#endif
        MatrixX3T<T> moving_transformed =
            computeTransformedMoving(moving_T, T_lines, line_sizes_);
#ifdef DEBUG
        toc = std::chrono::high_resolution_clock::now();
        double time_E =
            std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                .count();

        tic = std::chrono::high_resolution_clock::now();
#endif
        MatrixX6T<T> lambdaGW = T(lambda_) * G_T * W_matrix;
#ifdef DEBUG
        toc = std::chrono::high_resolution_clock::now();
        double time_F =
            std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                .count();

        tic = std::chrono::high_resolution_clock::now();
#endif
        const Eigen::Map<const VectorT<T>> lambdaGW_vec(lambdaGW.data(),
                                                        lambdaGW.size());
#ifdef DEBUG
        toc = std::chrono::high_resolution_clock::now();
        double time_G =
            std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                .count();

        tic = std::chrono::high_resolution_clock::now();
#endif
        const MatrixT<T> dP1Yr =
            P1_.cast<T>().asDiagonal() * moving_transformed;
#ifdef DEBUG
        toc = std::chrono::high_resolution_clock::now();
        double time_H =
            std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                .count();

        tic = std::chrono::high_resolution_clock::now();
#endif
        const Eigen::Map<const VectorT<T>> dP1Yr_vec(dP1Yr.data(),
                                                     dP1Yr.size());
#ifdef DEBUG
        toc = std::chrono::high_resolution_clock::now();
        double time_I =
            std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                .count();

        tic = std::chrono::high_resolution_clock::now();
#endif
        // This matrix multiplication is the bottleneck (70% - 90%)
        SparseT<T> aux0 = (F * jacobianG_.cast<T>()).transpose();
#ifdef DEBUG
        toc = std::chrono::high_resolution_clock::now();
        double time_J =
            std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                .count();

        tic = std::chrono::high_resolution_clock::now();
#endif
        VectorT<T> aux1 = (-PX_vec_.cast<T>() + dP1Yr_vec);
#ifdef DEBUG
        toc = std::chrono::high_resolution_clock::now();
        double time_K =
            std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                .count();

        tic = std::chrono::high_resolution_clock::now();
#endif
        const VectorT<T> res = aux0 * aux1 + lambdaGW_vec;
#ifdef DEBUG
        toc = std::chrono::high_resolution_clock::now();
        double time_L =
            std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
                .count();

        double time_total = time_A + time_B + time_C + time_D + time_E +
                            time_F + time_G + time_H + time_I + time_J +
                            time_K + time_L;
#if true
        std::cout << "****************" << std::endl;
        std::cout << " Times cost_W.h " << std::endl;
        std::cout << "****************" << std::endl;
        std::cout << "time_A: " << time_A << " (" << time_A * 100 / time_total
                  << "%)" << std::endl;
        std::cout << "time_B: " << time_B << " (" << time_B * 100 / time_total
                  << "%)" << std::endl;
        std::cout << "time_C: " << time_C << " (" << time_C * 100 / time_total
                  << "%)" << std::endl;
        std::cout << "time_D: " << time_D << " (" << time_D * 100 / time_total
                  << "%)" << std::endl;
        std::cout << "time_E: " << time_E << " (" << time_E * 100 / time_total
                  << "%)" << std::endl;
        std::cout << "time_F: " << time_F << " (" << time_F * 100 / time_total
                  << "%)" << std::endl;
        std::cout << "time_G: " << time_G << " (" << time_G * 100 / time_total
                  << "%)" << std::endl;
        std::cout << "time_H: " << time_H << " (" << time_H * 100 / time_total
                  << "%)" << std::endl;
        std::cout << "time_I: " << time_I << " (" << time_I * 100 / time_total
                  << "%)" << std::endl;
        std::cout << "time_J: " << time_J << " (" << time_J * 100 / time_total
                  << "%)" << std::endl;
        std::cout << "time_K: " << time_K << " (" << time_K * 100 / time_total
                  << "%)" << std::endl;
        std::cout << "time_L: " << time_L << " (" << time_L * 100 / time_total
                  << "%)" << std::endl;
        std::cout << "     ........    " << std::endl;
#endif
#endif

        for (Eigen::Index i = 0; i < res.rows(); i++)
            residual[i] = res(i);
        return true;
    }

private:
    MatrixX3 moving_;
    Matrix G_;
    Sparse jacobianG_;
    Vector PX_vec_;
    Vector P1_;
    VectorInt line_sizes_;
    double lambda_;
};
} // namespace lnrr