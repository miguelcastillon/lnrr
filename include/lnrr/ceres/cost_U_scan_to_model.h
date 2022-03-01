#pragma once

#include <ceres/ceres.h>
#include <lnrr/scan_to_model.h>
#include <lnrr/utils/conversions.h>
#include <lnrr/utils/types.h>

namespace lnrr {
class CostFunctionScanToModelRot {
public:
    CostFunctionScanToModelRot(const Matrix& G, const Matrix& C, //
                               const Matrix& D, const double& lambda)
        : G_(G), C_(C), D_(D), lambda_(lambda) {}

    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, 3> computeRotationMatrices(
        const Eigen::Matrix<T, Eigen::Dynamic, 3>& GU) const {
        Eigen::Matrix<T, Eigen::Dynamic, 3> RT =
            Eigen::Matrix<T, Eigen::Dynamic, 3>::Zero(3 * GU.rows(), 3);
        for (Eigen::Index l = 0; l < GU.rows(); l++) {
            /* RPY */
            RT.block(3 * l, 0, 3, 3) =
                conversions::rpy2RotationMatrix(GU(l, 0), GU(l, 1), GU(l, 2))
                    .transpose();
        }
        return RT;
    }

    template <typename T>
    Eigen::Matrix<T, 3, 9> jacobianRotationMatrix(const T& roll, const T& pitch,
                                                  const T& yaw) const {
        T cx = cos(roll);
        T sx = sin(roll);
        T cy = cos(pitch);
        T sy = sin(pitch);
        T cz = cos(yaw);
        T sz = sin(yaw);

        // Column-major order
        Eigen::Matrix<T, 3, 9> output;
        output << T(0), T(0), T(0),                                  // wrt roll
            cz * sy * cx + sz * sx, sz * sy * cx - cz * sx, cy * cx, //
            -cz * sy * sx + sz * cx, -sz * sy * sx - cz * cx, -cy * sx, //
                                                                        //
            -cz * sy, -sz * sy, -cy,                               // wrt pitch
            cz * cy * sx, sz * cy * sx, -sy * sx,                  //
            cz * cy * cx, sz * cy * cx, -sy * cx,                  //
                                                                   //
            -sz * cy, cz * cy, T(0),                               // wrt yaw
            -sz * sy * sx - cz * cx, cz * sy * sx - sz * cx, T(0), //
            -sz * sy * cx + cz * sx, cz * sy * cx + sz * sx, T(0);

        return output;
    }

    template <typename T>
    Eigen::SparseMatrix<T> jacobianRotationMatrices(
        const Eigen::Matrix<T, Eigen::Dynamic, 3>& GU) const {
        int number_rows = GU.rows();
        Eigen::SparseMatrix<T> output(3 * number_rows, 9 * number_rows);
        std::vector<Eigen::Triplet<T>> tripletList;
        tripletList.reserve(3 * 9 * number_rows);

        for (int line = 0; line < number_rows; line++) {
            Eigen::Matrix<T, 3, 9> JR =
                jacobianRotationMatrix(GU(line, 0), GU(line, 1), GU(line, 2));
            Eigen::Matrix<T, 1, 9> JRx = JR.row(0);
            Eigen::Matrix<T, 1, 9> JRy = JR.row(1);
            Eigen::Matrix<T, 1, 9> JRz = JR.row(2);
            for (Eigen::Index j = 0; j < JR.cols(); j++) {
                tripletList.push_back(
                    Eigen::Triplet<T>(line, 9 * line + j, JRx(j)));
                tripletList.push_back(Eigen::Triplet<T>(number_rows + line,
                                                        9 * line + j, JRy(j)));
                tripletList.push_back(Eigen::Triplet<T>(2 * number_rows + line,
                                                        9 * line + j, JRz(j)));
            }
        }
        output.setFromTriplets(tripletList.begin(), tripletList.end());
        return output;
    }

    template <typename T>
    bool operator()(T const* const* parameters, T* residual) const {
        const int number_lines_scan = G_.rows();
        const Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 3>> U_matrix(
            parameters[0], number_lines_scan, 3);

        Eigen::Matrix<T, Eigen::Dynamic, 3> GU = G_.cast<T>() * U_matrix;
        Eigen::Matrix<T, Eigen::Dynamic, 3> lambdaGU = lambda_ * GU;
        const Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>>
            lambdaGU_vec(lambdaGU.data(), lambdaGU.size());

        Eigen::Matrix<T, 3, Eigen::Dynamic> R =
            computeRotationMatrices(GU).transpose();
        Eigen::Matrix<T, 3, Eigen::Dynamic> RCD =
            R * C_.cast<T>() + D_.cast<T>();
        const Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> RCD_vec(
            RCD.data(), RCD.size());

        Eigen::SparseMatrix<T> JR = jacobianRotationMatrices(GU);
        Eigen::Matrix<T, Eigen::Dynamic, 1> res = JR * RCD_vec + lambdaGU_vec;

        // residual = res.data();
        for (Eigen::Index i = 0; i < res.rows(); i++)
            residual[i] = res(i);

        return true;
    }

private:
    Matrix G_;
    Matrix C_;
    Matrix D_;
    double lambda_;
};
} // namespace lnrr
