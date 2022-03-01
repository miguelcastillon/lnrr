#pragma once

#include <lnrr/utils/conversions.h>

namespace lnrr {

inline MatrixX3 computeRotationMatrices(const Matrix& G, const MatrixX3& U) {
    MatrixX3 RT = MatrixX3::Zero(3 * G.rows(), 3);
    MatrixX3 GU = G * U;
    for (Eigen::Index l = 0; l < G.rows(); l++)
        RT.block(3 * l, 0, 3, 3) =
            conversions::rpy2RotationMatrix(GU(l, 0), GU(l, 1), GU(l, 2))
                .transpose();
    return RT;
}

inline SparseMatrix matrixAsSparseBlockDiag(const Matrix& input) {
    SparseMatrix output(input.rows(), input.rows() * input.cols());
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(input.rows() * input.cols());
    for (Eigen::Index i = 0; i < input.rows(); i++)
        for (Eigen::Index j = 0; j < input.cols(); j++)
            tripletList.push_back(
                Eigen::Triplet<double>(i, input.cols() * i + j, input(i, j)));
    output.setFromTriplets(tripletList.begin(), tripletList.end());
    return output;
}

inline Matrix computeG(const double& beta, const int& number_lines) {
    double k = -2.0 * beta * beta;
    Matrix G = Matrix(number_lines, number_lines);
    for (int i = 0; i < number_lines; i++)
        for (int j = 0; j < number_lines; j++)
            G(i, j) = std::abs(i - j);
    G = (G.array().pow(2) / k).exp();
    return G;
}

inline SparseMatrix computeF(const int& M, const Vector& line_sizes) {
    SparseMatrix F(M, line_sizes.rows());
    int m = 0;
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(M);
    for (Eigen::Index col = 0; col < line_sizes.rows(); col++)
        for (size_t i = 0; i < line_sizes[col]; i++) {
            tripletList.push_back(Eigen::Triplet<double>(m, col, 1.0));
            m++;
        }
    F.setFromTriplets(tripletList.begin(), tripletList.end());
    return F;
}

inline SparseMatrix computeH(const int& M, const Vector& line_sizes) {
    SparseMatrix H(3 * M, 3 * line_sizes.rows());
    int m = 0;
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(3 * M);
    for (Eigen::Index col = 0; col < line_sizes.rows(); col++)
        for (size_t i = 0; i < line_sizes[col]; i++) {
            tripletList.push_back(Eigen::Triplet<double>(m, 3 * col, 1.0));
            m++;
            tripletList.push_back(Eigen::Triplet<double>(m, 3 * col + 1, 1.0));
            m++;
            tripletList.push_back(Eigen::Triplet<double>(m, 3 * col + 2, 1.0));
            m++;
        }
    H.setFromTriplets(tripletList.begin(), tripletList.end());
    return H;
}

} // namespace lnrr