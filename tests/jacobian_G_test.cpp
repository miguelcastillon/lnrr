#include <gtest/gtest.h>
#include <lnrr/utils/operations.h>
#include <lnrr/utils/types.h>

using namespace lnrr;

Matrix computeT(const Matrix& G, const MatrixX6& W, const double& epsilon) {
    return G * (W + epsilon * Matrix::Ones(W.rows(), W.cols()));
}

const int L = 10;
const double BETA = 1.0;
const double EPSILON = 1e-3;

struct Data {
    MatrixX6 W0;
};

Data data0 = {MatrixX6::Random(L, 6)};
Data data1 = {MatrixX6::Random(L, 6)};
Data data2 = {MatrixX6::Random(L, 6)};
std::vector<Data> data = {data0, data1, data2};

TEST(Tests_Jacobian_G, ComputeJacobianG) {
    for (auto d : data) {
        Matrix G = computeG(BETA, L);
        Matrix T0 = computeT(G, d.W0, 0.0);
        Matrix T1 = computeT(G, d.W0, EPSILON);
        Matrix jTW = computeJacobianG(G);
        Vector epsilon_vector = Vector::Ones(6 * L) * EPSILON;
        Vector jTW_eps_vec = jTW * epsilon_vector;
        Eigen::Map<const Matrix> jTW_eps(&jTW_eps_vec[0], L, 6);
        Matrix T1_approx = T0 + jTW_eps;
        EXPECT_TRUE(T1.isApprox(T1_approx, 1e-6));
    }
}
