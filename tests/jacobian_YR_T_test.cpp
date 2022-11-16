#include <gtest/gtest.h>
#include <lnrr_se3/utils/operations.h>
#include <lnrr_se3/utils/types.h>

using namespace lnrr;

Matrix computeT(const Matrix& G, const MatrixX6& W, const double& epsilon) {
    return G * (W + epsilon * Matrix::Ones(W.rows(), W.cols()));
}

const int L = 10;
const int NUMBER_POINTS_LINE = 10;
const int M = L * NUMBER_POINTS_LINE;
const VectorInt LINE_SIZES = VectorInt::Ones(L) * NUMBER_POINTS_LINE;
const Matrix G = Matrix::Identity(L, L); // If G is identity, then T = W
const double BETA = 1.0;
const double EPSILON = 1e-5;

struct Data {
    MatrixX3 moving;
    MatrixX6 W0;
};

Data data0 = {MatrixX3::Random(M, 3), MatrixX6::Random(L, 6)};
Data data1 = {MatrixX3::Random(M, 3), MatrixX6::Random(L, 6)};
Data data2 = {MatrixX3::Random(M, 3), MatrixX6::Random(L, 6)};
std::vector<Data> data = {data0, data1, data2};

TEST(Tests_Jacobian_YR_T, ComputeJacobianYR_T) {
    for (auto d : data) {
        MatrixX6 W1 = d.W0 + EPSILON * Matrix::Ones(d.W0.rows(), d.W0.cols());

        MatrixX3 YR0 = computeTransformedMoving(d.moving, G, d.W0, LINE_SIZES);
        MatrixX3 YR1 = computeTransformedMoving(d.moving, G, W1, LINE_SIZES);
        std::vector<RigidTransform<double>> T_lines =
            computeTransformations(G, d.W0);

        Matrix jYRT =
            computeJacobianPointComposition(T_lines, d.moving, LINE_SIZES);

        Vector epsilon_vector = Vector::Ones(6 * L) * EPSILON;
        Vector jYRT_eps_vec = jYRT * epsilon_vector;

        Eigen::Map<const Matrix> jYRT_eps(&jYRT_eps_vec[0], M, 3);

        MatrixX3 YR1_approx = YR0 + jYRT_eps;
        EXPECT_TRUE(YR1.isApprox(YR1_approx, 1e-4))
            << "YR1 - YR1_approx:\n"
            << YR1 - YR1_approx << std::endl;
    }
}
