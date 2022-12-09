#include <gtest/gtest.h>
#include <lnrr/scan_to_model.h>
#include <lnrr/utils/test_utils.h>

using namespace lnrr;

struct Data {
    std::vector<RigidTransform<double>> T;
    MatrixX3 X;
    MatrixX3 Y;
    Probabilities P;
};

const double LAMBDA = 1.0e-30;
const double BETA = 1e-6;
const int L = 3;
const Matrix G = Matrix::Identity(L, L); // If G is identity, then T = W
const int NUMBER_POINTS_LINE = 3;
const int M = L * NUMBER_POINTS_LINE;
const VectorInt LINE_SIZES = VectorInt::Ones(L) * NUMBER_POINTS_LINE;

// No rotation
Data data_translation0 = {
    generateRandomTransformations(L, 1.0e-3, 0.0),
    MatrixX3::Random(M, 3),
    computeTransformedMoving(data_translation0.X, data_translation0.T,
                             LINE_SIZES),
    {Vector::Ones(M), Vector::Ones(M), data_translation0.X, 0.0}};
Data data_translation1 = {
    generateRandomTransformations(L, 1.0, 0.0),
    MatrixX3::Random(M, 3),
    computeTransformedMoving(data_translation1.X, data_translation1.T,
                             LINE_SIZES),
    {Vector::Ones(M), Vector::Ones(M), data_translation1.X, 0.0}};
Data data_translation2 = {
    generateRandomTransformations(L, 10.0, 0.0),
    MatrixX3::Random(M, 3),
    computeTransformedMoving(data_translation2.X, data_translation2.T,
                             LINE_SIZES),
    {Vector::Ones(M), Vector::Ones(M), data_translation2.X, 0.0}};

// No translation
Data data_rotation0 = {
    generateRandomTransformations(L, 0.0, 0.1),
    MatrixX3::Random(M, 3),
    computeTransformedMoving(data_rotation0.X, data_rotation0.T, LINE_SIZES),
    {Vector::Ones(M), Vector::Ones(M), data_rotation0.X, 0.0}};
Data data_rotation1 = {
    generateRandomTransformations(L, 0.0, 1.0),
    MatrixX3::Random(M, 3),
    computeTransformedMoving(data_rotation1.X, data_rotation1.T, LINE_SIZES),
    {Vector::Ones(M), Vector::Ones(M), data_rotation1.X, 0.0}};
Data data_rotation2 = {
    generateRandomTransformations(L, 0.0, 50.0),
    MatrixX3::Random(M, 3),
    computeTransformedMoving(data_rotation2.X, data_rotation2.T, LINE_SIZES),
    {Vector::Ones(M), Vector::Ones(M), data_rotation2.X, 0.0}};

// Both
Data data_both0 = {
    generateRandomTransformations(L, 1e-3, 0.1),
    MatrixX3::Random(M, 3),
    computeTransformedMoving(data_both0.X, data_both0.T, LINE_SIZES),
    {Vector::Ones(M), Vector::Ones(M), data_both0.X, 0.0}};
Data data_both1 = {
    generateRandomTransformations(L, 0.10, 1.0),
    MatrixX3::Random(M, 3),
    computeTransformedMoving(data_both1.X, data_both1.T, LINE_SIZES),
    {Vector::Ones(M), Vector::Ones(M), data_both1.X, 0.0}};
Data data_both2 = {
    generateRandomTransformations(L, 1.0, 5.0),
    MatrixX3::Random(M, 3),
    computeTransformedMoving(data_both2.X, data_both2.T, LINE_SIZES),
    {Vector::Ones(M), Vector::Ones(M), data_both2.X, 0.0}};

std::vector<Data> data_translation = {data_translation0, data_translation1,
                                      data_translation2};
std::vector<Data> data_rotation = {data_rotation0, data_rotation1,
                                   data_rotation2};
std::vector<Data> data_both = {data_both0, data_both1, data_both2};

///////////////////
////// TESTS //////
///////////////////

void test_computeW(const std::vector<Data>& data,
                   const bool& zero_translation) {
    for (auto d : data) {
        ScanToModel scan_to_model(d.X, d.Y, BETA, LAMBDA, LINE_SIZES);
        scan_to_model.initialize();
        scan_to_model.setP(d.P);
        Probabilities P = scan_to_model.getP();
        try {
            scan_to_model.computeW();
        } catch (const lnrr_error& e) {
            std::cerr << e.what() << '\n';
        }

        std::vector<lnrr::RigidTransform<double>> T_computed =
            computeTransformations(G, scan_to_model.getW());
        std::cout << "W:\n" << scan_to_model.getW() << std::endl;
        std::vector<lnrr::RigidTransform<double>> T_inverse =
            invertTransformations(T_computed);
        for (size_t i = 0; i < d.T.size(); i++) {
            if (zero_translation) {
                EXPECT_TRUE(d.T[i].xyz.isMuchSmallerThan(1e-3))
                    << "i:" << i << std::endl
                    << "Expected translation:\n"
                    << d.T[i].xyz << std::endl
                    << "Returned translation:\n"
                    << T_inverse[i].xyz << std::endl;
            } else {
                EXPECT_TRUE(d.T[i].xyz.isApprox(T_inverse[i].xyz, 1e-3))
                    << "i:" << i << std::endl
                    << "Expected translation:\n"
                    << d.T[i].xyz << std::endl
                    << "Returned translation:\n"
                    << T_inverse[i].xyz << std::endl;
            }

            EXPECT_TRUE(d.T[i].rot_mat.isApprox(T_inverse[i].rot_mat, 1e-3))
                << "i:" << i << std::endl
                << "Expected rot_mat:\n"
                << d.T[i].rot_mat << std::endl
                << "Returned rot_mat:\n"
                << T_inverse[i].rot_mat << std::endl;
        }
    }
}

TEST(Tests_M_Step, computeW_translation) {
    test_computeW(data_translation, false);
}
TEST(Tests_M_Step, computeW_rotation) { test_computeW(data_rotation, true); }
TEST(Tests_M_Step, computeW_both) { test_computeW(data_both, false); }
