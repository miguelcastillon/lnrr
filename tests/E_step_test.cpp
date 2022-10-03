#include <gtest/gtest.h>
#include <lnrr/scan_to_model.h>

using namespace lnrr;

struct Data {
    MatrixX3 X;
    MatrixX3 Y;
    Matrix P;
    double w;
    double sigma;
    double distance_threshold;
};

Data data0 = {MatrixX3::Random(100, 3),
              data0.X,
              Matrix::Identity(data0.Y.rows(), data0.X.rows()),
              0.0,
              1.0,
              0.0001};
Data data1 = {(MatrixX3(10, 3) << 0, 0.1, 0.2, 1, 0.1, 0.2, 2, 0.1, 0.2, 3, 0.1,
               0.2, 4, 0.1, 0.2, 5, 0.1, 0.2, 6, 0.1, 0.2, 7, 0.1, 0.2, 8, 0.1,
               0.2, 9, 0.1, 0.2)
                  .finished(),
              (MatrixX3(8, 3) << 0, 1, 0, 1, 1, 1, 2, 1, 2, 3, 1, 3, 4, 1, 3, 5,
               1, 2, 6, 1, 1, 7, 1, 0)
                  .finished(),
              (Matrix(8, 10) << 2.76376098e-01, 1.67818281e-01, 4.59473521e-02,
               4.54614881e-03, 1.37281799e-04, 0.00000000e+00, 0.00000000e+00,
               0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.24183786e-01,
               2.04973711e-01, 1.52550581e-01, 4.10290544e-02, 3.36786987e-03,
               8.43733422e-05, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
               0.00000000e+00, 7.55162380e-03, 3.38819265e-02, 6.85453947e-02,
               5.01130002e-02, 1.11817218e-02, 7.61470553e-04, 1.87395640e-05,
               0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 6.21479537e-05,
               7.57964848e-04, 4.16824974e-03, 8.28362322e-03, 5.02427145e-03,
               9.30062233e-04, 6.22175435e-05, 1.87670559e-06, 0.00000000e+00,
               0.00000000e+00, 1.87670559e-06, 6.22175435e-05, 9.30062233e-04,
               5.02427145e-03, 8.28362322e-03, 4.16824974e-03, 7.57964848e-04,
               6.21479537e-05, 2.38225320e-06, 0.00000000e+00, 0.00000000e+00,
               1.87395640e-05, 7.61470553e-04, 1.11817218e-02, 5.01130002e-02,
               6.85453947e-02, 3.38819265e-02, 7.55162380e-03, 7.86857183e-04,
               2.96405030e-05, 0.00000000e+00, 0.00000000e+00, 8.43733422e-05,
               3.36786987e-03, 4.10290544e-02, 1.52550581e-01, 2.04973711e-01,
               1.24183786e-01, 3.51734481e-02, 3.60162990e-03, 0.00000000e+00,
               0.00000000e+00, 0.00000000e+00, 1.37281799e-04, 4.54614881e-03,
               4.59473521e-02, 1.67818281e-01, 2.76376098e-01, 2.12786961e-01,
               5.92275314e-02)
                  .finished(),
              0.1,
              1.0,
              5};

std::vector<Data> data = {data0, data1};

const double BETA = 1.0;
const double LAMBDA = 1.0;

///////////////////
////// TESTS //////
///////////////////

TEST(Tests_E_Step, computeP) {
    for (auto d : data) {
        Vector line_sizes = Vector::Ones(d.Y.rows());
        ScanToModel scan_to_model(d.X, d.Y, BETA, LAMBDA, line_sizes);
        scan_to_model.setThresholdTruncate(d.distance_threshold);
        scan_to_model.setOutlierRate(d.w);
        scan_to_model.setSigma2(pow(d.sigma, 2));
        scan_to_model.initialize();
        scan_to_model.computeP();
        EXPECT_TRUE(scan_to_model.P_.p1.isApprox(d.P.rowwise().sum(), 1e-4))
            << "Expected P\n:" << d.P << std::endl
            << "Expected P1\n:" << d.P.rowwise().sum() << std::endl
            << "Returned P1\n:" << scan_to_model.P_.p1 << std::endl;
        EXPECT_TRUE(scan_to_model.P_.pt1.isApprox(
            d.P.colwise().sum().transpose(), 1e-4))
            << "Expected P\n:" << d.P << std::endl
            << "Expected Pt1\n:" << d.P.colwise().sum().transpose() << std::endl
            << "Returned Pt1\n:" << scan_to_model.P_.pt1 << std::endl;
        EXPECT_TRUE(scan_to_model.P_.px.isApprox(d.P * d.X, 1e-4))
            << "Expected P\n:" << d.P << std::endl
            << "Expected PX\n:" << d.P * d.X << std::endl
            << "Returned PX\n:" << scan_to_model.P_.px << std::endl;
    }
}

TEST(Tests_E_Step, computeP_FGT) {
    for (auto d : data) {
        Vector line_sizes = Vector::Ones(d.Y.rows());
        ScanToModel scan_to_model(d.X, d.Y, BETA, LAMBDA, line_sizes);
        scan_to_model.setThresholdTruncate(d.distance_threshold);
        scan_to_model.setOutlierRate(d.w);
        scan_to_model.setSigma2(pow(d.sigma, 2));
        scan_to_model.initialize();
        scan_to_model.computeP_FGT();
        EXPECT_TRUE(scan_to_model.P_.p1.isApprox(d.P.rowwise().sum(), 1e-4))
            << "Expected P\n:" << d.P << std::endl
            << "Expected P1\n:" << d.P.rowwise().sum() << std::endl
            << "Returned P1\n:" << scan_to_model.P_.p1 << std::endl;
        EXPECT_TRUE(scan_to_model.P_.pt1.isApprox(
            d.P.colwise().sum().transpose(), 1e-4))
            << "Expected P\n:" << d.P << std::endl
            << "Expected Pt1\n:" << d.P.colwise().sum().transpose() << std::endl
            << "Returned Pt1\n:" << scan_to_model.P_.pt1 << std::endl;
        EXPECT_TRUE(scan_to_model.P_.px.isApprox(d.P * d.X, 1e-4))
            << "Expected P\n:" << d.P << std::endl
            << "Expected PX\n:" << d.P * d.X << std::endl
            << "Returned PX\n:" << scan_to_model.P_.px << std::endl;
    }
}