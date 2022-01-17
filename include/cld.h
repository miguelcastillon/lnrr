#pragma once

#include <ceres/cost_cld_U.h>
#include <utils/conversions.h>
#include <utils/types.h>
#include <utils/utils.h>

#include <fgt.hpp>

namespace cld {
// const size_t DEFAULT_MAX_ITERATIONS = 150;
// const double DEFAULT_OUTLIERS = 0.1;
// const double DEFAULT_THRESHOLD_TRUNCATE = 1e6;
// const double DEFAULT_TOLERANCE = 1e-5;
// const double DEFAULT_SIGMA2 = 0.0;
// const double DEFAULT_BETA = 3.0;
// const double DEFAULT_LAMBDA = 3.0;
const double FGT_EPSILON = 1e-4;
const double FGT_THRESHOLD = 0.2;

class CLD {
private:
    MatrixX3 fixed_;
    MatrixX3 moving_;
    MatrixX3 moving_transformed_;
    int M_;
    Matrix G_;
    MatrixX3 C_;
    SparseMatrix D_;
    MatrixX3 U_; // rotation
    MatrixX3 V_; // translation
    SparseMatrix F_;
    Matrix FG_;
    SparseMatrix FT_;
    SparseMatrix H_;
    SparseMatrix YD_;
    Probabilities P_;
    MatrixX3 RT_; // tranpose of the rotation matrices

    MatrixX3 A_;
    Matrix B_;
    Matrix GB_;

    double fgt_epsilon_ = FGT_EPSILON;
    double fgt_threshold_ = FGT_THRESHOLD;

    size_t max_iterations_;
    double outliers_;
    double threshold_truncate_;
    double sigma2_;
    double tolerance_;
    double lambda_;
    double beta_;
    Vector line_sizes_;
    int number_lines_;

    double defaultSigma2();
    void computeSigma2();
    // void computeG();
    // void computeF();
    // void computeH();
    void computeP();
    void computeP_FGT();
    void computeU();
    double computeOptimalRotationCeres(const Matrix& S, const Matrix& T);

    // void computeRotationMatrices();
    SparseMatrix matrixAsSparseBlockDiag(const Matrix& input);

public:
    CLD(const size_t& max_iterations, const double& outliers,
        const double& threshold_truncate, const double& sigma2,
        const double& tolerance, const double& lambda, const double& beta,
        const Vector& line_sizes)
        : max_iterations_(max_iterations),
          outliers_(outliers),
          threshold_truncate_(threshold_truncate),
          sigma2_(sigma2),
          tolerance_(tolerance),
          lambda_(lambda),
          beta_(beta),
          line_sizes_(line_sizes),
          number_lines_(line_sizes.rows()) {}
    ~CLD() {}

    void initialize();
    void computeOne();
    Matrix getTransformedMoving();
    Result run();

    // Set functions

    void setFixed(const MatrixX3& fixed) { fixed_ = fixed; }
    void setMoving(const MatrixX3& moving) { moving_ = moving; }
    void setMaxIterations(const size_t& max_iterations) {
        max_iterations_ = max_iterations;
    }
    void setOutlierRate(const double& outliers) { outliers_ = outliers; }
    void setThresholdTruncate(const double& threshold) {
        threshold_truncate_ = threshold;
    }
    void setSigma2(const double& sigma2) { sigma2_ = sigma2; }
    void setTolerance(const double& tolerance) { tolerance_ = tolerance; }
    void setLambda(const double& lambda) { lambda_ = lambda; }
    void setLineSizes(const Vector& line_sizes) {
        line_sizes_ = line_sizes;
        number_lines_ = line_sizes.rows();
    }

    // Get functions

    MatrixX3 getFixed() { return fixed_; }
    MatrixX3 getMoving() { return moving_; }
    size_t getMaxIterations() { return max_iterations_; }
    double getOutlierRate() { return outliers_; }
    double getThresholdTruncate() { return threshold_truncate_; }
    double getSigma2() { return sigma2_; }
    double getTolerance() { return tolerance_; }
    double getLambda() { return lambda_; }
    Vector getLineSizes() { return line_sizes_; }
    int getNumberLines() { return number_lines_; }
};
} // namespace cld