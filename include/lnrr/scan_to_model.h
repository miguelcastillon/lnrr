#pragma once

#include <lnrr/ceres/cost_W.h>
#include <lnrr/utils/conversions.h>
#include <lnrr/utils/operations.h>
#include <lnrr/utils/types.h>

#include <fgt.hpp>

namespace lnrr {
const size_t DEFAULT_MAX_ITERATIONS = 150;
const double DEFAULT_OUTLIERS = 0.1;
const double DEFAULT_THRESHOLD_TRUNCATE = 1e6;
const double DEFAULT_TOLERANCE = 1e-5;
const double DEFAULT_SIGMA2 = 0.0;
const double FGT_EPSILON = 1e-4;
const double FGT_THRESHOLD = 0.2;

class ScanToModel {
private:
    MatrixX3 fixed_;
    MatrixX3 moving_;
    MatrixX3 moving_transformed_;

    Probabilities P_;

    Matrix G_;
    MatrixX6 W_;

    double fgt_epsilon_ = FGT_EPSILON;
    double fgt_threshold_ = FGT_THRESHOLD;

    double beta_;
    double lambda_;
    VectorInt line_sizes_;
    int number_lines_;
    size_t max_iterations_;
    double outliers_;
    double threshold_truncate_2_;
    double sigma2_;
    double tolerance_;

    double defaultSigma2();
    void computeSigma2();

public:
    ScanToModel(const MatrixX3& fixed, const MatrixX3& moving,
                const double& beta, const double& lambda,
                const VectorInt& line_sizes)
        : fixed_(fixed),
          moving_(moving),
          beta_(beta),
          lambda_(lambda),
          line_sizes_(line_sizes),
          number_lines_(line_sizes.rows()),
          max_iterations_(DEFAULT_MAX_ITERATIONS),
          outliers_(DEFAULT_OUTLIERS),
          threshold_truncate_2_(DEFAULT_THRESHOLD_TRUNCATE),
          sigma2_(DEFAULT_SIGMA2),
          tolerance_(DEFAULT_TOLERANCE) {
        assert(line_sizes.sum() == moving.rows() &&
               "The vector line_sizes must agree with the number of points of "
               "each line in the scan");
    }

    ~ScanToModel() {}

    void initialize();
    void computeP();
    void computeP_FGT();
    void computeW();
    void computeOne();
    MatrixX3 getTransformedMoving() {
        return computeTransformedMoving(moving_, G_, W_, line_sizes_);
    }
    Result run();

    // Set functions
    void setP(const Probabilities& P) { P_ = P; }
    void setFixed(const MatrixX3& fixed) { fixed_ = fixed; }
    void setMoving(const MatrixX3& moving) { moving_ = moving; }
    void setLineSizes(const VectorInt& line_sizes) {
        line_sizes_ = line_sizes;
        number_lines_ = line_sizes.rows();
    }
    void setBeta(const double& beta) { beta_ = beta; }
    void setLambda(const double& lambda) { lambda_ = lambda; }

    void setMaxIterations(const size_t& max_iterations) {
        max_iterations_ = max_iterations;
    }
    void setOutlierRate(const double& outliers) { outliers_ = outliers; }
    void setThresholdTruncate(const double& threshold) {
        threshold_truncate_2_ = pow(threshold, 2);
    }
    void setSigma2(const double& sigma2) { sigma2_ = sigma2; }
    void setTolerance(const double& tolerance) { tolerance_ = tolerance; }

    // Get functions
    Probabilities getP() { return P_; }
    MatrixX6 getW() { return W_; }
    MatrixX3 getFixed() { return fixed_; }
    MatrixX3 getMoving() { return moving_; }
    VectorInt getLineSizes() { return line_sizes_; }
    int getNumberLines() { return number_lines_; }
    double getBeta() { return beta_; }
    double getLambda() { return lambda_; }
    size_t getMaxIterations() { return max_iterations_; }
    double getOutlierRate() { return outliers_; }
    double getThresholdTruncate() { return threshold_truncate_2_; }
    double getSigma2() { return sigma2_; }
    double getTolerance() { return tolerance_; }
};
} // namespace lnrr