#include <lnrr_se3/scan_to_model.h>

namespace lnrr {

void ScanToModel::computeP() {
    double ksig = -2.0 * sigma2_;
    size_t cols = fixed_.cols();
    double outlier_tmp = (outliers_ * moving_transformed_.rows() *
                          std::pow(-ksig * M_PI, 0.5 * cols)) /
                         ((1 - outliers_) * fixed_.rows());
    Vector p = Vector::Zero(moving_transformed_.rows());
    Vector p1 = Vector::Zero(moving_transformed_.rows());
    Vector p1_max = Vector::Zero(moving_transformed_.rows());
    Vector pt1 = Vector::Zero(fixed_.rows());
    Matrix px = Matrix::Zero(moving_transformed_.rows(), cols);
    double l = 0.0;

    for (Matrix::Index i = 0; i < fixed_.rows(); ++i) {
        double sp = 0;
        for (Matrix::Index j = 0; j < moving_transformed_.rows(); ++j) {
            double razn = (fixed_.row(i) - moving_transformed_.row(j))
                              .array()
                              .pow(2)
                              .sum();
            if (razn >= threshold_truncate_2_) {
                p(j) = 0.0;
            } else {
                p(j) = std::exp(razn / ksig);
                sp += p(j);
            }
        }
        sp += outlier_tmp;
        pt1(i) = 1 - outlier_tmp / sp;
        for (Matrix::Index j = 0; j < moving_transformed_.rows(); ++j) {
            p1(j) += p(j) / sp;
            px.row(j) += fixed_.row(i) * p(j) / sp;
            if (p(j) / sp > p1_max(j)) {
                p1_max(j) = p(j) / sp;
            }
        }
        l += -std::log(sp);
    }
    l += cols * fixed_.rows() * std::log(sigma2_) / 2;
    P_ = {p1, pt1, px, l};
}

void ScanToModel::computeP_FGT() {
    double bandwidth = std::sqrt(2.0 * sigma2_);
    size_t cols = fixed_.cols();
    std::unique_ptr<fgt::Transform> transform;
    transform = std::unique_ptr<fgt::Transform>(
        new fgt::DirectTree(moving_transformed_, bandwidth, fgt_epsilon_));
    auto kt1 = transform->compute(fixed_, threshold_truncate_2_);
    double ndi = outliers_ / (1.0 - outliers_) * moving_transformed_.rows() /
                 fixed_.rows() * std::pow(2.0 * M_PI * sigma2_, 0.5 * cols);
    Eigen::ArrayXd denom_p = kt1.array() + ndi;
    Vector pt1 = 1 - ndi / denom_p;

    transform = std::unique_ptr<fgt::Transform>(
        new fgt::DirectTree(fixed_, bandwidth, fgt_epsilon_));
    Vector p1 = transform->compute(moving_transformed_, 1 / denom_p,
                                   threshold_truncate_2_);
    Matrix px(moving_transformed_.rows(), cols);
    for (size_t i = 0; i < cols; ++i) {
        px.col(i) = transform->compute(moving_transformed_,
                                       fixed_.col(i).array() / denom_p,
                                       threshold_truncate_2_);
    }
    double l =
        -denom_p.log().sum() + cols * fixed_.rows() * std::log(sigma2_) / 2;
    P_ = {p1, pt1, px, l};
}

/// TODO: Are we sure this sigma equals to equation in Fig. 2 of [Myronenko
/// 2010] ??
double ScanToModel::defaultSigma2() {
    return ((moving_.rows() * (fixed_.transpose() * fixed_).trace()) +
            (fixed_.rows() * (moving_.transpose() * moving_).trace()) -
            2 * fixed_.colwise().sum() * moving_.colwise().sum().transpose()) /
           (fixed_.rows() * moving_.rows() * fixed_.cols()) /
           1000.; // 1000 here is a magic number
}

void ScanToModel::computeW() {
    ceres::Problem problem;
    ceres::LossFunction* loss = NULL;

    const Eigen::Map<const Vector> PX_vec(P_.px.data(), P_.px.size());

    Matrix jacobianG = computeJacobianG(G_);

    ceres::DynamicAutoDiffCostFunction<CostFunctionW>* cost =
        new ceres::DynamicAutoDiffCostFunction<CostFunctionW>(
            new CostFunctionW(moving_, G_, jacobianG / sigma2_, PX_vec, P_.p1,
                              line_sizes_, lambda_));
    std::vector<double*> parameter_blocks;
    parameter_blocks.push_back(W_.data());
    cost->AddParameterBlock(W_.size());
    cost->SetNumResiduals(W_.size());
    problem.AddResidualBlock(cost, loss, parameter_blocks);
    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.check_gradients = false;
    // options.minimizer_type = ceres::LINE_SEARCH;
    // options.logging_type = ceres::SILENT;
    // options.function_tolerance = 1e-4;

    ceres::Solver::Summary summary;
    Solve(options, &problem, &summary);

    return;
}

void ScanToModel::computeSigma2() {
    sigma2_ = 0.0;
    sigma2_ += (fixed_.transpose() * P_.pt1.asDiagonal() * fixed_).trace();
    sigma2_ -= 2 * (P_.px.transpose() * moving_transformed_).trace();
    sigma2_ += (moving_transformed_.transpose() * P_.p1.asDiagonal() *
                moving_transformed_)
                   .trace();
    sigma2_ /= (P_.p1.sum() * 3.0);
}

void ScanToModel::initialize() {
    assert(fixed_.rows() > 0);
    assert(moving_.rows() > 0);

    W_ = MatrixX6::Zero(number_lines_, 6);
    G_ = computeG(beta_, number_lines_);

    if (sigma2_ == 0.0)
        sigma2_ = defaultSigma2();
    moving_transformed_ = getTransformedMoving();
}

void ScanToModel::computeOne() {
    // #ifdef MODE_DEBUG
    auto tic = std::chrono::high_resolution_clock::now();
    std::cout << "Computing P..." << std::endl;
    // #endif
    computeP_FGT();
    // #ifdef MODE_DEBUG
    std::cout << "P computed." << std::endl;
    auto toc = std::chrono::high_resolution_clock::now();
    double time_E =
        std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
            .count();
    // #endif

    // #ifdef MODE_DEBUG
    tic = std::chrono::high_resolution_clock::now();
    std::cout << "Computing M..." << std::endl;
    // #endif
    computeW();
    // #ifdef MODE_DEBUG
    std::cout << "M computed." << std::endl;
    toc = std::chrono::high_resolution_clock::now();
    double time_M =
        std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
            .count();
    // #endif

    // #ifdef MODE_DEBUG
    tic = std::chrono::high_resolution_clock::now();
    std::cout << "Computing other..." << std::endl;
    // #endif
    moving_transformed_ = getTransformedMoving();
    computeSigma2();
    // #ifdef MODE_DEBUG
    std::cout << "Other computed." << std::endl;
    toc = std::chrono::high_resolution_clock::now();
    double time_other =
        std::chrono::duration_cast<std::chrono::microseconds>(toc - tic)
            .count();

    double time_total = time_E + time_M + time_other;
    std::cout << "-----" << std::endl;
    std::cout << "time_E: " << time_E << " (" << time_E * 100 / time_total
              << "%)" << std::endl;
    std::cout << "time_M: " << time_M << " (" << time_M * 100 / time_total
              << "%)" << std::endl;
    std::cout << "time_other: " << time_other << " ("
              << time_other * 100 / time_total << "%)" << std::endl;
    std::cout << std::endl;
    // #endif
}

Result ScanToModel::run() {
    auto tic = std::chrono::high_resolution_clock::now();

    // #ifdef MODE_DEBUG
    std::cout << "Initializing non-rigid registration in DEBUG mode"
              << std::endl;
    // #endif
    initialize();
    // #ifdef MODE_DEBUG
    std::cout << "Initialization complete" << std::endl;
    // #endif

    size_t iter = 0;

    double l = 0.0;
    double ntol = tolerance_ + 10.0;

    while (iter < max_iterations_ && ntol > tolerance_ &&
           sigma2_ > 10 * std::numeric_limits<double>::epsilon()) {
        ScanToModel::computeOne();
        ntol = std::abs((P_.l - l) / P_.l);
        l = P_.l;
        ++iter;
        // #ifdef MODE_DEBUG
        std::cout << "Iteration " << iter << " complete. " << std::endl;
        std::cout << "  --> Tolerance = " << ntol << ". ";
        std::cout << "Convergence criteria: " << tolerance_ << std::endl;
        std::cout << "  --> sigma2 = " << sigma2_ << ". ";
        std::cout << "Convergence criteria: "
                  << 10 * std::numeric_limits<double>::epsilon() << std::endl;
        // #endif
    }

    auto toc = std::chrono::high_resolution_clock::now();
    Result result;
    result.runtime =
        std::chrono::duration_cast<std::chrono::microseconds>(toc - tic);
    result.points = moving_transformed_;
    result.sigma2 = sigma2_;
    result.iterations = iter;
    result.line_transforms = computeTransformations(G_, W_);

    // #ifdef MODE_DEBUG
    std::cout << "Non-rigid registration complete: ";
    std::cout << result.iterations << " iterations, ";
    std::cout << result.runtime.count() / 1e6 << " seconds. ";
    std::cout << "Final sigma2 = " << result.sigma2 << std::endl;
    // #endif

    return result;
}

} // namespace lnrr