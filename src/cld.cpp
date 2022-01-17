#include <cld.h>

namespace cld {

Matrix CLD::getTransformedMoving() {
    RT_ = computeRotationMatrices(G_, U_);
    return C_ + D_ * RT_;
}

void CLD::computeP() {
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
            if (razn >= threshold_truncate_) {
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

void CLD::computeP_FGT() {
    /// TODO: how to use threshold_truncate?
    double bandwidth = std::sqrt(2.0 * sigma2_);
    size_t cols = fixed_.cols();
    std::unique_ptr<fgt::Transform> transform;
    // if (bandwidth > fgt_threshold_)
    // {
    //   std::cout << "IFGT" << std::endl;
    //   transform = std::unique_ptr<fgt::Transform>(new
    //   fgt::Ifgt(moving_transformed_, bandwidth, fgt_epsilon_));
    // }
    // else
    transform = std::unique_ptr<fgt::Transform>(
        new fgt::DirectTree(moving_transformed_, bandwidth, fgt_epsilon_));
    auto kt1 = transform->compute(fixed_, threshold_truncate_);
    double ndi = outliers_ / (1.0 - outliers_) * moving_transformed_.rows() /
                 fixed_.rows() * std::pow(2.0 * M_PI * sigma2_, 0.5 * cols);
    Eigen::ArrayXd denom_p = kt1.array() + ndi;
    Vector pt1 = 1 - ndi / denom_p;

    // if (bandwidth > threshold_truncate)
    //   bandwidth = threshold_truncate;

    // if (bandwidth > fgt_threshold_)
    //   transform = std::unique_ptr<fgt::Transform>(new fgt::Ifgt(fixed_,
    //   bandwidth, fgt_epsilon_));
    // else
    transform = std::unique_ptr<fgt::Transform>(
        new fgt::DirectTree(fixed_, bandwidth, fgt_epsilon_));
    Vector p1 = transform->compute(moving_transformed_, 1 / denom_p,
                                   threshold_truncate_);
    Matrix px(moving_transformed_.rows(), cols);
    for (size_t i = 0; i < cols; ++i) {
        px.col(i) = transform->compute(moving_transformed_,
                                       fixed_.col(i).array() / denom_p,
                                       threshold_truncate_);
    }
    double l =
        -denom_p.log().sum() + cols * fixed_.rows() * std::log(sigma2_) / 2;
    P_ = {p1, pt1, px, l};
}

/// TODO: Are we sure this sigma equals to equation in Fig. 2 of [Myronenko
/// 2010] ??
double CLD::defaultSigma2() {
    return ((moving_.rows() * (fixed_.transpose() * fixed_).trace()) +
            (fixed_.rows() * (moving_.transpose() * moving_).trace()) -
            2 * fixed_.colwise().sum() * moving_.colwise().sum().transpose()) /
           (fixed_.rows() * moving_.rows() * fixed_.cols()) /
           1000.; // 1000 here is a magic number
}

double CLD::computeOptimalRotationCeres(const Matrix& S, const Matrix& T) {
    ceres::Problem problem;
    ceres::LossFunction* loss = NULL;

    /// TODO: stride size? =4?
    ceres::DynamicAutoDiffCostFunction<CostFunctionCLDRot, 4>* cost =
        new ceres::DynamicAutoDiffCostFunction<CostFunctionCLDRot, 4>(
            new CostFunctionCLDRot(G_, S, T, lambda_));
    std::vector<double*> parameter_blocks;
    parameter_blocks.push_back(U_.data());
    cost->AddParameterBlock(U_.size());
    cost->SetNumResiduals(U_.size());
    problem.AddResidualBlock(cost, loss, parameter_blocks);

    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    // options.max_num_iterations = 500;
    // options.linear_solver_type = ceres::DENSE_QR;
    // options.linear_solver_type = ceres::DENSE_SCHUR;
    options.check_gradients = false;
    // options.use_nonmonotonic_steps = true;
    // options.max_num_consecutive_invalid_steps = 100;
    options.minimizer_type = ceres::LINE_SEARCH;
    // options.line_search_direction_type = ceres::BFGS;
    options.logging_type = ceres::SILENT;
    // options.function_tolerance = 1e-10;
    options.function_tolerance = 1e-4;

    ceres::Solver::Summary summary;
    Solve(options, &problem, &summary);
    // std::cout << summary.BriefReport() << "\n";

    return summary.final_cost;
}

void CLD::computeU() {
    Matrix lhs =
        sigma2_ * lambda_ * Matrix::Identity(number_lines_, number_lines_) +
        (FT_ * P_.p1).asDiagonal() * G_;

    A_ = lhs.partialPivLu().solve(FT_ * P_.px);
    B_ = lhs.partialPivLu().solve(
        Matrix(matrixAsSparseBlockDiag(FT_ * P_.p1.asDiagonal() * moving_)));
    GB_ = G_ * B_;

    /// TODO: Why does a sparse D give better results?
    D_ = YD_ - F_ * GB_;
    C_ = FG_ * A_;
    Matrix S = lambda_ * B_.transpose() * GB_ +
               D_.transpose() * P_.p1.asDiagonal() * D_ / sigma2_;
    Matrix T = -lambda_ * A_.transpose() * GB_ +
               (P_.p1.asDiagonal() * C_ - P_.px).transpose() * D_ / sigma2_;
    double cost_aux = computeOptimalRotationCeres(S, T);
}

void CLD::computeSigma2() {
    sigma2_ = 0.0;
    sigma2_ += (fixed_.transpose() * P_.pt1.asDiagonal() * fixed_).trace();
    sigma2_ -= 2 * (P_.px.transpose() * moving_transformed_).trace();
    sigma2_ += (moving_transformed_.transpose() * P_.p1.asDiagonal() *
                moving_transformed_)
                   .trace();
    sigma2_ /= (P_.p1.sum() * 3.0);
}

void CLD::initialize() {
    assert(fixed_.rows() > 0);
    assert(moving_.rows() > 0);
    M_ = moving_.rows();
    G_ = computeG(beta_, number_lines_);
    F_ = computeF(M_, line_sizes_);
    FT_ = F_.transpose();
    FG_ = F_ * G_;
    H_ = computeH(M_, line_sizes_);
    YD_ = matrixAsSparseBlockDiag(moving_) * H_;

    C_ = Matrix::Zero(moving_.rows(), 3);
    D_ = YD_;

    U_ = Matrix::Zero(number_lines_, 3);
    V_ = Matrix::Zero(number_lines_, 3);

    if (sigma2_ == 0.0)
        sigma2_ = defaultSigma2();
    moving_transformed_ = getTransformedMoving();
}

void CLD::computeOne() {
    computeP_FGT();
    computeU();
    moving_transformed_ = getTransformedMoving();
    computeSigma2();
}

Result CLD::run() {
    auto tic = std::chrono::high_resolution_clock::now();

    initialize();

    size_t iter = 0;

    double l = 0.0;
    double ntol = tolerance_ + 10.0;

    while (iter < max_iterations_ && //
           ntol > tolerance_ &&      //
           sigma2_ >
               10 *
                   std::numeric_limits<double>::epsilon()) // Maybe also if
                                                           // sigma2_ increases?
    {
        CLD::computeOne();
        ntol = std::abs((P_.l - l) / P_.l);
        l = P_.l;
        ++iter;
    }

    auto toc = std::chrono::high_resolution_clock::now();
    Result result;
    result.runtime =
        std::chrono::duration_cast<std::chrono::microseconds>(toc - tic);
    result.points = moving_transformed_;
    result.sigma2 = sigma2_;
    result.iterations = iter;

    MatrixX3 GV = G_ * A_ + GB_ * RT_;
    MatrixX3 GU = G_ * U_;
    std::vector<RigidTransformRPY> line_transforms(number_lines_);
    for (size_t i = 0; i < number_lines_; i++) {
        line_transforms[number_lines_].xyz = GV.row(i);
        line_transforms[number_lines_].rpy = GU.row(i);
    }
    result.line_transforms = line_transforms;

    return result;
}

} // namespace cld