#include <lnrr_se3/utils/test_utils.h>

int main(int argc, char const* argv[]) {
    std::string filename_model, filename_scan, filename_scan_size,
        filename_output;
    double beta, lambda, threshold_truncate;
    if (argc == 8) {
        filename_model = argv[1];
        filename_scan = argv[2];
        filename_scan_size = argv[3];
        filename_output = argv[4];
        beta = atof(argv[5]);
        lambda = atof(argv[6]);
        threshold_truncate = atof(argv[7]);
    } else {
        std::cout << "Usage: " << argv[0]
                  << " <filename_model> <filename_scan> "
                     "<filename_scan_size> <filename_output> <beta> <lambda> "
                     "<threshold_truncate>"
                  << std::endl;
        return -1;
    }
    std::cout << "Read input parameters:" << std::endl;
    std::cout << "   -- filename_model: " << filename_model << std::endl;
    std::cout << "   -- filename_scan: " << filename_scan << std::endl;
    std::cout << "   -- filename_scan_size: " << filename_scan_size
              << std::endl;
    std::cout << "   -- filename_output: " << filename_output << std::endl;
    std::cout << "   -- beta: " << beta << std::endl;
    std::cout << "   -- lambda: " << lambda << std::endl;
    std::cout << "   -- threshold_truncate: " << threshold_truncate
              << std::endl;

    MatrixX3 model = readFromTxtFile(filename_model);
    MatrixX3 scan = readFromTxtFile(filename_scan);
    VectorInt line_sizes = readLineSizesFromTxtFile(filename_scan_size);

    std::cout << "Number of points in model: " << model.rows() << std::endl;
    std::cout << "Number of points in scan: " << scan.rows() << std::endl;
    std::cout << "Number of lines in scan: " << line_sizes.rows() << std::endl;
    std::cout << "Number of points in scan: " << line_sizes.sum() << std::endl;

    std::cout << "Initializing..." << std::endl;
    ScanToModel scan_to_model(model, scan, beta, lambda, line_sizes);
    scan_to_model.setMaxIterations(30);
    scan_to_model.setThresholdTruncate(threshold_truncate);
    std::cout << "Starting registration (compile LNRR in DEBUG mode to get "
                 "information at each iteration)"
              << std::endl;
    Result result = scan_to_model.run();
    std::cout
        << "Registration finished, writing registered scan in output txt file"
        << std::endl;
    writeResultToTxtFile(filename_output, result.points);
    std::cout << "Registration finished after " << result.iterations
              << " iterations in " << (double)result.runtime.count() / 1e6
              << " seconds." << std::endl;
    return 0;
}
