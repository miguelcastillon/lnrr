#include <lnrr/utils/test_utils.h>
#include <lnrr/utils/types.h>
#include <pcl/common/transforms.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/visualization/pcl_visualizer.h>

const double VOXEL_SIZE = 0.01;
const int NUMBER_ITERATIONS = 15;

using namespace lnrr;

int main(int argc, char const* argv[]) {
    std::string filename_model, folder_scans, filename_output;
    double beta, lambda, threshold_truncate;
    if (argc == 7) {
        filename_model = argv[1];
        folder_scans = argv[2];
        filename_output = argv[3];
        beta = atof(argv[4]);
        lambda = atof(argv[5]);
        threshold_truncate = atof(argv[6]);
    } else {
        std::cout << "Usage: " << argv[0]
                  << " <filename_model>"
                     " <folder_scans>"
                     " <filename_output>"
                     " <beta>"
                     " <lambda>"
                     " <threshold_truncate>"
                  << std::endl;
        return -1;
    }
    std::cout << "Read input parameters:" << std::endl;
    std::cout << "   -- filename_model: " << filename_model << std::endl;
    std::cout << "   -- folder_scans: " << folder_scans << std::endl;
    std::cout << "   -- filename_output: " << filename_output << std::endl;
    std::cout << "   -- beta: " << beta << std::endl;
    std::cout << "   -- lambda: " << lambda << std::endl;
    std::cout << "   -- threshold_truncate: " << threshold_truncate
              << std::endl;

    PointCloudPtr model(new PointCloud);
    pcl::io::loadPCDFile(filename_model, *model);

    std::vector<PointCloudPtr> scan = loadPCDFilesFromPath(folder_scans);

    // subsample model
    PointCloudPtr model_subsampled(new PointCloud);
    pcl::VoxelGrid<pcl::PointXYZ> voxel_grid;
    voxel_grid.setInputCloud(model);
    voxel_grid.setLeafSize(VOXEL_SIZE / 2.0, VOXEL_SIZE / 2.0,
                           VOXEL_SIZE / 2.0);
    voxel_grid.filter(*model_subsampled);

    // subsample scans (take only half of the lines)
    std::vector<PointCloudPtr> scan_subsampled;
    for (size_t i = 0; i < scan.size(); i += 3) {
        PointCloudPtr scan_subsampled_i(new PointCloud);
        voxel_grid.setInputCloud(scan[i]);
        voxel_grid.setLeafSize(VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE);
        voxel_grid.filter(*scan_subsampled_i);
        scan_subsampled.push_back(scan_subsampled_i);
    }

    int number_lines = scan_subsampled.size();
    int number_points = 0;
    for (size_t i = 0; i < number_lines; i++)
        number_points += scan_subsampled[i]->size();

    std::cout << "Number of points in model: " << model_subsampled->size()
              << std::endl;
    std::cout << "Number of lines in scan: " << number_lines << std::endl;
    std::cout << "Number of points in scan: " << number_points << std::endl;

    std::cout << "Initializing..." << std::endl;
    ScanToModel scan_to_model(model_subsampled, scan_subsampled, beta, lambda);
    scan_to_model.setMaxIterations(NUMBER_ITERATIONS);
    scan_to_model.setThresholdTruncate(threshold_truncate);
    std::cout << "Starting registration (compile LNRR in DEBUG mode to get "
                 "information at each iteration)"
              << std::endl;
    scan_to_model.initialize();
    // Result result = scan_to_model.run();

    // compute one step of the algorithm and visualize the resulting point
    // clouds
    pcl::transformPointCloud(*model_subsampled, *model_subsampled,
                             Eigen::Vector3d(0.0, 0.0, 0.0),
                             conversions::rpy2Quaternion(M_PI_2, 0.0, 0.0));
    pcl::visualization::PCLVisualizer::Ptr viewer(
        new pcl::visualization::PCLVisualizer("LNRR"));
    // viewer->setCameraPosition(0, 0, 0.5, 0, 0, -1, 0);
    viewer->setBackgroundColor(0.9, 0.9, 0.9);
    viewer->addPointCloud(model_subsampled, "model");
    viewer->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 1.0, "model");
    PointCloudPtr scan_aggregated(new PointCloud);
    for (size_t i = 0; i < number_lines; i++)
        *scan_aggregated += *scan_subsampled[i];
    viewer->addPointCloud(scan_aggregated, "scan");
    viewer->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "scan");

    viewer->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "model");
    viewer->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "scan");

    // viewer->spin();

    for (int i = 0; i < NUMBER_ITERATIONS; i++) {
        scan_to_model.computeOne();
        PointCloudPtr scan_transformed(new PointCloud);
        eigenCloudToPCL(scan_to_model.getTransformedMoving(), scan_transformed);
        // visualize the resulting point clouds
        viewer->updatePointCloud(scan_transformed, "scan");
        // viewer->addCoordinateSystem(1.0);
        viewer->spinOnce();
        std::cout << "Finished teration " << i + 1 << " / " << NUMBER_ITERATIONS
                  << std::endl;
    }

    std::cout
        << "Registration finished, writing registered scan in output file."
        << std::endl;

    writeResultToPCDFile(filename_output, scan_to_model.getTransformedMoving());
    // std::cout << "Registration finished after " << result.iterations
    //           << " iterations in " << (double)result.runtime.count() / 1e6
    //           << " seconds." << std::endl;
    viewer->spin();
    std::cout << "Exiting." << std::endl;
    return 0;
}
