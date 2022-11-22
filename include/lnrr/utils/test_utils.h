#include <fstream>
#include <lnrr/scan_to_model.h>
#include <random>

using namespace lnrr;

MatrixX3 readFromTxtFile(const std::string& file_name) {
    std::ifstream file(file_name);
    std::string line;
    std::vector<double> values;
    while (std::getline(file, line, ',')) {
        std::stringstream line_stream(line);
        double value;
        while (line_stream >> value)
            values.push_back(value);
    }
    size_t rows = values.size() / 3;

    MatrixX3 result(rows, 3);
    for (size_t i = 0; i < rows; i++) {
        result(i, 0) = values[i * 3];
        result(i, 1) = values[i * 3 + 1];
        result(i, 2) = values[i * 3 + 2];
    }
    file.close();
    return result;
}

VectorInt readLineSizesFromTxtFile(const std::string& file_name) {
    std::ifstream file(file_name);
    std::string line;
    std::vector<double> values;
    while (std::getline(file, line)) {
        std::stringstream line_stream(line);
        double value;
        while (line_stream >> value)
            values.push_back(value);
    }
    file.close();
    VectorInt vector(values.size());
    for (size_t i = 0; i < values.size(); i++)
        vector[i] = values[i];
    return vector;
}

void writeResultToTxtFile(const std::string& file_name,
                          const MatrixX3& points) {
    std::ofstream file(file_name);
    for (size_t i = 0; i < points.rows(); i++)
        file << points(i, 0) << " " << points(i, 1) << " " << points(i, 2)
             << std::endl;
    file.close();
}

std::vector<PointCloudPtr> loadPCDFilesFromPath(const std::string& path) {
    std::vector<PointCloudPtr> clouds;
    std::vector<std::string> files;
    boost::filesystem::path dir(path);
    boost::filesystem::directory_iterator end_iter;
    if (boost::filesystem::exists(dir) && boost::filesystem::is_directory(dir))
        for (boost::filesystem::directory_iterator dir_itr(dir);
             dir_itr != end_iter; ++dir_itr)
            if (boost::filesystem::is_regular_file(dir_itr->status()))
                files.push_back(dir_itr->path().string());
    std::sort(files.begin(), files.end());
    for (auto& file : files) {
        PointCloudPtr cloud(new PointCloud);
        pcl::io::loadPCDFile(file, *cloud);
        clouds.push_back(cloud);
    }
    return clouds;
}

// wrote result to pcd file
void eigenCloudToPCL(const MatrixX3& result, PointCloudPtr cloud) {
    if (cloud == nullptr)
        cloud.reset(new PointCloud);
    cloud->width = result.rows();
    cloud->height = 1;
    cloud->points.resize(cloud->width * cloud->height);
    for (size_t i = 0; i < result.rows(); i++) {
        cloud->points[i].x = result(i, 0);
        cloud->points[i].y = result(i, 1);
        cloud->points[i].z = result(i, 2);
    }
}

void writeResultToPCDFile(const std::string& file_name,
                          const MatrixX3& result) {
    PointCloudPtr cloud(new PointCloud);
    eigenCloudToPCL(result, cloud);
    pcl::io::savePCDFile(file_name, *cloud);
}

//////////////////////////////

double generateRandomNumber(const double& min, const double& max) {
    // std::random_device r;
    // std::default_random_engine re{r()};
    std::default_random_engine re{static_cast<long unsigned int>(time(0))};
    std::uniform_real_distribution<double> unif(min, max);
    return unif(re);
}

RigidTransform<double>
generateRandomTransformation(const double& max_translation,
                             const double& max_rotation_deg) {
    RigidTransform<double> t;
    t.xyz = max_translation * Vector3::Random();
    double angle = M_PI / 180.0 *
                   generateRandomNumber(-max_rotation_deg, max_rotation_deg);
    Vector3 axis = Vector3(generateRandomNumber(-1.0, 1.0),
                           generateRandomNumber(-1.0, 1.0),
                           generateRandomNumber(-1.0, 1.0))
                       .normalized();
    assert(abs(axis.norm() - 1.0) < 1e-8);
    // just in case random vector generation
    // does something strange
    if (abs(angle) > 0.0)
        t.rot_mat = Matrix3(Eigen::AngleAxisd(angle, axis));
    else
        t.rot_mat.setIdentity();

    return t;
}

std::vector<RigidTransform<double>>
generateRandomTransformations(const int& number_transformations,
                              const double& max_translation,
                              const double& max_rotation_deg) {
    std::vector<RigidTransform<double>> t(number_transformations);
    for (size_t i = 0; i < number_transformations; i++)
        t[i] = generateRandomTransformation(max_translation, max_rotation_deg);
    return t;
}

RigidTransform<double> invertTransformation(const RigidTransform<double>& T) {
    RigidTransform<double> T_inv;
    T_inv.rot_mat = T.rot_mat.transpose();
    T_inv.xyz = -T_inv.rot_mat * T.xyz;
    return T_inv;
}

std::vector<RigidTransform<double>>
invertTransformations(const std::vector<RigidTransform<double>>& T) {
    std::vector<RigidTransform<double>> T_inv(T.size());
    for (size_t i = 0; i < T.size(); i++)
        T_inv[i] = invertTransformation(T[i]);
    return T_inv;
}