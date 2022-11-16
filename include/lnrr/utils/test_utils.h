#include <fstream>
#include <lnrr/scan_to_model.h>

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