# Linewise Non-Rigid Registration (LNRR)



This repository provides a basic implementation of the non-rigid point cloud registration method described in [1].
For a detailed explanation of the method and its motivation, please refer to the [original paper](https://doi.org/10.1109/LRA.2022.3180038).
In a nutshell, the method is able to find the set of rigid transformations that need to be applied to *each line* in the scan in order to match an existing model:

<img src="docs/images/lnrr.png" width="400">

You can also watch the video that summarizes this work:

[![Watch the video](https://img.youtube.com/vi/4QDZ7z1WER8/sddefault.jpg)](https://youtu.be/4QDZ7z1WER8)


Thank you for citing the original publication [1] if you use our method in academic work:
```
@ARTICLE{9788023,  
  author={Castillón, Miguel and Ridao, Pere and Siegwart, Roland and Cadena, César},
  journal={IEEE Robotics and Automation Letters},
  title={Linewise Non-Rigid Point Cloud Registration},
  year={2022},
  volume={7},
  number={3},
  pages={7044-7051},
  doi={10.1109/LRA.2022.3180038}}
```

If you want to know more about our underwater 3D scanner, check out [the paper](https://doi.org/10.1109/TMECH.2022.3170504) [2] and [this blog post](https://miguelcastillon.github.io/project/underwater-3d-scanner/).

## Publications

[1] M. Castillón, P. Ridao, R. Siegwart and C. Cadena, "Linewise Non-Rigid Point Cloud Registration," in IEEE Robotics and Automation Letters, doi: 10.1109/LRA.2022.3180038. [[pdf](https://doi.org/10.1109/LRA.2022.3180038)]

[2] M. Castillón, J. Forest and P. Ridao, "Underwater 3D Scanner to Counteract Refraction: Calibration and Experimental Results," in IEEE/ASME Transactions on Mechatronics, doi: 10.1109/TMECH.2022.3170504. [[pdf](https://doi.org/10.1109/TMECH.2022.3170504)]


## Installation

### Dependencies

Our methods depends on [CMake](https://cmake.org/), [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), and [Ceres](http://ceres-solver.org/index.html).
Please note that so far, it has only been tested on Ubuntu 20.04 + Eigen 3.3.7 + Ceres 2.0.

Moreover, our method uses Fast Gauss Transforms to compute the correspondence probability between each pair of points.
Therefore, our method depends on
[fgt](https://github.com/miguelcastillon/fgt_threshold), which is a fork of [this repository](https://github.com/gadomski/fgt).

### Compilation
As usual, just download and unzip this repository in your preferred location and `cd` into it.
Then:
```bash
mkdir build
cd build
cmake ..
make
sudo make install
```

If you want Debug messages to be printed for each iteration, compile using
```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
```



## Usage

```cpp
#include <lnrr/scan_to_model.h>

int main(int argc, char** argv) {
    lnrr::Matrix fixed = loadModel();
    lnrr::Matrix moving = loadScan();
    lnrr::Vector line_sizes;  // Vector containing the number of points in each line
    double beta = ...;
    double lambda = ...;

    lnrr::ScanToModel lnrr(fixed, moving, beta, lambda, line_sizes);
    lnrr::Result result = lnrr.run();
    return 0;
}
```

And your `CMakeLists.txt` should include:
```cmake
find_package(OpenMP REQUIRED)
find_package(Ceres REQUIRED)
find_package(Fgt REQUIRED)
find_package(Lnrr REQUIRED)

add_library(my-new-library
    my_program.cpp
    )
target_link_libraries(my-new-library
    PUBLIC
    Lnrr::Library-C++
    ${CERES_LIBRARIES}
    )
```

## Example

To run the code with the example model and scan in the folder `data`, you can just run the test:
```bash
./test_lnrr data/stanford-bunny_dense_occluded.txt data/scan.txt data/scan_linesizes.txt data/scan_registered.txt 15 100 0.005
```
Converting between `.pcd` and `.txt` files is easy using PCL, but the code is not added here to limit the number of dependencies.

## Contributing

Please feel free to create [issues](https://github.com/miguelcastillon/lnrr/issues) and [pull requests](https://github.com/miguelcastillon/lnrr/pulls), they will be much appreciated.

## Documentation

Please be aware that this repository is only a simple implementation of the method and may therefore unfortunately not always show a robust behaviour.
There is no documentation yet but we hope the code is self-explanatory.

## License

This library is GPL2, copyright 2022 Miguel Castillón. See LICENSE.txt for the full license text.

In the creation of this library we have drawn inspiration from [this cpd implementation](https://github.com/gadomski/cpd) by Gadomski.
