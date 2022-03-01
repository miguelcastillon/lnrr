# Linewise Non-Rigid Registration (LNRR)



This repository provides a basic implementation of the non-rigid point cloud registration method described in [1].
For a detailed explanation of the method and its motivation, please refer to the [original paper](https://doi.org/10.1109/ACCESS.2021.3069189).
You can also watch the [video](https://youtu.be/ZsPw2voKi10) that summarizes this work.
In a nutshell, the method is able to find the set of rigid transformations that need to be applied to *each line* in the scan in order to match an existing model:

<img src="docs/images/lnrr.png" width="400">


Thank you for citing the original publication [1] if you use our method in academic work:
```
@ARTICLE{castillon2021,  
author={Castillón, Miguel and Palomer, Albert and Forest, Josep and Ridao, Pere},
journal={IEEE Access},
title={Underwater 3D Scanner Model Using a Biaxial MEMS Mirror},
year={2021},
volume={9},
number={},
pages={50231-50243},
doi={10.1109/ACCESS.2021.3069189}}
```

## Publications

[1] M. Castillón, A. Palomer, J. Forest and P. Ridao, "Underwater 3D Scanner Model Using a Biaxial MEMS Mirror," in IEEE Access, vol. 9, pp. 50231-50243, 2021, doi: 10.1109/ACCESS.2021.3069189 

## Installation

### Dependencies

Our methods depends on [CMake](https://cmake.org/), [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), and [Ceres](http://ceres-solver.org/index.html).
Please note that so far, it has only been tested on Ubuntu 20.04 + Eigen 3.3.7 + Ceres 2.0.

Moreover, our method uses Fast Gauss Transforms to compute the correspondence probability between each pair of points.
Therefore, our method depends on
[fgt](https://github.com/miguelcastillon/fgt), which is a fork of [this repository](https://github.com/gadomski/fgt).

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
    double lambda = ...

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

## Contributing

Please feel free to create [issues](https://github.com/miguelcastillon/lnrr/issues) and [pull requests](https://github.com/miguelcastillon/lnrr/pulls), they will be much appreciated.

## Documentation

Please be aware that this repository is only a simple implementation of the method and may therefore unfortunately not always show a robust behaviour.
There is no documentation yet but we hope the code is self-explanatory.

## License

This library is GPL2, copyright 2022 Miguel Castillón. See LICENSE.txt for the full license text.

In the creation of this library we have drawn inspiration from [this cpd implementation](https://github.com/gadomski/cpd) by Gadomski.