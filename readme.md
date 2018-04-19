
# LiDAR-Camera Calibration Program

## Description
This is a LiDAR and camera calibration program using motion.

Paper(https://arxiv.org/abs/1804.05178)

## Requirement
Visual studio 2017
Eigen
OpenCV 3.x (https://opencv.org/) (We use 3.30.ver)
OpenGV (https://laurentkneip.github.io/opengv/)
Ceres solver (http://ceres-solver.org/)
## Install Dependencies

### Eigen

    $ wget http://bitbucket.org/eigen/eigen/get/3.3.4.zip
    $ unzip 3.3.4.zip

### OpenCV 3.x

    open bash in c:\\repo
    $ git clone https://github.com/opencv/opencv.git
    $ mkdir opencv_build

#### cmake OpenCV

    Open cmake-gui
    Set source directory "c:\\path\to\opencv", build directory "c:\\path\to\opencv_build"
    Push "Configure" and select your Visual Studio Version (Win64)

### Ceres(minimum config)

    $ git clone https://ceres-solver.googlesource.com/ceres-solver
    $ mkdir ceres_build
    Open cmake-gui
    Set source directory "c:\\path\to\ceres-solver", build directory "c:\\path\to\ceres_build"
    Push "Configure" and select your Visual Studio Version (Win64)
    If you meet error "EIGEN_INCLUDE_DIR not found", please set EIGEN_INCLUDE_DIR "c:\\repo\\eigen-eigen-xxxxxxxxxx" (Please input unziped file name)
    Check "Miniglog"
    Push "Generate"

#### Build on Visual Studio

    Open Visual Studio as administrator
    Open Ceres.sln in "ceres_build" folder
    Go to project property > C/C++ > Optimization > Optimization and select "/Od"
    right click "INSTALL" and build`
    (alternative on cmd as administrator)
    $ "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvars64.bat" && cd C:\path\to\ceres_build && msbuild /p:Configuration=Release INSTALL.vcxproj

### install opengv

    open bash in c:\\repo
    $ git clone https://github.com/laurentkneip/opengv
    $ mkdir opengv_build  

#### cmake opengv

    Open cmake-gui
    Set source directory "c:\\repo\opengv", build directory "c:\\repo\opengv_build"
    Push "Configure" and select your Visual Studio Version (Win64)
    If you meet error "EIGEN_INCLUDE_DIR not found", please set EIGEN_INCLUDE_DIR "c:\\repo\\eigen-eigen-xxxxxxxxxx" (Please input unziped file name)
    Push "Generate"

#### Build on Visual Studio

    Open Visual Studio as administrator
    Open c:\\repo\opengv\opengv.sln
    Remove -Wextra -Wwrite-strings -Wno-unused-parameter option from project "opengv" and "random_generator" at preference page
    if you meet compiler error at about line 773 in main.cpp, change as following
    Eigen::Matrix<double, 1, 1> tt; tt << 2.0;
    Eigen::Matrix<double,1,1> valueM = s.transpose() * M * s + tt * C * s;
    also add "#define _USE_MATH_DEFINES" before #include <math.h> in random_generator.cpp
    right click "INSTALL" and build
    #often compile stops and visual studio becomes busy in the case of using VS 2015 or 2017 but sometime it succeed. You may compile easily on cmd as following.
    (alternative on cmd as administrator)$ "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvars64.bat" && cd C:\repo\opengv_build && msbuild /p:Configuration=Release INSTALL.vcxproj

## Install

Please check the additional include directories, additional library directories, additional dependencies. (Configure: Release, Platform x64)

## Dataset

Indoor Velodyne HDL-64E - Ladybug 3 dataset is here
http://www.cvl.iis.u-tokyo.ac.jp/~ishikawa/data/dataset.zip
This dataset consists of 3D ply files, 2D image files, motion text, and ground truth data.  
The same number jpg file and ply file are scanned in the same place.

### 3D scans

Format : ply
these scans are already aligned and translation matrix is included in header

### motion files

position and orientation changing between two positions in one row is computed and used for calibration. First 26 motions are vertical rotation motion and remains are horizontal rotation.

## Usage
    $ path\to\exe path\to\dataset <number of motion> <separation index>

example

    $ path\to\exe path\to\dataset 8 26
    $ path\to\exe path\to\dataset 0

The number of motion is the number of vertical motion and horizontal motion used for calibration. If this number is 8, 16 motions of 8 vertical motions and 8 horizontal motions are used. If number of motion is zero, all motions described in motion.txt are used.
<separation index> is the separated motions index between vertical and horizontal motion.

## License

## Author
Ryoichi Ishikawa, Computer Vision Laboratory, The University of Tokyo
http://www.cvl.iis.u-tokyo.ac.jp/~ishikawa/
http://www.cvl.iis.u-tokyo.ac.jp/
