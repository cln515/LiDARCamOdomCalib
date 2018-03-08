#pragma once

#include <opencv2\opencv.hpp>
#include <opengv\relative_pose\CentralRelativeAdapter.hpp>
#include <opengv\relative_pose\methods.hpp>
#include <opengv\sac\Ransac.hpp>
#include <opengv\sac_problems\relative_pose\CentralRelativePoseSacProblem.hpp>

#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac_problems/absolute_pose/AbsolutePoseSacProblem.hpp>

#include <vector>

#include <ceres\ceres.h>

#include "Eigen\Core"
#include "Eigen\Dense"

#include "utility.h"

using namespace Eigen;
using namespace std;

#ifndef MOTION_COMPUTE
#define MOTION_COMPUTE

void cameraPoseEst_lin(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof& mot,double threshold= (1.0 - cos(atan(sqrt(2.0)*0.5 / 400.0))),int max_itr=1000);
void cameraPoseAbsEst_lin(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof& mot, double threshold = (1.0 - cos(atan(sqrt(2.0)*0.5 / 400.0))), int max_itr = 1000);
void cameraPoseRelEst_non_lin(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof& mot, double loss_scale = 1.0e-7);
void cameraPoseAbsEst_non_lin(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof& mot, double loss_scale = 1.0e-7);
#endif