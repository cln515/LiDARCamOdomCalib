#pragma once

#include <vector>

#include "Eigen\Core"
#include "Eigen\Dense"

#include "utility.h"

void Reject2D2Dcorrespondence(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof mot, double threshold, vector<int>& inlier);
void Reject2D3Dcorrespondence(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof mot, double threshold, vector<int>& inlier);
void Reject2D3Dcorrespondence(vector<Vector3d> bvs1, vector<Vector3d> bvs2, double threshold, vector<int>& inlier);
