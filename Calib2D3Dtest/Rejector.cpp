#include "Rejector.h"

void Reject2D2Dcorrespondence(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof mot, double threshold,vector<int>& inlier) {
	Matrix3d R, E;
	Vector3d c1, c2, c21;
	c1 << 0, 0, 0;
	c2 << mot.x, mot.y, mot.z;
	R = axisRot2R(mot.rx, mot.ry, mot.rz).inverse();
	c21 = R*(c1 - c2);
	Matrix3d tx;
	tx << 0, -c21(2), c21(1),
		c21(2), 0, -c21(0),
		-c21(1), c21(0), 0;
	E = tx*R;
	for (int i = 0;i < bvs1.size();i++) {	
		Vector3d v1, v2;
		v1 = bvs1.at(i).normalized();
		v2 = bvs2.at(i).normalized();
		Vector3d Ev = E*v1;
		double l = Ev.norm();
		double err = (v2.transpose()*Ev);
		err = err / l;//sin(theta)
		if (fabs(err) < threshold)inlier.push_back(i);
	}
	return;
}

void Reject2D3Dcorrespondence(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof mot, double threshold, vector<int>& inlier) {
	Matrix3d Ra = axisRot2R(mot.rx, mot.ry, mot.rz);
	Vector3d ta; ta << mot.x, mot.y, mot.z;
	inlier.clear();
	for (int i = 0;i < bvs1.size();i++) {
		Vector3d pointc1=bvs1.at(i);
		Vector3d pointc2 = Ra.inverse()*(pointc1-ta);
		Vector3d vis = bvs2.at(i);
		
		double err = pointc2.cross(vis).norm() / pointc2.norm();
		if (err < threshold)inlier.push_back(i);
	}	
	return;
}

void Reject2D3Dcorrespondence(vector<Vector3d> bvs1, vector<Vector3d> bvs2, double threshold, vector<int>& inlier) {
	
	inlier.clear();
	for (int i = 0;i < bvs1.size();i++) {
		Vector3d pointc1 = bvs1.at(i);
		Vector3d vis = bvs2.at(i);

		double err = (pointc1.cross(vis).norm()) / pointc1.norm();
		if (err < threshold)inlier.push_back(i);
	}
	return;
}