#include "MotionCompute.h"

void cameraPoseEst_lin(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof& mot, double threshold, int max_itr){
	opengv::bearingVectors_t bvs1_, bvs2_;
	bvs1_ = opengv::bearingVectors_t(bvs1.begin(), bvs1.end());
	bvs2_ = opengv::bearingVectors_t(bvs2.begin(), bvs2.end());
	opengv::relative_pose::CentralRelativeAdapter adapter(bvs1_, bvs2_);
	opengv::sac::Ransac<opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem> ransac;
	std::shared_ptr<opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem> absposeproblem_ptr(
		new opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem(
			adapter,
			opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem::NISTER));
	ransac.sac_model_ = absposeproblem_ptr;
	ransac.threshold_ = threshold;
	ransac.max_iterations_ = max_itr;
	ransac.computeModel();
	double optPara[6];
	R2axisRot(ransac.model_coefficients_.block(0, 0, 3, 3), mot.rx, mot.ry, mot.rz);
	mot.x = ransac.model_coefficients_(0, 3);
	mot.y = ransac.model_coefficients_(1, 3);
	mot.z = ransac.model_coefficients_(2, 3);
}

void cameraPoseAbsEst_lin(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof& mot, double threshold, int max_itr) {
	opengv::points_t points = opengv::points_t(bvs1.begin(), bvs1.end());
	opengv::bearingVectors_t bearingVectors(bvs2.begin(), bvs2.end());
	opengv::absolute_pose::CentralAbsoluteAdapter adapter(bearingVectors, points);
	opengv::sac::Ransac<opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem> ransac;
	std::shared_ptr<opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem> absposeproblem_ptr(
		new opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem(
			adapter,
			opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem::GP3P));
	ransac.sac_model_ = absposeproblem_ptr;
	ransac.threshold_ = threshold;
	ransac.max_iterations_ = max_itr;
	ransac.computeModel();
	_6dof optimizedPara;
	Matrix3d cpRot = ransac.model_coefficients_.block(0, 0, 3, 3);
	R2axisRot(cpRot, optimizedPara.rx, optimizedPara.ry, optimizedPara.rz);
	optimizedPara.x = ransac.model_coefficients_(0, 3);
	optimizedPara.y = ransac.model_coefficients_(1, 3);
	optimizedPara.z = ransac.model_coefficients_(2, 3);
	mot = optimizedPara;
}

struct projectionCostFunc {
public:
	projectionCostFunc(Vector3d& p3d_, Vector3d& eyeDirec_)
	{
		p3d = p3d_;
		eyeDirec = eyeDirec_;
	};
	bool operator()(const double* parameters, double* residual) const {
		double rx = parameters[0];
		double ry = parameters[1];
		double rz = parameters[2];
		double x = parameters[3];
		double y = parameters[4];
		double z = parameters[5];

		Matrix3d R = axisRot2R(rx, ry, rz);
		Vector3d T;T << x, y, z;
		Vector3d plot = (R.transpose()*(p3d - T)).normalized();

		residual[0] = eyeDirec.cross(plot).norm();
		return true;
	}
private:
	Vector3d p3d;
	double x2d, y2d;
	Vector3d eyeDirec;

};


struct angleFromEpiLine {
public:
	angleFromEpiLine(Vector3d& v1_, Vector3d& v2_) {
		v1 = v1_;
		v2 = v2_;
	}
	bool operator()(const double* parameters, double* residual)const {
		double rx = parameters[0];
		double ry = parameters[1];
		double rz = parameters[2];
		double x = parameters[3];
		double y = parameters[4];
		double z = parameters[5];

		Matrix3d R, E;
		Vector3d c1, c2, c21;
		c1 << 0, 0, 0;
		R = axisRot2R(rx, ry, rz).inverse();
		c2 << x, y, z;
		c21 = c1 - c2;
		c21 = R*c21;
		Matrix3d tx;
		tx << 0, -c21(2), c21(1),
			c21(2), 0, -c21(0),
			-c21(1), c21(0), 0;
		E = tx*R;
		Vector3d Ev = E*v1;
		double l = Ev.norm();
		residual[0] = (v2.transpose()*Ev);
		residual[0] = residual[0] / l;

		return true;
	}

private:
	Vector3d v1, v2;

};


void cameraPoseRelEst_non_lin(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof& mot, double loss_scale) {
	double optPara[6];
	optPara[0] = mot.rx;
	optPara[1] = mot.ry;
	optPara[2] = mot.rz;
	optPara[3] = mot.x;
	optPara[4] = mot.y;
	optPara[5] = mot.z;//redundant

	ceres::Problem problem;
	ceres::CauchyLoss* loss;
	if (loss_scale >= 0)loss = NULL;
	else loss = new ceres::CauchyLoss(loss_scale);
	for (int i = 0;i<bvs1.size();i++) {
		angleFromEpiLine* p = new angleFromEpiLine(bvs1.at(i), bvs2.at(i));
		ceres::CostFunction* c = new ceres::NumericDiffCostFunction<angleFromEpiLine, ceres::CENTRAL, 1, 6>(
			p);
		problem.AddResidualBlock(c, loss, optPara);
	}
	ceres::Solver::Options options;
	options.max_num_iterations = 1e5;
	options.function_tolerance = 1e-9;
	options.parameter_tolerance = 1e-9;
	options.linear_solver_type = ceres::DENSE_QR;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	mot.rx = optPara[0];
	mot.ry = optPara[1];
	mot.rz = optPara[2];
	mot.x = optPara[3];
	mot.y = optPara[4];
	mot.z = optPara[5];
}



void cameraPoseAbsEst_non_lin(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof& mot, double loss_scale) {
	double optPara[6];
	optPara[0] = mot.rx;
	optPara[1] = mot.ry;
	optPara[2] = mot.rz;
	optPara[3] = mot.x;
	optPara[4] = mot.y;
	optPara[5] = mot.z;

	ceres::Problem problem;
	ceres::CauchyLoss* loss;
	if (loss_scale <= 0)loss=NULL;
	else loss=new ceres::CauchyLoss(loss_scale);
	for (int i = 0;i<bvs1.size();i++) {
		projectionCostFunc* p = new projectionCostFunc(bvs1.at(i), bvs2.at(i));
		ceres::CostFunction* c = new ceres::NumericDiffCostFunction<projectionCostFunc, ceres::CENTRAL, 1, 6>(
			p);
		problem.AddResidualBlock(c, loss, optPara);
	}
	ceres::Solver::Options options;
	options.max_num_iterations = 1e4;
	options.function_tolerance = 1e-6;
	options.parameter_tolerance = 1e-6;
	options.linear_solver_type = ceres::DENSE_QR;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	mot.rx = optPara[0];
	mot.ry = optPara[1];
	mot.rz = optPara[2];
	mot.x = optPara[3];
	mot.y = optPara[4];
	mot.z = optPara[5];
}














