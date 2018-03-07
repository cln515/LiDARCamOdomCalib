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
	ransac.threshold_ = threshold;//reprojection error(sqrt(2.0)+0.5)pixel / focal length;
															   //ransac.threshold_ = 7e-7;
	ransac.max_iterations_ = max_itr;
	//cout<<ransac.probability_<<endl;
	ransac.computeModel();

	cout << ransac.model_coefficients_ << endl;
	cout << ransac.inliers_.size() << endl;

	//nonlinear optimization

	double optPara[6];
	R2axisRot(ransac.model_coefficients_.block(0, 0, 3, 3), mot.rx, mot.ry, mot.rz);
	mot.x = ransac.model_coefficients_(0, 3);
	mot.y = ransac.model_coefficients_(1, 3);
	mot.z = ransac.model_coefficients_(2, 3);//redundant
}

void cameraPoseAbsEst_lin(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof& mot, double threshold, int max_itr) {
	opengv::points_t points = opengv::points_t(bvs1.begin(), bvs1.end());
	opengv::bearingVectors_t bearingVectors(bvs2.begin(), bvs2.end());
	
	opengv::absolute_pose::CentralAbsoluteAdapter adapter(bearingVectors, points);
	//ransac
	opengv::sac::Ransac<opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem> ransac;
	std::shared_ptr<opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem> absposeproblem_ptr(
		new opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem(
			adapter,
			opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem::GP3P));
	ransac.sac_model_ = absposeproblem_ptr;
	ransac.threshold_ = threshold;//reprojection error(sqrt(2.0)+0.5)pixel / focal length;
															   //ransac.threshold_ = 7e-7;
	ransac.max_iterations_ = max_itr;
	//cout<<ransac.probability_<<endl;
	ransac.computeModel();
	_6dof optimizedPara;
	Matrix3d cpRot = ransac.model_coefficients_.block(0, 0, 3, 3);
	R2axisRot(cpRot, optimizedPara.rx, optimizedPara.ry, optimizedPara.rz);
	optimizedPara.x = ransac.model_coefficients_(0, 3);
	optimizedPara.y = ransac.model_coefficients_(1, 3);
	optimizedPara.z = ransac.model_coefficients_(2, 3);
	mot = optimizedPara;



}


void cameraPoseAbsTransEst_lin(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof& mot, double threshold, int max_itr) {
	opengv::points_t points = opengv::points_t(bvs1.begin(), bvs1.end());
	opengv::bearingVectors_t bearingVectors(bvs2.begin(), bvs2.end());
	
	opengv::absolute_pose::CentralAbsoluteAdapter adapter(bearingVectors, points);
	Matrix3d priorR = axisRot2R(mot.rx, mot.ry, mot.rz);

	adapter.setR(priorR);
	//ransac
	opengv::sac::Ransac<opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem> ransac;
	std::shared_ptr<opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem> absposeproblem_ptr(
		new opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem(
			adapter,
			opengv::sac_problems::absolute_pose::AbsolutePoseSacProblem::TWOPT));
	ransac.sac_model_ = absposeproblem_ptr;
	ransac.threshold_ = threshold;//reprojection error(sqrt(2.0)+0.5)pixel / focal length;
								  //ransac.threshold_ = 7e-7;
	ransac.max_iterations_ = max_itr;
	//cout<<ransac.probability_<<endl;
	ransac.computeModel();
	//_6dof optimizedPara;
	mot.x = ransac.model_coefficients_(0, 3);
	mot.y = ransac.model_coefficients_(1, 3);
	mot.z = ransac.model_coefficients_(2, 3);
	//mot = optimizedPara;
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

		residual[0] = eyeDirec.cross(plot).norm();//eyeDirec.cross(plot).norm();
		return true;
	}
private:
	Vector3d p3d;
	double x2d, y2d;
	Vector3d eyeDirec;

};

struct projectionCostFunc_Trans {
public:
	projectionCostFunc_Trans(Vector3d& p3d_, Vector3d& eyeDirec_, Matrix3d& R_)
	{
		p3d = p3d_;
		eyeDirec = eyeDirec_;
		R=R_;
	};
	bool operator()(const double* parameters, double* residual) const {
		double x = parameters[0];
		double y = parameters[1];
		double z = parameters[2];

		Vector3d T;T << x, y, z;
		Vector3d plot = (R.transpose()*(p3d - T)).normalized();

		residual[0] = (eyeDirec.cross(plot)).norm();//eyeDirec.cross(plot).norm();
		return true;
	}
private:
	Vector3d p3d;
	Matrix3d R;
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
		//c2=c2.normalized();
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

struct angleFromEpiLine_Rot {
public:
	angleFromEpiLine_Rot(Vector3d& v1_, Vector3d& v2_, Vector3d& c2_) {
		v1 = v1_;
		v2 = v2_;
		c2 = c2_;
	}
	bool operator()(const double* parameters, double* residual)const {
		double rx = parameters[0];
		double ry = parameters[1];
		double rz = parameters[2];

		Matrix3d R, E;
		Vector3d c1, c21;
		c1 << 0, 0, 0;
		R = axisRot2R(rx, ry, rz).inverse();
		c21 = c1 - c2;

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
	Vector3d v1, v2, c2;

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
	optPara[5] = mot.z;//redundant

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

void cameraPoseAbsEstHybrid_non_lin(vector<Vector3d> bvs1, vector<Vector3d> bvs2, vector<Vector3d> bvsi1, vector<Vector3d> bvsi2, _6dof& mot, double loss_scale) {
	double optParar[3],optParat[3];
	optParar[0] = mot.rx;
	optParar[1] = mot.ry;
	optParar[2] = mot.rz;
	optParat[0] = mot.x;
	optParat[1] = mot.y;
	optParat[2] = mot.z;//redundant



	ceres::Problem problem,problem2;
	ceres::CauchyLoss* loss,*loss2;
	if (loss_scale <= 0) {
		loss = NULL;
		loss2=NULL;
	}
	else {
		loss = new ceres::CauchyLoss(loss_scale);
		loss2 = new ceres::CauchyLoss(loss_scale);
	}
	Vector3d t;t << mot.x, mot.y, mot.z;
	for (int i = 0;i<bvsi1.size();i++) {
		angleFromEpiLine_Rot* p = new angleFromEpiLine_Rot(bvsi1.at(i), bvsi2.at(i), t);
		ceres::CostFunction* c = new ceres::NumericDiffCostFunction<angleFromEpiLine_Rot, ceres::CENTRAL, 1, 3>(
			p);
		problem.AddResidualBlock(c, loss, optParar);
	}
	ceres::Solver::Options options;
	options.max_num_iterations = 1e4;
	options.function_tolerance = 1e-6;
	options.parameter_tolerance = 1e-6;
	options.linear_solver_type = ceres::DENSE_QR;
	ceres::Solver::Summary summary,summary2;
	ceres::Solve(options, &problem, &summary);
	Matrix3d R = axisRot2R(optParar[0], optParar[1], optParar[2]);
	for (int i = 0;i<bvs1.size();i++) {
		projectionCostFunc_Trans* p = new projectionCostFunc_Trans(bvs1.at(i), bvs2.at(i), R);
		ceres::CostFunction* c = new ceres::NumericDiffCostFunction<projectionCostFunc_Trans, ceres::CENTRAL, 1, 3>(
			p);
		problem2.AddResidualBlock(c, loss2, optParat);
	}
	ceres::Solve(options, &problem2, &summary2);
	mot.rx = optParar[0];
	mot.ry = optParar[1];
	mot.rz = optParar[2];
	mot.x = optParat[0];
	mot.y = optParat[1];
	mot.z = optParat[2];
}


void cameraPoseAbsTransEst_non_lin(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof& mot, double loss_scale) {
	double optPara[3];
	optPara[0] = mot.x;
	optPara[1] = mot.y;
	optPara[2] = mot.z;
	Matrix3d R = axisRot2R(mot.rx, mot.ry, mot.rz);
	ceres::Problem problem;
	ceres::CauchyLoss* loss;
	if (loss_scale <= 0)loss = NULL;
	else loss = new ceres::CauchyLoss(loss_scale);
	for (int i = 0;i<bvs1.size();i++) {
		projectionCostFunc_Trans* p = new projectionCostFunc_Trans(bvs1.at(i), bvs2.at(i),R);
		ceres::CostFunction* c = new ceres::NumericDiffCostFunction<projectionCostFunc_Trans, ceres::CENTRAL, 1, 3>(
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
	mot.x = optPara[0];
	mot.y = optPara[1];
	mot.z = optPara[2];
}

void cameraPoseRotEst_non_lin(vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof& mot, double loss_scale) {
	double optPara[3];
	optPara[0] = mot.rx;
	optPara[1] = mot.ry;
	optPara[2] = mot.rz;
	Vector3d t; t << mot.x, mot.y, mot.z;
	ceres::Problem problem;
	ceres::CauchyLoss* loss;
	if (loss_scale <= 0)loss = NULL;
	else loss = new ceres::CauchyLoss(loss_scale);
	for (int i = 0;i<bvs1.size();i++) {
		angleFromEpiLine_Rot* p = new angleFromEpiLine_Rot(bvs1.at(i), bvs2.at(i), t);
		ceres::CostFunction* c = new ceres::NumericDiffCostFunction<angleFromEpiLine_Rot, ceres::CENTRAL, 1, 3>(
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
}


struct projectionICPCostFunc {
public:
	projectionICPCostFunc(Vector3d& p3d_, Vector3d& bv1_, Vector3d& bv2_,_6dof mot)
	{
		p3d = p3d_;
		bv1 = bv1_;
		bv2 = bv2_;
		Ra = axisRot2R(mot.rx, mot.ry, mot.rz);
		ta << mot.x, mot.y, mot.z;
	};
	bool operator()(const double* parameters, double* residual) const {
		double rx = parameters[0];
		double ry = parameters[1];
		double rz = parameters[2];
		double x = parameters[3];
		double y = parameters[4];
		double z = parameters[5];

		Matrix3d R = axisRot2R(rx, ry, rz);
		Vector3d t;t << x, y, z;

		Vector3d pt_c1 = R*p3d + t;
		Vector3d pt_c1n = pt_c1.normalized();
		Vector3d pt_c2n = (Ra.transpose()*(pt_c1 - ta)).normalized();

		residual[0] = 1 - fabs(bv1.dot(pt_c1n));
		residual[1] = 1 - fabs(bv2.dot(pt_c2n));
		//eyeDirec.cross(plot).norm();
		return true;
	}
private:
	Vector3d p3d;
	Matrix3d Ra;
	Vector3d ta;
	Vector3d bv1,bv2;

};
//pts: world coordinates
//bvs1: camera 1 coordinates
//bvs2: camera 2 coordinates

void relativePoseEst_non_lin_oneMotion(vector<Vector3d> pts, vector<Vector3d> bvs1, vector<Vector3d> bvs2, _6dof& param, _6dof& mot, double loss_scale) {
	
	double optPara[6];
	optPara[0] = param.rx;
	optPara[1] = param.ry;
	optPara[2] = param.rz;
	optPara[3] = param.x;
	optPara[4] = param.y;
	optPara[5] = param.z;//redundant

	ceres::Problem problem;
	ceres::CauchyLoss* loss;
	if (loss_scale <= 0)loss = NULL;
	else loss = new ceres::CauchyLoss(loss_scale);
	for (int i = 0;i<pts.size();i++) {
		projectionICPCostFunc* p = new projectionICPCostFunc(pts.at(i),bvs1.at(i), bvs2.at(i),mot);
		ceres::CostFunction* c = new ceres::NumericDiffCostFunction<projectionICPCostFunc, ceres::CENTRAL, 2, 6>(
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
	param.rx = optPara[0];
	param.ry = optPara[1];
	param.rz = optPara[2];
	param.x = optPara[3];
	param.y = optPara[4];
	param.z = optPara[5];



}

struct projectionICPCostFunc_withMotion {
public:
	projectionICPCostFunc_withMotion(Vector3d& p3d_, Vector3d& bv1_, Vector3d& bv2_,double* weight)
	{
		p3d = p3d_;
		bv1 = bv1_;
		bv2 = bv2_;
		lambda = weight;
	};
	bool operator()(const double* parameters, const double* parameters2, double* residual) const {
		double rx = parameters[0];
		double ry = parameters[1];
		double rz = parameters[2];
		double x = parameters[3];
		double y = parameters[4];
		double z = parameters[5];

		double rxa = parameters2[0];
		double rya = parameters2[1];
		double rza = parameters2[2];
		double xa = parameters2[3];
		double ya = parameters2[4];
		double za = parameters2[5];


		Matrix3d R = axisRot2R(rx, ry, rz);
		Vector3d t;t << x, y, z;

		Matrix3d Ra = axisRot2R(rxa, rya, rza);
		Vector3d ta;ta << xa, ya, za;

		Vector3d pt_c1 = R*p3d + t;//sensor --> camera 1
		Vector3d pt_c1n = pt_c1.normalized();
		Vector3d pt_c2 = (Ra.transpose()*(pt_c1 - ta));//.normalized();//camera 1 --> camera 2
		Vector3d pt_c2n = pt_c2.normalized();

		residual[0] = lambda[0] * (bv1.cross(pt_c1n)).norm();
		residual[1] = lambda[0] * (bv2.cross(pt_c2n)).norm();
		//eyeDirec.cross(plot).norm();
		return true;
	}
private:
	Vector3d p3d;
	Vector3d bv1, bv2;
	double* lambda;
};

struct projectionICPCostFuncTrans_withMotion {
public:
	projectionICPCostFuncTrans_withMotion(Vector3d& p3d_, Vector3d& bv1_, Vector3d& bv2_,Matrix3d& R_ , double* weight)
	{
		R = R_;
		p3d = p3d_;
		bv1 = bv1_;
		bv2 = bv2_;
		lambda = weight;
	};
	bool operator()(const double* parameters, const double* parameters2, double* residual) const {
		//double rx = parameters[0];
		//double ry = parameters[1];
		//double rz = parameters[2];
		double x = parameters[0];
		double y = parameters[1];
		double z = parameters[2];

		double rxa = parameters2[0];
		double rya = parameters2[1];
		double rza = parameters2[2];
		double xa = parameters2[3];
		double ya = parameters2[4];
		double za = parameters2[5];


		// = axisRot2R(rx, ry, rz);
		Vector3d t;t << x, y, z;

		Matrix3d Ra = axisRot2R(rxa, rya, rza);
		Vector3d ta;ta << xa, ya, za;

		Vector3d pt_c1 = R*p3d + t;//sensor --> camera 1
		Vector3d pt_c1n = pt_c1.normalized();
		Vector3d pt_c2n = (Ra.transpose()*(pt_c1 - ta)).normalized();//camera 1 --> camera 2

		residual[0] = lambda[0] * (1 - fabs(bv1.dot(pt_c1n)));
		residual[1] = lambda[0] * (1 - fabs(bv2.dot(pt_c2n)));
		//eyeDirec.cross(plot).norm();
		return true;
	}
private:
	Matrix3d R;
	Vector3d p3d;
	Vector3d bv1, bv2;
	double* lambda;
};

void mat2axis_angle_loc(Matrix3d m, Vector3d& retv, double& angle) {
	double x, y, z;
	double r = sqrt((m(2, 1) - m(1, 2))*(m(2, 1) - m(1, 2)) + (m(0, 2) - m(2, 0))*(m(0, 2) - m(2, 0)) + (m(1, 0) - m(0, 1))*(m(1, 0) - m(0, 1)));
	x = (m(2, 1) - m(1, 2)) / r;
	y = (m(0, 2) - m(2, 0)) / r;
	z = (m(1, 0) - m(0, 1)) / r;
	Vector3d t;
	t << x, y, z;
	retv = t;
	angle = acos((m(0, 0) + m(1, 1) + m(2, 2) - 1) / 2);
}

struct motionConstraint {
public:
	motionConstraint(_6dof motb,double* weight) {
		
		Rb = axisRot2R(motb.rx, motb.ry, motb.rz);
		tb;tb << motb.x, motb.y, motb.z;
		lambda = weight;

	};
	bool operator()(const double* parameters, const double* parameters2, double* residual) const {
		double rx = parameters[0];
		double ry = parameters[1];
		double rz = parameters[2];
		double x = parameters[3];
		double y = parameters[4];
		double z = parameters[5];

		double rxa = parameters2[0];
		double rya = parameters2[1];
		double rza = parameters2[2];
		double xa = parameters2[3];
		double ya = parameters2[4];
		double za = parameters2[5];

		Matrix3d R = axisRot2R(rx, ry, rz);
		Vector3d t;t << x, y, z;

		Matrix3d Ra = axisRot2R(rxa, rya, rza);
		Vector3d ta;ta << xa, ya, za;

		Vector3d ka, kb;
		double angle1, angle2;
		mat2axis_angle_loc(Ra,ka,angle1);
		mat2axis_angle_loc(Rb, kb, angle2);

		Vector3d diffax = ka - R*kb;
		/*Matrix3d diffR = Ra*R - R*Rb;*/
		Vector3d difft = (Ra*t+ ta)-(R*tb+t);

		residual[0] = lambda[0] * diffax.norm();
		residual[1] = lambda[1] * difft.norm();

		return true;
	}
private:
	Matrix3d Rb;
	Vector3d tb;
	double* lambda;
};

struct motionConstraint2 {
public:
	motionConstraint2(_6dof mota, _6dof motb, double* weight) {

		Ra = axisRot2R(mota.rx, mota.ry, mota.rz);
		ta;ta << mota.x, mota.y, mota.z;

		Rb = axisRot2R(motb.rx, motb.ry, motb.rz);
		tb;tb << motb.x, motb.y, motb.z;
		lambda = weight;

	};
	bool operator()(const double* parameters, double* residual) const {
		double rx = parameters[0];
		double ry = parameters[1];
		double rz = parameters[2];
		double x = parameters[3];
		double y = parameters[4];
		double z = parameters[5];

		Matrix3d R = axisRot2R(rx, ry, rz);
		Vector3d t;t << x, y, z;

		
		Vector3d ka, kb;
		double angle1, angle2;
		mat2axis_angle_loc(Ra, ka, angle1);
		mat2axis_angle_loc(Rb, kb, angle2);

		Vector3d diffax = ka - R*kb;
		/*Matrix3d diffR = Ra*R - R*Rb;*/
		Vector3d difft = (Ra*t + ta) - (R*tb + t);

		residual[0] = lambda[0] * diffax.norm();
		residual[1] = lambda[1] * difft.norm();

		return true;
	}
private:
	Matrix3d Ra,Rb;
	Vector3d ta,tb;
	double* lambda;
};

struct motionConstraint_Trans {
public:
	motionConstraint_Trans(_6dof motb,Matrix3d R_ , double* weight) {

		Rb = axisRot2R(motb.rx, motb.ry, motb.rz);
		tb;tb << motb.x, motb.y, motb.z;
		lambda = weight;
		R = R_ ;
	};
	bool operator()(const double* parameters, const double* parameters2, double* residual) const {
		double x = parameters[0];
		double y = parameters[1];
		double z = parameters[2];

		double rxa = parameters2[0];
		double rya = parameters2[1];
		double rza = parameters2[2];
		double xa = parameters2[3];
		double ya = parameters2[4];
		double za = parameters2[5];

		Vector3d t;t << x, y, z;

		Matrix3d Ra = axisRot2R(rxa, rya, rza);
		Vector3d ta;ta << xa, ya, za;

		/*Matrix3d diffR = Ra*R - R*Rb;*/
		Vector3d difft = (Ra*t + ta) - (R*tb + t);

		residual[0] = lambda[0] * difft.norm();

		return true;
	}
private:
	Matrix3d R;
	Matrix3d Rb;
	Vector3d tb;
	double* lambda;
};
void allOptimization(vector<vector<Vector3d>> pts, vector<vector<Vector3d>> bvs1, vector<vector<Vector3d>> bvs2, _6dof& param, vector<_6dof>& mota, vector<Matrix4d> motb, double loss_scale) {

	double optPara[6];
	optPara[0] = param.rx;
	optPara[1] = param.ry;
	optPara[2] = param.rz;
	optPara[3] = param.x;
	optPara[4] = param.y;
	optPara[5] = param.z;//redundant

	double** optPara_A=(double**)malloc(sizeof(double*)*mota.size());
	
	for (int i = 0;i < mota.size();i++) {
		optPara_A[i] = (double*)malloc(sizeof(double) * 6);
		optPara_A[i][0] = mota.at(i).rx;
		optPara_A[i][1] = mota.at(i).ry;
		optPara_A[i][2] = mota.at(i).rz;
		optPara_A[i][3] = mota.at(i).x;
		optPara_A[i][4] = mota.at(i).y;
		optPara_A[i][5] = mota.at(i).z;
	}

	ceres::Problem problem;
	ceres::CauchyLoss* loss;
	if (loss_scale <= 0)loss = NULL;
	else loss = new ceres::CauchyLoss(loss_scale);

	double weight[1];weight[0] = 1;
	int cnt = 0;
	for (int transitionIdx = 0;transitionIdx<pts.size();transitionIdx++) {
		for (int i = 0;i<pts.at(transitionIdx).size();i++) {//Icp func
			projectionICPCostFunc_withMotion* p = new projectionICPCostFunc_withMotion(pts.at(transitionIdx).at(i), bvs1.at(transitionIdx).at(i), bvs2.at(transitionIdx).at(i),weight);
			ceres::CostFunction* c = new ceres::NumericDiffCostFunction<projectionICPCostFunc_withMotion, ceres::CENTRAL, 2, 6, 6>(
				p);
			problem.AddResidualBlock(c, loss, optPara,optPara_A[transitionIdx]);
			cnt++;
		}
	}
	double weight_[2];
	weight_[0] = sqrt(cnt*loss_scale)*50/mota.size();
	weight_[1] = weight_[0]*0.5;


	for (int transitionIdx = 0;transitionIdx<mota.size();transitionIdx++) {
		_6dof motb_ = m2_6dof(motb.at(transitionIdx));
		motionConstraint* p = new motionConstraint(motb_, weight_);
		ceres::CostFunction* c = new ceres::NumericDiffCostFunction<motionConstraint, ceres::CENTRAL, 2, 6, 6>(
			p);
		problem.AddResidualBlock(c, NULL, optPara, optPara_A[transitionIdx]);
	}


	ceres::Solver::Options options;
	options.max_num_iterations = 20;
	options.function_tolerance = 1e-4;
	options.parameter_tolerance = 1e-4;
	options.linear_solver_type = ceres::DENSE_QR;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	param.rx = optPara[0];
	param.ry = optPara[1];
	param.rz = optPara[2];
	param.x = optPara[3];
	param.y = optPara[4];
	param.z = optPara[5];

	for (int i = 0;i < mota.size();i++) {
		mota.at(i).rx = optPara_A[i][0];
		mota.at(i).ry = optPara_A[i][1];
		mota.at(i).rz = optPara_A[i][2];
		mota.at(i).x = optPara_A[i][3];
		mota.at(i).y = optPara_A[i][4];
		mota.at(i).z = optPara_A[i][5];
		free(optPara_A[i]);
	}
	free(optPara_A);

}

void allTransOptimization(vector<vector<Vector3d>> pts, vector<vector<Vector3d>> bvs1, vector<vector<Vector3d>> bvs2, _6dof& param, vector<_6dof>& mota, vector<Matrix4d> motb, double loss_scale) {

	double optPara[3];
	optPara[0] = param.x;
	optPara[1] = param.y;
	optPara[2] = param.z;//redundant
	Matrix3d R = axisRot2R(param.rx, param.ry, param.rz);


	double** optPara_A = (double**)malloc(sizeof(double*)*mota.size());

	for (int i = 0;i < mota.size();i++) {
		optPara_A[i] = (double*)malloc(sizeof(double) * 6);
		optPara_A[i][0] = mota.at(i).rx;
		optPara_A[i][1] = mota.at(i).ry;
		optPara_A[i][2] = mota.at(i).rz;
		optPara_A[i][3] = mota.at(i).x;
		optPara_A[i][4] = mota.at(i).y;
		optPara_A[i][5] = mota.at(i).z;
	}

	ceres::Problem problem;
	ceres::CauchyLoss* loss;
	if (loss_scale <= 0)loss = NULL;
	else loss = new ceres::CauchyLoss(loss_scale);

	double weight[1];weight[0] = 1;
	int cnt = 0;
	for (int transitionIdx = 0;transitionIdx<pts.size();transitionIdx++) {
		for (int i = 0;i<pts.at(transitionIdx).size();i++) {//Icp func
			projectionICPCostFuncTrans_withMotion* p = new projectionICPCostFuncTrans_withMotion(pts.at(transitionIdx).at(i), bvs1.at(transitionIdx).at(i), bvs2.at(transitionIdx).at(i), R, weight);
			ceres::CostFunction* c = new ceres::NumericDiffCostFunction<projectionICPCostFuncTrans_withMotion, ceres::CENTRAL, 2, 3, 6>(
				p);
			problem.AddResidualBlock(c, loss, optPara, optPara_A[transitionIdx]);
			cnt++;
		}
	}
	double weight_[1];
	weight_[0] = sqrt(cnt)*0.2 / mota.size();
	//weight_[1] = weight_[0] * 0.01;


	for (int transitionIdx = 0;transitionIdx<mota.size();transitionIdx++) {
		_6dof motb_ = m2_6dof(motb.at(transitionIdx));
		motionConstraint_Trans* p = new motionConstraint_Trans(motb_, R, weight_);
		ceres::CostFunction* c = new ceres::NumericDiffCostFunction<motionConstraint_Trans, ceres::CENTRAL, 1, 3, 6>(
			p);
		//problem.AddResidualBlock(c, NULL, optPara, optPara_A[transitionIdx]);
	}


	ceres::Solver::Options options;
	options.max_num_iterations = 20;
	options.function_tolerance = 1e-6;
	options.parameter_tolerance = 1e-6;
	options.linear_solver_type = ceres::DENSE_QR;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);

	param.x = optPara[0];
	param.y = optPara[1];
	param.z = optPara[2];

	for (int i = 0;i < mota.size();i++) {
		mota.at(i).rx = optPara_A[i][0];
		mota.at(i).ry = optPara_A[i][1];
		mota.at(i).rz = optPara_A[i][2];
		mota.at(i).x = optPara_A[i][3];
		mota.at(i).y = optPara_A[i][4];
		mota.at(i).z = optPara_A[i][5];
		free(optPara_A[i]);
	}
	free(optPara_A);

}