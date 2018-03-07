// Calib2D3Dtest.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include "stdafx.h"
#include <iostream>
#include <fstream>

#include <opencv2\opencv.hpp>
#include <opengv\relative_pose\CentralRelativeAdapter.hpp>
#include <opengv\relative_pose\methods.hpp>
#include <opengv\sac\Ransac.hpp>
#include <opengv\sac_problems\relative_pose\CentralRelativePoseSacProblem.hpp>

#include <vector>

#include <algorithm>

#include <ceres\ceres.h>

#include "Eigen\Core"
#include "Eigen\Dense"

#include "utility.h"
#include "BasicPly.h"
#include "image_utility.h"
#include "MotionCompute.h"
#include "Rejector.h"
#include "FileFinder.h"

using namespace Eigen;
using namespace std;



void solvelin(vector<Matrix4d> pepdMat, vector<Matrix4d> holdMat,Matrix3d R, Vector3d& t) {
	//solve least square problem
	MatrixXd A(pepdMat.size() * 3, 3);
	VectorXd B(pepdMat.size() * 3);
	for (int i = 0;i<pepdMat.size();i++) {
		Vector3d tpep, thol;
		tpep << pepdMat.at(i)(0, 3), pepdMat.at(i)(1, 3), pepdMat.at(i)(2, 3);
		thol << holdMat.at(i)(0, 3), holdMat.at(i)(1, 3), holdMat.at(i)(2, 3);
		Vector3d rightt = tpep - R*thol;
		Matrix3d leftm = Matrix3d::Identity() - pepdMat.at(i).block(0, 0, 3, 3);
		A.block(i * 3, 0, 3, 3) = leftm;
		B(i * 3) = rightt(0);
		B(i * 3 + 1) = rightt(1);
		B(i * 3 + 2) = rightt(2);
	}
	t = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(B);
	//std::cout << "Rotation (Initial)" << std::endl << R << std::endl;
	//std::cout << "Translation (Initial)" << std::endl << T << std::endl;
}

void solvet_non_lin(vector<Matrix4d> pepdMat, vector<Matrix4d> holdMat, Matrix3d R, Vector3d& t) {
	struct motionAlignCostFunc {
	public:
		motionAlignCostFunc(Matrix4d& ma_, Matrix4d& mb_, Matrix3d& R_)
		{
			ma = ma_;
			mb = mb_;
			R = R_;

		};
		bool operator()(const double* parameters, double* residual) const {
			double x = parameters[0];
			double y = parameters[1];
			double z = parameters[2];

			Vector3d t;t << x, y, z;

			Vector3d err = (ma.block(0,0,3,3)*t+ma.block(0,3,3,1)) - (R*mb.block(0, 3, 3, 1) + t);

			residual[0] = err.norm();//eyeDirec.cross(plot).norm();
			return true;
		}
	private:
		Matrix4d ma,mb;
		Matrix3d R;

	};
	double opt[3];
	ceres::Problem problem;
	opt[0] = t(0);
	opt[1] = t(1);
	opt[2] = t(2);
	for (int i = 0;i < pepdMat.size();i++) {
		;
		motionAlignCostFunc* p = new motionAlignCostFunc(pepdMat.at(i), holdMat.at(i),R);
		ceres::CostFunction* c = new ceres::NumericDiffCostFunction<motionAlignCostFunc, ceres::CENTRAL, 1, 3>(
			p);
		problem.AddResidualBlock(c, new ceres::CauchyLoss(5e-3), opt);
	}
	ceres::Solver::Options options;
	options.max_num_iterations = 1e4;
	options.function_tolerance = 1e-6;
	options.parameter_tolerance = 1e-6;
	options.linear_solver_type = ceres::DENSE_QR;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	t << opt[0], opt[1], opt[2];
}

void solvelin2(vector<Matrix4d>& pepdMat, vector<Matrix4d> holdMat, Matrix3d R, Vector3d& t) {
	//solve least square problem
	MatrixXd A(pepdMat.size() * 3, 3+ pepdMat.size());
	A.setZero();
	VectorXd B(pepdMat.size() * 3);
	for (int i = 0;i<pepdMat.size();i++) {
		Vector3d tpep, thol;
		tpep << pepdMat.at(i)(0, 3), pepdMat.at(i)(1, 3), pepdMat.at(i)(2, 3);
		tpep = tpep.normalized();
		thol << holdMat.at(i)(0, 3), holdMat.at(i)(1, 3), holdMat.at(i)(2, 3);
		Vector3d rightt = -R*thol;
		Matrix3d leftm = Matrix3d::Identity() - pepdMat.at(i).block(0, 0, 3, 3);
		A.block(i * 3, 0, 3, 3) = leftm;
		A.block(i * 3, 3+i, 3, 1) = -tpep;
		B(i * 3) = rightt(0);
		B(i * 3 + 1) = rightt(1);
		B(i * 3 + 2) = rightt(2);
	}
	VectorXd ans(3 + pepdMat.size());
	cout << A << endl;

	ans = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(B);
	cout << "x,y,z,t1,t2,t3..." << endl;
	cout << ans.transpose() << endl;
	for (int i = 0;i < pepdMat.size();i++) {
		pepdMat.at(i).block(0, 3, 3, 1) = ans(i+3) *pepdMat.at(i).block(0, 3, 3, 1);
	}
	t << ans(0), ans(1), ans(2);

}

struct energyFunc {
public:
	energyFunc(Vector3d& v1_,Vector3d& v2_) {
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

		Matrix3d R,E;
		Vector3d c1,v2, c2, c21;
		c1 << 0, 0, 0;
		R = axisRot2R(rx, ry, rz).inverse();
		c2 << x, y, z;
		c21 = c1 - c2;

		Matrix3d tx;
		tx << 0, -c21(2), c21(1),
			c21(2), 0, -c21(0),
			-c21(1),c21(0),0;
		E=tx*R;
		Vector3d Ev = E*v1;
		double l = Ev.norm();
		residual[0] = (v2.transpose()*Ev);
		residual[0] = residual[0] / l;

		return true;
	}

private:
	Vector3d v1,v2;


};


void mat2axis_angle(Matrix3d m,  Vector3d& retv, double& angle) {
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

void solverot(vector<Matrix4d> pepdMat, vector<Matrix4d> holdMat, Matrix3d& R) {
	MatrixXd KA(3, pepdMat.size()), KB(3, pepdMat.size());
	std::cout << "axis of rotation matrix" << std::endl;
	struct axisAlignCostFunc {
	public:
		axisAlignCostFunc(Vector3d& p3d_, Vector3d& eyeDirec_)
		{
			p3d = p3d_;
			eyeDirec = eyeDirec_;
		};
		bool operator()(const double* parameters, double* residual) const {
			double rx = parameters[0];
			double ry = parameters[1];
			double rz = parameters[2];

			Matrix3d R = axisRot2R(rx, ry, rz);
			Vector3d plot = (R.transpose()*p3d).normalized();

			residual[0] = eyeDirec.cross(plot).norm();//eyeDirec.cross(plot).norm();
			return true;
		}
	private:
		Vector3d p3d;
		Vector3d eyeDirec;

	};

	struct rotErrCostFunc {
	public:
		rotErrCostFunc(Matrix3d& RA_, Matrix3d& RB_)
		{
			RA = RA_;
			RB = RB_;
		};
		bool operator()(const double* parameters, double* residual) const {
			double rx = parameters[0];
			double ry = parameters[1];
			double rz = parameters[2];

			Matrix3d R = axisRot2R(rx, ry, rz);
			Matrix3d err = RA*R-R*RB;

			residual[0] = err.norm();//eyeDirec.cross(plot).norm();
			return true;
		}
	private:
		Matrix3d RA;
		Matrix3d RB;

	};


	double opt[3];
	ceres::Problem problem;
	for (int i = 0;i<pepdMat.size();i++) {
		double angle1,angle2;
		cout << pepdMat.at(i) << endl;
		cout << holdMat.at(i) << endl;
		Vector3d pepax; mat2axis_angle(pepdMat.at(i).block(0,0,3,3), pepax,angle1);
		Vector3d holax; mat2axis_angle(holdMat.at(i).block(0, 0, 3, 3), holax, angle2);
		KA.col(i) = pepax;
		KB.col(i) = holax;
		std::cout << pepax << std::endl;
		std::cout << angle1 << endl;
		std::cout << holax << std::endl << std::endl;
		std::cout << angle2 << endl<<endl;
		//axisAlignCostFunc* p = new axisAlignCostFunc(pepax, holax);
		Matrix3d RA = pepdMat.at(i).block(0, 0, 3, 3);
		Matrix3d RB = holdMat.at(i).block(0, 0, 3, 3);
		rotErrCostFunc* p= new rotErrCostFunc(RA, RB);
		ceres::CostFunction* c = new ceres::NumericDiffCostFunction<rotErrCostFunc, ceres::CENTRAL, 1, 3>(
			p);
		problem.AddResidualBlock(c,new ceres::HuberLoss(1.0e-2), opt);
	}
	std::cout << "KA,KB" << std::endl;
	std::cout << KA << std::endl;
	std::cout << KB << std::endl;
	MatrixXd KBKA = KB*KA.transpose();
	//calc rotation
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(KBKA, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::Matrix3d Hm;
	Eigen::Matrix3d uvt = svd.matrixU()*svd.matrixV().transpose();
	Hm << 1, 0, 0,
		0, 1, 0,
		0, 0, uvt.determinant();
	R = svd.matrixV()*Hm*svd.matrixU().transpose();
	R2axisRot(R, opt[0], opt[1], opt[2]);
	ceres::Solver::Options options;
	options.max_num_iterations = 1e4;
	options.function_tolerance = 1e-6;
	options.parameter_tolerance = 1e-6;
	options.linear_solver_type = ceres::DENSE_QR;
	ceres::Solver::Summary summary;
	R = axisRot2R(opt[0], opt[1], opt[2]);
	cout << KA - R*KB << endl;
	ceres::Solve(options, &problem, &summary);
	


	R = axisRot2R(opt[0], opt[1], opt[2]);
	cout << KA - R*KB << endl;
	Matrix3d tR = axisRot2R(0.004563, 0.000525, 2.18017);
	cout << KA - tR*KB << endl;
}

void getTransition(Matrix3d R, Vector3d t, Matrix3d Ra, Vector3d ta, Matrix3d& Rb, Vector3d& tb) {
	tb = R.transpose()*(Ra*t+ta-t);
	Vector3d axis;
	double angle;
	mat2axis_angle(Ra,axis,angle);
	axis = R*axis;
	Matrix3d raxis;
	raxis << 0, -axis(2), axis(1),
		axis(2), 0, -axis(0),
		-axis(1), axis(0), 0;

	Rb = Matrix3d::Identity() + sin(angle)*raxis + (1-cos(angle))*raxis*raxis;

}



int main(int argc,char* argv[])
{
	Matrix4d mgt = readCPara("C:\\data\\calib\\outdooreval\\manual.cpara");
	Matrix4d mgtinv = mgt.inverse();
	cout << m2_6dof(mgtinv) << endl;
	//config parameteres
	//strategy
	bool rotOnly = false;
	bool trackerStrategy = false;//initial position of the tracked point
	//outlier rejection
	double thres_init = 5e-2;
	double thresDecreaseRate = 0.5;
	double trackerror = 2.0;
	double cornerThres = 0.0;
	int outloopInlier = 500;
	//initial parameter
	bool giveInitPara = false;

	if (argc < 2) {
	
		return 0;
	}
	
	string targetFolder=argv[1];
	ofstream csvlog(targetFolder + "\\log"+getTimeStamp()+".csv");
	ofstream txtlog(targetFolder + "\\log" + getTimeStamp() + ".txt");
	//
	ifstream motionfile(targetFolder+"\\motion.txt");
	//
	int m_cnt = stoi(argv[2]);
	std::random_device rdev{}; //乱数生成器を生成
	std::mt19937 mt(rdev()); //メルセンヌ・ツイスタ生成器を使用
	vector<int> vec;
	vector<int> vec2;




	//
	string str;
	vector<int> motionList_,motionList;

	while (getline(motionfile, str)) {
		int mbase=stoi(str.substr(0, str.find_first_of("-")))-1;
		str.erase(0, str.find_first_of("-") + 1);
		int mobj= stoi(str)-1;
		if (m_cnt > 0) {
			motionList_.push_back(mbase);
			motionList_.push_back(mobj);
		}
		else {
			motionList.push_back(mbase);
			motionList.push_back(mobj);
		}
	}

	if (m_cnt > 0) {
		int sep = stoi(argv[3]);
		int cnt = 0;
		while (cnt <= sep) {
			vec.push_back(cnt);
			cnt++;
		}
		while (cnt <= motionList_.size()/2-1) {
			vec2.push_back(cnt);
			cnt++;
		}
		std::shuffle(vec.begin(), vec.end(), mt);
		std::shuffle(vec2.begin(), vec2.end(), mt);
	}


	if (m_cnt > 0) {
		for (int i = 0;i < m_cnt;i++) {
			motionList.push_back(motionList_.at(vec.at(i) * 2));
			motionList.push_back(motionList_.at(vec.at(i) * 2 + 1));
		}
		for (int i = 0;i < m_cnt;i++) {
			motionList.push_back(motionList_.at(vec2.at(i) * 2));
			motionList.push_back(motionList_.at(vec2.at(i) * 2 + 1));
		}
	}

	csvlog << ",rx,ry,rz,x,y,z" << endl;
	Matrix4d gt;
	if (true) {
		//gt calculation
		Matrix4d gtal = getMatrixFlomPly(targetFolder + "/groundtruth.ply");

		Matrix4d faro2im;
		faro2im << 1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1;
		gt = faro2im*gtal.inverse();
		_6dof gtpos;
		R2axisRot(gt.block(0, 0, 3, 3),gtpos.rx , gtpos.ry, gtpos.rz);
		gtpos.x = gt(0,3);
		gtpos.y = gt(1, 3);
		gtpos.z = gt(2, 3);

		//Matrix4d m1 = getMatrixFlomPly("C:\\data\\calib\\velo2\\tozf.ply");
		//Matrix4d m2 = readCPara("C:\\data\\calib\\velo2\\tozf.cpara");
		//gt = m2.inverse()*m1;
		//cout << m2_6dof(m2) << endl;

		csvlog << "gt," << m2_6dof(gt) << endl;

	}
	char szFullPath[MAX_PATH] = { '\0' };
	char *szFilePart;

	DWORD dwRet = GetFullPathNameA(
		targetFolder.c_str(),
		sizeof(szFullPath) / sizeof(szFullPath[0]),
		szFullPath,
		&szFilePart);
	targetFolder = szFullPath;
	//data read
	vector<cv::Mat> imgs,descriptors;
	vector<string>plyFileList,imgFileList;
	plyFileList=getFileWithNumber(targetFolder, "position", "ply");
	imgFileList=getFileWithNumber(targetFolder, "position", "jpg");
//	gtaFileList= getFileWithNumber(targetFolder, "gta", "ply");

	vector<bool> included(plyFileList.size(),false);
	for (int i = 0;i < motionList.size();i++) {
		included.at(motionList.at(i)) = true;
	}
	

	for (int i = 0;i < imgFileList.size();i++) {
		cout<< imgFileList.at(i) <<endl;
		cv::Mat img = cv::imread(imgFileList.at(i), cv::IMREAD_COLOR);
		cv::resize(img, img, cv::Size(5000, 2500));
	
		imgs.push_back(img);
	}
	cv::Mat img1g, img2g;

	int motionNumber = motionList.size()/2;
	vector<Matrix4d> lidarPos,gtAPos;
	
	for (int i = 0;i < plyFileList.size();i++) {
		Matrix4d m1 = getMatrixFlomPly(plyFileList.at(i));
		lidarPos.push_back(m1);
	}
	Matrix3d R;
	Vector3d t;
	vector<Matrix4d> Ma_array, Mb_array;
	for (int i = 0;i < motionNumber;i++) {
		int mbase = motionList.at(i * 2);
		int mdst = motionList.at(i * 2 + 1);
		Mb_array.push_back(lidarPos.at(mbase).inverse()*lidarPos.at(mdst));
	}
	if (!giveInitPara) {
		auto algorithm = cv::AKAZE::create();

		vector<vector<cv::KeyPoint>> keypoints;
		for (int i = 0;i < imgs.size();i++) {
			cout << "image " << i << "keypoint detection" << endl;
			vector<cv::KeyPoint> keypoint;
			cv::Mat descriptor;
			if (included.at(i)) {
				algorithm->detect(imgs.at(i), keypoint);
				algorithm->compute(imgs.at(i), keypoint, descriptor);
			}
			keypoints.push_back(keypoint);
			descriptors.push_back(descriptor);
		}
		//motionNumber = imgs.size() - 1;
		cv::Ptr<cv::DescriptorMatcher> matcher = cv::DescriptorMatcher::create("BruteForce");
		std::vector<cv::DMatch> match12, match13, match34, match12_, match13_, match_;

		cout << "camera pose estimation" << endl;
		//camera pose estimation
		vector<_6dof> amots;
		for (int i = 0;i < motionNumber;i++) {
			int mbase = motionList.at(i * 2);
			int mdst = motionList.at(i * 2 + 1);
			matcher->match(descriptors.at(mbase), descriptors.at(mdst), match12);
			vector<Vector3d> bearingVectors1, bearingVectors2;
			opengv::bearingVector_t bv;
			for (int j = 0;j < match12.size();j++) {
				double ix = keypoints.at(mbase).at(match12.at(j).queryIdx).pt.x;
				double iy = keypoints.at(mbase).at(match12.at(j).queryIdx).pt.y;
				rev_omniTrans(ix, iy, imgs.at(mbase).cols, imgs.at(mbase).rows, bv);
				bearingVectors1.push_back(bv);
				ix = keypoints.at(mdst).at(match12.at(j).trainIdx).pt.x;
				iy = keypoints.at(mdst).at(match12.at(j).trainIdx).pt.y;
				rev_omniTrans(ix, iy, imgs.at(mdst).cols, imgs.at(mdst).rows, bv);
				bearingVectors2.push_back(bv);
			}
			_6dof mot;
			cameraPoseEst_lin(bearingVectors1, bearingVectors2, mot, sin(1.0e-2), 10000);
			txtlog << "position " << mbase << " - " << (mdst) << endl;
			double thres = 1.0e-2;
			vector<int>inlier;
			//Reject2D2Dcorrespondence(bearingVectors1, bearingVectors2, mot, sin(3.0e-3), inlier);
			vector<Vector3d> in_bvs1, in_bvs2;
			for (int loop = 0;loop < 60;loop++) {
				inlier.clear();
				Reject2D2Dcorrespondence(bearingVectors1, bearingVectors2, mot, thres, inlier);
				txtlog << "initial inlier" << inlier.size() << endl;
				txtlog << mot << endl;
				if (inlier.size() < 100)break;
				in_bvs1.clear();
				in_bvs2.clear();
				for (int j = 0;j < inlier.size();j++) {
					in_bvs1.push_back(bearingVectors1.at(inlier.at(j)));
					in_bvs2.push_back(bearingVectors2.at(inlier.at(j)));
				}
				cameraPoseRelEst_non_lin(in_bvs1, in_bvs2, mot, thres / 2);
				thres = thres*0.8;
			}

			txtlog << "initial: " << mot << endl;


			txtlog << mot << endl;
			//re estimation
			if (false) {
				cv::Mat resynthimg1, resynthimg2;
				Matrix3d R0, R1;
				Matrix3d Ra = axisRot2R(mot.rx, mot.ry, mot.rz);;
				Vector3d ta; ta << mot.x, mot.y, mot.z;
				panoramaRectification(imgs.at(mbase), imgs.at(mdst), resynthimg1, resynthimg2, ta, Ra, R0, R1);
				Matrix3d R0inv = R0.inverse();
				Matrix3d R1inv = R1.inverse();
				vector<cv::Point2f> src, dst;

				//cv::cvtColor(resynthimg1, img1g, cv::COLOR_BGR2GRAY);
				//cv::cvtColor(resynthimg2, img2g, cv::COLOR_BGR2GRAY);
				//cv::equalizeHist(img1g, img1g);
				//cv::equalizeHist(img2g, img2g);

				//cv::imwrite(targetFolder + "/haaa.jpg",img1g);
				//cv::goodFeaturesToTrack(img1g, src, 0, 1e-5, 2.0);
				//vector<uchar> status;
				//vector<float> err;
				//cv::calcOpticalFlowPyrLK(img1g, img2g, src, dst, status, err, cv::Size(21, 21));

				vector<cv::KeyPoint> key1, key2;
				cv::Mat desc1, desc2;
				algorithm->detect(resynthimg1, key1);
				algorithm->detect(resynthimg2, key2);
				algorithm->compute(resynthimg1, key1, desc1);
				algorithm->compute(resynthimg2, key2, desc2);
				matcher->match(desc1, desc2, match12);

				bearingVectors1.clear();
				bearingVectors2.clear();


				for (int j = 0;j < match12.size();j++) {
					double ix = key1.at(match12.at(j).queryIdx).pt.x;
					double iy = key1.at(match12.at(j).queryIdx).pt.y;
					rev_omniTrans(ix, iy, imgs.at(mbase).cols, imgs.at(mbase).rows, bv);
					bv = R0inv*bv;

					bearingVectors1.push_back(bv);
					ix = key2.at(match12.at(j).trainIdx).pt.x;
					iy = key2.at(match12.at(j).trainIdx).pt.y;
					rev_omniTrans(ix, iy, imgs.at(mdst).cols, imgs.at(mdst).rows, bv);

					bv = R1inv*bv;
					bearingVectors2.push_back(bv);
				}
				//for (int j = 0;j < src.size();j++) {
				//	if (status[j] != '\1' || err[j] > trackerror)continue;
				//	double ix = src.at(j).x;
				//	double iy = src.at(j).y;
				//	rev_omniTrans(ix, iy, imgs.at(mbase).cols, imgs.at(mbase).rows, bv);
				//	bv = R0inv*bv;
				//	bearingVectors1.push_back(bv);
				//	ix = dst.at(j).x;
				//	iy = dst.at(j).y;
				//	rev_omniTrans(ix, iy, imgs.at(mdst).cols, imgs.at(mdst).rows, bv);
				//	bv = R1inv*bv;
				//	bearingVectors2.push_back(bv);
				//}
				//cameraPoseEst_lin(bearingVectors1, bearingVectors2, mot, sin(1.0e-2),10000);
				thres = 5.0e-2;
				for (int loop = 0;loop < 60;loop++) {
					inlier.clear();
					Reject2D2Dcorrespondence(bearingVectors1, bearingVectors2, mot, thres, inlier);
					txtlog << "refinement inlier" << inlier.size() << endl;
					txtlog << mot << endl;
					if (inlier.size() < 500)break;
					in_bvs1.clear();
					in_bvs2.clear();
					for (int j = 0;j < inlier.size();j++) {
						in_bvs1.push_back(bearingVectors1.at(inlier.at(j)));
						in_bvs2.push_back(bearingVectors2.at(inlier.at(j)));
					}
					cameraPoseRelEst_non_lin(in_bvs1, in_bvs2, mot, thres / 2);
					thres = thres*0.8;
				}

				for (int j = 0;j < in_bvs1.size();j++) {
					double phi, theta;
					Vector3d pt = R0 * in_bvs1.at(j);
					Vector3d pt2 = R1 * in_bvs2.at(j);
					omniTrans(pt(0), pt(1), pt(2), phi, theta, resynthimg1.rows);
					cv::Point2f p1(theta, phi);
					omniTrans(pt2(0), pt2(1), pt2(2), phi, theta, resynthimg1.rows);
					cv::Point2f p2(theta, phi);
					cv::circle(resynthimg1, p1, 3, cv::Scalar(255, 0, 0));
					cv::line(resynthimg1, p1, p2, cv::Scalar(0, 255, 0));
					cv::circle(resynthimg2, p2, 3, cv::Scalar(0, 0, 255));
				}
				stringstream ss, ss2;
				ss << "imotion" << i << "_src.jpg";
				cv::imwrite(targetFolder + "/" + ss.str(), resynthimg1);
				ss2 << "imotion" << i << "_dst.jpg";
				cv::imwrite(targetFolder + "/" + ss2.str(), resynthimg2);
			}
			txtlog << mot << endl;

			//gt
			Matrix4d gtA;
			Matrix4d estA = _6dof2m(mot);
			gtA = gt * lidarPos.at(mbase).inverse() * lidarPos.at(mdst) * gt.inverse();
			_6dof gtAmot = m2_6dof(gtA);
			txtlog << gtAmot << endl;

			double anglea, angleb;
			Vector3d axisa,axisb;
			mat2axis_angle(gtA.block(0, 0, 3, 3), axisa, anglea);
			mat2axis_angle(estA.block(0, 0, 3, 3), axisb, angleb);

			txtlog << "axis gt-est: " << axisa << endl << axisb << endl;

			txtlog << "angle gt-est: " << anglea << "," << angleb << endl;
			amots.push_back(mot);
		}
		cout << "initial calibration parameter computation" << endl;


		for (int i = 0;i < amots.size();i++) {
			Matrix4d ma;
			ma.block(0, 0, 3, 3) = axisRot2R(amots.at(i).rx, amots.at(i).ry, amots.at(i).rz);
			ma.block(0, 3, 3, 1) << amots.at(i).x, amots.at(i).y, amots.at(i).z;
			ma.block(3, 0, 1, 4) << 0, 0, 0, 1;
			Ma_array.push_back(ma);
		}


		//solving
		solverot(Ma_array, Mb_array, R);


		solvelin2(Ma_array, Mb_array, R, t);
	}
	cout << R << endl;
	Matrix3d Rinv = R.inverse();
	//R = R.transpose();
	Vector3d tinv =- Rinv*t;
	cout << Rinv << endl;
	Vector4d q = dcm2q(Rinv);
	
	//log calc
	_6dof param_;
	{
		_6dof param;
		R2axisRot(R, param.rx, param.ry, param.rz);
		param.x = t(0);
		param.y = t(1);
		param.z = t(2);
		
		double errt,errr;
		_6dof gt6mot= m2_6dof(gt);
		_6dof errParam = gt6mot - param;
		errt = errParam.x*errParam.x+ errParam.y*errParam.y+ errParam.z*errParam.z;
		errt = sqrt(errt);
		errr = errParam.rx*errParam.rx + errParam.ry*errParam.ry + errParam.rz*errParam.rz;
		errr = sqrt(errr);
		csvlog << "init," << param.rx << "," << param.ry << "," << param.rz << "," << param.x << "," << param.y << "," << param.z <<",,"<<errr<<","<<errt<< endl;

		

		param_ = param;
	}
	ofstream ofs(targetFolder+"/checkinit.cpara");
	ofs << "Camera Parameter, Trans" << endl;
	ofs << tinv(0)<<" "<<tinv(1) << " " << tinv(2) << endl;
	ofs << "Rotation" << endl;
	ofs << q(0)<<" " << q(1) << " " << q(2) << " " << q(3)  << endl;
	
	


	double thresholds[] = {1e-2,1e-3,8e-4,6e-4,4e-4};
	double losses[] = { 3e-3,3e-4,2.4e-4,1.8e-4,1.2e-4 };
	
	int paramIdx = 0;
	int paramCnt = 0;

	vector<BasicPly> bps;
	for (int i = 0;i < lidarPos.size();i++) {
		BasicPly bp;
		vector<string> fns;
		fns.push_back(plyFileList.at(i));
		bp.readPlyFile(fns, 1);
		bps.push_back(bp);
	}

	int prevMaxInlier = -1;
	

	for (int itr = 0;itr < 500;itr++) {
		vector<vector<Vector3d>> vin_bvs1, vin_bvs_n1, vin_bvs2, vin_bvsi1, vin_bvsi2, vin_pt_s1;
		//vector<Vector3d> 
		int maxInlier = -1;
		for (int imgid = 0;imgid < motionNumber;imgid++) {
			int mbase = motionList.at(imgid * 2);
			int mdst = motionList.at(imgid * 2 + 1);
			cv::Mat resynthimg1, resynthimg2;
			Matrix3d R0, R1;
			Matrix3d Ra = Ma_array.at(imgid).block(0, 0, 3, 3);
			Vector3d ta = Ma_array.at(imgid).block(0, 3, 3, 1);
			panoramaRectification(imgs.at(mbase), imgs.at(mdst),resynthimg1,resynthimg2,ta, Ra,R0,R1);
			Matrix3d R1inv = R1.inverse();

			//cv::normalize(resynthimg1, resynthimg1);
			//cv::normalize(resynthimg2, resynthimg2);

			cv::cvtColor(resynthimg1, img1g, cv::COLOR_BGR2GRAY);
			cv::cvtColor(resynthimg2, img2g, cv::COLOR_BGR2GRAY);
			cv::equalizeHist(img1g, img1g);
			cv::equalizeHist(img2g, img2g);

			Matrix4d mrev = getMatrixFlomPly(plyFileList.at(mbase)).inverse();
			Matrix4d mrev2 = getMatrixFlomPly(plyFileList.at(mdst)).inverse();

			BasicPly bp = bps.at(mbase);
			float*vps = bp.getVertecesPointer();


			vector<Vector3d> bearingVectors1, bearingVectors2, bearingVectors3, bearingVectors4;
			vector<Vector3d> bearingVectors1c, bearingVectors2c, bearingVectorsBase, bearingVectorsBaseback;
			vector<cv::Point2f> src, dst, dst2;
			src.clear();
			vector<cv::Point2f> srcc, dstc;
			for (int i = 0;i < bp.getVertexNumber();i++) {
				Vector3d pt, pt2, pt3;pt << vps[i * 3], vps[i * 3 + 1], vps[i * 3 + 2];
				//world->lidar->camera2
				pt3 = R*(mrev2.block(0, 0, 3, 3)*pt + mrev2.block(0, 3, 3, 1)) + t;
				pt = R*(mrev.block(0, 0, 3, 3)*pt + mrev.block(0, 3, 3, 1)) + t; //camera cooridate
				pt2 = R0 * pt; //camera rect coordinates
				if (pt2.norm() < 0.5)continue;
				double ix, iy;
				omniTrans(pt2(0), pt2(1), pt2(2), iy, ix, resynthimg1.rows);
				cv::Point2f p(ix, iy);
				if (iy <= 0 || ix <= 0 || iy > resynthimg1.rows - 1 || ix > resynthimg1.cols - 1)continue;
				bearingVectors1c.push_back(pt);
				//camera1 --> camera2 rect
				pt2 = R1*(Ra.inverse()*(pt - ta));

				bearingVectors2c.push_back(pt3);
				srcc.push_back(p);

				if (trackerStrategy) {
					omniTrans(pt2(0), pt2(1), pt2(2), iy, ix, resynthimg1.rows);
					cv::Point2f p2(ix, iy);
					dstc.push_back(p2);//initial lidar1->camera1->c2
				}
				else {
					pt3 = R1*pt3;
					omniTrans(pt3(0), pt3(1), pt3(2), iy, ix, resynthimg1.rows);
					cv::Point2f p3(ix, iy);
					dstc.push_back(p3);//initial l1->l2->c2 
				}
			}
			vector<int>selectedIdx;
			goodFeatureToTrack_onProjection(img1g, srcc, selectedIdx, 3.0, 0, cornerThres);
			for (int sIdx = 0;sIdx < selectedIdx.size();sIdx++) {
				src.push_back(srcc.at(selectedIdx.at(sIdx)));
				bearingVectors1.push_back(bearingVectors1c.at(selectedIdx.at(sIdx)));
				bearingVectors3.push_back(bearingVectors2c.at(selectedIdx.at(sIdx)));
				dst.push_back(dstc.at(selectedIdx.at(sIdx)));
			}
			//dst = vector<cv::Point2f>(src);
			vector<uchar> status;
			vector<float> err;
			cv::calcOpticalFlowPyrLK(img1g, img2g, src, dst, status, err, cv::Size(21, 21),0);
			for (int idx = 0;idx < dst.size();idx++) {
				Vector3d bv2;
				rev_omniTrans(dst.at(idx).x, dst.at(idx).y,img2g.size().width, img2g.size().height,bv2);
				bv2 = R1inv*bv2;
				bearingVectors2.push_back(bv2);
			}
			//initial
			_6dof mot;
			if (itr == 0) {
				R2axisRot(Ra, mot.rx, mot.ry, mot.rz);
				cameraPoseAbsTransEst_lin(bearingVectors1, bearingVectors2, mot, 1.0e-2, 10000);
				//txtlog << mot << endl;
			}
			else {
				Matrix4d Ma = Ma_array.at(imgid);
				mot=m2_6dof(Ma);
			}
			vector<int>inlier, inlierImg;
			int itcnt = -1;
			double thres = thres_init;
			vector<Vector3d> in_bvs1, in_bvs_n1, in_bvs2, in_bvsi1, in_bvsi2, in_pt_s1;
			for (int loop = 0;loop < 60;loop++) {
				
				Reject2D2Dcorrespondence(bearingVectors1, bearingVectors2, mot, thres, inlierImg);
				Reject2D3Dcorrespondence(bearingVectors1, bearingVectors2, mot, thres, inlier);
				if (inlier.size() < outloopInlier) {
					itcnt = loop;
					break;
				}
				//if (inlier.size() < prevMaxInlier / 10 && prevMaxInlier>0) {
				//	cameraPoseAbsTransEst_lin(bearingVectors1, bearingVectors2, mot, 2.0e-2, 10000);
				//	//AX=XB: A=XBX-1
				//	//Matrix4d X;
				//	//X.block(0, 0, 3, 3) = R;X.block(0, 3, 3, 1) = t;X.block(3, 0, 1, 4) << 0, 0, 0, 1;
				//	//Matrix4d Ma = X*Mb_array.at(imgid)*X.inverse();
				//	//mot = m2_6dof(Ma);
				//	Reject2D3Dcorrespondence(bearingVectors1, bearingVectors2, mot, thresholds[paramIdx], inlier);
				//}
				/*cout << maxInlier << "," << inlier.size() << endl;
				int inlierSize = inlier.size();
				if (inlierSize > maxInlier) {
					maxInlier = inlier.size();
				}*/
				//Reject2D3Dcorrespondence(bearingVectors3, bearingVectors2, thresholds[paramIdx], inlier);
				cout << maxInlier << "," << prevMaxInlier << endl;
				//txtlog << "position " << mbase << " - " << mdst << endl;
				//txtlog << "inlier " << inlier.size() << endl;
				//txtlog << mot << endl;

				in_bvs1.clear(); in_bvs_n1.clear(); in_bvs2.clear(); in_bvsi1.clear(); in_bvsi2.clear();in_pt_s1.clear();
				for (int j = 0;j < inlierImg.size();j++) {
					if (status.at(inlierImg.at(j)) != '\1' || err.at(inlierImg.at(j)) > trackerror)continue;
					in_bvsi1.push_back(bearingVectors1.at(inlierImg.at(j)).normalized());
					in_bvsi2.push_back(bearingVectors2.at(inlierImg.at(j)).normalized());
				}
				for (int j = 0;j < inlier.size();j++) {
					if (status.at(inlier.at(j)) != '\1' || err.at(inlier.at(j)) > trackerror)continue;
					in_bvs1.push_back(bearingVectors1.at(inlier.at(j)));//point camera 1 coordinates
					in_pt_s1.push_back(R.inverse()*(bearingVectors1.at(inlier.at(j)) - t));//point sensor 1 coordinates
					in_bvs2.push_back(bearingVectors2.at(inlier.at(j)));//normalized, camera 2 coordinates
					in_bvs_n1.push_back(bearingVectors1.at(inlier.at(j)).normalized());//normalized sensor 1 coordiantes
				}
				if (rotOnly) {
					cameraPoseAbsTransEst_non_lin(in_bvs1, in_bvs2, mot, thres / 2);
				}
				else {
					cameraPoseAbsEst_non_lin(in_bvs1, in_bvs2, mot, thres / 2);
				}
				thres = thres*0.9;
			}
			if (itr == 0) {
				for (int j = 0;j < in_bvs1.size();j++) {
					double phi, theta;
					Vector3d pt = R0 * in_bvs1.at(j);
					Vector3d pt2 = R1 * in_bvs2.at(j);
					omniTrans(pt(0), pt(1), pt(2), phi, theta, resynthimg1.rows);
					cv::Point2f p1(theta, phi);
					omniTrans(pt2(0), pt2(1), pt2(2), phi, theta, resynthimg1.rows);
					cv::Point2f p2(theta, phi);
					cv::circle(resynthimg1, p1, 3, cv::Scalar(255, 0, 0));
					cv::line(resynthimg1, p1, p2, cv::Scalar(0, 255, 0));
					cv::circle(resynthimg2, p2, 3, cv::Scalar(0, 0, 255));
				}
				stringstream ss, ss2;
				ss << "motion" << imgid << "_src.jpg";
				cv::imwrite(targetFolder + "/" + ss.str(), resynthimg1);
				ss2 << "motion" << imgid << "_dst.jpg";
				cv::imwrite(targetFolder + "/" + ss2.str(), resynthimg2);
			}
			if (itcnt == 0 && itr<3) {
				cameraPoseAbsEst_lin(bearingVectors1, bearingVectors2, mot, thres_init, 10000);
			}
			//txtlog << "position " << mbase << " - " << mdst << endl;
			//txtlog << "loop " << itcnt << ":inlier " << inlier.size() << endl;
			txtlog << mot <<",,";
			//gt
			Matrix4d gtA;
			
			gtA = gt * lidarPos.at(mbase).inverse() * lidarPos.at(mdst) * gt.inverse();
			_6dof gtAmot = m2_6dof(gtA);
			//txtlog << gtAmot << endl;
			Ma_array.at(imgid) = _6dof2m(mot);
			Vector3d axis;
			double anglea,angleb;
			mat2axis_angle(Ma_array.at(imgid).block(0,0,3,3),axis,anglea);
			mat2axis_angle(Mb_array.at(imgid).block(0, 0, 3, 3), axis, angleb);
			//txtlog << "angle A-B: " << anglea << "," << angleb << endl;
		}
		//txtlog << endl;

		//for (int i = 0;i < vin_bvs1.size();i++) {
		//	_6dof mot = m2_6dof(Ma_array.at(i));
		//	//cameraPoseAbsEst_non_lin(vin_bvs1.at(i), vin_bvs2.at(i), mot, losses[paramIdx]);
		//	cameraPoseAbsTransEst_non_lin(vin_bvs1.at(i), vin_bvs2.at(i), mot, losses[paramIdx]);
		//	//cameraPoseAbsEstHybrid_non_lin(vin_bvs1.at(i), vin_bvs2.at(i),vin_bvsi1.at(i),vin_bvsi2.at(i), mot, losses[paramIdx]);
		//	Ma_array.at(i) = _6dof2m(mot);
		//}

	
		solverot(Ma_array, Mb_array, R);
		solvelin(Ma_array, Mb_array, R, t);
		solvet_non_lin(Ma_array, Mb_array, R, t);
		txtlog <<  endl;
			_6dof param;
			R2axisRot(R, param.rx, param.ry, param.rz);
			param.x = t(0);
			param.y = t(1);
			param.z = t(2);
			csvlog << "itr"<<itr<<"," << param.rx << "," << param.ry << "," << param.rz << "," << param.x << "," << param.y << "," << param.z <<",,"<<paramIdx<< endl;
			
			
			_6dof mot=m2_6dof(Ma_array.at(Ma_array.size()-1));
			//relativePoseEst_non_lin_oneMotion(in_pt_s1, in_bvs_n1, in_bvs2, param, mot);
			//vector<_6dof> mota;
			//for (int i = 0;i < Ma_array.size();i++) {
			//	_6dof mota_ = m2_6dof(Ma_array.at(i));
			//	mota.push_back(mota_);
			//}
			//

			//allOptimization(vin_pt_s1, vin_bvs_n1, vin_bvs2, param, mota, Mb_array, losses[paramIdx]);
			//
			//for (int i = 0;i < Ma_array.size();i++) {
			//	Ma_array.at(i) = _6dof2m(mota.at(i));
			//}
			_6dof paramd = param - param_;
			double diffr = (paramd.rx*paramd.rx)+ (paramd.ry*paramd.ry)+ (paramd.rz*paramd.rz);
			double difft = (paramd.x*paramd.x) + (paramd.y*paramd.y) + (paramd.z*paramd.z);
			prevMaxInlier = maxInlier;
			if ((diffr < (5e-4- 1e-4*paramIdx)*(5e-4 - 1e-4*paramIdx) && difft < (5e-4 - 1e-4*paramIdx)*(5e-4 - 1e-4*paramIdx)) || paramCnt>20) {

				paramIdx+=5;
				paramCnt = 0;
				prevMaxInlier = -1;
			}


			param_ = param;
			paramCnt++;
			R = axisRot2R(param.rx, param.ry, param.rz);
			t << param.x, param.y, param.z;
			thres_init = thres_init*thresDecreaseRate;
		if (paramIdx == 5)break;
	}
	_6dof param;
	R2axisRot(R, param.rx, param.ry, param.rz);
	param.x = t(0);
	param.y = t(1);
	param.z = t(2);
	double errt, errr;
	_6dof gt6mot = m2_6dof(gt);
	_6dof errParam = gt6mot - param;
	errt = errParam.x*errParam.x + errParam.y*errParam.y + errParam.z*errParam.z;
	errt = sqrt(errt);
	errr = errParam.rx*errParam.rx + errParam.ry*errParam.ry + errParam.rz*errParam.rz;
	errr = sqrt(errr);
	csvlog << "result," << param.rx << "," << param.ry << "," << param.rz << "," << param.x << "," << param.y << "," << param.z << ",," << errr<<","<<errt << endl;
	ofs.close();
	csvlog.close();




	for (int i = 0;i < bps.size();i++) {
		bps.at(i).release();	
	}


	Rinv = R.inverse();
	tinv = -Rinv*t;
	q = dcm2q(Rinv);

	ofstream ofs2(targetFolder + "/result.cpara");
	ofs2 << "Camera Parameter, Trans" << endl;
	ofs2 << tinv(0) << " " << tinv(1) << " " << tinv(2) << endl;
	ofs2 << "Rotation" << endl;
	ofs2 << q(0) << " " << q(1) << " " << q(2) << " " << q(3) << endl;




    return 0;
}

