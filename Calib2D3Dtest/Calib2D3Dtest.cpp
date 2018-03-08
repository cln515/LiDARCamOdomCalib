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
			residual[0] = err.norm();
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

void solvetlin(vector<Matrix4d>& pepdMat, vector<Matrix4d> holdMat, Matrix3d R, Vector3d& t) {
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
	ans = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(B);
	for (int i = 0;i < pepdMat.size();i++) {
		pepdMat.at(i).block(0, 3, 3, 1) = ans(i+3) *pepdMat.at(i).block(0, 3, 3, 1);
	}
	t << ans(0), ans(1), ans(2);

}


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

			residual[0] = eyeDirec.cross(plot).norm();
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

			residual[0] = err.norm();
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
		Matrix3d RA = pepdMat.at(i).block(0, 0, 3, 3);
		Matrix3d RB = holdMat.at(i).block(0, 0, 3, 3);
		rotErrCostFunc* p= new rotErrCostFunc(RA, RB);
		ceres::CostFunction* c = new ceres::NumericDiffCostFunction<rotErrCostFunc, ceres::CENTRAL, 1, 3>(
			p);
		problem.AddResidualBlock(c,new ceres::HuberLoss(1.0e-2), opt);
	}
	MatrixXd KBKA = KB*KA.transpose();
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
	ceres::Solve(options, &problem, &summary);
	R = axisRot2R(opt[0], opt[1], opt[2]);
}



//Main code=========================================================================================================================
int main(int argc,char* argv[])
{
	//:bool rotOnly = false;
	//bool trackerStrategy = false;//initial position of the tracked point
	//outlier rejection
	double thres_init = 1e-2;
	double thres_iter = 5e-2;
	double thresDecreaseRate = 0.5;
	double trackerror = 2.0;
	double cornerThres = 0.0;
	int outloopInlier = 500;
	//initial parameter
	//bool giveInitPara = false;

	if (argc < 2) {
		return 0;
	}
	
	string targetFolder=argv[1];
	ofstream csvlog(targetFolder + "\\log" + getTimeStamp() + ".csv");
	ifstream motionfile(targetFolder+"\\motion.txt");
	int m_cnt = stoi(argv[2]);
	std::random_device rdev{};
	std::mt19937 mt(rdev());
	vector<int> vec;
	vector<int> vec2;
	string str;
	vector<int> motionList_,motionList;

	//get motions
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
	
	//read ground truth
	Matrix4d gt;
	if (true) {
		Matrix4d m1 = getMatrixFlomPly("C:\\data\\calib\\velo2\\tozf.ply");//velodyne  - panorama range sensor
		Matrix4d m2 = readCPara("C:\\data\\calib\\velo2\\tozf.cpara");//ladybug  - panorama range sensor
		gt = m2.inverse()*m1;
		cout << m2_6dof(m2) << endl;
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
	//load 2D and 3D data
	vector<cv::Mat> imgs,descriptors;
	vector<string>plyFileList,imgFileList;
	plyFileList=getFileWithNumber(targetFolder, "position", "ply");
	imgFileList=getFileWithNumber(targetFolder, "position", "jpg");

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
	//get LiDAR position	
	for (int i = 0;i < plyFileList.size();i++) {
		Matrix4d m1 = getMatrixFlomPly(plyFileList.at(i));
		lidarPos.push_back(m1);
	}
	Matrix3d R;
	Vector3d t;
	vector<Matrix4d> Ma_array, Mb_array;
	//initial motion estimation
	for (int i = 0;i < motionNumber;i++) {
		int mbase = motionList.at(i * 2);
		int mdst = motionList.at(i * 2 + 1);
		Mb_array.push_back(lidarPos.at(mbase).inverse()*lidarPos.at(mdst));
	}
	auto algorithm = cv::AKAZE::create();
		vector<vector<cv::KeyPoint>> keypoints;
		for (int i = 0;i < imgs.size();i++) {
			vector<cv::KeyPoint> keypoint;
			cv::Mat descriptor;
			if (included.at(i)) {
				cout << "image " << i << "keypoint detection" << endl;
				algorithm->detect(imgs.at(i), keypoint);
				algorithm->compute(imgs.at(i), keypoint, descriptor);
			}
			keypoints.push_back(keypoint);
			descriptors.push_back(descriptor);
		}
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
			cameraPoseEst_lin(bearingVectors1, bearingVectors2, mot, thres_init, 10000);
			double thres = thres_init;
			vector<int>inlier;
			vector<Vector3d> in_bvs1, in_bvs2;
			for (int loop = 0;loop < 60;loop++) {
				inlier.clear();
				Reject2D2Dcorrespondence(bearingVectors1, bearingVectors2, mot, thres, inlier);
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
		solvetlin(Ma_array, Mb_array, R, t);
		cout << "initial rotation" << endl;
	cout << R << endl;

	cout << "initial translation" << endl;
	cout << t << endl;

	//log output
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
	//motion estimation
	cout << "refine phase" << endl;
	for (int itr = 0;itr < 30;itr++) {
		cout << "Iteration "<< itr << endl;
		for (int imgid = 0;imgid < motionNumber;imgid++) {
			int mbase = motionList.at(imgid * 2);
			int mdst = motionList.at(imgid * 2 + 1);
			cv::Mat resynthimg1, resynthimg2;
			Matrix3d R0, R1;
			Matrix3d Ra = Ma_array.at(imgid).block(0, 0, 3, 3);
			Vector3d ta = Ma_array.at(imgid).block(0, 3, 3, 1);

			//panorama image is rectified as
			// epipolar-> z axis
			// (0,0,1)^T=R0*ta
			// R1=R*R0
			panoramaRectification(imgs.at(mbase), imgs.at(mdst),resynthimg1,resynthimg2,ta, Ra,R0,R1);
			
			Matrix3d R1inv = R1.inverse();

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
				Vector3d pt, pt2, pt3;
				pt << vps[i * 3], vps[i * 3 + 1], vps[i * 3 + 2];
				pt3 = R*(mrev2.block(0, 0, 3, 3)*pt + mrev2.block(0, 3, 3, 1)) + t;//lidar2->camera2
				pt = R*(mrev.block(0, 0, 3, 3)*pt + mrev.block(0, 3, 3, 1)) + t; //camera1 cooridate
				pt2 = R0 * pt; //camera rect coordinates
				if (pt2.norm() < 0.5)continue;
				double ix, iy;
				omniTrans(pt2(0), pt2(1), pt2(2), iy, ix, resynthimg1.rows);
				cv::Point2f p(ix, iy);
				if (iy <= 0 || ix <= 0 || iy > resynthimg1.rows - 1 || ix > resynthimg1.cols - 1)continue;
				bearingVectors1c.push_back(pt);
				pt2 = R1*(Ra.inverse()*(pt - ta));
				bearingVectors2c.push_back(pt3);
				srcc.push_back(p);
				//initial point
				pt3 = R1*pt3;
				omniTrans(pt3(0), pt3(1), pt3(2), iy, ix, resynthimg1.rows);
				cv::Point2f p3(ix, iy);
				dstc.push_back(p3);//lidar1->lidar2->camera2 
				
			}
			vector<int>selectedIdx;
			goodFeatureToTrack_onProjection(img1g, srcc, selectedIdx, 3.0, 0, cornerThres);
			for (int sIdx = 0;sIdx < selectedIdx.size();sIdx++) {
				src.push_back(srcc.at(selectedIdx.at(sIdx)));
				bearingVectors1.push_back(bearingVectors1c.at(selectedIdx.at(sIdx)));
				bearingVectors3.push_back(bearingVectors2c.at(selectedIdx.at(sIdx)));
				dst.push_back(dstc.at(selectedIdx.at(sIdx)));
			}
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
				cameraPoseAbsEst_lin(bearingVectors1, bearingVectors2, mot, 1.0e-2, 10000);
				//cameraPoseAbsTransEst_lin(bearingVectors1, bearingVectors2, mot, 1.0e-2, 10000);
			}
			else {
				Matrix4d Ma = Ma_array.at(imgid);
				mot=m2_6dof(Ma);
			}
			vector<int>inlier, inlierImg;
			double thres = thres_iter;
			vector<Vector3d> in_bvs1, in_bvs2;
			for (int loop = 0;loop < 60;loop++) {
				Reject2D3Dcorrespondence(bearingVectors1, bearingVectors2, mot, thres, inlier);
				if (inlier.size() < outloopInlier) {
					break;
				}
				in_bvs1.clear();  in_bvs2.clear();
				for (int j = 0;j < inlier.size();j++) {
					if (status.at(inlier.at(j)) != '\1' || err.at(inlier.at(j)) > trackerror)continue;
					in_bvs1.push_back(bearingVectors1.at(inlier.at(j)));//point camera 1 coordinates
					in_bvs2.push_back(bearingVectors2.at(inlier.at(j)));//normalized, camera 2 coordinates
				}
				cameraPoseAbsEst_non_lin(in_bvs1, in_bvs2, mot, thres / 2);
				thres = thres*0.9;
			}
			if (false) {//Image Output
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
			}
			Ma_array.at(imgid) = _6dof2m(mot);
		}
	
		solverot(Ma_array, Mb_array, R);
		solvelin(Ma_array, Mb_array, R, t);
		solvet_non_lin(Ma_array, Mb_array, R, t);
		_6dof param;
		R2axisRot(R, param.rx, param.ry, param.rz);
		param.x = t(0);
		param.y = t(1);
		param.z = t(2);
		csvlog << "itr"<<itr<<"," << param.rx << "," << param.ry << "," << param.rz << "," << param.x << "," << param.y << "," << param.z <<",,"<<paramIdx<< endl;
			
		_6dof mot=m2_6dof(Ma_array.at(Ma_array.size()-1));
		_6dof paramd = param - param_;
		double diffr = (paramd.rx*paramd.rx)+ (paramd.ry*paramd.ry)+ (paramd.rz*paramd.rz);
		double difft = (paramd.x*paramd.x) + (paramd.y*paramd.y) + (paramd.z*paramd.z);
		
		param_ = param;
		R = axisRot2R(param.rx, param.ry, param.rz);
		t << param.x, param.y, param.z;
		thres_init = thres_init*thresDecreaseRate;
		if ((diffr < (5e-4)*(5e-4) && difft < (5e-4)*(5e-4)) || paramCnt>20) {
			break;
		}
		paramCnt++;

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
	csvlog.close();
	for (int i = 0;i < bps.size();i++) {
		bps.at(i).release();	
	}
	cout << "optimized rotation" << endl;
	cout << R << endl;

	cout << "optimized translation" << endl;
	cout << t << endl;

	return 0;
}

