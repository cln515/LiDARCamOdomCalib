#define _USE_MATH_DEFINES
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen\Eigen>
#include <Eigen\Core>
#include <windows.system.h>
#include <time.h> 
#include <imagehlp.h>
#pragma comment(lib, "imagehlp.lib")

using namespace std;
using namespace Eigen;

#ifndef MYUTILITY
#define MYUTILITY

struct POINT_3D{
	float x,y,z;
};
struct POINT3D_N :POINT_3D{
	float normx,normy;
	//
	int meshindex;
};
struct MESH{
	int index1,index2,index3;
};

struct _6dof{
	double rx,ry,rz;
	double x,y,z;
	friend ostream& operator<<(ostream& os, const _6dof& dt);
};




inline _6dof operator+(_6dof a,_6dof b){
	_6dof c={
		a.rx+b.rx,
		a.ry+b.ry,
		a.rz+b.rz,
		a.x+b.x,
		a.y+b.y,
		a.z+b.z,
	};
	return c;
}

inline _6dof operator-(_6dof a,_6dof b){
	_6dof c={
		a.rx-b.rx,
		a.ry-b.ry,
		a.rz-b.rz,
		a.x-b.x,
		a.y-b.y,
		a.z-b.z,
	};
	return c;
}

inline _6dof operator*(double a,_6dof b){
	_6dof c={
		a*b.rx,
		a*b.ry,
		a*b.rz,
		a*b.x,
		a*b.y,
		a*b.z,
	};
	return c;
}


bool headString(string line,string chara);

double getDist2(POINT_3D a,POINT_3D b);
double getDist(POINT_3D a, POINT_3D b);


struct matching{
	POINT_3D objectPoint;
	POINT_3D correspondPoint;
};

Vector4d dcm2q(Matrix3d& dcm);
Matrix3d q2dcm(Vector4d& q);
Matrix3d axisRot2R(double rx,double ry,double rz);
Matrix3d ladybug_rot2xyz (double rph[3]);
void R2axisRot(Matrix3d R,double& rx,double& ry,double& rz);

Matrix4d _6dof2m(_6dof dof);
_6dof m2_6dof(Matrix4d& m);

Matrix4d getMatrixFlomPly(string fn);

void writePlyHeader(ofstream& ofs,int vertSize,int faceSize);
void writePlyHeaderRGB(ofstream& ofs,int vertSize,int faceSize);
void writePlyOnePointRGB(ofstream& ofs,Vector3f& p,unsigned char* rgba);


void timeSequencePtx2ply(string in_ptxFn,string out_plyFn);

void omniTrans(double x,double y, double z,double& phi,double& theta);
void omniTrans(double x,double y, double z,double& phi,double& theta,int height);
void rev_omniTrans(double ix,double iy,int width,int height,Vector3d& ret);

void int2rgba(int color,unsigned char& r,unsigned char& g,unsigned char& b,unsigned char& a);
void rgba2int(int& color,unsigned char r,unsigned char g,unsigned char b,unsigned char a);
int getSubpixelColor(int topLeftColor,int topRightColor,int bottomLeftColor,int bottomRightColor,double dx,double dy);
void getSubpixelColor(unsigned char* topLeftColor,unsigned char* topRightColor,unsigned char* bottomLeftColor,unsigned char* bottomRightColor,double dx,double dy,unsigned char* rgb);

//double miComputing(double* histogram,int width,int height);
//double miComputing(double* histogram,int width,int height,int start);

//
//Matrix4d readVannoPara(string fileName);
Matrix4d readCPara(string fileName);

bool makeFolder(string folderPath);

Vector3d cubic_spline(Vector3d& p0,Vector3d& p1,Vector3d& p2,Vector3d& p3,double t0,double t1,double t2,double t3,double t);
Vector3d cubic_spline(Vector3d& p0,Vector3d& p1,Vector3d& p2,Vector3d& p3,double t);

double get2Line_Distance(Vector3d& p1,Vector3d& v1,Vector3d& p2,Vector3d& v2);
double get2Line_Distance(Vector3d& p1,Vector3d& v1,Vector3d& p2,Vector3d& v2,Vector3d& r1,Vector3d& r2);
double get_point2lineDistance(Vector3d& p1,Vector3d& v1,Vector3d& p2);

string getTimeStamp();

#endif


