
#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <string>
#include "utility.h"
#include <vector>

#include "Eigen/Core"
#include "Eigen/Eigen"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

#ifndef BASICPLY
#define BASICPLY



class BasicPly{
private:
	typedef struct {
		unsigned char	nindex_;
		unsigned int			index0_, index1_, index2_;
	}ply_index;


protected:
	vector<POINT_3D> dataPtsVector;

	Matrix4d GlobalPose;
	Vector3d g;
	float* verteces;
	float* norm;
	float* reflectance;
	unsigned int* faces;
	unsigned char* rgba;//{r,g,b,a,r,g,b,a,...}
	__int64 facenum;
	__int64 vertexnum;
	double alpha;
	bool isANN;
	bool bRead;

public:
	BasicPly(){
		vertexnum=0;
		facenum=0;
		bRead=false;
	}

	bool readPlyFile(vector<string> fileName,int dataNum);
	bool readPlyFileRGB(vector<string> fileName,int dataNum);

	void panoramaTexture(unsigned char* rgbArray,int width,int height,Matrix4d& transMat);

	int getVertexNumber(){return vertexnum;};
	int getFaceNumber(){return facenum;};
	void release();
	float* getVertecesPointer(){return verteces;};
	void setVertecesPointer(float* vertices,int vtnum){verteces=vertices;vertexnum=vtnum;};
	float* getNormPointer(){return norm;};
	float* getReflectancePointer(){return reflectance;};
	void setReflectancePointer(float* reflectance_,int vtnum){reflectance=reflectance_;vertexnum=vtnum;};
	unsigned int* getFaces(){return faces;};
	unsigned char* getRgbaPointer(){return rgba;};
	Vector3d getCentroid(){return g;};

	void writePlyFile(string fileName);
	void writePlyFileRGB(string fileName);
	void writePlyFileRGBForMeshlab(string fileName);

};

#endif // !REFDATA