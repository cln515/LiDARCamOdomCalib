#include "utility.h"

ostream& operator<<(ostream& os, const _6dof& dof)
{
    os << dof.rx << ',' << dof.ry << ',' << dof.rz << ','<< dof.x << ','<< dof.y << ','<< dof.z;
    return os;
}

bool headString(string line,string chara){
	if(line.find(chara)==0)return true;
	else return false;
}
double getDist2(POINT_3D a,POINT_3D b){
		return (float)((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
	}
double getDist(POINT_3D a, POINT_3D b) {
	return sqrt(getDist2(a,b));
}

Vector4d dcm2q(Matrix3d& dcm){
	Vector4d q;
	if(dcm.trace()>0){
		double sr=sqrt(1+dcm.trace());
		double sr2=sr*2;
		
		q(0)=(dcm(1,2)-dcm(2,1))/sr2;
		q(1)=(dcm(2,0)-dcm(0,2))/sr2;
		q(2)=(dcm(0,1)-dcm(1,0))/sr2;
		q(3)=0.5*sr;	
	}else{
		if(dcm(0,0)>dcm(1,1) && dcm(0,0)>dcm(2,2)){
			double sr=sqrt(1+(dcm(0,0)-(dcm(1,1)+dcm(2,2))));
			double sr2=sr*2;
			q(3)=(dcm(1,2)-dcm(2,1))/sr2;
			q(2)=(dcm(2,0)+dcm(0,2))/sr2;
			q(1)=(dcm(0,1)+dcm(1,0))/sr2;
			q(0)=0.5*sr;
		}else if(dcm(1,1)>dcm(2,2)){
		    double  sr  = sqrt( 1 + (dcm(1,1) - ( dcm(2,2) + dcm(0,0) )) );
			double  sr2 = 2*sr;
			  q(0) = ( dcm(1,0) + dcm(0,1) ) / sr2;
			  q(1) = 0.5 * sr;
			  q(2) = ( dcm(1,2) + dcm(2,1) ) / sr2;
			  q(3) = ( dcm(2,0) - dcm(0,2) ) / sr2;

		}else{
		    double  sr  = sqrt( 1 + (dcm(2,2) - ( dcm(0,0) + dcm(1,1) )) );
		    double  sr2 = 2*sr;
			  q(0) = ( dcm(2,0) + dcm(0,2) ) / sr2;
			  q(1) = ( dcm(1,2) + dcm(2,1) ) / sr2;
			  q(2) = 0.5 * sr;
			  q(3) = ( dcm(0,1) - dcm(1,0) ) / sr2;
		}
	}
	return q;

}

Matrix3d q2dcm(Vector4d& q){
	Matrix3d R;

	// Build quaternion element products
double q1q1=q(0)*q(0);
double q1q2=q(0)*q(1);
double q1q3=q(0)*q(2);
double q1q4=q(0)*q(3);

double q2q2=q(1)*q(1);
double q2q3=q(1)*q(2);
double q2q4=q(1)*q(3);

double q3q3=q(2)*q(2);
double q3q4=q(2)*q(3);
  
double q4q4=q(3)*q(3);

	// Build DCM
R(0,0) =  q1q1 - q2q2 - q3q3 + q4q4;
R(0,1) = 2*(q1q2 + q3q4);
R(0,2) = 2*(q1q3 - q2q4);
  
R(1,0) = 2*(q1q2 - q3q4);
R(1,1) = -q1q1 + q2q2 - q3q3 + q4q4;
R(1,2) = 2*(q2q3 + q1q4);
  
R(2,0) = 2*(q1q3 + q2q4);
R(2,1) = 2*(q2q3 - q1q4);
R(2,2) = -q1q1 - q2q2 + q3q3 + q4q4;

return R;

}

Matrix3d axisRot2R(double rx,double ry,double rz){
	Matrix4d R,rotx,roty,rotz;
	double sinv,cosv;
	sinv=sin(rx),cosv=cos(rx);



			rotx<<1,0,0,0,0,cosv,-sinv,0,0,sinv,cosv,0,0,0,0,1;
			sinv=sin(ry);cosv=cos(ry);
			roty<<cosv,0,sinv,0,0,1,0,0,-sinv,0,cosv,0,0,0,0,1;
			sinv=sin(rz);cosv=cos(rz);
			rotz<<cosv,-sinv,0,0,sinv,cosv,0,0,0,0,1,0,0,0,0,1;
			R=rotx*roty*rotz;
			Matrix3d retMat=R.block(0,0,3,3);
			return retMat;
}

Matrix3d ladybug_rot2xyz (double rph[3])
{
    double sr, sp, sh, cr, cp, ch;
	sr=sin(rph[0]);cr=cos(rph[0]);
	sp=sin(rph[1]);cp=cos(rph[1]);
	sh=sin(rph[2]);ch=cos(rph[2]);

	Matrix3d R;
	R<<  ch*cp, -sh*cr + ch*sp*sr,  sh*sr + ch*sp*cr,
                             sh*cp,  ch*cr + sh*sp*sr, -ch*sr + sh*sp*cr,
                            -sp,     cp*sr,             cp*cr             ;
	return R;
}

void R2axisRot(Matrix3d R,double& rx,double& ry,double& rz){
	ry=asin(R(0,2));
		rx=-atan2(R(1,2),R(2,2));				
		rz=-atan2(R(0,1),R(0,0));


}

Matrix4d _6dof2m(_6dof dof){
	Matrix4d ret;
	ret.block(0,0,3,3)=axisRot2R(dof.rx,dof.ry,dof.rz);
	ret(0,3)=dof.x;
	ret(1,3)=dof.y;
	ret(2,3)=dof.z;
	ret.block(3,0,1,4)<<0,0,0,1;
	return ret;
};
_6dof m2_6dof(Matrix4d& m){
	Matrix3d r=m.block(0,0,3,3);
	_6dof dof;
	R2axisRot(r,dof.rx,dof.ry,dof.rz);
	dof.x=m(0,3);
	dof.y=m(1,3);
	dof.z=m(2,3);
	return dof;

};


void int2rgba(int color,unsigned char& r,unsigned char& g,unsigned char& b,unsigned char& a){
		r=(color&0x0000ff);
		g=(color&0x00ff00)>>8;
		b=(color&0xff0000)>>16;
		a=(color&0xff000000)>>24;
}

void rgba2int(int& color,unsigned char r,unsigned char g,unsigned char b,unsigned char a){
		color=((a&0x0000ff)<<24)|((b&0x0000ff)<<16)|((g&0x0000ff)<<8)|((r&0x0000ff));
}

int getSubpixelColor(int topLeftColor,int topRightColor,int bottomLeftColor,int bottomRightColor,double dx,double dy){
	unsigned char rs[4],gs[4],bs[4],as[4];
	unsigned char r,g,b,a;
	int2rgba(topLeftColor,rs[0],gs[0],bs[0],as[0]);
	int2rgba(topRightColor,rs[1],gs[1],bs[1],as[1]);
	int2rgba(bottomLeftColor,rs[2],gs[2],bs[2],as[2]);
	int2rgba(bottomRightColor,rs[3],gs[3],bs[3],as[3]);
	r=rs[0]*(1-dx)*(1-dy)+rs[1]*dx*(1-dy)+rs[2]*(1-dx)*(dy)+rs[3]*dx*dy;
	g=gs[0]*(1-dx)*(1-dy)+gs[1]*dx*(1-dy)+gs[2]*(1-dx)*(dy)+gs[3]*dx*dy;
	b=bs[0]*(1-dx)*(1-dy)+bs[1]*dx*(1-dy)+bs[2]*(1-dx)*(dy)+bs[3]*dx*dy;
	a=as[0]*(1-dx)*(1-dy)+as[1]*dx*(1-dy)+as[2]*(1-dx)*(dy)+as[3]*dx*dy;
	int color;
	rgba2int(color,r,g,b,a);
	return color;
}


void getSubpixelColor(unsigned char* topLeftColor,unsigned char* topRightColor,unsigned char* bottomLeftColor,unsigned char* bottomRightColor,double dx,double dy,unsigned char* rgb){
	unsigned char rs[4],gs[4],bs[4],as[4];
//	unsigned char r,g,b,a;
	rgb[0]=topLeftColor[0]*(1-dx)*(1-dy)+topRightColor[0]*dx*(1-dy)+bottomLeftColor[0]*(1-dx)*(dy)+bottomRightColor[0]*dx*dy;
	rgb[1]=topLeftColor[1]*(1-dx)*(1-dy)+topRightColor[1]*dx*(1-dy)+bottomLeftColor[1]*(1-dx)*(dy)+bottomRightColor[1]*dx*dy;
	rgb[2]=topLeftColor[2]*(1-dx)*(1-dy)+topRightColor[2]*dx*(1-dy)+bottomLeftColor[2]*(1-dx)*(dy)+bottomRightColor[2]*dx*dy;
	int color;

}

void FAROColorPtx2ply(string in_ptxFn, string out_plyFn) {
	ifstream ifs(in_ptxFn, ios::binary);
	string line_str;
	unsigned int* lineStartIdx = (unsigned int*)malloc(sizeof(int) * 3000);
	unsigned int* idxpointer = lineStartIdx;
	getline(ifs, line_str);
	int line = stoi(line_str);
	getline(ifs, line_str);
	int oneLine = stoi(line_str);
	for (int i = 0;i<8;i++) {
		getline(ifs, line_str);
	}
	vector<float> vert;


	int lineNum = 0, vertCnt = 0;

	lineStartIdx[lineNum] = vert.size() / 3;
	float thresMax = 120 * 120;
	float thresMin = 2 * 2;
	POINT_3D* arrayPoint = (POINT_3D*)malloc(sizeof(POINT_3D)*line*oneLine);
	float* reflectanceArray = (float*)malloc(sizeof(float)*line*oneLine);
	unsigned char* colorArray = (unsigned char*)malloc(sizeof(unsigned char)*3*line*oneLine);

	while (getline(ifs, line_str)) {
		float x, y, z;
		float f[5];
		unsigned char rgb[4];
		for (int i = 0;i<7;i++) {
			if (i < 4)f[i] = stof(line_str.substr(0, line_str.find_first_of(" ")));
			else if (i < 6) {
				rgb[i - 4] = (unsigned char)(stoi(line_str.substr(0, line_str.find_first_of(" "))));
			}
			else rgb[i - 4] = (unsigned char)(stoi(line_str));
			line_str.erase(0, line_str.find_first_of(" ") + 1);
		}

		arrayPoint[vertCnt].x = f[0];
		arrayPoint[vertCnt].y = f[1];
		arrayPoint[vertCnt].z = f[2];
		reflectanceArray[vertCnt] = f[4];//scan48:4
		vertCnt++;
	}
	vector<float> refVect;
	vector<unsigned int>vertIdx;
	for (int i = 0;i<line;i++) {
		lineStartIdx[i] = vert.size() / 3;
		cout << "line:" << i << " start:" << lineStartIdx[i] << endl;
		for (int j = 0;j<oneLine;j++) {
			int idx = i*oneLine + j;
			float r2 = arrayPoint[idx].x*arrayPoint[idx].x + arrayPoint[idx].y*arrayPoint[idx].y + arrayPoint[idx].z*arrayPoint[idx].z;
			if (r2>thresMax || r2<thresMin || reflectanceArray[idx]<0.01)continue;
			vert.push_back(arrayPoint[idx].x);
			vert.push_back(arrayPoint[idx].y);
			vert.push_back(arrayPoint[idx].z);
			vertIdx.push_back(idx);//多分ここバグ
			refVect.push_back(reflectanceArray[idx]);
		}
	}
	int vertexNum = vert.size() / 3;
	float* vertex = (float*)malloc(sizeof(float)*vert.size());
	for (int i = 0;i<vert.size();i++) {
		vertex[i] = vert.at(i);
	}
	unsigned int* vertexIdx = (unsigned int*)malloc(sizeof(unsigned int)*vertIdx.size());

	vert.clear();
	free(reflectanceArray);
	free(arrayPoint);

}

void timeSequencePtx2ply(string in_ptxFn,string out_plyFn){
	ifstream ifs(in_ptxFn,ios::binary);
	string line_str;
	unsigned int* lineStartIdx=(unsigned int*)malloc(sizeof(int)*3000);
	unsigned int* idxpointer=lineStartIdx;
	getline(ifs,line_str);
	int line=stoi(line_str);
	getline(ifs,line_str);
	int oneLine=stoi(line_str);
	for(int i=0;i<8;i++){
	getline(ifs,line_str);
	}
	vector<float> vert;
	

	int lineNum=0,vertCnt=0;

	lineStartIdx[lineNum]=vert.size()/3;
	float thresMax=120*120;
	float thresMin=2*2;
	POINT_3D* arrayPoint=(POINT_3D*)malloc(sizeof(POINT_3D)*line*oneLine);
	float* reflectanceArray=(float*)malloc(sizeof(float)*line*oneLine);

	while(getline(ifs,line_str)){
		float x,y,z;
		float f[5];
		for(int i=0;i<5;i++){
			if(i!=4)f[i]=stof(line_str.substr(0,line_str.find_first_of(" ")));
			else f[i]=stof(line_str);
			line_str.erase(0,line_str.find_first_of(" ")+1);
		}

		arrayPoint[vertCnt].x=f[0];
		arrayPoint[vertCnt].y=f[1];
		arrayPoint[vertCnt].z=f[2];
		reflectanceArray[vertCnt]=f[4];//scan48:4
		vertCnt++;
	}
	vector<float> refVect;
	vector<unsigned int>vertIdx;
	for(int i=0;i<line;i++){
			lineStartIdx[i]=vert.size()/3;
			cout<<"line:"<<i<<" start:"<<lineStartIdx[i]<<endl;
		for(int j=0;j<oneLine;j++){
			int idx=i*oneLine+j;
			float r2=arrayPoint[idx].x*arrayPoint[idx].x+arrayPoint[idx].y*arrayPoint[idx].y+arrayPoint[idx].z*arrayPoint[idx].z;
			if(r2>thresMax||r2<thresMin||reflectanceArray[idx]<0.01)continue;
			vert.push_back(arrayPoint[idx].x);
			vert.push_back(arrayPoint[idx].y);
			vert.push_back(arrayPoint[idx].z);
			vertIdx.push_back(idx);//多分ここバグ
			refVect.push_back(reflectanceArray[idx]);
		}
	}
	int vertexNum=vert.size()/3;
	float* vertex=(float*)malloc(sizeof(float)*vert.size());
	for(int i=0;i<vert.size();i++){
		vertex[i]=vert.at(i);
	}
	unsigned int* vertexIdx=(unsigned int*)malloc(sizeof(unsigned int)*vertIdx.size());
	
	vert.clear();
	free(reflectanceArray);
	free(arrayPoint);

}

Matrix4d getMatrixFlomPly(string fn){
	ifstream ifs(fn,ios::binary);
	string line;
	Matrix4d globalPose=Matrix4d::Identity();
	int n=0;
	while(getline(ifs,line)){
//		ofs<<line<<endl;
		cout<<line<<endl;
		if(headString(line,"matrix")){
			float f[4];
			int i;
			for(i=0;i<5;i++){
				if(i!=0)f[i-1]=stof(line.substr(0,line.find_first_of(" ")));
				line.erase(0,line.find_first_of(" ")+1);
			}
			globalPose(0,n)=f[0];
			globalPose(1,n)=f[1];
			globalPose(2,n)=f[2];
			globalPose(3,n)=f[3];
			n++;
		}
		if(headString(line,"end_header"))break;		

	}
	ifs.close();
	return globalPose;
}



void writePlyHeader(ofstream& ofs,int vertSize,int faceSize){
	ofs<<"ply"<<endl
	<<"format binary_little_endian 1.0"<<endl
	<<"element vertex "<<vertSize<<endl
	<<"property float x"<<endl
	<<"property float y"<<endl
	<<"property float z"<<endl
	<<"property float confidence"<<endl
	<<"property float intensity"<<endl
	<<"element face "<<faceSize<<endl
	<<"property list uchar int vertex_index"<<endl
	<<"end_header"<<endl;
}

void writePlyHeaderRGB(ofstream& ofs,int vertSize,int faceSize){
	ofs<<"ply"<<endl
	<<"format binary_little_endian 1.0"<<endl
	<<"element vertex "<<vertSize<<endl
	<<"property float x"<<endl
	<<"property float y"<<endl
	<<"property float z"<<endl
	<<"property float confidence"<<endl
	<<"property uchar red"<<endl
	<<"property uchar green"<<endl
	<<"property uchar blue"<<endl
	<<"property uchar alpha"<<endl
	<<"element face "<<faceSize<<endl
	<<"property list uchar int vertex_index"<<endl
	<<"end_header"<<endl;
}

void writePlyOnePointRGB(ofstream& ofs,Vector3f& p,unsigned char* rgba){
	float xyz[]={p(0),p(1),p(2)};
	ofs.write((char*) xyz,sizeof(float)*3);
	float conf=1.0;
	ofs.write((char*)&conf,sizeof(float));
	ofs.write((char*)rgba,sizeof(unsigned char)*4);

}

//p_3d=M*p_2d: cam->3D scanner translation
Matrix4d readVannoPara(string fileName){
					ifstream ifs(fileName,ios::binary);
					string line;	
					int n=0;
					getline(ifs,line);//"Camera position"
					getline(ifs,line);
					float position[3];
					for(int i=0;i<3;i++){
						if(i!=2)position[i]=stof(line.substr(0,line.find_first_of(" ")));
						else position[i]=stof(line);
						line.erase(0,line.find_first_of(" ")+1);
					}
					getline(ifs,line);//"Camera pose"
					getline(ifs,line);
					float pose[4];
					for(int i=0;i<4;i++){
						if(i!=3)pose[i]=stof(line.substr(0,line.find_first_of(" ")));
						else pose[i]=stof(line);
						line.erase(0,line.find_first_of(" ")+1);
					}
					ifs.close();
					Matrix4d P;

					float q3=pose[0],
					q0=pose[1],
					q1=pose[2],
					q2=pose[3];
					//阪野さんのcparaファイルはおそらくladybugの座標が正しい定義よりz軸回りに180°回転した定義になってる
					P<<-(1-2.0*(q1*q1+q2*q2)),
						-2*(q0*q1+q2*q3),
						2*(q2*q0-q1*q3),
						position[0],

					-2*(q0*q1-q2*q3),
					-(1.0-2.0*(q2*q2+q0*q0)),
					2*(q1*q2+q0*q3),
					position[1],

					-2*(q2*q0+q1*q3),
					-2*(q1*q2-q0*q3),
					1.0-2.0*(q1*q1+q0*q0),
					position[2],

					0,0,0,1;


//					Matrix4d rot180;
//					rot180<<-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1;

	return P;
}

//p_3d=M*p_2d: cam->3D scanner translation
Matrix4d readCPara(string fileName){
					ifstream ifs(fileName,ios::binary);
					string line;	
					int n=0;
					getline(ifs,line);//"Camera position"
					getline(ifs,line);
					float position[3];
					for(int i=0;i<3;i++){
						if(i!=2)position[i]=stof(line.substr(0,line.find_first_of(" ")));
						else position[i]=stof(line);
						line.erase(0,line.find_first_of(" ")+1);
					}
					getline(ifs,line);//"Camera pose"
					getline(ifs,line);
					float pose[4];
					for(int i=0;i<4;i++){
						if(i!=3)pose[i]=stof(line.substr(0,line.find_first_of(" ")));
						else pose[i]=stof(line);
						line.erase(0,line.find_first_of(" ")+1);
					}
					ifs.close();
					Vector4d q;
					q<<pose[0],pose[1],pose[2],pose[3];
					Matrix4d P;
					P.block(0,0,3,3)=q2dcm(q);
					P.block(0,3,3,1)<<position[0],position[1],position[2];
					P.block(3,0,1,4)<<0,0,0,1;

			/*		float q3=pose[0],
					q0=pose[1],
					q1=pose[2],
					q2=pose[3];
					P<<-(1-2.0*(q1*q1+q2*q2)),
						-2*(q0*q1+q2*q3),
						2*(q2*q0-q1*q3),
						position[0],

					-2*(q0*q1-q2*q3),
					-(1.0-2.0*(q2*q2+q0*q0)),
					2*(q1*q2+q0*q3),
					position[1],

					-2*(q2*q0+q1*q3),
					-2*(q1*q2-q0*q3),
					1.0-2.0*(q1*q1+q0*q0),
					position[2],

					0,0,0,1;


					Matrix4d rot180;
					rot180<<-1,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,1;*/

	return P;
}


void omniTrans(double x,double y, double z,double& phi,double& theta){
	double r=sqrt(x*x+y*y+z*z);
	theta=atan2(y,x);
	phi=acos(z/r);
};

void omniTrans(double x,double y, double z,double& phi_pix,double& theta_pix,int height){
	double r=sqrt(x*x+y*y+z*z);
	double theta=atan2(y,x);
	double phi=acos(z/r);
	int width=height*2;
	theta_pix=(-theta+M_PI)*width/(M_PI*2);
	phi_pix=phi/M_PI*height;
};

void rev_omniTrans(double x2d,double y2d,int width,int height,Vector3d& ret){
		double theta=-(x2d*(M_PI*2)/width-M_PI);
		double phi=y2d/height*M_PI;
		ret<<sin(phi)*cos(theta),sin(phi)*sin(theta),cos(phi);

}

bool makeFolder(string folderPath){
	
	if (MakeSureDirectoryPathExists(folderPath.c_str())) {
	    cout<<"succeeded."<<endl;
		return true;
	  } else {
		cout<<"failed"<<endl;
		return false;
	}
	
}

Vector3d cubic_spline(Vector3d& p0,Vector3d& p1,Vector3d& p2,Vector3d& p3,double t0,double t1,double t2,double t3,double t){
	Vector3d a1=(t1-t)/(t1-t0)*p0+(t-t0)/(t1-t0)*p1;
	Vector3d a2=(t2-t)/(t2-t1)*p1+(t-t1)/(t2-t1)*p2;
	Vector3d a3=(t3-t)/(t3-t2)*p2+(t-t2)/(t3-t2)*p3;
	Vector3d b1=(t2-t)/(t2-t0)*a1+(t-t0)/(t2-t0)*a2;
	Vector3d b2=(t3-t)/(t3-t1)*a2+(t-t1)/(t3-t1)*a3;
	return 	(t2-t)/(t2-t1)*b1+(t-t1)/(t2-t1)*b2;

};
Vector3d cubic_spline(Vector3d& p0,Vector3d& p1,Vector3d& p2,Vector3d& p3,double t){return cubic_spline(p0,p1,p2,p3,0,1/3.0,2/3.0,1,t);}

double get2Line_Distance(Vector3d& p1,Vector3d& v1,Vector3d& p2,Vector3d& v2){
	Vector3d prod=v1.cross(v2);
	double prod_d=prod.norm();
	Vector3d p1p2=p2-p1;
	if(prod_d==0){
		return v1.cross(p1p2).norm()/v1.norm();
	}else{
		return abs(prod.dot(p1p2)/prod_d);
	}
}

//
double get2Line_Distance(Vector3d& p1,Vector3d& v1,Vector3d& p2,Vector3d& v2,Vector3d& r1,Vector3d& r2){
	Vector3d nv1=v1.normalized();
	Vector3d nv2=v2.normalized();
	Vector3d p1p2=p1-p2;
	double d=1-nv1.dot(nv2)*nv1.dot(nv2);
	double alpha=(nv1.dot(nv2)*(nv2.dot(p1p2))-nv1.dot(p1p2))/d;
	double beta=(-nv1.dot(nv2)*(nv1.dot(p1p2))+nv2.dot(p1p2))/d;
	r1=p1+alpha*nv1;
	r2=p2+beta*nv2;
	return (r2-r1).norm();


}
double get_point2lineDistance(Vector3d& p1,Vector3d& v1,Vector3d& p2){
	return (p2-p1).cross(v1).norm()/v1.norm();
}

string getTimeStamp() {
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 80, "%y%m%d_%H%M%S", timeinfo);
	string ts(buffer);
	return ts;
}