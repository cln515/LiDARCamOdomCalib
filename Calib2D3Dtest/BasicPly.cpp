
#include "BasicPly.h"

using namespace std;

POINT_3D getNorm(POINT_3D p1,POINT_3D p2,POINT_3D p3);




bool BasicPly::readPlyFile(vector<string> fileName,int dataNum){
	
	int rf=0;
//	vector<float> vec;
	float* vec=(float*)malloc(1);
	float* norm_=(float*)malloc(1);
	float* reflectance_=(float*)malloc(1);
	faces=(unsigned int*)malloc(1);
	int validVertex=0;
	int facesNum=0;
	double gx=0,gy=0,gz=0;
	for(rf=0;rf<dataNum;rf++){
	ifstream ifs(fileName[rf],ios::binary);
	string line;
	string format;
	long n=0;
	int xi=-1,yi=-1,zi=-1,reflecti=-1;
	
	int paranum=0;
	int vertex=0;
	int face=0;
	int matIdx=0;
	while(getline(ifs,line)){
		cout<<line<<endl;
		if(headString(line,"format")){
			format=line.erase(0,7);
		}
		if(headString(line,"end_header"))break;		
		if(headString(line,"property")){
			
			line.erase(0,9);
			if(headString(line,"float")){
				line.erase(0,6);
				if(headString(line,"x"))xi=paranum;
				else if(headString(line,"y"))yi=paranum;
				else if(headString(line,"z"))zi=paranum;
				else if(headString(line,"intensity"));
				else if(headString(line,"confidence"))reflecti=paranum;
				else return false;
			}
			paranum++;
		}
		if(headString(line,"element")){
			line.erase(0,8);
			if(headString(line,"vertex")){
				line.erase(0,7);
				vertex=stoi(line);
			
			}
			if(headString(line,"face")){
				line.erase(0,5);
				face=stoi(line);			
			}
		}
		if(headString(line,"matrix")){
		float f[4];

		for(int i=0;i<5;i++){
			line.erase(0,line.find_first_not_of(" "));
			if(i>0)f[i-1]=stof(line.substr(0,line.find_first_of(" ")));
			line.erase(0,line.find_first_of(" "));
		}
			GlobalPose(0,matIdx)=f[0];
			GlobalPose(1,matIdx)=f[1];
			GlobalPose(2,matIdx)=f[2];
			GlobalPose(3,matIdx)=f[3];
			matIdx++;
		}
	}
	float a=0;
	int time=0;
	int num=0;
//	paranum--;
	int offset=validVertex*3;
	int offseti=facesNum*3;
	float* tmp;
	unsigned int*tmp_i;
	tmp=(float*)realloc(norm_,sizeof(float)*(offset+vertex*3));
	if(tmp){norm_=tmp;}
	else exit(EXIT_FAILURE);
	tmp=(float*)realloc(vec,sizeof(float)*(offset+vertex*3));
	if(tmp){vec=tmp;}
	else exit(EXIT_FAILURE);

	tmp=(float*)realloc(reflectance_,sizeof(float)*(offset/3+vertex));
	if(tmp){reflectance_=tmp;}
	else exit(EXIT_FAILURE);
	if(face>0){
		tmp_i=(unsigned int*)realloc(faces,sizeof(unsigned int)*(offseti+face*3));
		if(tmp_i){faces=tmp_i;}
		else exit(EXIT_FAILURE);
	}
	paranum--;
	for(int i=offset;i<offset+vertex*3;i++){norm_[i]=0;}
	
	if(headString(format,"binary_little_endian")){
		//vertex読み込み
		
		while(n<vertex&&!ifs.eof()){
			float f;
			float fa[10];
			for(time=0;time<paranum;time++){
			ifs.read((char *) &f,sizeof(float));
			fa[time]=f;
			}	
			vec[offset+n*3]=fa[xi];
			vec[offset+n*3+1]=fa[yi];
			vec[offset+n*3+2]=fa[zi];
			gx=gx*((double)(n)/(n+1))+vec[offset+n*3]*((double)(1)/(n+1));
			gy=gy*((double)(n)/(n+1))+vec[offset+n*3+1]*((double)(1)/(n+1));
			gz=gz*((double)(n)/(n+1))+vec[offset+n*3+2]*((double)(1)/(n+1));
			reflectance_[offset/3+n]=fa[reflecti];
			validVertex++;
			n++;
		}
		n=0;
		//norm読み込み
		while(n<face&&!ifs.eof()){
			int f;
			int fa[10];
			for(time=0;time<4;time++){
				if(time!=0){
					ifs.read((char *) &f,sizeof(int));
					fa[time-1]=f;
				}else ifs.read((char *) &f,sizeof(char));	
				
			}
			faces[offseti+n*3]=fa[0];
			faces[offseti+n*3+1]=fa[1];
			faces[offseti+n*3+2]=fa[2];

			POINT_3D p1,p2,p3;
			p1.x=vec[offset+fa[0]*3];			p1.y=vec[offset+fa[0]*3+1];			p1.z=vec[offset+fa[0]*3+2];
			p2.x=vec[offset+fa[1]*3];			p2.y=vec[offset+fa[1]*3+1];			p2.z=vec[offset+fa[1]*3+2];
			p3.x=vec[offset+fa[2]*3];			p3.y=vec[offset+fa[2]*3+1];			p3.z=vec[offset+fa[2]*3+2];

			POINT_3D nc=getNorm(p1,p2,p3);
			norm_[offset+fa[0]*3]+=nc.x;norm_[offset+fa[0]*3+1]+=nc.y;norm_[offset+fa[0]*3+2]+=nc.z;
			norm_[offset+fa[1]*3]+=nc.x;norm_[offset+fa[1]*3+1]+=nc.y;norm_[offset+fa[1]*3+2]+=nc.z;
			norm_[offset+fa[2]*3]+=nc.x;norm_[offset+fa[2]*3+1]+=nc.y;norm_[offset+fa[2]*3+2]+=nc.z;
			facesNum++;
			n++;
		}
	}
	if(headString(format,"ascii")){//永遠に未完成, Incompleting Eternally 
		while(n<vertex&&!ifs.eof()){
			float f[5];
				getline(ifs,line);
				int i;
				for(i=0;i<5;i++){
					if(i!=4)f[i]=stof(line.substr(0,line.find_first_of(" ")));
					else f[i]=stof(line);
					line.erase(0,line.find_first_of(" ")+1);
				}

			vec[offset+n*3]=f[xi];
			vec[offset+n*3+1]=f[yi];
			vec[offset+n*3+2]=f[zi];
			
			n++;
		}
	}
	ifs.close();
	}
	int allVnum=validVertex;
	validVertex=0;
	//normのチェック
//	for (int i=0;i<allVnum;i++){
//		if(norm_[i*3]!=0||norm_[i*3+1]!=0||norm_[i*3+2]!=0)validVertex++;
//	}

	norm=(float*)malloc(sizeof(float)*allVnum*3);
//	dataPts=annAllocPts(validVertex,3);
	reflectance=(float*)malloc(sizeof(float)*allVnum);
	verteces=(float*)malloc(sizeof(float)*allVnum*3);
	cout<<allVnum<<endl;
	int idx=0;
	for(int i=0;i<allVnum;i++){
		verteces[idx*3]=vec[i*3];
		verteces[idx*3+1]=vec[i*3+1];
		verteces[idx*3+2]=vec[i*3+2];
		norm[idx*3]=norm_[i*3];
		norm[idx*3+1]=norm_[i*3+1];
		norm[idx*3+2]=norm_[i*3+2];
		reflectance[idx]=reflectance_[i];
		idx++;
	}

	idx=0;
	free(vec);
	vertexnum=allVnum;
	facenum=facesNum;
	free(norm_);
	free(reflectance_);
		g<<gx,gy,gz;
	bRead=true;
	return true;
	

}


bool BasicPly::readPlyFileRGB(vector<string> fileName,int dataNum){
	
	int rf=0;
	float* vec=(float*)malloc(1);
	float* norm_=(float*)malloc(1);
	//float* reflectance_=(float*)malloc(1);
	char* rgba_=(char*)malloc(1);
	faces=(unsigned int*)malloc(1);
	int validVertex=0;
	int facesNum=0;
	for(rf=0;rf<dataNum;rf++){
	ifstream ifs(fileName[rf],ios::binary);
	string line;
	string format;
	long n=0;
	int xi=-1,yi=-1,zi=-1,reflecti=-1,redi=-1,greeni=-1,bluei=-1,alphai=3;
	int paranum=0;
	int paranumb=0;
	int vertex=0;
	int face=0;
	int matIdx=0;
	while(getline(ifs,line)){
		cout<<line<<endl;
		if(headString(line,"format")){
			format=line.erase(0,7);
		}
		if(headString(line,"end_header"))break;		
		if(headString(line,"property")){
			
			line.erase(0,9);
			if(headString(line,"float")){
				line.erase(0,6);
				if(headString(line,"x"))xi=paranum;
				else if(headString(line,"y"))yi=paranum;
				else if(headString(line,"z"))zi=paranum;
				else if(headString(line,"confidence"));
				else return false;
				paranum++;
			}
			else if(headString(line,"uchar")){
				line.erase(0,6);
				if(headString(line,"red"))redi=paranumb;
				else if(headString(line,"green"))greeni=paranumb;
				else if(headString(line,"blue"))bluei=paranumb;
				else if(headString(line,"alpha"));
				else return false;
				paranumb++;
			}
			
		}
		if(headString(line,"element")){
			line.erase(0,8);
			if(headString(line,"vertex")){
				line.erase(0,7);
				vertex=stoi(line);
			
			}
			if(headString(line,"face")){
				line.erase(0,5);
				face=stoi(line);			
			}
		}
		if(headString(line,"matrix")){
		float f[4];

		for(int i=0;i<5;i++){
			line.erase(0,line.find_first_not_of(" "));
			if(i>0)f[i-1]=stof(line.substr(0,line.find_first_of(" ")));
			line.erase(0,line.find_first_of(" "));
		}
			GlobalPose(0,matIdx)=f[0];
			GlobalPose(1,matIdx)=f[1];
			GlobalPose(2,matIdx)=f[2];
			GlobalPose(3,matIdx)=f[3];
			matIdx++;
		}
	}
	float a=0;
	int time=0;
	int num=0;
	int offset=validVertex*3;
	int offseti=facesNum*3;
	float* tmp;
	char* tmpc;
	unsigned int*tmp_i;
	tmp=(float*)realloc(norm_,sizeof(float)*(offset+vertex*3));
	if(tmp){norm_=tmp;}
	else exit(EXIT_FAILURE);
	tmp=(float*)realloc(vec,sizeof(float)*(offset+vertex*3));
	if(tmp){vec=tmp;}
	else exit(EXIT_FAILURE);

	tmpc=(char*)realloc(rgba_,sizeof(char)*(offset/3*4+vertex*4));
	if(tmpc){rgba_=tmpc;}
	else exit(EXIT_FAILURE);
	tmp_i=(unsigned int*)realloc(faces,sizeof(unsigned int)*(offseti+face*3));
	if(tmp_i){faces=tmp_i;}
	else exit(EXIT_FAILURE);
	//paranum--;
	for(int i=offset;i<offset+vertex*3;i++){norm_[i]=0;}
	
	if(headString(format,"binary_little_endian")){
		//vertex読み込み
		
		while(n<vertex&&!ifs.eof()){
			float f;
			float fa[10];
			for(time=0;time<paranum;time++){
			ifs.read((char *) &f,sizeof(float));
			fa[time]=f;
			}	
			char c;
			char ca[10];
			for(time=0;time<paranumb;time++){
			ifs.read((char *) &c,sizeof(char));
			ca[time]=c;
			}


			vec[offset+n*3]=fa[xi];
			vec[offset+n*3+1]=fa[yi];
			vec[offset+n*3+2]=fa[zi];
			rgba_[offset/3*4+n*4]=ca[redi];
			rgba_[offset/3*4+n*4+1]=ca[greeni];
			rgba_[offset/3*4+n*4+2]=ca[bluei];
			rgba_[offset/3*4+n*4+3]=ca[alphai];
			//reflectance_[offset/3+n]=fa[reflecti];
			validVertex++;
			n++;
		}
		n=0;
		//norm読み込み
		while(n<face&&!ifs.eof()){
			int f;
			int fa[10];
			for(time=0;time<4;time++){
				if(time!=0){
					ifs.read((char *) &f,sizeof(int));
					fa[time-1]=f;
				}else ifs.read((char *) &f,sizeof(char));	
				
			}
			if(fa[0]>=vertex||fa[1]>=vertex||fa[2]>=vertex)continue;
			faces[offseti+n*3]=fa[0]+offset/3;
			faces[offseti+n*3+1]=fa[1]+offset/3;
			faces[offseti+n*3+2]=fa[2]+offset/3;
			
			POINT_3D p1,p2,p3;
			p1.x=vec[offset+fa[0]*3];			p1.y=vec[offset+fa[0]*3+1];			p1.z=vec[offset+fa[0]*3+2];
			p2.x=vec[offset+fa[1]*3];			p2.y=vec[offset+fa[1]*3+1];			p2.z=vec[offset+fa[1]*3+2];
			p3.x=vec[offset+fa[2]*3];			p3.y=vec[offset+fa[2]*3+1];			p3.z=vec[offset+fa[2]*3+2];

			POINT_3D nc=getNorm(p1,p2,p3);
			norm_[offset+fa[0]*3]+=nc.x;norm_[offset+fa[0]*3+1]+=nc.y;norm_[offset+fa[0]*3+2]+=nc.z;
			norm_[offset+fa[1]*3]+=nc.x;norm_[offset+fa[1]*3+1]+=nc.y;norm_[offset+fa[1]*3+2]+=nc.z;
			norm_[offset+fa[2]*3]+=nc.x;norm_[offset+fa[2]*3+1]+=nc.y;norm_[offset+fa[2]*3+2]+=nc.z;
			facesNum++;
			n++;
		}
	}
	if(headString(format,"ascii")){//永遠に未完成, Incompleting Eternally 
		while(n<vertex&&!ifs.eof()){
			float f[5];
				getline(ifs,line);
				int i;
				for(i=0;i<5;i++){
					if(i!=4)f[i]=stof(line.substr(0,line.find_first_of(" ")));
					else f[i]=stof(line);
					line.erase(0,line.find_first_of(" ")+1);
				}

			vec[offset+n*3]=f[xi];
			vec[offset+n*3+1]=f[yi];
			vec[offset+n*3+2]=f[zi];
			
			n++;
		}
	}
	ifs.close();
	}
	int allVnum=validVertex;
	validVertex=0;
	norm=(float*)malloc(sizeof(float)*allVnum*3);
	rgba=(unsigned char*)malloc(sizeof(char)*allVnum*4);
	verteces=(float*)malloc(sizeof(float)*allVnum*3);
	cout<<allVnum<<endl;
	int idx=0;
	for(int i=0;i<allVnum;i++){
		verteces[idx*3]=vec[i*3];
		verteces[idx*3+1]=vec[i*3+1];
		verteces[idx*3+2]=vec[i*3+2];
		norm[idx*3]=norm_[i*3];
		norm[idx*3+1]=norm_[i*3+1];
		norm[idx*3+2]=norm_[i*3+2];
		rgba[idx*4]=rgba_[i*4];
		rgba[idx*4+1]=rgba_[i*4+1];
		rgba[idx*4+2]=rgba_[i*4+2];
		rgba[idx*4+3]=rgba_[i*4+3];
		idx++;
	}
	free(vec);
	vertexnum=allVnum;
	facenum=facesNum;
	free(norm_);
	free(rgba_);
	bRead=true;
	return true;
	

}

void BasicPly::writePlyFile(string fileName){
			ofstream ofs(fileName,ios::out|ios::binary);
	ofs<<"ply"<<endl;
	ofs<<"format binary_little_endian 1.0"<<endl;
	ofs<<"element vertex "<<vertexnum<<endl;
	ofs<<"property float x"<<endl;
	ofs<<"property float y"<<endl;
	ofs<<"property float z"<<endl;
	ofs<<"property float confidence"<<endl;
	ofs<<"property float intensity"<<endl;
	ofs<<"element face "<<facenum<<endl;
	ofs<<"property list uchar int vertex_index"<<endl;
	ofs<<"end_header"<<endl;
	int time;
	int i;
	for(i=0;i<vertexnum;i++){
	
		float fa[5];
		fa[0]=verteces[i*3];
		fa[1]=verteces[i*3+1];
		fa[2]=verteces[i*3+2];
		fa[3]=1.0f;
		fa[4]=reflectance[i];
		for(time=0;time<5;time++){
			ofs.write((char *) &fa[time],sizeof(float));	
		}
	}
	unsigned char ftri=3;
	for(i=0;i<facenum;i++){
		ofs.write((char *) &ftri,sizeof(char));
		ofs.write((char *) &faces[i*3],sizeof(int));
		ofs.write((char *) &faces[i*3+1],sizeof(int));
		ofs.write((char *) &faces[i*3+2],sizeof(int));
		
	}
	ofs.close();
}

void BasicPly::writePlyFileRGB(string fileName){
			ofstream ofs(fileName,ios::out|ios::binary);
	ofs<<"ply"<<endl;
	ofs<<"format binary_little_endian 1.0"<<endl;
	ofs<<"element vertex "<<vertexnum<<endl;
	ofs<<"property float x"<<endl;
	ofs<<"property float y"<<endl;
	ofs<<"property float z"<<endl;
	ofs<<"property float confidence"<<endl;
	ofs<<"property uchar red"<<endl;
	ofs<<"property uchar green"<<endl;
	ofs<<"property uchar blue"<<endl;
	ofs<<"property uchar alpha"<<endl;
	ofs<<"element face "<<facenum<<endl;
	ofs<<"property list uchar int vertex_index"<<endl;
	ofs<<"end_header"<<endl;
	int time;
	int i;
	for(i=0;i<vertexnum;i++){
	
		float fa[5];
		fa[0]=verteces[i*3];
		fa[1]=verteces[i*3+1];
		fa[2]=verteces[i*3+2];
		fa[3]=1.0f;
		for(time=0;time<4;time++){
			ofs.write((char *) &fa[time],sizeof(float));	
		}
		ofs.write((char *) &rgba[i*4],sizeof(char));	
		ofs.write((char *) &rgba[i*4+1],sizeof(char));		
		ofs.write((char *) &rgba[i*4+2],sizeof(char));
		ofs.write((char *) &rgba[i*4+3],sizeof(char));
	}
	unsigned char ftri=3;
	for(i=0;i<facenum;i++){
		ofs.write((char *) &ftri,sizeof(char));
		ofs.write((char *) &faces[i*3],sizeof(int));
		ofs.write((char *) &faces[i*3+1],sizeof(int));
		ofs.write((char *) &faces[i*3+2],sizeof(int));
		
	}
	ofs.close();
}

void BasicPly::writePlyFileRGBForMeshlab(string fileName) {
	ofstream ofs(fileName, ios::out | ios::binary);
	ofs << "ply" << endl;
	ofs << "format binary_little_endian 1.0" << endl;
	ofs << "element vertex " << vertexnum << endl;
	ofs << "property float x" << endl;
	ofs << "property float y" << endl;
	ofs << "property float z" << endl;
	ofs << "property uchar red" << endl;
	ofs << "property uchar green" << endl;
	ofs << "property uchar blue" << endl;
	ofs << "element face " << facenum << endl;
	ofs << "property list uchar int vertex_index" << endl;
	ofs << "end_header" << endl;
	int time;
	int i;
	for (i = 0;i<vertexnum;i++) {

		float fa[5];
		fa[0] = verteces[i * 3];
		fa[1] = verteces[i * 3 + 1];
		fa[2] = verteces[i * 3 + 2];
		//fa[3] = 1.0f;
		for (time = 0;time<3;time++) {
			ofs.write((char *)&fa[time], sizeof(float));
		}
		ofs.write((char *)&rgba[i * 4], sizeof(char));
		ofs.write((char *)&rgba[i * 4 + 1], sizeof(char));
		ofs.write((char *)&rgba[i * 4 + 2], sizeof(char));
	}
	unsigned char ftri = 3;
	for (i = 0;i<facenum;i++) {
		ofs.write((char *)&ftri, sizeof(char));
		ofs.write((char *)&faces[i * 3], sizeof(int));
		ofs.write((char *)&faces[i * 3 + 1], sizeof(int));
		ofs.write((char *)&faces[i * 3 + 2], sizeof(int));

	}
	ofs.close();
}

void BasicPly::panoramaTexture(unsigned char* rgbArray,int width,int height,Matrix4d& transMat){
	if(rgba==NULL){
		rgba=(unsigned char*)malloc(sizeof(unsigned char)*vertexnum*4);
	}
	for(int i=0;i<vertexnum;i++){
		Vector4d p;p<<verteces[3*i],verteces[3*i+1],verteces[3*i+2],1;
		
		p=transMat*p;
		float pitch,yaw;
		pitch=-(atan2(p(1),p(0))/M_PI)*180;
		if(pitch<0)pitch+=360;
		yaw=acos(p(2)/sqrt(p(0)*p(0)+p(1)*p(1)+p(2)*p(2)))/M_PI*180;
		if(pitch<360&&pitch>=0&&yaw>=0&&yaw<180){
			int x=(int)(width*pitch/360);
			int y=(int)(height*(yaw/180));
		
			int idx=y*width+x;
			rgba[i*4]=(unsigned char)rgbArray[idx*3+2];
			rgba[i*4+1]=(unsigned char)rgbArray[idx*3+1];
			rgba[i*4+2]=(unsigned char)rgbArray[idx*3];
			rgba[i*4+3]=255;
		}else{
			rgba[i*4]=0;
			rgba[i*4+1]=0;
			rgba[i*4+2]=0;
			rgba[i*4+3]=0;
		}
	}

}


void BasicPly::release(){
	if(!bRead)return;
	free(verteces);
	free(norm);
	free(reflectance);
	bRead=false;	

}

void BasicPly::convertFAROColoredPtx2Ply(string inFileName, double threshold) {
	//read ptx
	ifstream ifs(inFileName, ios::binary);
	string line_str;
	//unsigned int* lineStartIdx = (unsigned int*)malloc(sizeof(int) * 3000);
	//unsigned int* idxpointer = lineStartIdx;
	getline(ifs, line_str);
	int column = stoi(line_str);
	getline(ifs, line_str);
	int row = stoi(line_str);
	for (int i = 0;i<8;i++) {
		getline(ifs, line_str);
	}
	//vector<float> vert;


	int lineNum = 0, vertCnt = 0;

	//lineStartIdx[lineNum] = vert.size() / 3;
	float thresMax = 120 * 120;
	float thresMin = 2 * 2;
	POINT_3D* arrayPoint = (POINT_3D*)malloc(sizeof(POINT_3D)*row*column);
	float* reflectanceArray = (float*)malloc(sizeof(float)*row*column);
	unsigned char* colorArray = (unsigned char*)malloc(sizeof(unsigned char) * 4 * row*column);
	cout << "row,column:(" << row << "," << column << ")" << endl;

	cout << "read data" << endl;
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
		reflectanceArray[vertCnt] = f[3];
		colorArray[vertCnt*4+0] = rgb[0];
		colorArray[vertCnt * 4 + 1] = rgb[1];
		colorArray[vertCnt * 4 + 2] = rgb[2];
		colorArray[vertCnt * 4 + 3] = 255;
		vertCnt++;
	}
	//vector<float> refVect;
	//vector<unsigned int>vertIdx;
	//for (int i = 0;i<line;i++) {
	//	lineStartIdx[i] = vert.size() / 3;
	//	cout << "line:" << i << " start:" << lineStartIdx[i] << endl;
	//	for (int j = 0;j<oneLine;j++) {
	//		int idx = i*oneLine + j;
	//		float r2 = arrayPoint[idx].x*arrayPoint[idx].x + arrayPoint[idx].y*arrayPoint[idx].y + arrayPoint[idx].z*arrayPoint[idx].z;
	//		if (r2>thresMax || r2<thresMin || reflectanceArray[idx]<0.01)continue;
	//		vert.push_back(arrayPoint[idx].x);
	//		vert.push_back(arrayPoint[idx].y);
	//		vert.push_back(arrayPoint[idx].z);
	//		vertIdx.push_back(idx);//多分ここバグ
	//		refVect.push_back(reflectanceArray[idx]);
	//	}
	//}
	//int vertexNum = vert.size() / 3;
	//float* vertex = (float*)malloc(sizeof(float)*vert.size());
	//for (int i = 0;i<vert.size();i++) {
	//	vertex[i] = vert.at(i);
	//}
	//unsigned int* vertexIdx = (unsigned int*)malloc(sizeof(unsigned int)*vertIdx.size());

	//vert.clear();
	//free(reflectanceArray);
	//free(arrayPoint);

	int			_column = column;
	int			_mcolumn = _column - 1;
	int			_row = row;
	int			_mrow = _row - 1;

	int			_avertex = _column * _row;

	ply_index*	_indices = new ply_index[_mcolumn * _mrow * 2];
	int			_idx = 0;

	int*		_flags = new int[_avertex];
	::memset((char*)_flags, 0, sizeof(int) * _avertex);

	float		_lenc0, _lenc1;
	float		_len0, _len1, _len2, _len3;
	float		_ref0, _ref1, _ref2, _ref3;
	//pts_vertex*	_ptsv = pts.vertex_;

	int			i;

	float		_threshold = threshold;
	
	_threshold *= _threshold;
	float reflectanceMin = 0.01;
	float reflectanceMax = 1.0;
	cout << "meshing" << endl;
	for (int i = 0; i < _mcolumn; ++i) {
		int		_cindex = i * _row;
		for (int j = 0; j < _mrow; ++j) {
			//	Calculate length of vertices.
			_len0 = getDist2(arrayPoint[_cindex], arrayPoint[_cindex + 1]);
			_len1 = getDist2(arrayPoint[_cindex], arrayPoint[_cindex + _row]);
			_len2 = getDist2(arrayPoint[_cindex + _row], arrayPoint[_cindex + _row + 1]);
			_len3 = getDist2(arrayPoint[_cindex + _row + 1], arrayPoint[_cindex + 1]);

			_ref0 = reflectanceArray[_cindex];
			_ref1 = reflectanceArray[_cindex + 1];
			_ref2 = reflectanceArray[_cindex + _row];
			_ref3 = reflectanceArray[_cindex + _row + 1];


			if (_len0 < _threshold && _len0) {
				if (_len1 < _threshold && _len1) {
					if (reflectanceMin <= _ref0 && reflectanceMax >= _ref0 &&
						reflectanceMin <= _ref1 && reflectanceMax >= _ref1 &&
						reflectanceMin <= _ref2 && reflectanceMax >= _ref2) {
						_indices[_idx].nindex_ = 3;
						_indices[_idx].index0_ = _cindex;
						_indices[_idx].index2_ = _cindex + _row;
						_indices[_idx].index1_ = _cindex + 1;
						++_flags[_cindex];
						++_flags[_cindex + _row];
						++_flags[_cindex + 1];
						++_idx;
					}
					if (_len2 < _threshold && _len2) {
						if (_len3 < _threshold && _len3) {
							if (reflectanceMin <= _ref0 && reflectanceMax >= _ref1 &&
								reflectanceMin <= _ref1 && reflectanceMax >= _ref2 &&
								reflectanceMin <= _ref2 && reflectanceMax >= _ref3) {
								_indices[_idx].nindex_ = 3;
								_indices[_idx].index0_ = _cindex + 1;
								_indices[_idx].index2_ = _cindex + _row;
								_indices[_idx].index1_ = _cindex + _row + 1;
								++_flags[_cindex + 1];
								++_flags[_cindex + _row];
								++_flags[_cindex + _row + 1];
								++_idx;
							}
						}
					}
				}
			}

			else if (_len1 < _threshold && _len1) {
				if (_len2 < _threshold && _len2) {
					if (reflectanceMin <= _ref0 && reflectanceMax >= _ref0 &&
						reflectanceMin <= _ref1 && reflectanceMax >= _ref1 &&
						reflectanceMin <= _ref2 && reflectanceMax >= _ref2) {

						_indices[_idx].nindex_ = 3;
						_indices[_idx].index0_ = _cindex;
						_indices[_idx].index2_ = _cindex + _row;
						_indices[_idx].index1_ = _cindex + _row + 1;
						++_flags[_cindex];
						++_flags[_cindex + _row];
						++_flags[_cindex + _row + 1];
						++_idx;
					}
				}
			}
			else if (_len2 < _threshold && _len2) {
				if (_len3 < _threshold && _len3) {
					if (reflectanceMin <= _ref1 && reflectanceMax >= _ref1 &&
						reflectanceMin <= _ref2 && reflectanceMax >= _ref2 &&
						reflectanceMin <= _ref3 && reflectanceMax >= _ref3) {
						_indices[_idx].nindex_ = 3;
						_indices[_idx].index0_ = _cindex + 1;
						_indices[_idx].index2_ = _cindex + _row;
						_indices[_idx].index1_ = _cindex + _row + 1;
						++_flags[_cindex + 1];
						++_flags[_cindex + _row];
						++_flags[_cindex + _row + 1];
						++_idx;
					}
				}
			}
			++_cindex;
		}
	}
	

	int		_rvnum = 0;
	for (i = 0; i < _avertex; ++i) {
		if (_flags[i]) {
			_flags[i] = _rvnum++;
		}
		else
			_flags[i] = -1;
	}

	verteces = (float*)malloc(sizeof(float)*_rvnum*3);
	reflectance = (float*)malloc(sizeof(float)*_rvnum);
	rgba = (unsigned char*)malloc(sizeof(unsigned char)*_rvnum*4);
	::memset(reinterpret_cast<char*>(verteces), 0, sizeof(float) * _rvnum*3);
	::memset(reinterpret_cast<char*>(reflectance), 0, sizeof(float) * _rvnum);
	::memset(reinterpret_cast<char*>(rgba), 0, sizeof(unsigned char) * _rvnum * 4);
	int			_vnum = 0;
	for (i = 0; i < _avertex; ++i) {
		if (_flags[i] >= 0) {
			(POINT_3D&)(verteces[_vnum * 3]) = arrayPoint[i];
			reflectance[_vnum] = reflectanceArray[i];
			rgba[_vnum * 4 + 0] = colorArray[i*4+0];
			rgba[_vnum * 4 + 1] = colorArray[i * 4 + 1];
			rgba[_vnum * 4 + 2] = colorArray[i * 4 + 2];
			rgba[_vnum * 4 + 3] = colorArray[i * 4 + 3];
			_vnum++;
		}
	}
	vertexnum = _rvnum;
	faces = (unsigned int*)malloc(sizeof(unsigned int) * 3 * _idx);
	for (i = 0; i < _idx; ++i) {
		ply_index&		_oindex = _indices[i];
		
		faces[i*3] = _flags[_oindex.index2_];
		faces[i * 3 + 1] = _flags[_oindex.index1_];
		faces[i * 3 + 2] = _flags[_oindex.index0_];
	}

	facenum = _idx;

	free(arrayPoint);
	free(reflectanceArray);
	free(colorArray);
	delete[]	_flags;
	delete[]	_indices;

	return;
}

//向きが逆な気がしてきたぞ・・・
POINT_3D getNorm(POINT_3D p1,POINT_3D p2,POINT_3D p3){
	double v1x=p2.x-p1.x,v1y=p2.y-p1.y,v1z=p2.z-p1.z,
		v2x=p3.x-p2.x,v2y=p3.y-p2.y,v2z=p3.z-p2.z;
	Vector3d nc_;
	POINT_3D nc;
	nc_(0)=v1y*v2z-v1z*v2y;
	nc_(1)=v1z*v2x-v1x*v2z;
	nc_(2)=v1x*v2y-v1y*v2x;
	nc_.normalize();
	nc.x=nc_(0);
	nc.y=nc_(1);
	nc.z=nc_(2);
	return nc;
};

