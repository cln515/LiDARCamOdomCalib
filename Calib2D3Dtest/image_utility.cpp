#include "image_utility.h";

void getGraySubPixel_uchar(cv::Mat image,cv::Point2f p,double *ret){
	int xl=(int)(p.x-0.5);
	int yt=(int)(p.y-0.5);
	int xr=xl+1;
	int yb=yt+1;

	if(xl<0)xl=0;
	if(yt<0)yt=0;
	if(xr>=image.cols)xr=image.cols-1;
	if(yb>=image.rows)yb=image.rows-1;
	uchar* imArray=image.data;
	double dx=(p.x-0.5)-xl;
	double dy=(p.y-0.5)-yt;

	ret[0]=(1-dx)*(1-dy);
	ret[1]=imArray[xl+yt*image.cols]/256.0;
	ret[2]=(dx)*(1-dy);
	ret[3]=imArray[xr+yt*image.cols]/256.0;
	ret[4]=(1-dx)*(dy);
	ret[5]=imArray[xl+yb*image.cols]/256.0;
	ret[6]=(dx)*(dy);
	ret[7]=imArray[xr+yb*image.cols]/256.0;
}

void getGraySubPixel(cv::Mat image,cv::Point2f p,double *ret){
	int xl=(int)(p.x-0.5);
	int yt=(int)(p.y-0.5);
	int xr=xl+1;
	int yb=yt+1;

	if(xl<0)xl=0;
	if(yt<0)yt=0;
	if(xr>=image.cols)xr=image.cols-1;
	if(yb>=image.rows)yb=image.rows-1;
	uchar* imArray=image.data;
	double dx=(p.x-0.5)-xl;
	double dy=(p.y-0.5)-yt;

	ret[0]=(1-dx)*(1-dy);
	ret[1]=imArray[xl+yt*image.cols]/256.0;
	ret[2]=(dx)*(1-dy);
	ret[3]=imArray[xr+yt*image.cols]/256.0;
	ret[4]=(1-dx)*(dy);
	ret[5]=imArray[xl+yb*image.cols]/256.0;
	ret[6]=(dx)*(dy);
	ret[7]=imArray[xr+yb*image.cols]/256.0;
}

void getGraySubPixel_float(cv::Mat image,cv::Point2f p,double *ret){
	int xl=(int)(p.x-0.5);
	int yt=(int)(p.y-0.5);
	int xr=xl+1;
	int yb=yt+1;

	if(xl<0)xl=0;
	if(yt<0)yt=0;
	if(xr>=image.cols)xr=image.cols-1;
	if(yb>=image.rows)yb=image.rows-1;
	float* imArray=(float*)image.data;
	double dx=(p.x-0.5)-xl;
	double dy=(p.y-0.5)-yt;

	ret[0]=(1-dx)*(1-dy);
	ret[1]=imArray[xl+yt*image.cols]/256.0;
	ret[2]=(dx)*(1-dy);
	ret[3]=imArray[xr+yt*image.cols]/256.0;
	ret[4]=(1-dx)*(dy);
	ret[5]=imArray[xl+yb*image.cols]/256.0;
	ret[6]=(dx)*(dy);
	ret[7]=imArray[xr+yb*image.cols]/256.0;
}

void getSubPixel_float(cv::Mat image,cv::Point2f p,double *ret){
	int xl=(int)(p.x-0.5);
	int yt=(int)(p.y-0.5);
	int xr=xl+1;
	int yb=yt+1;

	if(xl<0)xl=0;
	if(yt<0)yt=0;
	if(xr>=image.cols)xr=image.cols-1;
	if(yb>=image.rows)yb=image.rows-1;
	float* imArray=(float*)image.data;
	double dx=(p.x-0.5)-xl;
	double dy=(p.y-0.5)-yt;

	ret[0]=(1-dx)*(1-dy);
	ret[1]=imArray[xl+yt*image.cols];
	ret[2]=(dx)*(1-dy);
	ret[3]=imArray[xr+yt*image.cols];
	ret[4]=(1-dx)*(dy);
	ret[5]=imArray[xl+yb*image.cols];
	ret[6]=(dx)*(dy);
	ret[7]=imArray[xr+yb*image.cols];
}

void getGraySubPixel(cv::Mat image,cv::Point2f p,double *ret, double *dret){
	int xl=(int)(p.x-0.5);
	int yt=(int)(p.y-0.5);
	int xr=xl+1;
	int yb=yt+1;

	if(xl<0)xl=0;
	if(yt<0)yt=0;
	if(xr>=image.cols)xr=image.cols-1;
	if(yb>=image.rows)yb=image.rows-1;
	uchar* imArray=image.data;
	double dx=(p.x-0.5)-xl;
	double dy=(p.y-0.5)-yt;

	ret[0]=(1-dx)*(1-dy);
	ret[1]=imArray[xl+yt*image.cols]/256.0;
	ret[2]=(dx)*(1-dy);
	ret[3]=imArray[xr+yt*image.cols]/256.0;
	ret[4]=(1-dx)*(dy);
	ret[5]=imArray[xl+yb*image.cols]/256.0;
	ret[6]=(dx)*(dy);
	ret[7]=imArray[xr+yb*image.cols]/256.0;

	dret[0]=dx;
	dret[1]=dy;
}

void getColorSubPixel(cv::Mat image, cv::Point2f p, uchar *ret) {
	int xl = (int)(p.x - 0.5);
	int yt = (int)(p.y - 0.5);
	int xr = xl + 1;
	int yb = yt + 1;

	if (xl<0)xl = 0;
	if (yt<0)yt = 0;
	if (xr >= image.cols)xr = image.cols - 1;
	if (yb >= image.rows)yb = image.rows - 1;
	uchar* imArray = image.data;//3 channel
	double dx = (p.x - 0.5) - xl;
	double dy = (p.y - 0.5) - yt;


	double rgb[3];

	rgb[0] = (1 - dx)*(1 - dy)*imArray[(xl + yt*image.cols) * 3]
		+ (dx)*(1 - dy) * imArray[(xr + yt*image.cols) * 3]
		+ (1 - dx)*(dy)* imArray[(xl + yb*image.cols) * 3]
		+ (dx)*(dy)* imArray[(xr + yb*image.cols) * 3];
	rgb[1] = (1 - dx)*(1 - dy)*imArray[(xl + yt*image.cols) * 3 + 1]
		+ (dx)*(1 - dy) * imArray[(xr + yt*image.cols) * 3 + 1]
		+ (1 - dx)*(dy)* imArray[(xl + yb*image.cols) * 3 + 1]
		+ (dx)*(dy)* imArray[(xr + yb*image.cols) * 3 + 1];
	rgb[2] = (1 - dx)*(1 - dy)*imArray[(xl + yt*image.cols) * 3 + 2]
		+ (dx)*(1 - dy) * imArray[(xr + yt*image.cols) * 3 + 2]
		+ (1 - dx)*(dy)* imArray[(xl + yb*image.cols) * 3 + 2]
		+ (dx)*(dy)* imArray[(xr + yb*image.cols) * 3 + 2];
	ret[0] = (uchar)rgb[0];
	ret[1] = (uchar)rgb[1];
	ret[2] = (uchar)rgb[2];
}


void panoramaRectification(cv::Mat image1, cv::Mat image2, cv::Mat& dstimage1, cv::Mat& dstimage2, Vector3d epi_, Matrix3d R, Matrix3d& R0, Matrix3d& R1) {
	Vector3d axz, ep, axr;
	axz << 0, 0, 1;
	ep = epi_.normalized();
	axr = axz.cross(ep);
	double angle = -asin(axr.norm());
	if (axz.dot(ep) < 0)angle = M_PI - angle;

	axr = axr.normalized();



	Matrix3d raxis;
	raxis << 0, -axr(2), axr(1),
		axr(2), 0, -axr(0),
		-axr(1), axr(0), 0;



	R0 = Matrix3d::Identity() + sin(angle)*raxis + (1 - cos(angle))*raxis*raxis;

	/*R0 << -cos(phie), 0, sin(phie),
	cos(thetae)*sin(phie), sin(thetae), cos(thetae)*cos(phie),
	-sin(thetae)*sin(phie), cos(thetae), -sin(thetae)*cos(phie);*/

	R1 = R0*R;
	cout << R0 << endl;
	cout << R0*ep << endl;
	//image resynthesis test

	cv::Mat resynthimg1(image1.rows, image1.cols, image1.type());
	cv::Mat resynthimg2(image2.rows, image2.cols, image2.type());
	Matrix3d R0inv, R1inv;
	R0inv = R0.inverse();
	R1inv = R1.inverse();


	for (int ix = 0;ix < resynthimg1.cols;ix++) {
		for (int iy = 0;iy < resynthimg1.rows;iy++) {
			Vector3d ret, transed;
			rev_omniTrans(ix, iy, resynthimg1.cols, resynthimg1.rows, ret);
			transed = R0inv*ret;
			double nix, niy;
			omniTrans(transed(0), transed(1), transed(2), niy, nix, resynthimg1.rows);
			cv::Point2f p(nix, niy);
			uchar rgb[3];
			if (image1.type() == CV_8UC1) {
				if (nix<0 || nix>image1.cols || niy<0 || niy>image1.rows) {
					rgb[0] = 0;
				}
				else {
					double ret[8];
					getGraySubPixel_uchar(image1, p, ret);
					rgb[0]=(uchar)((ret[0] * ret[1] + ret[2] * ret[3] + ret[4] * ret[5] + ret[6] * ret[7])*256);
				}
				resynthimg1.data[(ix + iy*resynthimg1.cols)] = rgb[0];
			}
			else {
				if (nix<0 || nix>image1.cols || niy<0 || niy>image1.rows) {
					rgb[0] = 0;rgb[1] = 0;rgb[2] = 0;
				}
				else {
					getColorSubPixel(image1, p, rgb);

				}
				resynthimg1.data[(ix + iy*resynthimg1.cols) * 3] = rgb[0];
				resynthimg1.data[(ix + iy*resynthimg1.cols) * 3 + 1] = rgb[1];
				resynthimg1.data[(ix + iy*resynthimg1.cols) * 3 + 2] = rgb[2];
			}
			
			
		}
	}
	for (int ix = 0;ix < resynthimg2.cols;ix++) {
		for (int iy = 0;iy < resynthimg2.rows;iy++) {
			Vector3d ret, transed;
			rev_omniTrans(ix, iy, resynthimg2.cols, resynthimg2.rows, ret);
			transed = R1inv*ret;
			double nix, niy;
			omniTrans(transed(0), transed(1), transed(2), niy, nix, resynthimg2.rows);
			cv::Point2f p(nix, niy);
			uchar rgb[3];
			if (image2.type() == CV_8UC1) {
				if (nix<0 || nix>image1.cols || niy<0 || niy>image1.rows) {
					rgb[0] = 0;
				}
				else {
					double ret[8];
					getGraySubPixel_uchar(image2, p, ret);
					rgb[0] = (uchar)((ret[0] * ret[1] + ret[2] * ret[3] + ret[4] * ret[5] + ret[6] * ret[7]) * 256);
				}
				resynthimg2.data[(ix + iy*resynthimg2.cols)] = rgb[0];
			}else{
				if (nix<0 || nix>image2.cols || niy<0 || niy>image2.rows) {
					rgb[0] = 0;rgb[1] = 0;rgb[2] = 0;
				}
				else {
					getColorSubPixel(image2, p, rgb);
				}
				resynthimg2.data[(ix + iy*resynthimg2.cols) * 3] = rgb[0];
				resynthimg2.data[(ix + iy*resynthimg2.cols) * 3 + 1] = rgb[1];
				resynthimg2.data[(ix + iy*resynthimg2.cols) * 3 + 2] = rgb[2];
			}
		}
	}
	dstimage1 = resynthimg1;
	dstimage2 = resynthimg2;

}

namespace cv
{

	template<typename T> struct greaterThanPtr
	{
		bool operator()(const T* a, const T* b) const { return *a > *b; }
	};

}

void goodFeatureToTrack_onProjection(cv::Mat image, vector<cv::Point2f> proj_points, vector<int>& selectedIdx, double minDistance, int maxCorners, double cornerThres) {
	cv::Mat eig, tmp;
	cv::cornerHarris(image, eig, 3, 5, 0.04);//corner calculation
	cv::Mat tmpMask = cv::Mat::zeros(image.rows, image.cols, CV_8U);
	map<unsigned int, int> point_idx;
	map<unsigned int, double> point_dist;
	for (int j = 0;j < proj_points.size();j++) {
		cv::Point2f p = proj_points.at(j);
		unsigned int pint = ((int)p.x) + ((int)p.y)*image.cols;;
		tmpMask.at<uchar>(p) = 255;//remove mask on projected point
		double cendx = p.x - (int)p.x - 0.5;
		double cendy = p.y - (int)p.y - 0.5;
		double distcenter = cendx*cendx + cendy*cendy;
		auto itr = point_dist.find(pint);
		if (itr == point_dist.end() || itr->second > distcenter) {
			point_dist.insert(map<unsigned int, double>::value_type(pint, distcenter));
			point_idx.insert(map<unsigned int, int>::value_type(pint, j));
		}
	}

	cv::dilate(tmpMask, tmpMask, cv::Mat());//expansion non-masked area
	double maxVal;
	//non-maxima suppression
	cv::minMaxLoc(eig, 0, &maxVal, 0, 0, tmpMask);
	//cv::threshold(eig, eig, maxVal*1e-12, 0, cv::THRESH_TOZERO);
	cv::dilate(eig, tmp, cv::Mat());

	cv::Size imgsize = image.size();

	vector<const float*> tmpCorners;

	// collect list of pointers to features - put them into temporary image
	for (int y = 1; y < imgsize.height - 1; y++)
	{
		const float* eig_data = (const float*)eig.ptr(y);
		const float* tmp_data = (const float*)tmp.ptr(y);
		const uchar* mask_data = tmpMask.data ? tmpMask.ptr(y) : 0;

		for (int x = 1; x < imgsize.width - 1; x++)
		{
			float val = eig_data[x];
			if (val > cornerThres /*&& val == tmp_data[x]*/ && (!mask_data || mask_data[x]))
				tmpCorners.push_back(eig_data + x);
		}
	}
	//cv::sort(tmpCorners,tmpCorners,CV_SORT_DESCENDING);
	std::sort(tmpCorners.begin(), tmpCorners.end(), [](const float*& a, const float*& b) {return (*a) >= (*b);});

	vector<cv::Point2f> corners;
	size_t i, j, total = tmpCorners.size(), ncorners = 0;
	if (minDistance >= 1)
	{
		// Partition the image into larger grids
		int w = image.cols;
		int h = image.rows;

		const int cell_size = cvRound(minDistance);
		const int grid_width = (w + cell_size - 1) / cell_size;
		const int grid_height = (h + cell_size - 1) / cell_size;

		std::vector<std::vector<cv::Point2f> > grid(grid_width*grid_height);

		minDistance *= minDistance;

		for (i = 0; i < total; i++)
		{
			int ofs = (int)((const uchar*)tmpCorners[i] - eig.data);
			int y = (int)(ofs / eig.step);
			int x = (int)((ofs - y*eig.step) / sizeof(float));

			bool good = true;

			int x_cell = x / cell_size;
			int y_cell = y / cell_size;

			int x1 = x_cell - 1;
			int y1 = y_cell - 1;
			int x2 = x_cell + 1;
			int y2 = y_cell + 1;

			// boundary check
			x1 = max(0, x1);
			y1 = max(0, y1);
			x2 = min(grid_width - 1, x2);
			y2 = min(grid_height - 1, y2);

			for (int yy = y1; yy <= y2; yy++)
			{
				for (int xx = x1; xx <= x2; xx++)
				{
					vector <cv::Point2f> &m = grid[yy*grid_width + xx];

					if (m.size())
					{
						for (j = 0; j < m.size(); j++)
						{
							float dx = x - m[j].x;
							float dy = y - m[j].y;

							if (dx*dx + dy*dy < minDistance)
							{
								good = false;
								goto break_out;
							}
						}
					}
				}
			}

		break_out:

			if (good)
			{
				// printf("%d: %d %d -> %d %d, %d, %d -- %d %d %d %d, %d %d, c=%d\n",
				//    i,x, y, x_cell, y_cell, (int)minDistance, cell_size,x1,y1,x2,y2, grid_width,grid_height,c);
				grid[y_cell*grid_width + x_cell].push_back(cv::Point2f((float)x, (float)y));

				corners.push_back(cv::Point2f((float)x, (float)y));
				++ncorners;

				if (maxCorners > 0 && (int)ncorners == maxCorners)
					break;
			}
		}
	}
	else
	{
		for (i = 0; i < total; i++)
		{
			int ofs = (int)((const uchar*)tmpCorners[i] - eig.data);
			int y = (int)(ofs / eig.step);
			int x = (int)((ofs - y*eig.step) / sizeof(float));

			corners.push_back(cv::Point2f((float)x, (float)y));
			++ncorners;
			if (maxCorners > 0 && (int)ncorners == maxCorners)
				break;
		}
	}
	selectedIdx.clear();
	for (int j = 0;j<corners.size();j++) {
		cv::Point2f p = corners.at(j);
		unsigned int pint = ((int)p.x) + ((int)p.y)*image.cols;
		auto itr2 = point_idx.find(pint);
		if (itr2 != point_idx.end()) {
			selectedIdx.push_back(point_idx.at(pint));
		}
	}


}