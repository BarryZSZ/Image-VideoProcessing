#include "head\ClassProjectHead.h"

void CourseProject1Main() {
	Mat I1 = imread("Images\\lena.pgm", 0);
	Mat I2 = imread("Images\\camera.pgm", 0);
	Mat I3 = imread("Images\\bridge.pgm", 0);

	//Question 1
	Mat I4 = LogTrans(I1, 1);
	imshow("Log Transformation", I4);
	I4 = GammaCorrection(I1, 0.25, 1);
	imshow("Gama = 0.25", I4);
	I4 = GammaCorrection(I1, 0.5, 1);
	imshow("Gama = 0.5", I4);
	I4 = GammaCorrection(I1, 1, 1);
	imshow("Gama = 1", I4);
	I4 = GammaCorrection(I1, 1.5, 1);
	imshow("Gama = 1.5", I4);
	I4 = GammaCorrection(I1, 2, 1);
	imshow("Gama = 2", I4);
	
	//Question 2
	Q2();

	//Question 3
	Mat Averagemask(3, 3, CV_8U, Scalar::all(1));
	Mat Gaussianmask = (Mat_<uchar>(3, 3) << 1, 2, 1, 2, 4, 2, 1, 2, 1);
	Mat I5 = SpatialFilter(I1, Averagemask);
	imshow("Lena Average Filter", I5);
	I5 = SpatialFilter(I2, Averagemask);
	imshow("Camera Average Filter", I5);
	I5 = SpatialFilter(I3, Averagemask);
	imshow("Bridge Average Filter", I5);
	Mat I6 = SpatialFilter(I1, Averagemask);
	imshow("Lena Gaussian Filter", I6);
	I6 = SpatialFilter(I2, Averagemask);
	imshow("Camera Gaussian Filter", I6);
	I6 = SpatialFilter(I3, Averagemask);
	imshow("Bridge Gaussian Filter", I6);

	//Question 4
	Mat histI;
	MatND hist;
	I5 = AddRandNoise(I1, 0.1, 0, 255);
	imshow("Lena Added Noise", I5);
	hist = Hist(I5);
	histI = HistImg(hist, 300, 500);
	imshow("Lena Add Noise Histogram", histI);
	I6 = MedianFilter(I5, 3, 3);
	imshow("Lena 3*3 Median filted", I6);
	hist = Hist(I6);
	histI = HistImg(hist, 300, 500);
	imshow("Lena 3*3 Median Histogram", histI);
	I6 = MedianFilter(I5, 5, 5);
	imshow("Lena 5*5 Median filted", I6);
	hist = Hist(I6);
	histI = HistImg(hist, 300, 500);
	imshow("Lena 5*5 Median Histogram", histI);

	I5 = AddRandNoise(I2, 0.1, 0, 255);
	imshow("Camera Added Noise", I5);
	hist = Hist(I5);
	histI = HistImg(hist, 300, 500);
	imshow("Camera Add Noise Histogram", histI);
	I6 = MedianFilter(I5, 3, 3);
	imshow("Camera 3*3 Median filted", I6);
	hist = Hist(I6);
	histI = HistImg(hist, 300, 500);
	imshow("Camera 3*3 Median Histogram", histI);
	I6 = MedianFilter(I5, 3, 3);
	imshow("Camera 5*5 Median filted", I6);
	hist = Hist(I6);
	histI = HistImg(hist, 300, 500);
	imshow("Camera 5*5 Median Histogram", histI);

	I5 = AddRandNoise(I3, 0.1, 0, 255);
	imshow("Bridge Added Noise", I5);
	hist = Hist(I5);
	histI = HistImg(hist, 300, 500);
	imshow("Bridge Add Noise Histogram", histI);
	I6 = MedianFilter(I5, 3, 3);
	imshow("Bridge 3*3 Median filted", I6);
	hist = Hist(I6);
	histI = HistImg(hist, 300, 500);
	imshow("Bridge 3*3 Median Histogram", histI);
	I6 = MedianFilter(I5, 3, 3);
	imshow("Bridge 5*5 Median filted", I6);
	hist = Hist(I6);
	histI = HistImg(hist, 300, 500);
	imshow("Bridge 5*5 Median Histogram", histI);

	//Question 5
	Mat I7 = imread("Images\\Plant.tif", 0);
	imshow("Plant", I7);
	Mat unsharpmask = (Mat_<uchar>(3, 3) << 1, 2, 1, 2, 4, 2, 1, 2, 1);
	Mat Masked = SpatialFilter(I7, unsharpmask);
	I7.convertTo(I7, CV_32F);
	Mat Mask = I7 - Masked;	
	Mat I8 = I7 + Mask;
	I8.convertTo(I8, CV_8U);
	Mat I9 = I7 + 2 * Mask;
	I9.convertTo(I9, CV_8U);
	Mask.convertTo(Mask, CV_8U);
	Masked.convertTo(Masked, CV_8U);
	normalize(Mask, Mask, 0, 255, NORM_MINMAX);
	imshow("Mask", Mask);
	imshow("Blurred Image", Masked);
	imshow("k = 1 unshark mask", I8);
	imshow("k = 2 unshark mask", I9);

	//Question 6
	Mat I10, I11, I12;
	I10 = IDHP(I1, 40);	
	imshow("Lena IDHP", I10);
	I10 = IDHP(I2, 40);
	imshow("Camera IDHP", I10);
	I10 = IDHP(I3, 40);
	imshow("Bridge IDFT", I10);

	I11 = BHP(I1, 40, 1);
	imshow("Lena BHP", I11);
	I11 = BHP(I2, 40, 1);
	imshow("Camera BHP", I11);
	I11 = BHP(I3, 40, 1);
	imshow("Bridge BFT", I11);

	I12 = GHP(I1, 40);
	imshow("Lena GHP", I12);
	I12 = GHP(I2, 40);
	imshow("Camera GHP", I12);
	I12 = IDHP(I3, 40);
	imshow("Bridge GFT", I12);



	waitKey(0);
}

/*@brief Image log transformation.*/
Mat LogTrans(Mat srcImage, float c) {
	srcImage.convertTo(srcImage, CV_32F);
	double minv = 0, maxv = 0;
	double *minp = &minv, *maxp = &maxv;
	minMaxIdx(srcImage, minp, maxp);
	float *psrcData;
	for (int i = 0; i < srcImage.rows; i++) {
		psrcData = srcImage.ptr<float>(i);
		for (int j = 0; j < srcImage.cols; j++) {
			psrcData[j] = c*log(psrcData[j] + 1);
		}
	}
	normalize(srcImage, srcImage, minv, maxv, NORM_MINMAX);
	srcImage.convertTo(srcImage, CV_8U);
	return srcImage;
}

/*@brief Image Gamma Correction.*/
Mat GammaCorrection(Mat srcImage, float Gama, float c) {
	srcImage.convertTo(srcImage, CV_32F);
	double minv = 0, maxv = 0;
	double *minp = &minv, *maxp = &maxv;
	minMaxIdx(srcImage, minp, maxp);
	float *psrcData;
	for (int i = 0; i < srcImage.rows; i++) {
		psrcData = srcImage.ptr<float>(i);
		for (int j = 0; j < srcImage.cols; j++) {
			psrcData[j] = c*pow(psrcData[j], Gama);
		}
	}
	normalize(srcImage, srcImage, minv, maxv, NORM_MINMAX);
	srcImage.convertTo(srcImage, CV_8U);
	return srcImage;
}

/*@brief Get the histogram of the srcImage.*/
MatND Hist(Mat srcImage) {
	MatND hist;
	const int channels[1] = { 0 };
	const int histSize[1] = { 256 };
	float pranges[2] = { 0,255 };
	const float* ranges[1] = { pranges };
	calcHist(&srcImage, 1, channels, Mat(), hist, 1, histSize, ranges);
	return hist;
}

/*@brief Generate the histogram Image of the srcImage.
@param srcImage is the source Image.
@param hist_w is the rows of the histogram.
@param hist_h is the cols of the histogran.*/
Mat HistImg(MatND hist, int hist_w, int hist_h) {
	/*MatND hist;
	const int channels[1] = { 0 };
	const int histSize[1] = { 256 };
	float pranges[2] = { 0,255 };
	const float* ranges[1] = { pranges };
	calcHist(&srcImage, 1, channels, Mat(), hist, 1, histSize, ranges);*/
	int zoom_size = round((float)hist_h / 256.0);
	hist_h = zoom_size * 256;
	int hist_rows = hist.rows;
	int i, intensity;
	float binVal;
	double maxv = 0;
	minMaxIdx(hist, 0, &maxv, 0, 0);
	Mat histImg(hist_w, hist_h, CV_8U, Scalar::all(255));
	for (i = 0; i < 256; i++) {
		binVal = hist.at<float>(i);
		intensity = static_cast<int>(binVal*hist_w / maxv);
		line(histImg, Point(zoom_size*i, hist_w), Point(zoom_size*i, hist_w - intensity), Scalar::all(0));
	}
	return histImg;
}

/*@brief Generate the random number with gaussian distribution.*/
float GaussianRand(float start, float end, float sigma, float mu) {
	if (start > end) {
		float count = start;
		start = end;
		end = count;
	}
	float u1, u2;
	float z0 = start - 1;
	while (z0<start || z0>end) {
		u1 = rand()*(1.0 / RAND_MAX);
		u2 = rand()*(1.0 / RAND_MAX);
		z0 = sqrt(-2.0*log(u1))*cos(2 * CV_PI*u2)*sigma + mu;
	}
	return z0;
}

/*@brief Generate the background with gaussian distribution.*/
Mat GenerateBackground(int rows, int cols, float start, float end, float sigma, float mu) {
	Mat dstImage(rows, cols, CV_32F, Scalar::all(0));
	int i, j;
	float *pdstData;
	for (i = 0; i < rows; i++) {
		pdstData = dstImage.ptr<float>(i);
		for (j = 0; j < cols; j++) {
			pdstData[j] = GaussianRand(start, end, sigma, mu);
		}
	}
	dstImage.convertTo(dstImage, CV_8U);
	return dstImage;
}

/*@brief Using the Image histogram imformation to get a threshold to remove the background.*/
int GetThreshold(MatND hist) {
	hist = hist / sum(hist)(0);
	MatND unithist;
	hist.copyTo(unithist);
	int i;
	int rows = hist.rows;
	for (i = 0; i < rows; i++) {
		unithist.at<float>(i) = (float)i*hist.at<float>(i);
	}
	MatND sighist(rows, hist.cols, hist.type(), Scalar::all(0));
	float w0, w1, u0, u1;
	for (i = 1; i < rows; i++) {
		w0 = sum(hist(Range(0, i), Range(0, 1)))(0);
		if (w0 > 0 && w0 != 1 && hist.at<float>(i) == 0) {
			sighist.at<float>(i) = sighist.at<float>(i - 1) + 1.0;
		}
	}
	Mat ss = HistImg(sighist,300,512);
	imshow("Sighist", ss);
	double maxv;
	int maxI[2];
	minMaxIdx(sighist, 0, &maxv, 0, maxI);	
	return maxI[0];
}

void Q2() {
	Mat BackG = GenerateBackground(256, 256, 200, 220, 10, 210);
	imshow("Background", BackG);
	Mat Object = GenerateBackground(100, 100, 80, 100, 10, 90);
	imshow("Object", Object);
	Mat srcImage;
	BackG.copyTo(srcImage);
	int i, j;
	int Brows = BackG.rows, Bcols = BackG.cols;
	int Orows = Object.rows, Ocols = Object.cols;
	int Start_row = (Brows - Orows) / 2, Start_col = (Bcols - Ocols) / 2;
	int End_row = Start_row + Orows, End_col = Start_col + Ocols;
	uchar *psrcData, *pOData, *pdstData;
	for (i = Start_row; i < End_row; i++) {
		psrcData = srcImage.ptr<uchar>(i);
		pOData = Object.ptr<uchar>(i - Start_row);
		for (j = Start_col; j < End_col; j++) {
			psrcData[j] = pOData[j - Start_col];
		}
	}
	imshow("Result", srcImage);
	MatND src_hist = Hist(srcImage);
	Mat s_h = HistImg(src_hist, 300, 512);
	imshow("Result Hist", s_h);
	int th = GetThreshold(src_hist);
	Mat dstImage(Brows, Bcols, CV_8U, Scalar::all(0));
	for (i = 0; i < Brows; i++) {
		psrcData = srcImage.ptr<uchar>(i);
		pdstData = dstImage.ptr<uchar>(i);
		for (j = 0; j < Bcols; j++) {
			if (psrcData[j] < th) {
				pdstData[j] = psrcData[j];
			}
		}
	}
	imshow("Removed Background", dstImage);
}

/*@brief Perform spatial filtering of an image.
@param srcImage is the source image.
@param mask is the filter mask.*/
Mat SpatialFilter(Mat srcImage, Mat mask) {
	int srows = srcImage.rows;
	int scols = srcImage.cols;
	int mrows = mask.rows;
	int mcols = mask.cols;
	int mrd = mrows / 2;
	int mcd = mcols / 2;
	int i, j, m, n;
	int S = sum(mask)(0);
	srcImage.convertTo(srcImage, CV_32F);
	mask.convertTo(mask, CV_32F);
	Mat dstImage;
	srcImage.copyTo(dstImage);
	float *psrcData, *pdstData, *pmData, counter;
	for (i = mrd; i < srows - mrd; i++) {
		pdstData = dstImage.ptr<float>(i);
		for (j = mcd; j < scols - mcd; j++) {
			counter = 0;
			for (m = -mrd; m < mrd + 1; m++) {
				psrcData = srcImage.ptr<float>(i + m);
				pmData = mask.ptr<float>(m + mrd);
				for (n = -mcd; n < mcd + 1; n++) {
					counter += psrcData[j + n] * pmData[n + mcd];
				}
			}
			pdstData[j] = counter / S;
		}
	}
	dstImage.convertTo(dstImage, CV_8U);
	return dstImage;
}

/*@brief Add a random noise to an image.
@param srcImage is the source Image.
@param p is the proporation of the noise.
@param start and end are the range of the noise pixels value.*/
Mat AddRandNoise(Mat srcImage, float p, uchar start, uchar end) {
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int channels = srcImage.channels();
	if (start > end) {
		uchar temp = end;
		end = start;
		start = temp;
	}
	uchar *psrcData, noise;
	uchar gap = end - start;
	int i, j;
	float r1, r2;
	for (i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<uchar>(i);
		for (j = 0; j < cols*channels; j++) {
			r1 = rand() * (1.0/RAND_MAX);
			if (r1 < p) {
				r2 = rand()*(1.0 / RAND_MAX);
				psrcData[j] = r2*gap + start;
			}
		}
	}
	return srcImage;
}