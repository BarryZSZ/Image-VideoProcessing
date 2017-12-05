#include "head\ClassProjectHead.h"

void Lab7Main() {
	Lab7Q1();
	Lab7Q2();
	Lab7Q3();
}

void Lab7Q1() {
	Mat I1 = imread("Images\\headCT-Vandy.pgm", 0);
	Mat I2 = imread("Images\\building_original.pgm", 0);
	Mat I3 = imread("Images\\noisy_fingerprint.tif", 0);


	imshow("Head CT-Vandy", I1);
	imshow("Buliding Original", I2);
	imshow("Noisy fingerprint", I3);

	Mat I4, I5, I6;
	//Roberts
	I4 = RobertsGradient(I1);
	imshow("Robert Head CT-Vandy", I4);
	I4 = GlobalThreshold(I4, 1);
	imshow("Thresholded Robert Head CT-Vandy", I4);

	I4 = RobertsGradient(I2);
	imshow("Robert Buliding Original", I4);
	I4 = GlobalThreshold(I4, 1);
	imshow("Thresholded Robert Buliding Original", I4);

	I4 = RobertsGradient(I3);
	imshow("Robert Noisy fingerprint", I4);
	I4 = GlobalThreshold(I4, 1);
	imshow("Thresholded Robert Noisy fingerprint", I4);

	//Prewitt
	I5 = PrewittGradient(I1);
	imshow("Prewitt Head CT-Vandy", I5);
	I5 = GlobalThreshold(I5, 1);
	imshow("Thresholded Prewitt Head CT-Vandy", I5);

	I5 = PrewittGradient(I2);
	imshow("Prewitt Buliding Original", I5);
	I5 = GlobalThreshold(I5, 1);
	imshow("Thresholded Prewitt Buliding Original", I5);

	I5 = PrewittGradient(I3);
	imshow("Prewitt Noisy fingerprint", I5);
	I5 = GlobalThreshold(I5, 1);
	imshow("Thresholded Prewitt Noisy fingerprint", I5);

	//Sobel
	I6 = SobelGradient(I1);
	imshow("Sobel Head CT-Vandy", I6);
	I6 = GlobalThreshold(I6, 1);
	imshow("Thresholded Sobel Head CT-Vandy", I6);

	I6 = SobelGradient(I2);
	imshow("Sobel Buliding Original", I6);
	I6 = GlobalThreshold(I6, 1);
	imshow("Thresholded Sobel Buliding Original", I6);

	I6 = SobelGradient(I3);
	imshow("Sobel Noisy fingerprint", I6);
	I6 = GlobalThreshold(I6, 1);
	imshow("Thresholded Sobel Noisy fingerprint", I6);
}

void Lab7Q2() {
	Mat I1 = imread("Images\\headCT-Vandy.pgm", 0);
	Mat I2 = imread("Images\\noisy_fingerprint.tif", 0);

	imshow("Head CT-Vandy", I1);
	imshow("Noisy Fingerprint", I2);

	Mat I3 = LaplacianOfGaussian(I1, 1);
	Mat I4 = LaplacianOfGaussian(I2, 1);
	I3.convertTo(I3, CV_8U);
	I4.convertTo(I4, CV_8U);
	imshow("LOG Head CT-Vandy", I3);
	imshow("LOG Noisy Fingerprint", I4);

	int edgeThresh1 = 40;
	int edgeThresh2 = 50;
	Mat I5 = CannyEdgeDetection(I1, edgeThresh1, edgeThresh1 * 2);
	Mat I6 = CannyEdgeDetection(I2, edgeThresh2, edgeThresh1 * 3);

	imshow("Canny Edge headCT-Vandy", I5);
	imshow("Canny Edge Noisy Fingerprint", I6);

}

void Lab7Q3() {
	Mat I1 = imread("Images\\polymersomes.pgm", 0);
	Mat I2 = imread("Images\\noisy_fingerprint.tif", 0);

	imshow("Polymersomes", I1);
	imshow("Noisy fingerprint", I2);


	Mat I3 = GlobalThreshold(I1, 1);
	Mat I4 = GlobalThreshold(I2, 1);
		
	imshow("Global thresholded Polymersomes", I3);
	imshow("Global thresholded Noisy fingerprint", I4);

	

}


/*@brief Global threshold image.*/
Mat GlobalThreshold(Mat InputArray, float deltaT) {
	int rows = InputArray.rows;
	int colch = InputArray.cols*InputArray.channels();
	double maxv, minv;
	minMaxIdx(InputArray, &minv, &maxv, 0, 0);
	float T0 = (maxv + minv) / 2;
	MatND hist = Hist(InputArray);
	MatND histM;
	hist.copyTo(histM);
	hist = hist / rows / colch;
	int i, j, T;
	for (i = 0; i < 256; i++) {
		histM.at<float>(i, 0) = histM.at<float>(i, 0)*i / rows / colch;
	}
	float temp = deltaT + 1;
	float m1, m2;
	while (temp > deltaT) {
		T = (int)(T0 + 0.5);
		if (sum(hist(Range(0, T), Range(0, 1)))(0) != 0) {
			m1 = sum(histM(Range(0, T), Range(0, 1)))(0) / sum(hist(Range(0, T), Range(0, 1)))(0);
		}
		else {
			m1 = 0;
		}
		if (sum(hist(Range(T, 256), Range(0, 1)))(0) != 0) {
			m2 = sum(histM(Range(T, 256), Range(0, 1)))(0) / sum(hist(Range(T, 256), Range(0, 1)))(0);
		}
		else {
			m2 = 0;
		}
		temp = abs(T0 - (m1 + m2) / 2);
		T0 = (m1 + m2) / 2;
	}
	T = (int)(T0 + 0.5);
	uchar *pInputData, *pOutputData;
	Mat OutputArray(InputArray.size(), InputArray.type(), Scalar::all(0));
	for (i = 0; i < rows; i++) {
		pInputData = InputArray.ptr<uchar>(i);
		pOutputData = OutputArray.ptr<uchar>(i);
		for (j = 0; j < colch; j++) {
			if (pInputData[j] > T) {
				pOutputData[j] = 255;
			}
		}
	}
	return OutputArray;
}

/*@brief Get Roberts gradian image.*/
Mat RobertsGradient(Mat inputArray) {
	int rows = inputArray.rows;
	int cols = inputArray.cols;
	int i, j;
	float gx, gy;
	Mat outputArray;
	inputArray.copyTo(outputArray);
	for (i = 0; i < rows - 1; i++) {
		for (j = 0; j < cols - 1; j++) {
			gx = inputArray.at<uchar>(i + 1, j + 1) - inputArray.at<uchar>(i, j);
			gy = inputArray.at<uchar>(i + 1, j) - inputArray.at<uchar>(i, j + 1);
			outputArray.at<uchar>(i, j) = (uchar)sqrt(gx*gx + gy*gy);
		}
	}
	return outputArray;
}

/*@brief Get Prewitt gradian image.*/
Mat PrewittGradient(Mat srcImage) {
	Mat xkern = (Mat_<float>(3, 3) << -1, -1, -1, 0, 0, 0, 1, 1, 1);
	Mat ykern = (Mat_<float>(3, 3) << -1, 0, 1, -1, 0, 1, -1, 0, 1);
	Mat xdstImage, ydstImage;
	filter2D(srcImage, xdstImage, CV_32F, xkern);
	filter2D(srcImage, ydstImage, CV_32F, ykern);
	Mat dstImage(srcImage.size(), CV_32F);
	multiply(xdstImage, xdstImage, xdstImage);
	multiply(ydstImage, ydstImage, ydstImage);
	dstImage = xdstImage + ydstImage;
	sqrt(dstImage, dstImage);
	dstImage.convertTo(dstImage, CV_8U);
	return dstImage;
}

/*@brief Get Sobel gradian image.*/
Mat SobelGradient(Mat srcImage) {
	Mat xkern = (Mat_<float>(3, 3) << -1, -2, -1, 0, 0, 0, 1, 2, 1);
	Mat ykern = (Mat_<float>(3, 3) << -1, 0, 1, -2, 0, 2, -1, 0, 1);
	Mat xdstImage, ydstImage;
	filter2D(srcImage, xdstImage, CV_32F, xkern);
	filter2D(srcImage, ydstImage, CV_32F, ykern);
	Mat dstImage(srcImage.size(), CV_32F);
	multiply(xdstImage, xdstImage, xdstImage);
	multiply(ydstImage, ydstImage, ydstImage);
	dstImage = xdstImage + ydstImage;
	sqrt(dstImage, dstImage);
	dstImage.convertTo(dstImage, CV_8U);
	return dstImage;
}

/*@brief Get Gaussian window from sigma.*/
Mat GetGaussianWindow(float sigma) {
	int wsize = 3 * sigma + 0.5;
	Mat window(2 * wsize + 1, 2 * wsize + 1, CV_32F);
	int i, j;
	for (i = -wsize; i <= wsize; i++) {
		for (j = -wsize; j <= wsize; j++) {
			window.at<float>(i + wsize, j + wsize) = (((float)(i*i + j*j) - 2 * sigma*sigma) / pow(sigma, 4))*exp(-((float)(i*i + j*j)) / (2 * sigma*sigma));
		}
	}
	window /= sum(window)(0);
	return window;
}

/*@brief Laplacian of Gaussian filter.*/
Mat LaplacianOfGaussian(Mat srcImage, float sigma) {
	Mat kernG = GetGaussianWindow(sigma);
	Mat kernL = (Mat_<float>(3, 3) << 1, 1, 1, 1, -8, 1, 1, 1, 1);
	double minv, maxv;
	Mat dstImage;
	srcImage.convertTo(srcImage, CV_32F);
	filter2D(srcImage, dstImage, CV_32F, kernG);
	minMaxIdx(dstImage, &minv, &maxv, 0, 0);
	filter2D(dstImage, dstImage, CV_32F, kernL);
	threshold(dstImage, dstImage, maxv*0.04, maxv, THRESH_TOZERO);
	normalize(dstImage, dstImage, 0, 255, NORM_MINMAX);
	return dstImage;
}

/*@brief Canny Edge Detection.*/
Mat CannyEdgeDetection(Mat srcImage, float thresholdLower, float thresholdUpper) {
	Mat Image(srcImage.size(), CV_32F);
	srcImage.convertTo(srcImage, CV_32F);
	//Gaussian filter
	//GaussianBlur(srcImage, Image, Size(7, 7), 1);
	Image = LaplacianOfGaussian(srcImage, 1.0);

	//Use Sobel operater to computer the gradient magnitude and direction
	Mat magX(srcImage.size(), CV_32F);
	Mat magY(srcImage.size(), CV_32F);
	Sobel(Image, magX, CV_32F, 1, 0, 3);
	Sobel(Image, magY, CV_32F, 0, 1, 3);

	//Compute the sloaps
	Mat directionI(srcImage.size(), CV_32F);
	divide(magY, magX, directionI);

	//Compute the gradient of each pixel
	Mat prodX(srcImage.size(), CV_32F);
	Mat prodY(srcImage.size(), CV_32F);
	Mat magnitudeI(srcImage.size(), CV_32F);

	multiply(magX, magX, prodX);
	multiply(magY, magY, prodY);
	magnitudeI = prodX + prodY;
	sqrt(magnitudeI, magnitudeI);
	
	//Apply nonmaxima suppression to the gradient magnitude image
	Mat checkImage(srcImage.size(), CV_8U);
	MatIterator_<float>itMag = magnitudeI.begin<float>();
	MatIterator_<float>itDirection=directionI.begin<float>();
	MatIterator_<float>itEnd = magnitudeI.end<float>();
	MatIterator_<uchar>itRet = checkImage.begin<uchar>();
	const Point pos;
	float currentDirection;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	for (; itMag != itEnd; itDirection++, itMag++, itRet++) {
		Point pos = itRet.pos();
		currentDirection = atan(*itDirection) * 180 / CV_PI;
		while (currentDirection < 0)
			currentDirection += 180;
		*itDirection = currentDirection;
		if (currentDirection > 22.5&&currentDirection <= 67.5) {
			if (pos.y > 0 && pos.x > 0 && pos.y < rows - 1 && pos.x < cols - 1) {
				if (*itMag <= magnitudeI.at<float>(pos.y - 1, pos.x - 1) || *itMag <= magnitudeI.at<float>(pos.y + 1, pos.x + 1)) {
					magnitudeI.at<float>(pos.y, pos.x) = 0;
				}
			}
			
		}
		else if (currentDirection > 67.5&&currentDirection <= 112.5) {
			if (pos.y > 0 && pos.x > 0 && pos.y < rows - 1 && pos.x < cols - 1) {
				if (*itMag <= magnitudeI.at<float>(pos.y - 1, pos.x) || *itMag <= magnitudeI.at<float>(pos.y + 1, pos.x)) {
					magnitudeI.at<float>(pos.y, pos.x) = 0;
				}
			}

		}
		else if (currentDirection > 112.5&&currentDirection <= 157.5) {
			if (pos.y > 0 && pos.x > 0 && pos.y < rows - 1 && pos.x < cols - 1) {
				if (*itMag <= magnitudeI.at<float>(pos.y + 1, pos.x - 1) || *itMag <= magnitudeI.at<float>(pos.y + 1, pos.x - 1)) {
					magnitudeI.at<float>(pos.y, pos.x) = 0;
				}
			}

		}
		else {
			if (pos.y > 0 && pos.x > 0 && pos.y < rows - 1 && pos.x < cols - 1) {
				if (*itMag <= magnitudeI.at<float>(pos.y, pos.x - 1) || *itMag <= magnitudeI.at<float>(pos.y, pos.x + 1)) {
					magnitudeI.at<float>(pos.y, pos.x) = 0;
				}
			}

		}
	}
	Mat edges;
	edgeDetect(magnitudeI, thresholdUpper, thresholdLower, edges);
	edges.convertTo(edges, CV_8U);
	return edges;
}

/*@brief Connect the edges.*/
void followEdges(int x, int y, Mat &magnitude, int tUpper, int tLower, Mat &edges){
	edges.at<float>(y, x) = 255;
	for (int i = -1; i < 2; i++) {
		for (int j = -1; j < 2; j++) {
			if ((i != 0) && (j != 0) && (x + i >= 0) &&
				(y + j >= 0) && (x + i < magnitude.cols) &&
				(y + j < magnitude.rows))
			{
				if ((magnitude.at<float>(y + j, x + i) > tLower)
					&& (edges.at<float>(y + j, x + i) != 255))
				{
					followEdges(x + i, y + j, magnitude, tUpper, tLower, edges);
				}
			}
		}
	}
}

/*@brief Edge detection*/
void edgeDetect(Mat &magnitude, int tUpper, int tLower, Mat &edges){
	int rows = magnitude.rows;
	int cols = magnitude.cols;
	edges = Mat(magnitude.size(), CV_32F, 0.0);
	for (int x = 0; x < cols; x++)
	{
		for (int y = 0; y < rows; y++)
		{
			if (magnitude.at<float>(y, x) >= tUpper)
			{
				followEdges(x, y, magnitude, tUpper, tLower, edges);
			}
		}
	}
}