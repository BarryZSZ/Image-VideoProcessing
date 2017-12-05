#include "head\ClassProjectHead.h"

void Lab8Main() {
	//Lab8Q1();
	//Lab8Q2();
	//Lab8Q3();
	//Lab8Q4();
}

void Lab8Q1() {
	Mat I1 = imread("Images\\large_septagon_gaussian_noise_mean_0_std_50_added.pgm", 0);

	Mat I2, I3, I4;
	Mat kern(5, 5, CV_32F, Scalar::all(0.04));
	I1.convertTo(I1, CV_32F);
	filter2D(I1, I2, CV_32F, kern);
	I1.convertTo(I1, CV_8U);
	I2.convertTo(I2, CV_8U);

	I3 = OtsuBinarization(I1);
	I4 = OtsuBinarization(I2);

	imshow("Original image", I1);
	imshow("Average smoothed 5*5", I2);
	imshow("Otsu original image", I3);
	imshow("Otsu smoothed image", I4);
}

void Lab8Q2() {
	Mat I1 = imread("Images\\septagon_noisy_shaded.pgm", 0);
	imshow("Original image", I1);
	Mat I2;
	Mat kern(5, 5, CV_32F, Scalar::all(0.04));
	I1.convertTo(I1, CV_32F);
	filter2D(I1, I2, CV_32F, kern);
	I1.convertTo(I1, CV_8U);
	I2.convertTo(I2, CV_8U);

	Mat I3 = PartitionOtsu(I2, 2, 3);
	imshow("Partition Otsu", I3);

	Mat I4 = OtsuBinarization(I1);
	imshow("Otsu original image", I4);
}

void Lab8Q3() {
	Mat I1 = imread("Images\\spot_shaded_text_image.pgm", 0);
	imshow("Original Image", I1);

	Mat I2 = MovingAverageThreshold(I1, 20, 0.1);
	imshow("Moving Average Threshold", I2);

	Mat I3 = OtsuBinarization(I1);
	imshow("Otsu Original Image", I3);
}

void Lab8Q4() {
	Mat I1 = imread("Images\\defective_weld.pgm", 0);
	Mat I2 = imread("Images\\noisy_region.pgm", 0);

	imshow("Defective weld", I1);
	imshow("Noisy Region", I2);
}

Mat OtsuBinarization(Mat srcImage) {
	int rows = srcImage.rows;
	int colch = srcImage.cols*srcImage.channels();
	MatND hist(256, 1, CV_64F, Scalar::all(0));
	uchar *psrcData;
	int i, k, T;
	for (i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<uchar>(i);
		for (k = 0; k < colch; k++) {
			hist.at<double>(psrcData[k], 0)++;
		}
	}
	hist = hist / ((double)(rows*colch));
	MatND histM(hist.size(), hist.type());
	hist.copyTo(histM);
	for (i = 0; i < 256; i++) {
		histM.at<double>(i, 0) = histM.at<double>(i, 0)*(double)i;
	}
	double mG = sum(histM)(0);
	double P1, m;
	double sig = 0.0, sigB = 0.0;
	for (k = 0; k < 255; k++) {
		P1 = (double)sum(hist(Range(0, k), Range(0, 1)))(0);
		m= (double)sum(histM(Range(0, k), Range(0, 1)))(0);
		if (P1*(1 - P1) != 0.0) {
			sigB = pow((mG*P1 - m), 2) / (P1*(1 - P1));
		}
		if (sigB > sig) {
			T = k;
		}
		sig = sigB;
	}
	Mat dstImage;
	threshold(srcImage, dstImage, T, 255, THRESH_BINARY);
	return dstImage;
}

/*@brief partition method first then Otsu's method to segement.
@param srcImage is source Image.
@param h is the partition height.
@param w is the partition width.*/
Mat PartitionOtsu(Mat srcImage, int h, int w) {
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int drows = rows / h;
	int dcols = cols / w;
	int i, j, rowend, colend;
	Mat part;
	Mat dstImage(rows, cols, CV_8U, Scalar::all(0));
	for (i = 0; i < h; i++) {
		for (j = 0; j < w; j++) {
			rowend = (i + 1)*drows;
			colend = (j + 1)*dcols;
			if (rowend > rows) {
				rowend = rows;
			}
			if (colend > cols) {
				colend = cols;
			}
			part = OtsuBinarization(srcImage(Range(i*drows, rowend), Range(j*dcols, colend)));
			part.copyTo(dstImage(Range(i*drows, rowend), Range(j*dcols, colend)));
		}
	}
	return dstImage;
}

/*@brief Moving Average Thresholding.*/
Mat MovingAverageThreshold(Mat srcImage, int n, float b) {
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	srcImage.convertTo(srcImage, CV_32F);
	Mat dstImage(rows, cols, CV_8U, Scalar::all(0));
	Mat Move(1, rows*cols, CV_32F, Scalar::all(0));
	float *psrcData, *pMoveData;
	uchar *pdstData;
	int i, j, k;
	int count = 0;
	float mk, Tk;
	pMoveData = Move.ptr<float>(0);
	for (i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<float>(i);
		for (j = 0; j < cols; j++) {
			if (i % 2 == 1) {
				k = cols - j - 1;
			}
			else
			{
				k = j;
			}
			pMoveData[i*rows + j] = psrcData[k];
		}
	}
	for (i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<float>(i);
		pdstData = dstImage.ptr<uchar>(i);
		for (j = 0; j < cols; j++) {
			if (count < 20) {
				mk = sum(Move(Range(0, 1), Range(0, count + 1)))(0) / n;
			}
			else {
				mk = sum(Move(Range(0, 1), Range(count - 20, count + 1)))(0) / n;
			}
			Tk = b*mk;
			if (psrcData[j] > Tk) {
				pdstData[j] = 255;
			}
			count++;
		}
	}
	return dstImage;
}

//Mat RegionGrowing()
