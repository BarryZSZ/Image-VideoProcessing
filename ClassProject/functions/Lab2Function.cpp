#include "head\ClassProjectHead.h"

void Lab2Main() {
	Mat im1 = imread("Images\\lena.tiff");
	Mat im2 = imread("Images\\Q2_3.tif");
	Mat im3 = imread("Images\\baboon.tiff");
	Mat im4 = imread("Images\\peppers(1).tiff");
	//Mat im  = imread("Images\\Plant.tif");
	Mat im  = imread("Images\\lena.pgm");

	Mat im_o1 = Translation(im, 50, -70);
	imshow("Translation1", im_o1);
	im_o1 = Translation(im, 50, 70);
	imshow("Translation2", im_o1);
	im_o1 = Translation(im, -50, -70);
	imshow("Translation3", im_o1);
	im_o1 = Translation(im, -50, 70);
	imshow("Translation4", im_o1);

	Mat im_o2 = Rotation(im, 30);
	imshow("Rotation1", im_o2);
	im_o2 = Rotation(im, -30);
	imshow("Rotation2", im_o2);

	Mat im_o3 = ShearVertical(im, 1);
	imshow("ShearVectical", im_o3);
	Mat im_o4 = ShearHorizontal(im, 1);
	imshow("ShearHorizonal", im_o4);

	Mat im_o5 = AveragingFilter(im, 3);
	imshow("AveragingFilter 3*3", im_o5);
	Mat im_o6 = AveragingFilter(im, 5);
	imshow("AveragingFilter 5*5", im_o6);
	Mat im_o17 = MedianFilter(im, 3, 3);
	imshow("MedianFilter 3*3", im_o17);
	im_o17 = MedianFilter(im, 5, 5);
	imshow("MedianFilter 5*5", im_o17);
	Mat im_o18 = OtsuBinarization(im);
	imshow("OtsuBinarization", im_o18);

	Mat im_o7 = LaplacianOperators(im, -1);
	imshow("Laplacian", im_o7);
	Mat im_o8 = SobelOperators(im, -0.15);
	imshow("Sobel1", im_o8);

	Mat im_o9 = GammaCorrection(im, 0.1);
	imshow("Gamma = 0.1", im_o9);
	Mat im_o10 = GammaCorrection(im, 0.4);
	imshow("Gamma = 0.4", im_o10);
	Mat im_o11 = GammaCorrection(im, 0.6);
	imshow("Gamma = 0.6", im_o11);
	Mat im_o12 = GammaCorrection(im, 0.8);
	imshow("Gamma = 0.8", im_o12);
	Mat im_o13 = GammaCorrection(im, 1);
	imshow("Gamma = 1", im_o13);

	Mat im_o14 = GlobalHistogram(im1);
	imshow("GlobalHistogram", im_o14);
	Mat im_o15 = LocalHistogram(im2, 3);
	imshow("LocalHistogram 3*3", im_o15);
	Mat im_o16 = LocalHistogram(im2, 5);
	imshow("LocalHistogram 5*5", im_o16);
	Mat im_o19 = MedianFilter(im, 5, 5);
	imshow("MedianFilter 3*3", im_o17);
	Mat im_o20 = OtsuBinarization(im);
	imshow("OtsuBinarization", im_o20);
	waitKey(0);
}
/*@para: srcImage is the Source Image
dx is the width translation
dy is the height translation*/
Mat Translation(Mat srcImage, int dx, int dy) {
	if (!srcImage.data) {
		cout << "Error: The image is empty..." << endl;
	}
	Mat dstImage = Mat(srcImage.size(), srcImage.type(), Scalar::all(0));
	int i, j;
	int iy = 0;
	int jx = 0;
	for (i = 0; i < srcImage.rows & iy < srcImage.rows - 1; i++) {
		iy = i - dy;
		if (iy < 0) {
			i = -iy;
			iy = i - dy;
		}
		jx = 0;
		for (j = 0; j < srcImage.cols & jx < srcImage.cols - 1; j++) {
			jx = j - dx;
			if (jx < 0) {
				j = -jx;
				jx = j - dx;
			}
			dstImage.at<Vec3b>(i, j) = srcImage.at<Vec3b>(iy, jx);
		}
	}
	return dstImage;
}

/*@para: srcImage is the source Image
angle is angle of rotation
Notice: angle is the angle, not the radian*/
Mat Rotation(Mat srcImage, double angle) {
	if (!srcImage.data) {
		cout << "Error: The image is empty..." << endl;
	}
	double rad = angle * CV_PI / 180;
	double sin_angle = sin(rad);
	double cos_angle = cos(rad);
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	double a1 = cols*cos_angle;
	double a2 = -rows*sin_angle;
	double a3 = cols*cos_angle - rows*sin_angle;
	double b1 = rows*cos_angle;
	double b2 = cols*sin_angle;
	double b3 = rows*cos_angle + cols*sin_angle;
	double xmin = min(min(min(0.0, a1), a2), a3);
	double ymin = min(min(min(0.0, b1), b2), b3);
	double xmax = max(max(max(0.0, a1), a2), a3);
	double ymax = max(max(max(0.0, b1), b2), b3);
	int nrows = ymax - ymin + 1;
	int ncols = xmax - xmin + 1;
	Mat dstImage = Mat(nrows, ncols, srcImage.type(), Scalar::all(0));
	int i, j, iy, jx;
	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {
			iy = (i + ymin)*cos_angle - (j + xmin)*sin_angle + 0.5;
			jx = (j + xmin)*cos_angle + (i + ymin)*sin_angle + 0.5;
			if (iy >= 0 && iy < rows && jx >= 0 && jx < cols) {
				dstImage.at<Vec3b>(i, j) = srcImage.at<Vec3b>(iy, jx);
			}
		}
	}
	return dstImage;
}

/*@para: srcImage is the source Image
s is the shearing value*/
Mat ShearVertical(Mat srcImage, double s) {
	if (!srcImage.data) {
		cout << "Error: The image is empty..." << endl;
	}
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int nrows = s*cols + rows;
	Mat dstImage(nrows, cols, srcImage.type(), Scalar::all(0));
	int i, j, iy;
	for (i = 0; i < nrows; i++) {
		for (j = 0; j < cols; j++) {
			iy = i - s*j;
			if (iy >= 0 && iy < rows) {
				dstImage.at<Vec3b>(i, j) = srcImage.at<Vec3b>(iy, j);
			}
		}
	}
	return dstImage;
}

/*@para: srcImage is the source Image
s is the shearing value*/
Mat ShearHorizontal(Mat srcImage, double s) {
	if (!srcImage.data) {
		cout << "Error: The image is empty..." << endl;
	}
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int ncols = cols + s*rows;
	Mat dstImage(rows, ncols, srcImage.type(), Scalar::all(0));
	int i, j, jx;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < ncols; j++) {
			jx = j - s*i;
			if (jx >= 0 && jx < cols) {
				dstImage.at<Vec3b>(i, j) = srcImage.at<Vec3b>(i, jx);
			}
		}
	}
	return dstImage;
}

/*@para: srcImage is the source Image
wsize is the window size*/
Mat AveragingFilter(Mat srcImage, int wsize) {
	if (!srcImage.data) {
		cout << "Error: The image is empty..." << endl;
	}
	Mat dstImage(srcImage.size(), srcImage.type());
	int i, j, m, n;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	Vec3f count;
	Vec3b a = srcImage.at<Vec3b>(1, 1);
	int d = wsize / 2;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			if (i < d || j < d || i >= rows - d || j >= cols - d) {
				dstImage.at<Vec3b>(i, j) = srcImage.at<Vec3b>(i, j);
			}
			else
			{
				count = { 0,0,0 };
				for (m = -d; m <= d; m++) {
					for (n = -d; n <= d; n++) {
						if (i + m >= 0 && i + m < rows && j + n >= 0 && j + n < cols) {
							count = (Vec3f)srcImage.at<Vec3b>(i + m, j + n) + count;
						}
					}
				}
				dstImage.at<Vec3b>(i, j) = (Vec3b)(count / wsize / wsize);
			}

		}
	}
	return dstImage;
}

/*@para: srcImage is the source Image
wsize is the window size*/
Mat MedianFilter(Mat srcImage, int wsize) {
	if (!srcImage.data) {
		cout << "Error: The image is empty..." << endl;
	}
	Mat dstImage = srcImage.clone();
	int i, j, m, n, k, p, q;
	int d = wsize / 2;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int channels = srcImage.channels();
	int halve = wsize*wsize / 2;
	int left_count = wsize*wsize;
	int right_count = wsize*wsize;
	int median_num;
	uchar* pdstData;
	uchar* psrcData;
	uchar* psrcData1;
	Mat mask(wsize, wsize, srcImage.type());
	for (i = d; i < rows - d; i++) {
		pdstData = dstImage.ptr<uchar>(i);
		for (j = d; j < cols - d; j++) {
			for (k = 0; k < channels; k++) {
				m = 0;
				while (m < wsize) {
					psrcData = srcImage.ptr<uchar>(i + m - d);
					n = 0;
					while (n < wsize) {
						median_num = psrcData[(j + n - d)*channels + k];
						left_count = 0;
						right_count = 0;
						for (p = 0; p < wsize; p++) {
							psrcData1 = srcImage.ptr<uchar>(i + p - d);
							for (q = 0; q < wsize; q++) {
								if (median_num > psrcData1[(j + q - d)*channels + k]) {
									left_count++;
								}
								else if (median_num < psrcData1[(j + q - d)*channels + k]) {
									right_count++;
								}
							}
						}
						if (left_count <= halve&&right_count <= halve) {
							pdstData[j*channels + k] = median_num;
							break;
						}
						n++;
					}
					if (left_count <= halve&&right_count <= halve) {
						break;
					}
					m++;
				}

			}
		}
	}
	return dstImage;
}

/*@para: srcImage is the source Image
k is the weight of the blurred image*/
Mat LaplacianOperators(Mat srcImage, float k) {
	if (!srcImage.data) {
		cout << "Error: The image is empty..." << endl;
	}
	Mat dstImage = srcImage.clone();
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	Vec3f gmask;
	int i, j;
	for (i = 1; i < rows - 1; i++) {
		for (j = 1; j < cols - 1; j++) {
			gmask = (Vec3f)srcImage.at<Vec3b>(i - 1, j) + (Vec3f)srcImage.at<Vec3b>(i, j - 1) +
				(Vec3f)srcImage.at<Vec3b>(i + 1, j) + (Vec3f)srcImage.at<Vec3b>(i, j + 1) - 4 * (Vec3f)srcImage.at<Vec3b>(i, j);
			dstImage.at<Vec3b>(i, j) = (Vec3b)(k*gmask + (Vec3f)dstImage.at<Vec3b>(i, j));
		}
	}
	return dstImage;
}

/*@para: srcImage is the source Image
k is the weight of the blurred image*/
Mat SobelOperators(Mat srcImage, float k) {
	if (!srcImage.data) {
		cout << "Error: The image is empty..." << endl;
	}
	Mat dstImage = srcImage.clone();
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	Vec3f gy, gx, mag;
	int i, j, m;
	for (i = 1; i < rows - 1; i++) {
		for (j = 1; j < cols - 1; j++) {
			gy = (Vec3f)srcImage.at<Vec3b>(i - 1, j - 1) + 2 * (Vec3f)srcImage.at<Vec3b>(i - 1, j) + (Vec3f)srcImage.at<Vec3b>(i - 1, j + 1)
				- (Vec3f)srcImage.at<Vec3b>(i + 1, j - 1) - 2 * (Vec3f)srcImage.at<Vec3b>(i + 1, j) - (Vec3f)srcImage.at<Vec3b>(i + 1, j + 1);
			gx = (Vec3f)srcImage.at<Vec3b>(i - 1, j - 1) + 2 * (Vec3f)srcImage.at<Vec3b>(i, j - 1) + (Vec3f)srcImage.at<Vec3b>(i + 1, j - 1)
				- (Vec3f)srcImage.at<Vec3b>(i - 1, j + 1) - 2 * (Vec3f)srcImage.at<Vec3b>(i, j + 1) - (Vec3f)srcImage.at<Vec3b>(i + 1, j + 1);
			for (m = 0; m < 3; m++) {
				gx[m] = gx[m] * gx[m];
				gy[m] = gy[m] * gy[m];
				mag[m] = k*sqrt(gx[m] + gy[m]);
			}
			dstImage.at<Vec3b>(i, j) = (Vec3b)((Vec3f)dstImage.at<Vec3b>(i, j) + mag);
		}
	}
	return dstImage;
}

void VarianceOfImage(Mat srcImage) {
	double average_value = 0;
	double variance = 0;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int channels = srcImage.channels();
	int i, j;
	uchar* psrcData;
	for (i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<uchar>(i);
		for (j = 0; j < cols*channels; j++) {
			average_value = (double)psrcData[j] + average_value;
		}
	}
	average_value = average_value / rows / cols / channels;
	for (i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<uchar>(i);
		for (j = 0; j < cols*channels; j++) {
			variance = variance + pow(((double)psrcData[j] - average_value), 2);
		}
	}
	variance = variance / rows / cols / channels;
	cout << "Averaging value: " << average_value << " Variance: " << variance << endl;
}
Mat GammaCorrection(Mat srcImage, float gamma) {
	if (!srcImage.data) {
		cout << "Error: The image is empty..." << endl;
	}
	Mat dstImage = srcImage.clone();
	int rows = dstImage.rows;
	int nStep = dstImage.cols*dstImage.channels();
	int i, j;
	float c = 255 / pow(255, gamma);
	for (i = 0; i < rows; i++) {
		uchar* pDstData = dstImage.ptr<uchar>(i);
		for (j = 0; j < nStep; j++) {
			pDstData[j] = c*pow(pDstData[j], gamma);
		}
	}
	VarianceOfImage(dstImage);
	return dstImage;
}

Mat GlobalHistogram(Mat srcImage) {
	Mat dstImage = srcImage.clone();
	int rows = dstImage.rows;
	int cols = dstImage.cols;
	int channels = dstImage.channels();
	Mat histogram(channels, 256, CV_64F, Scalar::all(0));
	int i, j, k, counter;
	uchar* psrcData;
	uchar* pdstData;
	for (i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<uchar>(i);
		for (j = 0; j < cols; j++) {
			for (k = 0; k < channels; k++) {
				counter = psrcData[channels*j + k];
				histogram.at<double>(k, counter) = histogram.at<double>(k, counter) + 1;
			}
		}
	}
	Mat equalize_histogram(histogram.size(), histogram.type());
	for (i = 0; i < channels; i++) {
		double* phisData = histogram.ptr<double>(i);
		double* pequData = equalize_histogram.ptr<double>(i);
		for (j = 0; j < 256; j++) {
			if (j>0) {
				phisData[j] = phisData[j] + phisData[j - 1];
			}
			pequData[j] = (int)(255 * phisData[j] / rows / cols + 0.5);
		}
	}
	for (i = 0; i < rows; i++) {
		pdstData = dstImage.ptr<uchar>(i);
		for (j = 0; j < cols; j++) {
			for (k = 0; k < channels; k++) {
				counter = pdstData[channels*j + k];
				pdstData[channels*j + k] = (uchar)equalize_histogram.at<double>(k, counter);
			}
		}
	}
	return dstImage;
}

Mat LocalHistogram(Mat srcImage, int wsize) {
	Mat dstImage = srcImage.clone();
	Mat local_Image(wsize, wsize, srcImage.type());
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int channels = srcImage.channels();
	int d = wsize / 2;
	int i, j, k, m, n, p;
	uchar* pdstData;
	uchar* plocData;
	uchar* psrcData;
	for (i = d; i < rows - d; i++) {
		pdstData = dstImage.ptr<uchar>(i);
		for (j = d; j < cols - d; j++) {
			for (m = 0; m < wsize; m++) {
				plocData = local_Image.ptr<uchar>(m);
				psrcData = srcImage.ptr<uchar>(i + m - d);
				for (n = 0; n < wsize; n++) {
					for (k = 0; k < channels; k++) {
						plocData[channels*n + k] = psrcData[channels*(j + n - d) + k];
					}
				}
			}
			local_Image = GlobalHistogram(local_Image);
			plocData = local_Image.ptr<uchar>(d + 1);
			for (p = 0; p < channels; p++) {
				pdstData[channels*j + p] = plocData[channels*(d + 1) + p];
			}
		}
	}
	return dstImage;
}