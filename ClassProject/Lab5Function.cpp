#include "ClassProjectHead.h"

void Lab5Main() {
	//LenaNotchFilter();

	/*Mat I = imread("Images\\cameraWithNoise.pgm", 0);
	MatND hist = Hist(I);
	Mat histIm = HistImg(hist, 300, 512);
	imshow("Camera histogram", histIm);
	Mat backG = I(Range(0, 35), Range(0, 256));
	imshow("Background1", backG);
	hist = Hist(backG);
	histIm = HistImg(hist, 300, 512);
	imshow("Camera histogram1", histIm);
	imshow("Original", I);
	Mat I_O = AdaptiveMedianFilter(I, 3, 3, 7);
	imshow("Adaptive Median", I_O);
	backG = I_O(Range(0, 35), Range(0, 256));
	imshow("Background2", backG);
	hist = Hist(backG);
	histIm = HistImg(hist, 300, 512);
	imshow("Camera histogram2", histIm);
	I_O = AdaptiveLocalNoiseReduction(I_O, 7, 7, 150);
	imshow("Adaptive Local Noise Reduction", I_O);*/

	Mat I1 = imread("Images\\lenaD3.pgm", 0);
	/*imshow("Original", I1);
	Mat I2 = ArithmeticMeanFilter(I1, 3, 3);
	imshow("Arithmetic Mean", I2);
	Mat I3 = GeometricMeanFilter(I1, 3, 3);
	imshow("Geometric Mean", I3);
	Mat I4 = AlphaTrimmedMeanFilter(I1, 5, 5, 20);
	imshow("AlphaTrimmed Mean1", I4);
	Mat I5 = AlphaTrimmedMeanFilter(I1, 5, 5, 10);
	imshow("AlphaTrimmed Mean2", I5);
	Mat I6 = HarmonicMeanFilter(I1, 3, 3);
	imshow("Harmonic Mean", I6);
	Mat I7 = ContraharmonicMeanFilter(I1, 3, 3, -2);
	imshow("Contraharmonic Mean1", I7);
	I7 = ContraharmonicMeanFilter(I1, 3, 3, 2);
	imshow("Contraharmonic Mean2", I7);
	Mat I8 = MedianFilter(I1, 3, 3);
	imshow("Median Mean", I8);
	Mat I9 = MinFilter(I1, 3, 3);
	imshow("Min Filter", I9);
	Mat I10 = MaxFilter(I1, 3, 3);
	imshow("Max Filter", I10);
	Mat I11 = MidpointFilter(I1, 3, 3);
	imshow("Midpoint Filter", I11);*/
	Mat I12 = AdaptiveMedianFilter(I1, 3, 3, 7);
	imshow("Adaptive Median", I12);
	Mat I13 = AdaptiveLocalNoiseReduction(I12, 7, 7, 200);
	imshow("Adaptive Local", I13);

	waitKey(0);
}

/*@param srcImage is the source Image.
@param m is the row size of the window.
@param n is the col size of the window.*/
Mat ArithmeticMeanFilter(Mat srcImage, int m, int n) {
	srcImage.convertTo(srcImage, CV_32F);
	float *pdstData;
	int dm = m / 2;
	int dn = n / 2;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	Scalar S = mean(srcImage);
	Mat dstImage;
	copyMakeBorder(srcImage, srcImage, dm, dm, dn, dn, BORDER_CONSTANT, S(0));
	srcImage.copyTo(dstImage);
	for (int i = dm; i < rows + dm; i++) {
		pdstData = dstImage.ptr<float>(i);
		for (int j = dn; j < rows + dn; j++) {
			S = mean(srcImage(Range(i - dm, i + dm + 1), Range(j - dn, j + dn + 1)));
			pdstData[j] = S(0);
		}
	}
	dstImage = dstImage(Range(dm, rows + dm), Range(dn, cols + dn));
	dstImage.convertTo(dstImage, CV_8U);
	return dstImage;
}

/*@param srcImage is the source Image.
@param m is the row size of the window.
@param n is the col size of the window.*/
Mat GeometricMeanFilter(Mat srcImage, int m, int n) {
	float *pdstData, *psrcData;
	int dm = m / 2;
	int dn = n / 2;
	float mn = 1 / (float)(m*n);
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	Scalar S = mean(srcImage);
	int i, j, p, q;
	copyMakeBorder(srcImage, srcImage, dm, dm, dn, dn, BORDER_CONSTANT, (uchar)S(0));
	srcImage += Scalar::all(1);
	srcImage.convertTo(srcImage, CV_32F);
	Mat dstImage(rows + 2 * dm, cols + 2 * dn, CV_32F, Scalar::all(1));
	for (i = dm; i < rows + dm; i++) {
		pdstData = dstImage.ptr<float>(i);
		for (j = dn; j < rows + dn; j++) {
			for (p = 0; p < m; p++) {
				psrcData = srcImage.ptr<float>(i - dm + p);
				for (q = 0; q < n; q++) {
					pdstData[j] = pdstData[j] * pow(psrcData[j - dn + q], mn);
				}
			}
		}
	}
	dstImage = dstImage(Range(dm, rows + dm), Range(dn, cols + dn));
	dstImage.convertTo(dstImage, CV_8U);
	dstImage -= Scalar::all(1);
	return dstImage;
}

/*@param srcImage is the source Image.
@param m is the row size of the window.
@param n is the col size of the window.
@param d is the trimmed parameter and d is even.*/
Mat AlphaTrimmedMeanFilter(Mat srcImage, int m, int n, int d) {
	uchar *pdstData;
	if (d % 2 == 1) {
		d = d - 1;
	}
	int dm = m / 2;
	int dn = n / 2;
	int hd = d / 2;
	int num = m*n - hd;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int i, j;
	Mat dstImage, mask, sequency, s;
	srcImage.copyTo(dstImage);
	for (i = dm; i < rows - dm; i++) {
		pdstData = dstImage.ptr<uchar>(i);
		for (j = dn; j < rows - dn; j++) {
			mask = srcImage(Range(i - dm, i + dm + 1), Range(j - dn, j + dn + 1));
			mask.copyTo(sequency);
			sequency.reshape(0, 1).copyTo(sequency);
			cv::sort(sequency, sequency, CV_SORT_EVERY_ROW + CV_SORT_ASCENDING);
			s = sequency(Range(0, 1), Range(hd, num));
			pdstData[j] = (uchar)mean(s)(0);
		}
	}
	return dstImage;
}

/*@param srcImage is the source Image.
@param m is the row size of the window.
@param n is the col size of the window.*/
Mat HarmonicMeanFilter(Mat srcImage, int m, int n) {
	srcImage.convertTo(srcImage, CV_32F);
	float *pdstData;
	int dm = m / 2;
	int dn = n / 2;
	int mn = m*n;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	Scalar S = mean(srcImage);
	Mat dstImage;
	copyMakeBorder(srcImage, srcImage, dm, dm, dn, dn, BORDER_CONSTANT, S(0));
	srcImage.copyTo(dstImage);
	for (int i = dm; i < rows + dm; i++) {
		pdstData = dstImage.ptr<float>(i);
		for (int j = dn; j < rows + dn; j++) {
			S = sum(1/srcImage(Range(i - dm, i + dm + 1), Range(j - dn, j + dn + 1)));
			pdstData[j] = mn / S(0);
		}
	}
	dstImage = dstImage(Range(dm, rows + dm), Range(dn, cols + dn));
	dstImage.convertTo(dstImage, CV_8U);
	return dstImage;
}

/*@param srcImage is the source Image.
@param m is the row size of the window.
@param n is the col size of the window.*/
Mat ContraharmonicMeanFilter(Mat srcImage, int m, int n, int Q) {
	srcImage.convertTo(srcImage, CV_32F);
	float *pdstData, *pmask1Data, *pmask2Data;
	int dm = m / 2;
	int dn = n / 2;
	int mn = m*n;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int i, j, p, q;
	Scalar S = mean(srcImage);
	Mat dstImage, mask1, mask2;
	copyMakeBorder(srcImage, srcImage, dm, dm, dn, dn, BORDER_CONSTANT, S(0));
	srcImage.copyTo(dstImage);
	for (i = dm; i < rows + dm; i++) {
		pdstData = dstImage.ptr<float>(i);
		for (j = dn; j < rows + dn; j++) {
			srcImage(Range(i - dm, i + dm + 1), Range(j - dn, j + dn + 1)).copyTo(mask1);
			mask1.copyTo(mask2);
			for (p = 0; p < m; p++) {
				pmask1Data = mask1.ptr<float>(p);
				pmask2Data = mask2.ptr<float>(p);
				for (q = 0; q < n; q++) {
					pmask1Data[q] = pow(pmask1Data[q], Q);
					pmask2Data[q] = pmask1Data[q] * pmask2Data[q];
				}
			}
			pdstData[j] = sum(mask2)(0) / sum(mask1)(0);
		}
	}
	dstImage = dstImage(Range(dm, rows + dm), Range(dn, cols + dn));
	dstImage.convertTo(dstImage, CV_8U);
	return dstImage;
}

/*@param srcImage is the source Image.
@param m is the row size of the window.
@param n is the col size of the window.*/
Mat MedianFilter(Mat srcImage, int m, int n) {
	uchar *pdstData;
	int dm = m / 2;
	int dn = n / 2;
	int med = (m*n - 1) / 2;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int i, j;
	Mat dstImage, mask, sequency;
	srcImage.copyTo(dstImage);
	for (i = dm; i < rows - dm; i++) {
		pdstData = dstImage.ptr<uchar>(i);
		for (j = dn; j < rows - dn; j++) {
			mask = srcImage(Range(i - dm, i + dm + 1), Range(j - dn, j + dn + 1));
			mask.copyTo(sequency);
			sequency.reshape(0, 1).copyTo(sequency);
			cv::sort(sequency, sequency, CV_SORT_EVERY_ROW + CV_SORT_ASCENDING);
			pdstData[j] = sequency.at<uchar>(0, med);
		}
	}
	return dstImage;
}

/*@param srcImage is the source Image.
@param m is the row size of the window.
@param n is the col size of the window.*/
Mat MinFilter(Mat srcImage, int m, int n) {
	uchar *pdstData;
	int dm = m / 2;
	int dn = n / 2;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int i, j;
	Mat dstImage, mask, sequency;
	srcImage.copyTo(dstImage);
	for (i = dm; i < rows - dm; i++) {
		pdstData = dstImage.ptr<uchar>(i);
		for (j = dn; j < rows - dn; j++) {
			mask = srcImage(Range(i - dm, i + dm + 1), Range(j - dn, j + dn + 1));
			mask.copyTo(sequency);
			sequency.reshape(0, 1).copyTo(sequency);
			cv::sort(sequency, sequency, CV_SORT_EVERY_ROW + CV_SORT_ASCENDING);
			pdstData[j] = sequency.at<uchar>(0, 0);
		}
	}
	return dstImage;
}

/*@param srcImage is the source Image.
@param m is the row size of the window.
@param n is the col size of the window.*/
Mat MaxFilter(Mat srcImage, int m, int n) {
	uchar *pdstData;
	int dm = m / 2;
	int dn = n / 2;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int i, j;
	Mat dstImage, mask, sequency;
	srcImage.copyTo(dstImage);
	for (i = dm; i < rows - dm; i++) {
		pdstData = dstImage.ptr<uchar>(i);
		for (j = dn; j < rows - dn; j++) {
			mask = srcImage(Range(i - dm, i + dm + 1), Range(j - dn, j + dn + 1));
			mask.copyTo(sequency);
			sequency.reshape(0, 1).copyTo(sequency);
			cv::sort(sequency, sequency, CV_SORT_EVERY_ROW + CV_SORT_DESCENDING);
			pdstData[j] = sequency.at<uchar>(0, 0);
		}
	}
	return dstImage;
}

/*@param srcImage is the source Image.
@param m is the row size of the window.
@param n is the col size of the window.*/
Mat MidpointFilter(Mat srcImage, int m, int n) {
	uchar *pdstData;
	int dm = m / 2;
	int dn = n / 2;
	int mn = m*n - 1;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int i, j;
	Mat dstImage, mask, sequency;
	srcImage.copyTo(dstImage);
	for (i = dm; i < rows - dm; i++) {
		pdstData = dstImage.ptr<uchar>(i);
		for (j = dn; j < rows - dn; j++) {
			mask = srcImage(Range(i - dm, i + dm + 1), Range(j - dn, j + dn + 1));
			mask.copyTo(sequency);
			sequency.reshape(0, 1).copyTo(sequency);
			cv::sort(sequency, sequency, CV_SORT_EVERY_ROW + CV_SORT_DESCENDING);
			pdstData[j] = ((float)sequency.at<uchar>(0, 0) + (float)sequency.at<uchar>(0, mn))/2.0;
		}
	}
	return dstImage;
}

/*@param srcImage is the source Image.
@param m is the row size of the window.
@param n is the col size of the window.
@param Smax is the max size of the window.*/
Mat AdaptiveMedianFilter(Mat srcImage, int m, int n, int Smax) {
	uchar *pdstData, *psrcData;
	int dm, dn;
	int dS = Smax / 2;
	int mn = m*n - 1;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int i, j;
	int zmin, zmax, zmed, zxy;
	int A1, A2, B1, B2;
	Mat dstImage, mask, sequency;
	Scalar S = mean(srcImage);
	copyMakeBorder(srcImage, srcImage, dS, dS, dS, dS, BORDER_CONSTANT, S(0));
	srcImage.copyTo(dstImage);
	for (i = dS; i < rows + dS; i++) {
		pdstData = dstImage.ptr<uchar>(i);
		psrcData = srcImage.ptr<uchar>(i);
		for (j = dS; j < rows + dS; j++) {
			dm = m;
			dn = n;
			zxy = psrcData[j];
			while (max(dm, dn) <= Smax) {
				mask = srcImage(Range(i - dm / 2, i + dm / 2 + 1), Range(j - dn / 2, j + dn / 2 + 1));
				mask.copyTo(sequency);
				sequency.reshape(0, 1).copyTo(sequency);
				cv::sort(sequency, sequency, CV_SORT_EVERY_ROW + CV_SORT_DESCENDING);
				zmin = sequency.at<uchar>(0, dm*dn - 1);
				zmax = sequency.at<uchar>(0, 0);
				zmed = sequency.at<uchar>(0, (dm*dn - 1) / 2);
				A1 = zmed - zmin;
				A2 = zmed - zmax;
				if (A1 > 0 && A2 < 0) {
					B1 = zxy - zmin;
					B2 = zxy - zmax;
					if (B1 > 0 && B2 < 0) {
						pdstData[j] = zxy;
					}
					else
					{
						pdstData[j] = zmed;
					}
					break;
				}
				else {
					dm += 1;
					dn += 1;
					if (max(dm, dn) > Smax) {
						pdstData[j] = zmed;
					}
				}
			}
		}
	}
	dstImage = dstImage(Range(dS, rows + dS), Range(dS, cols + dS));
	return dstImage;
}

/*@param srcImage is the source Image.
@param m is the row size of the window.
@param n is the col size of the window.
@param nvar is the noise variance.*/
Mat AdaptiveLocalNoiseReduction(Mat srcImage, int m, int n, float nvar) {
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int i, j;
	int md = m / 2, nd = n / 2;
	Mat mask;
	Scalar stddev, ml;
	float lvar, lmean;
	uchar temp;
	Mat dstImage;
	srcImage.convertTo(srcImage, CV_32F);
	srcImage.copyTo(dstImage);
	float *pdstData, *psrcData;
	for (i = md; i < rows - md; i++) {
		pdstData = dstImage.ptr<float>(i);
		psrcData = srcImage.ptr<float>(i);
		for (j = nd; j < cols - nd; j++) {
			mask = srcImage(Rect(j - nd, i - md, n, m));
			//cout << mask << endl;
			meanStdDev(mask, ml, stddev);
			lvar = stddev(0)*stddev(0);
			//nvar = lvar;
			lmean = ml(0);
			//cout << stddev << " " << ml << endl;
			temp = psrcData[j];
			pdstData[j] = temp - (nvar / lvar)*(temp - lmean);
		}
	}
	normalize(dstImage, dstImage, 0, 255, NORM_MINMAX);
	dstImage.convertTo(dstImage, CV_8U);
	return dstImage;
}

/*@brief Creat the Horizontal Notch Filter.
@param rows is the row of the filter.
@param cols is the col of the filter.
@param Cloc is the Vertical Notch start col.
@param CW is the thickness of the notch.
@param Rloc is the split start row.
@param RW is the thickness of the split.*/
Mat VerticalNotch(Mat srcImage, int Cloc, int CW, int Rloc, int RW) {
	Mat dstImage;
	srcImage.copyTo(dstImage);
	int i, j, k;
	int rows = srcImage.rows;
	int channels = srcImage.channels();
	float *pdstData;
	double minv = 0.0, maxv = 0.0;
	double *minp = &minv, *maxp = &maxv;
	minMaxIdx(srcImage, minp, maxp);
	for (i = 0; i < Rloc - 1; i++) {
		pdstData = dstImage.ptr<float>(i);
		for (j = Cloc; j < Cloc + CW; j++) {
			for (k = 0; k < channels; k++) {
				pdstData[j*channels + k] = 0;
			}
		}
	}
	for (i = Rloc + RW - 1; i < rows; i++) {
		pdstData = dstImage.ptr<float>(i);
		for (j = Cloc; j < Cloc + CW; j++) {
			for (k = 0; k < channels; k++) {
				pdstData[j*channels + k] = 0;
			}
		}
	}
	return dstImage;
}

void LenaNotchFilter() {
	Mat srcImage = imread("Images\\LenaWithNoise.pgm", 0);
	imshow("Original", srcImage);
	Mat complexI = DFT(srcImage);
	Mat comI = DFTshow(complexI);
	imshow("Frequency Domain1", comI);
	Mat complexIm = VerticalNotch(complexI, 243, 3, 235, 20);

	Mat complexIm1;
	complexI.copyTo(complexIm1);
	complexIm1 = complexIm1 - complexIm;

	comI = DFTshow(complexIm);
	imshow("Frequency Domain2", comI);
	complexIm = IDFT(complexIm);
	Mat planes[2];
	split(complexIm, planes);
	srcImage = planes[0];
	srcImage = Move2Center(srcImage);
	normalize(srcImage, srcImage, 0, 255, NORM_MINMAX);
	srcImage.convertTo(srcImage, CV_8U);
	imshow("Notch Lena1", srcImage);

	Mat Filter(srcImage.rows, srcImage.cols, CV_32F, Scalar::all(1));
	Filter = VerticalNotch(Filter, 243, 3, 235, 20);
	imshow("Filter", Filter);

	
	comI = DFTshow(complexIm1);
	imshow("Frequency Domain3", comI);
	complexIm1 = IDFT(complexIm1);
	Mat planes1[2];
	split(complexIm1, planes1);
	Mat Image = planes1[0];
	Image = Move2Center(Image);
	normalize(Image, Image, 0, 255, NORM_MINMAX);
	Image.convertTo(Image, CV_8U);
	imshow("Notch Lena2", Image);

	waitKey(0);
}
