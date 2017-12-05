#include "head\ClassProjectHead.h"


void Lab4Main() {
	Mat I1 = imread("Images\\fingerprint1.pgm", 0);
	Mat I2 = imread("Images\\fingerprint2.pgm", 0);
	Mat I3 = imread("Images\\bridge.pgm", 0);
	Mat I4 = imread("Images\\goldhill.pgm", 0);
	Mat I5 = imread("Images\\lena.pgm", 0);
	Mat I6 = imread("Images\\template.png", 0);
	Mat I7 = imread("Images\\image.png", 0);

	//The first question
	imshow("fingerprint1", I1);
	Mat I_O1 = IDHP(I1, 40);
	imshow("IDHP D0 = 40 (1)", I_O1);
	I_O1 = OtsuBinarization(I_O1);
	imshow("Thresholding IDHP D0 = 40 (1)", I_O1);

	I_O1 = BHP(I1, 40, 1);
	imshow("BHP D0 = 40 (1)", I_O1);
	I_O1 = OtsuBinarization(I_O1);
	imshow("Thresholding BHP D0 = 40 (1)", I_O1);

	I_O1 = GHP(I1, 40);
	imshow("GHP D0 = 40 (1)", I_O1);
	I_O1 = OtsuBinarization(I_O1);
	imshow("Thresholding GHP D0 = 40 (1)", I_O1);

	imshow("fingerprint2", I2);
	Mat I_O2 = IDHP(I2, 40);
	imshow("IDHP D0 = 40 (2)", I_O2);
	I_O2 = OtsuBinarization(I_O2);
	imshow("Thresholding IDHP D0 = 40 (2)", I_O2);

	I_O2 = BHP(I2, 40, 4);
	imshow("BHP D0 = 40 ()", I_O2);
	I_O2 = OtsuBinarization(I_O2);
	imshow("Thresholding BHP D0 = 40 (2)", I_O2);

	I_O2 = GHP(I2, 40);
	imshow("GHP D0 = 40 (2)", I_O2);
	I_O2 = OtsuBinarization(I_O2);
	imshow("Thresholding GHP D0 = 40 (2)", I_O2);

	//The second question
	imshow("Bridge", I3);
	I_O2 = Homomorphic(I3, 2, 0.25, 1, 200);
	imshow("Bridge Homomorphic", I_O2);
	I_O2 = Homomorphic(I3, 1.5, 0.5, 1, 80);
	imshow("Bridge Homomorphic1", I_O2);
	I_O2 = Homomorphic(I3, 2, 0.25, 1, 100);
	imshow("Bridge Homomorphic2", I_O2);
	I_O2 = Homomorphic(I3, 1.5, 0.5, 1, 100);
	imshow("Bridge Homomorphic3", I_O2);

	imshow("Goldhill", I4);
	Mat I_O3 = Homomorphic(I4, 2, 0.25, 1, 80);
	imshow("Goldhill Homomorphic", I_O3);
	I_O3 = Homomorphic(I4, 1.5, 0.5, 1, 80);
	imshow("Goldhill Homomorphic1", I_O3);
	I_O3 = Homomorphic(I4, 2, 0.25, 1, 100);
	imshow("Goldhill Homomorphic2", I_O3);
	I_O3 = Homomorphic(I4, 1.5, 0.5, 1, 100);
	imshow("Goldhill Homomorphic3", I_O3);

	//The third question
	int D0 = 120;
	float W = 5;
	Mat I_O4 = AddSinN(I5, D0);
	imshow("Added Sin Noise Lena", I_O4);

	Mat I_O5 = IDBandReject(I_O4, D0, 4);
	imshow("IDBandReject", I_O5);
	I_O5 = BBandReject(I_O4, D0, 2, 6);
	imshow("BBandReject", I_O5);
	I_O5 = GBandReject(I_O4, D0, 10);
	imshow("GBandReject", I_O5);

	//The fourth question
	imshow("template", I6);
	imshow("image", I7);
	Mat I8(I7.size(), I6.type(), Scalar::all(100));
	imshow("compare", I8);
	Mat I_O6 = CorrelationFrequency(I6, I7);
	imshow("Corrlation1", I_O6);

	I_O6 = CorrelationFrequency(I6, I8);
	imshow("Corrlation2", I_O6);
	waitKey(0);
}

Mat IDHPF(Mat srcImage, int D0) {
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int channels = srcImage.channels();
	float hrows = ((float)rows - 1) / 2;
	float hcols = ((float)cols - 1) / 2;
	float *psrcData;
	float d;
	float D = D0*D0;
	for (int i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<float>(i);
		for (int j = 0; j < cols; j++) {
			d = pow((float)i - hrows, 2) + pow((float)j - hcols, 2);
			if (d <= D) {
				for (int k = 0; k < channels; k++) {
					psrcData[j*channels + k] = 0;
				}
			}
		}
	}
	return srcImage;
}

Mat IDHP(Mat srcImage, int D0) {
	Mat complexImg = DFT(srcImage);
	complexImg = IDHPF(complexImg, D0);
	complexImg = IDFT(complexImg);
	Mat planes[2];
	split(complexImg, planes);
	srcImage = planes[0];
	srcImage = Move2Center(srcImage);
	normalize(srcImage, srcImage, 0, 255, NORM_MINMAX);
	srcImage.convertTo(srcImage, CV_8U);
	return srcImage;
}

Mat BHPF(Mat srcImage, int D0, int n) {
	float *psrcData;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int channels = srcImage.channels();
	float hrows = ((float)rows - 1) / 2;
	float hcols = ((float)cols - 1) / 2;
	float H, d;
	float D = D0*D0;
	for (int i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<float>(i);
		for (int j = 0; j < cols; j++) {
			d = pow((float)i - hrows, 2) + pow((float)j - hcols, 2);
			H = 1 / (1 + pow(D / d, n));
			for (int k = 0; k < channels; k++) {
				psrcData[j*channels + k] = psrcData[j*channels + k] * H;
			}
		}
	}
	return srcImage;
}

Mat BHP(Mat srcImage, int D0, int n) {
	Mat complexImg = DFT(srcImage);
	complexImg = BHPF(complexImg, D0, n);
	complexImg = IDFT(complexImg);
	Mat planes[2];
	split(complexImg, planes);
	srcImage = planes[0];
	srcImage = Move2Center(srcImage);
	normalize(srcImage, srcImage, 0, 255, NORM_MINMAX);
	srcImage.convertTo(srcImage, CV_8U);
	return srcImage;
}

Mat GHPF(Mat srcImage, int D0) {
	float *psrcData;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int channels = srcImage.channels();
	float hrows = (float)(rows - 1) / 2;
	float hcols = (float)(cols - 1) / 2;
	float H, d;
	float D = D0*D0;
	for (int i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<float>(i);
		for (int j = 0; j < cols; j++) {
			d = pow((float)i - hrows, 2) + pow((float)j - hcols, 2);
			H = 1 - exp(-(d / 2 / D));
			for (int k = 0; k < channels; k++) {
				psrcData[j*channels + k] = psrcData[j*channels + k] * H;
			}
		}
	}
	return srcImage;
}

Mat GHP(Mat srcImage, int D0) {
	Mat complexImg = DFT(srcImage);
	complexImg = GHPF(complexImg, D0);
	complexImg = IDFT(complexImg);
	Mat planes[2];
	split(complexImg, planes);
	srcImage = planes[0];
	srcImage = Move2Center(srcImage);
	normalize(srcImage, srcImage, 0, 255, NORM_MINMAX);
	srcImage.convertTo(srcImage, CV_8U);
	return srcImage;
}

/*@param srcImage is the sourse Image.
@param H is the high frequency, where H > 1.
@param L is the low frequency, where L < 1.
@param c is the rate of exp function.
@param D0 is the distance threshold.*/
Mat HomomorphicF(Mat srcImage, float GH, float GL, float c, int D0) {
	//srcImage.convertTo(srcImage, CV_32F);
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int channels = srcImage.channels();
	float hrows = ((float)rows - 1) / 2;
	float hcols = ((float)cols - 1) / 2;
	float D = (float)(D0*D0) / c;
	float HL = GH - GL;
	float *psrcData;
	float d, H;
	for (int i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<float>(i);
		for (int j = 0; j < cols; j++) {
			d = pow((float)i - hrows, 2) + pow((float)j - hcols, 2);
			H = HL*(1 - exp(-d / D)) + GL;
			for (int k = 0; k < channels; k++) {
				psrcData[j*channels + k] = psrcData[j*channels + k] * H;
			}
		}
	}
	return srcImage;
}

Mat Homomorphic(Mat srcImage, float GH, float GL, float c, int D0) {
	srcImage.convertTo(srcImage, CV_32F);
	srcImage += Scalar::all(1.0);
	log(srcImage, srcImage);
	Mat complexImg = DFT(srcImage);
	complexImg = HomomorphicF(complexImg, GH, GL, c, D0);
	complexImg = IDFT(complexImg);
	Mat planes[2];
	split(complexImg, planes);
	srcImage = planes[0];
	srcImage = Move2Center(srcImage);
	float *psrcData;
	for (int i = 0; i < srcImage.rows; i++) {
		psrcData = srcImage.ptr<float>(i);
		for (int j = 0; j < srcImage.cols; j++) {
			psrcData[j] = exp(psrcData[j]);
		}
	}
	normalize(srcImage, srcImage, 0, 255, NORM_MINMAX);
	srcImage.convertTo(srcImage, CV_8U);
	return srcImage;
}

/*@param srcImage is the source Image.
@param d is the distance of the noise in the spectrum.*/
Mat AddSinN(Mat srcImage, int d) {
	Mat complexImg = DFT(srcImage);
	float center;
	int row = (srcImage.rows + 0.5) / 2;
	int col = (srcImage.cols + 0.5) / 2;
	int d0 = sqrt((float)(d*d / 2) + 0.5);

	Mat planes[2];
	split(complexImg, planes);
	for (int i = 0; i < 2; i++) {
		center = planes[i].at<float>(row, col);
		planes[i].at<float>(row, col + d) = center;
		planes[i].at<float>(row, col - d) = center;
		planes[i].at<float>(row + d, col) = center;
		planes[i].at<float>(row - d, col) = center;
		planes[i].at<float>(row + d0, col + d0) = center;
		planes[i].at<float>(row + d0, col - d0) = center;
		planes[i].at<float>(row - d0, col - d0) = center;
		planes[i].at<float>(row - d0, col + d0) = center;
	}

	merge(planes, 2, complexImg);
	Mat spectrum = DFTshow(complexImg);
	imshow("Spectrum", spectrum);

	complexImg = IDFT(complexImg);
	split(complexImg, planes);
	srcImage = planes[0];
	srcImage = Move2Center(srcImage);
	normalize(srcImage, srcImage, 0, 255, NORM_MINMAX);
	srcImage.convertTo(srcImage, CV_8U);
	return srcImage;
}

Mat IDBandRejectF(Mat srcImage, float D0, float W) {
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int channels = srcImage.channels();
	float hrows = ((float)rows - 1) / 2;
	float hcols = ((float)cols - 1) / 2;
	float *psrcData;
	float d;
	for (int i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<float>(i);
		for (int j = 0; j < cols; j++) {
			d = sqrt(pow((float)i - hrows, 2) + pow((float)j - hcols, 2));
			if (d <= D0 + W / 2 && d >= D0 - W / 2) {
				for (int k = 0; k < channels; k++) {
					psrcData[j*channels + k] = 0;
				}
			}
		}
	}
	return srcImage;
}

Mat IDBandReject(Mat srcImage, float D0, float W) {
	Mat complexImg = DFT(srcImage);
	complexImg = IDBandRejectF(complexImg, D0, W);
	complexImg = IDFT(complexImg);
	Mat planes[2];
	split(complexImg, planes);
	srcImage = planes[0];
	srcImage = Move2Center(srcImage);
	normalize(srcImage, srcImage, 0, 255, NORM_MINMAX);
	srcImage.convertTo(srcImage, CV_8U);

	Mat spectrum(srcImage.rows, srcImage.cols, CV_32F, Scalar::all(1));
	spectrum = IDBandRejectF(spectrum, D0, W);
	normalize(spectrum, spectrum, 0, 255, NORM_MINMAX);
	spectrum.convertTo(spectrum, CV_8U);
	imshow("Ideal BandReject F", spectrum);
	return srcImage;
}

Mat BBandRejectF(Mat srcImage, int D0, int n, float W) {
	float *psrcData;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int channels = srcImage.channels();
	float hrows = ((float)rows - 1) / 2;
	float hcols = ((float)cols - 1) / 2;
	float H, d;
	float D = D0*D0;
	for (int i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<float>(i);
		for (int j = 0; j < cols; j++) {
			d = pow((float)i - hrows, 2) + pow((float)j - hcols, 2);
			H = 1 / (1 + pow(d*W*W / pow(d - D, 2), n));
			for (int k = 0; k < channels; k++) {
				psrcData[j*channels + k] = psrcData[j*channels + k] * H;
			}
		}
	}
	return srcImage;
}

Mat BBandReject(Mat srcImage, int D0, int n, float W) {
	Mat complexImg = DFT(srcImage);
	complexImg = BBandRejectF(complexImg, D0, n, W);
	complexImg = IDFT(complexImg);
	Mat planes[2];
	split(complexImg, planes);
	srcImage = planes[0];
	srcImage = Move2Center(srcImage);
	normalize(srcImage, srcImage, 0, 255, NORM_MINMAX);
	srcImage.convertTo(srcImage, CV_8U);
	Mat spectrum(srcImage.rows, srcImage.cols, CV_32F, Scalar::all(1));
	spectrum = BBandRejectF(spectrum, D0, n, W);
	normalize(spectrum, spectrum, 0, 255, NORM_MINMAX);
	spectrum.convertTo(spectrum, CV_8U);
	imshow("Butter BandReject F", spectrum);
	return srcImage;
}

Mat GBandRejectF(Mat srcImage, int D0, float W) {
	float *psrcData;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int channels = srcImage.channels();
	float hrows = (float)(rows - 1) / 2;
	float hcols = (float)(cols - 1) / 2;
	float H, d;
	float D = D0*D0;
	for (int i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<float>(i);
		for (int j = 0; j < cols; j++) {
			d = pow((float)i - hrows, 2) + pow((float)j - hcols, 2);
			H = 1 - exp(-(pow(d - D, 2) / d / W / W));
			for (int k = 0; k < channels; k++) {
				psrcData[j*channels + k] = psrcData[j*channels + k] * H;
			}
		}
	}
	return srcImage;
}

Mat GBandReject(Mat srcImage, int D0, float W) {
	Mat complexImg = DFT(srcImage);
	complexImg = GBandRejectF(complexImg, D0, W);
	complexImg = IDFT(complexImg);
	Mat planes[2];
	split(complexImg, planes);
	srcImage = planes[0];
	srcImage = Move2Center(srcImage);
	normalize(srcImage, srcImage, 0, 255, NORM_MINMAX);
	srcImage.convertTo(srcImage, CV_8U);

	Mat spectrum(srcImage.rows, srcImage.cols, CV_32F, Scalar::all(1));
	spectrum = GBandRejectF(spectrum, D0, W);
	normalize(spectrum, spectrum, 0, 255, NORM_MINMAX);
	spectrum.convertTo(spectrum, CV_8U);
	imshow("Gaussian BandReject F", spectrum);
	return srcImage;
}

Mat CorrelationFrequency(Mat templateI, Mat srcImage) {
	int srows = srcImage.rows;
	int scols = srcImage.cols;
	int trows = templateI.rows;
	int tcols = templateI.cols;
	int M = 2 * srows;
	int N = 2 * scols;
	if (trows > srows || tcols > scols) {
		cout << " Error: The size of template image is larger than the source image' ..." << endl;
	}
	Mat srcPadded;
	Mat temPadded;
	copyMakeBorder(srcImage, srcPadded, 0, M - srows, 0, N - scols, BORDER_CONSTANT, Scalar::all(0));
	copyMakeBorder(templateI, temPadded, 0, M - trows, 0, N - tcols, BORDER_CONSTANT, Scalar::all(0));
	Mat srcComplex = DFT(srcPadded);
	Mat temComplex = DFT(temPadded);
	Mat srcPlanes[2];
	Mat temPlanes[2];
	split(srcComplex, srcPlanes);
	split(temComplex, temPlanes);
	Mat Planes[2];
	split(srcComplex, Planes);
	float *srcP0, *srcP1;
	float *temP0, *temP1;
	float *PData0, *PData1;
	for (int i = 0; i < M; i++) {
		srcP0 = srcPlanes[0].ptr<float>(i);
		srcP1 = srcPlanes[1].ptr<float>(i);
		temP0 = temPlanes[0].ptr<float>(i);
		temP1 = temPlanes[1].ptr<float>(i);
		PData0 = Planes[0].ptr<float>(i);
		PData1 = Planes[1].ptr<float>(i);
		for (int j = 0; j < N; j++) {
			PData0[j] = srcP0[j] * temP0[j] + srcP1[j] * temP1[j];
			PData1[j] = srcP0[j] * temP0[j] - srcP1[j] * temP1[j];
		}
	}

	merge(Planes, 2, srcComplex);
	srcComplex = IDFT(srcComplex);
	split(srcComplex, Planes);
	srcPadded = Planes[0];
	srcPadded = Move2Center(srcPadded);
	normalize(srcPadded, srcPadded, 0, 255, NORM_MINMAX);
	srcPadded.convertTo(srcPadded, CV_8U);
	Mat dstImage = srcPadded(Rect(0, 0, scols, srows));
	return dstImage;
}