#include "ClassProjectHead.h"


void Lab3Main() {
	Mat I = imread("Images\\lena.pgm", 0);
	Mat I1 = imread("Images\\camera.pgm", 0);
	Mat I2 = imread("Images\\horiz.pgm", 0);

	imshow("Lena", I);
	imshow("Camera", I1);
	imshow("Horiz", I2);

	Mat I_O = DFT(I);
	I_O = DFTshow(I_O);
	imshow("Lena Spectrum", I_O);
	I_O = DFT(I1);
	I_O = DFTshow(I_O);
	imshow("Camera Spectrum", I_O);
	I_O = DFT(I2);
	I_O = DFTshow(I_O);
	imshow("Horiz Spectrum", I_O);

	Mat I_O1 = PhaseAngleReconstruct(I);
	imshow("Phase Reconstruct Lena", I_O1);
	I_O1 = PhaseAngleReconstruct(I1);
	imshow("Phase Reconstruct Camera", I_O1);
	I_O1 = PhaseAngleReconstruct(I2);
	imshow("Phase Reconstruct Horiz", I_O1);

	Mat I_O2 = IDLP(I, 5);
	Mat I_O3 = IDLP(I, 15);
	Mat I_O4 = IDLP(I, 30);
	Mat I_O5 = IDLP(I, 80);
	Mat I_O6 = IDLP(I, 120);
	imshow("IDLP D0 = 5", I_O2);
	imshow("IDLP D0 = 15", I_O3);
	imshow("IDLP D0 = 30", I_O4);
	imshow("IDLP D0 = 80", I_O5);
	imshow("IDLP D0 = 120", I_O6);

	I_O2 = BLP(I, 5, 2);
	I_O3 = BLP(I, 15, 2);
	I_O4 = BLP(I, 30, 2);
	I_O5 = BLP(I, 80, 2);
	I_O6 = BLP(I, 120, 2);
	imshow("BLP D0 = 5", I_O2);
	imshow("BLP D0 = 15", I_O3);
	imshow("BLP D0 = 30", I_O4);
	imshow("BLP D0 = 80", I_O5);
	imshow("BLP D0 = 120", I_O6);

	I_O2 = GLP(I, 5);
	I_O3 = GLP(I, 15);
	I_O4 = GLP(I, 30);
	I_O5 = GLP(I, 80);
	I_O6 = GLP(I, 120);
	imshow("GLP D0 = 5", I_O2);
	imshow("GLP D0 = 15", I_O3);
	imshow("GLP D0 = 30", I_O4);
	imshow("GLP D0 = 80", I_O5);
	imshow("GLP D0 = 120", I_O6);
	waitKey(0);
}

Mat DFT(Mat srcImage) {
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int M = getOptimalDFTSize(rows);
	int N = getOptimalDFTSize(cols);
	Mat padded;
	copyMakeBorder(srcImage, padded, 0, M - rows, 0, N - cols, BORDER_CONSTANT, Scalar::all(0));
	padded = Move2Center(padded);
	Mat planes[2] = { Mat_<float>(padded), Mat::zeros(padded.size(), CV_32F) };
	Mat complexImg;
	merge(planes, 2, complexImg);
	dft(complexImg, complexImg);
	return complexImg;
}

Mat IDFT(Mat srcImage) {
	Mat dstImage;
	idft(srcImage, dstImage, DFT_SCALE | DFT_COMPLEX_OUTPUT);
	return dstImage;
}

Mat DFTshow(Mat complexImg) {
	Mat mag;
	if (complexImg.channels() > 1) {
		Mat planes[2];
		split(complexImg, planes);
		magnitude(planes[0], planes[1], planes[0]);
		mag = planes[0];
	}
	else if (complexImg.channels() == 1) {
		mag = complexImg;
	}
	mag += Scalar::all(1);
	log(mag, mag);
	mag = mag(Rect(0, 0, mag.cols&-2, mag.rows&-2));
	normalize(mag, mag, 0, 255, NORM_MINMAX);
	mag.convertTo(mag, CV_8U);
	return mag;
}

Mat Move2Center(Mat srcImage) {
	Mat dstImage = Mat_<float>(srcImage);
	float *pdstData;
	for (int i = 0; i < srcImage.rows; i++) {
		pdstData = dstImage.ptr<float>(i);
		for (int j = 0; j < srcImage.cols; j++) {
			pdstData[j] = pow(-1.0, i + j)*pdstData[j];
		}
	}
	return dstImage;
}

Mat PhaseAngleReconstruct(Mat srcImage) {
	Mat complexImg = DFT(srcImage);
	Mat planes[2];
	split(complexImg, planes);
	float mag;
	for (int i = 0; i < complexImg.rows; i++) {
		for (int j = 0; j < complexImg.cols; j++) {
			mag = sqrt(pow(planes[0].at<float>(i, j), 2) + pow(planes[1].at<float>(i, j), 2));
			planes[0].at<float>(i, j) /= mag;
			planes[1].at<float>(i, j) /= mag;
		}
	}
	merge(planes, 2, complexImg);
	complexImg = IDFT(complexImg);
	complexImg = DFTshow(complexImg);
	return complexImg;
}

Mat IDLPF(Mat srcImage, int D0) {
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
			if (d > D) {
				for (int k = 0; k < channels; k++) {
					psrcData[j*channels + k] = 0;
				}
			}
		}
	}
	return srcImage;
}

Mat IDLP(Mat srcImage, int D0) {
	Mat complexImg = DFT(srcImage);
	complexImg = IDLPF(complexImg, D0);
	complexImg = IDFT(complexImg);
	Mat planes[2];
	split(complexImg, planes);
	srcImage = planes[0];
	Move2Center(srcImage);
	normalize(srcImage, srcImage, 0, 255, NORM_MINMAX);
	srcImage.convertTo(srcImage, CV_8U);
	return srcImage;
}

Mat BLPF(Mat srcImage, int D0, int n) {
	float *psrcData;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int channels = srcImage.channels();
	cout << channels << " B" << endl;
	float hrows = ((float)rows - 1) / 2;
	float hcols = ((float)cols - 1) / 2;
	float H, d;
	float D = D0*D0;
	n = 2 * n;
	for (int i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<float>(i);
		for (int j = 0; j < cols; j++) {
			d = pow((float)i - hrows, 2) + pow((float)j - hcols, 2);
			H = 1 / (1 + pow(d / D, n));
			for (int k = 0; k < channels; k++) {
				psrcData[j*channels + k] = psrcData[j*channels + k] * H;
			}
		}
	}
	return srcImage;
}

Mat BLP(Mat srcImage, int D0, int n) {
	Mat complexImg = DFT(srcImage);
	complexImg = BLPF(complexImg, D0, n);
	complexImg = IDFT(complexImg);
	Mat planes[2];
	split(complexImg, planes);
	srcImage = planes[0];
	Move2Center(srcImage);
	normalize(srcImage, srcImage, 0, 255, NORM_MINMAX);
	srcImage.convertTo(srcImage, CV_8U);
	return srcImage;
}

Mat GLPF(Mat srcImage, int D0) {
	float *psrcData;
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int channels = srcImage.channels();
	cout << channels << " G" << endl;
	float hrows = (float)(rows - 1) / 2;
	float hcols = (float)(cols - 1) / 2;
	float H, d;
	float D = D0*D0;
	for (int i = 0; i < rows; i++) {
		psrcData = srcImage.ptr<float>(i);
		for (int j = 0; j < cols; j++) {
			d = pow((float)i - hrows, 2) + pow((float)j - hcols, 2);
			H = exp(-(d / 2 / D));
			for (int k = 0; k < channels; k++) {
				psrcData[j*channels + k] = psrcData[j*channels + k] * H;
			}
		}
	}
	return srcImage;
}

Mat GLP(Mat srcImage, int D0) {
	Mat complexImg = DFT(srcImage);
	complexImg = GLPF(complexImg, D0);
	complexImg = IDFT(complexImg);
	Mat planes[2];
	split(complexImg, planes);
	srcImage = planes[0];
	Move2Center(srcImage);
	normalize(srcImage, srcImage, 0, 255, NORM_MINMAX);
	srcImage.convertTo(srcImage, CV_8U);
	return srcImage;
}