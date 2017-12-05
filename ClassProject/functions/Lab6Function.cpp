#include "head\ClassProjectHead.h"


void Lab6Main() {
	Lab6Q1();
	Lab6Q2();
	Lab6Q3();
	Lab6Q4();
}

void Lab6Q1() {
	Mat I1 = imread("Images\\noisy_fingerprint.pgm", 0);
	Mat I2 = imread("Images\\noisy_rectangle.pgm", 0);
	resize(I2, I2, Size(0,0), 0.25, 0.25, INTER_CUBIC);

	Mat element1 = getStructuringElement(MORPH_ELLIPSE, Size(5, 5));
	Mat element2 = getStructuringElement(MORPH_RECT, Size(5, 5));
	Mat element3 = getStructuringElement(MORPH_CROSS, Size(5, 5));
	cout << element1 << endl;
	cout << element2 << endl;
	cout << element3 << endl;

	imshow("Ellpse", 255*element1);
	imshow("Rect", 255*element2);
	imshow("Cross",255*element3);

	Mat I3, I4, I5, I6, I7, I8, I9, I10, I11, I12, I13, I14;
	dilate(I1, I3, element1);
	dilate(I1, I4, element2);
	dilate(I1, I5, element3);
	erode(I1, I6, element1);
	erode(I1, I7, element2);
	erode(I1, I8, element3);
	morphologyEx(I1, I9, MORPH_OPEN, element1);
	morphologyEx(I1, I10, MORPH_OPEN, element2);
	morphologyEx(I1, I11, MORPH_OPEN, element3);
	morphologyEx(I1, I12, MORPH_CLOSE, element1);
	morphologyEx(I1, I13, MORPH_CLOSE, element2);
	morphologyEx(I1, I14, MORPH_CLOSE, element3);
	imshow("Noisy Fingerprint", I1);
	imshow("Dilate Ellipse1", I3);
	imshow("Dilate Rect1", I4);
	imshow("Dilate Cross1", I5);
	imshow("Erode Ellipse1", I6);
	imshow("Erode Rect1", I7);
	imshow("Erode Cross1", I8);
	imshow("Opening Ellipse1", I9);
	imshow("Opening Rect1", I10);
	imshow("Opening Cross1", I11);
	imshow("Closing Ellipse1", I12);
	imshow("Closing Rect1", I13);
	imshow("Closing Cross1", I14);

	dilate(I2, I3, element1);
	dilate(I2, I4, element2);
	dilate(I2, I5, element3);
	erode(I2, I6, element1);
	erode(I2, I7, element2);
	erode(I2, I8, element3);
	morphologyEx(I2, I9, MORPH_OPEN, element1);
	morphologyEx(I2, I10, MORPH_OPEN, element2);
	morphologyEx(I2, I11, MORPH_OPEN, element3);
	morphologyEx(I2, I12, MORPH_CLOSE, element1);
	morphologyEx(I2, I13, MORPH_CLOSE, element2);
	morphologyEx(I2, I14, MORPH_CLOSE, element3);
	imshow("Noisy Rectangle", I2);
	imshow("Dilate Ellipse2", I3);
	imshow("Dilate Rect2", I4);
	imshow("Dilate Cross2", I5);
	imshow("Erode Ellipse2", I6);
	imshow("Erode Rect2", I7);
	imshow("Erode Cross2", I8);
	imshow("Opening Ellipse2", I9);
	imshow("Opening Rect2", I10);
	imshow("Opening Cross2", I11);
	imshow("Closing Ellipse2", I12);
	imshow("Closing Rect2", I13);
	imshow("Closing Cross2", I14);

}

void Lab6Q2() {
	Mat I1 = imread("Images\\licoln.pgm", 0);
	Mat I2 = imread("Images\\U.pgm", 0);
	resize(I2, I2, Size(0, 0), 0.5, 0.5, INTER_CUBIC);

	Mat element1 = getStructuringElement(MORPH_ELLIPSE, Size(5, 5));
	Mat element2 = getStructuringElement(MORPH_RECT, Size(5, 5));
	Mat element3 = getStructuringElement(MORPH_CROSS, Size(5, 5));

	imshow("Licoln", I1);
	imshow("U", I2);

	Mat I3, I4, I5;
	I3 = ExtractBoundary(I1, element1);
	I4 = ExtractBoundary(I1, element2);
	I5 = ExtractBoundary(I1, element3);
	imshow("Boundary Ellipse1", I3);
	imshow("Boundary Rect1", I4);
	imshow("Boundary Cross1", I5);

	I3 = ExtractBoundary(I2, element1);
	I4 = ExtractBoundary(I2, element2);
	I5 = ExtractBoundary(I2, element3);
	imshow("Boundary Ellipse2", I3);
	imshow("Boundary Rect2", I4);
	imshow("Boundary Cross2", I5);

}

void Lab6Q3() {
	Mat I1 = imread("Images\\connected.pgm", 0);
	imshow("Connected", I1);
	Mat element = getStructuringElement(MORPH_RECT, Size(5, 5));
	StatisticsConnectedComponent(I1, element);
}

void Lab6Q4() {
	Mat I1 = imread("Images\\bubbles_on_black_background.pgm", 0);
	imshow("Original", I1);

	Mat element = getStructuringElement(MORPH_ELLIPSE, Size(3, 3));

	Mat initArray(I1.size(), I1.type(), Scalar::all(0));
	initArray.at<uchar>(0, 0) = 255;

	Mat I2, I3, I4, I5;
	I2 = ExtractConnectedComponent(I1, initArray, element);
	imshow("Boundary", I2);

	I3 = I1 - I2;
	imshow("Exclude Boundary", I3);

	StatisticsConnectedComponent(I3, element);

	I4 = GetOverlappingPart(I3, element, 400);
	imshow("Overlapping Particles", I4);

	I5 = I3 - I4;
	imshow("Nonoverlapping Particles", I5);
	
}

/*@brief Using erode function to get the Boundary of srcImage.
@param srcImage is the source Image.
@param element is the erode structuring element.*/
Mat ExtractBoundary(Mat srcImage, Mat element) {
	Mat dstImage;
	erode(srcImage, dstImage, element);
	dstImage = srcImage - dstImage;
	return dstImage;
}

/*@brief Decide Image1 is whether equal to Image2.*/
bool ImagesEqual(Mat Image1, Mat Image2) {
	Mat temp = Image1 - Image2;
	Scalar diff = sum(temp);
	if (diff == Scalar(0, 0, 0, 0)) {
		temp = Image2 - Image1;
		diff = sum(temp);
		return diff == Scalar(0, 0, 0, 0);
	}
	else {
		return false;
	}
}

/*@brief Compute the Intersection of Image1 and Image2.*/
Mat ImageIntersection(Mat Image1, Mat Image2) {
	if (Image1.size() == Image2.size()) {
		uchar *pIm1Data, *pIm2Data;
		int rows = Image1.rows;
		int colch = Image1.cols*Image1.channels();
		int i, j;
		for (i = 0; i < rows; i++) {
			pIm1Data = Image1.ptr<uchar>(i);
			pIm2Data = Image2.ptr<uchar>(i);
			for (j = 0; j < colch; j++) {
				if (pIm1Data[j] != pIm2Data[j]) {
					pIm1Data[j] = 0;
				}
			}
		}
		return Image1;
	}
	else {
		cout << "Error: Image1 and Image2 have different size..." << endl;
	}
}

/*@brief Extract the connected component of srcImage.*/
Mat ExtractConnectedComponent(Mat srcImage, Mat initArray, Mat element) {
	Mat temp(srcImage.size(), srcImage.type(), Scalar::all(0));
	while (ImagesEqual(temp, initArray) == false) {
		initArray.copyTo(temp);
		dilate(initArray, initArray, element);
		initArray = ImageIntersection(initArray, srcImage);
	}
	return initArray;
}

/*@brief Compute the Statist of connected component.*/
void StatisticsConnectedComponent(Mat srcImage, Mat element) {
	Mat temp;
	srcImage.copyTo(temp);
	Mat InitArray(srcImage.size(), srcImage.type(), Scalar::all(0));
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int i, j, num;
	cout << "\tNumber\t" << "Location\t" << "Numbers" << endl;
	int counter = 0;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			if (temp.at<uchar>(i, j) == 255) {
				counter += 1;
				InitArray.at<uchar>(i, j) = 255;
				InitArray = ExtractConnectedComponent(srcImage, InitArray, element);
				num = sum(InitArray / 255)(0);
				if ((i < 10 && j < 100) || (i < 100 && j < 10)) {
					cout << " \t" << counter << "\t(" << i << ", " << j << ")\t\t" << num << endl;
				}
				else {
					cout << " \t" << counter << "\t(" << i << ", " << j << ")\t" << num << endl;
				}
				temp -= InitArray;
				InitArray = (Scalar::all(0));
			}
		}
	}
}

Mat GetOverlappingPart(Mat srcImage, Mat element, int threshold) {	
	Mat InitArray(srcImage.size(), srcImage.type(), Scalar::all(0));
	Mat temp, dstImage;
	srcImage.copyTo(temp);
	InitArray.copyTo(dstImage);
	int rows = srcImage.rows;
	int cols = srcImage.cols;
	int i, j, num;
	int counter = 0;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			if (temp.at<uchar>(i, j) == 255) {
				InitArray.at<uchar>(i, j) = 255;
				InitArray = ExtractConnectedComponent(srcImage, InitArray, element);
				num = sum(InitArray / 255)(0);
				if (num > threshold) {
					dstImage += InitArray;
				}
				temp -= InitArray;
				InitArray = (Scalar::all(0));
			}
		}
	}
	return dstImage;
}