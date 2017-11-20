#include "functions\head\ClassProjectHead.h"

void main() {
	//Lab1Main();
	//Lab2Main();
	//Lab3Main();
	//Lab4Main();
	//Lab5Main();
	CourseProject1Main();



	/*Mat I1 = imread("Images\\LenaWithNoise.pgm", 0);
	Mat I2 = imread("Images\\cameraWithNoise.pgm", 0);
	Mat I3 = DFT(I1);
	Mat I4 = DFTshow(I3);
	imshow("Lena", I4);
	Mat I5 = DFT(I2);
	Mat I6 = DFTshow(I5);
	imshow("cameramen", I6);
	Scalar S = mean(I6);
	Mat a = I6( Range(0, 2), Range(1, 3));
	imshow("test", a);
	cout << a;*/
	/*Mat a(4, 4, CV_32FC1, Scalar::all(0));
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			a.at<float>(i, j) = 4 * i + j;
		}
	}
	Scalar s = sum(a);
	cout << s << endl;*/
	//Mat I1 = imread("Images\\lena.pgm");
	//imshow(" ", I1);
	//Mat b = a.reshape(0, 1);
	//cout << "a: " << a << endl;
	//cout << "b: " << b << endl;
	//Mat c;
	//a.reshape(1, 1).copyTo(c);
	//cv::sort(c, c, CV_SORT_EVERY_ROW+CV_SORT_ASCENDING);
	//cout << "c: " << c << endl;
	
	waitKey(0);
	waitKey(0);

}