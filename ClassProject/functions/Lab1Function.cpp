#include "head\ClassProjectHead.h"

void Lab1Main() {
	Mat im2 = imread("Images\\camera.pgm");
	Mat im4 = imread("Images\\lena.tiff");
	Mat im  = imread("Images\\lena.pgm");
	Mat im1 = imread("Images\\peppers.tiff");
	Mat im3 = imread("Images\\boats.tiff");


	Mat im_o = Mat(im.size(), im.type());
	resize(im1, im_o, im.size(), 0, 0, INTER_CUBIC);
	Mat Pixel_dstim1 = PixelReplication(im, 0.6, 0.8);
	Mat Pixel_dstim2 = PixelReplication(im, 1.2, 1.8);

	Mat Nearest_dstim1 = NearestInterpolation(im, 0.6, 0.8);
	Mat Nearest_dstim2 = NearestInterpolation(im, 1.2, 1.8);

	Mat Bilinear_dstim1 = BilinearInterpolation(im, 0.6, 0.8);
	Mat Bilinear_dstim2 = BilinearInterpolation(im, 1.2, 1.8);

	Mat Bicubic_dstim1 = BicubicInterpolation(im, 0.6, 0.8);
	Mat Bicubic_dstim2 = BicubicInterpolation(im, 1.2, 1.8);


	Mat Negative_dstim1 = NegativeImage(im);
	Mat Negative_dstim2 = NegativeImage(im_o);

	imshow("OriginGray", im);
	imshow("OriginColor", im_o);
	imshow("Origin1", im);
	imshow("Origin2", im4);

	imshow("PixelRepliication1", Pixel_dstim1);
	imshow("PixelRepliication2", Pixel_dstim2);

	imshow("NearestInterpolation1", Nearest_dstim1);
	imshow("NearestInterpolation2", Nearest_dstim2);

	imshow("BilinearInterpolation1", Bilinear_dstim1);
	imshow("BilinearInterpolation2", Bilinear_dstim2);

	imshow("BicubicInterpolation1", Bicubic_dstim1);
	imshow("BicubicInterpolation2", Bicubic_dstim2);


	imshow("NeigetiveGrayImage", Negative_dstim1);
	imshow("NeigetiveColorImage", Negative_dstim2);



	waitKey(0);
}


Mat PixelReplication(Mat srcImage, double fx, double fy)
{
	if (!srcImage.data) {
		cout << "Error: The image is empty..." << endl;
	}
	int srcrow = srcImage.rows;
	int srccol = srcImage.cols;
	int dstrow = srcrow*fy + 0.5;
	int dstcol = srccol*fx + 0.5;
	Mat dstImage = Mat(Size(dstcol, dstrow), srcImage.type(), Scalar::all(0));
	for (int i = 0; i < dstrow; i++)
	{
		int ix = i / fy;
		if (ix > srcrow - 1)
			ix = srcrow - 1;
		for (int j = 0; j < dstcol; j++)
		{
			int jy = j / fx;
			if (jy > srccol - 1)
				jy = srccol - 1;
			dstImage.at<Vec3b>(i, j) = srcImage.at<Vec3b>(ix, jy);
		}
	}
	return dstImage;
}

Mat NearestInterpolation(Mat srcImage, double fx, double fy)
{
	if (!srcImage.data) {
		cout << "Error: The image is empty..." << endl;
	}
	int srcrow = srcImage.rows;
	int srccol = srcImage.cols;
	int dstrow = srcrow*fy + 0.5;
	int dstcol = srccol*fx + 0.5;
	Mat dstImage = Mat(Size(dstcol, dstrow), srcImage.type());
	for (int i = 0; i < dstrow; i++)
	{
		int ix = i / fy + 0.5;
		if (ix > srcrow - 1)
			ix = srcrow - 1;
		for (int j = 0; j < dstcol; j++)
		{
			int jy = j / fx + 0.5;
			if (jy > srccol - 1)
				jy = srccol - 1;
			dstImage.at<Vec3b>(i, j) = srcImage.at<Vec3b>(ix, jy);
		}
	}
	return dstImage;
}

Mat BilinearInterpolation(Mat srcImage, double fx, double fy)
{
	if (!srcImage.data) {
		cout << "Error: The image is empty..." << endl;
	}
	int srcrow = srcImage.rows;
	int srccol = srcImage.cols;
	int dstrow = srcrow*fy + 0.5;
	int dstcol = srccol*fx + 0.5;
	Mat dstImage = Mat(Size(dstcol, dstrow), srcImage.type());
	for (int i = 0; i < dstrow; i++)
	{
		double ix = i / fy;
		double u = ix - (int)ix;
		int ix_l = ix;
		if (ix_l > srcrow - 2)
			ix_l = srcrow - 2;
		int ix_r = ix_l + 1;
		for (int j = 0; j < dstcol; j++)
		{
			double jy = j / fx;
			double v = jy - (int)jy;
			int jy_u = jy;
			if (jy_u > srccol - 2)
				jy_u = srccol - 2;
			int jy_d = jy_u + 1;
			dstImage.at<Vec3b>(i, j) = (1 - v)*((1 - u)*srcImage.at<Vec3b>(ix_l, jy_u) + u*srcImage.at<Vec3b>(ix_r, jy_u)) + v*((1 - u)*srcImage.at<Vec3b>(ix_l, jy_d) + u*srcImage.at<Vec3b>(ix_r, jy_d));
		}
	}
	return dstImage;
}


double Spline(double w)
{
	if (w > -1.5 && w < -0.5)
	{
		return(0.5 * pow(w + 1.5, 2.0));
	}
	else if (w >= -0.5 && w <= 0.5)
	{
		return 3.0 / 4.0 - (w * w);
	}
	else if ((w > 0.5 && w < 1.5))
	{
		return(0.5 * pow(w - 1.5, 2.0));
	}
	return 0.0;
}

Mat BicubicInterpolation(Mat srcImage, double fx, double fy)
{
	if (!srcImage.data) {
		cout << "Error: The image is empty..." << endl;
	}
	int srcrow = srcImage.rows;
	int srccol = srcImage.cols;
	int dstrow = srcrow*fy + 0.5;
	int dstcol = srccol*fx + 0.5;
	Mat dstImage = Mat(Size(dstcol, dstrow), srcImage.type(), Scalar::all(0));
	Vec3b vec;
	double u, v, Sm, Sn;
	int i, j, ix, jy, m, n;
	for (i = 0; i < dstrow; i++)
	{
		ix = i / fy;
		u = (double)(i / fy) - (int)ix;
		if (ix < 1) {
			ix = 1;
		}
		else if (ix > srcrow - 3) {
			ix = srcrow - 3;
		}
		for (j = 0; j < dstcol; j++)
		{
			jy = j / fx;
			v = (double)(j / fx) - (int)jy;
			if (jy < 1) {
				jy = 1;
			}
			else if (jy > srccol - 3) {
				jy = srccol - 3;
			}
			for (m = -1; m < 3; m++)
			{
				Sm = Spline(m - u);
				for (n = -1; n < 3; n++)
				{
					Sn = Spline(n - v);
					dstImage.at<Vec3b>(i, j) += Sm*Sn*srcImage.at<Vec3b>(ix + m, jy + n);
				}
			}
		}
	}
	return dstImage;
}


double Spline1(double w)
{
	if (abs(w) < 1) {
		return 1 - (2 * pow(w, 2)) + pow(abs(w), 3);
	}
	else if (abs(w) >= 1 && abs(w) <= 2) {
		return 4 - 8 * abs(w) + 5 * pow(w, 2) - pow(abs(w), 3);
	}
	else {
		return 0;
	}
}
double Spline2(double w)
{
	if (abs(w) < 1) {
		return 2 / 3 - pow(w, 2) + 0.5*pow(abs(w), 3);
	}
	else if (abs(w) >= 1 && abs(w) <= 2) {
		return pow((2 - abs(w)), 3) / 6;
	}
	else {
		return 0;
	}
}

Mat NegativeImage(Mat srcImage)
{
	if (!srcImage.data) {
		cout << "Error: The image is empty..." << endl;
	}
	int row = srcImage.rows;
	int col = srcImage.cols;
	int i, j;
	Mat dstImage = Mat(srcImage.size(), srcImage.type());
	Vec3b vec = { 255,255,255 };
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			dstImage.at<Vec3b>(i, j) = vec - srcImage.at<Vec3b>(i, j);
		}
	}
	return dstImage;
}