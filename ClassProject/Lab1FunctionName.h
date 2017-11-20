#include <opencv2/opencv.hpp>

using namespace cv;

void Lab1Main();
Mat PixelReplication(Mat srcImage, double fx, double fy);
Mat NearestInterpolation(Mat srcImage, double fx, double fy);
Mat BilinearInterpolation(Mat srcImage, double fx, double fy);
Mat BicubicInterpolation(Mat srcImage, double fx, double fy);
double Spline(double w);
double Spline1(double w);
double Spline2(double w);
Mat NegativeImage(Mat srcImage);
