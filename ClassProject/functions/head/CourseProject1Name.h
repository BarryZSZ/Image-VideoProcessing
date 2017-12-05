#include <opencv2\opencv.hpp>

using namespace cv;

void CourseProject1Main();
Mat LogTrans(Mat srcImage, float c);
Mat GammaCorrection(Mat srcImage, float Gama, float c);
MatND Hist(Mat srcImage);
Mat HistImg(MatND hist, int hist_w, int hist_h);
float GaussianRand(float start, float end, float sigma, float mu);
int GetThreshold(MatND hist);
Mat GenerateBackground(int rows, int cols, float start, float end, float sigma, float mu);
void CourseProject1Q2();
Mat SpatialFilter(Mat srcImage, Mat mask);
Mat AddRandNoise(Mat srcImage, float p, uchar start, uchar end);
