#include <opencv2\opencv.hpp>

using namespace cv;

void Lab5Main();
Mat ArithmeticMeanFilter(Mat srcImage, int m, int n);
Mat GeometricMeanFilter(Mat srcImage, int m, int n);
Mat AlphaTrimmedMeanFilter(Mat srcImage, int m, int n, int d);
Mat HarmonicMeanFilter(Mat srcImage, int m, int n);
Mat ContraharmonicMeanFilter(Mat srcImage, int m, int n, int Q);
Mat MedianFilter(Mat srcImage, int m, int n);
Mat MinFilter(Mat srcImage, int m, int n);
Mat MaxFilter(Mat srcImage, int m, int n);
Mat MidpointFilter(Mat srcImage, int m, int n);
Mat AdaptiveMedianFilter(Mat srcImage, int m, int n, int Smax);
Mat AdaptiveLocalNoiseReduction(Mat srcImage, int m, int n, float nvar);
Mat VerticalNotch(Mat srcImage, int Rloc, int RW, int Cloc, int CW);
void LenaNotchFilter();