#include <opencv2/opencv.hpp>

using namespace cv;

void Lab2Main();
Mat Translation(Mat srcImage, int dx, int dy);
Mat Rotation(Mat srcImage, double angle);
Mat ShearVertical(Mat srcImage, double s);
Mat ShearHorizontal(Mat srcImage, double s);
Mat AveragingFilter(Mat srcImage, int wsize);
Mat LaplacianOperators(Mat srcImage, float k);
Mat SobelOperators(Mat srcImage, float k);
Mat GammaCorrection(Mat srcImage, float gamma);
Mat GlobalHistogram(Mat srcImage);
Mat LocalHistogram(Mat srcImage, int wsize);
Mat OtusBinarization(Mat srcImage);
