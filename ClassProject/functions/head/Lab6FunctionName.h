#include <opencv2\opencv.hpp>

using namespace cv;

void Lab6Main();
void Lab6Q1();
void Lab6Q2();
void Lab6Q3();
void Lab6Q4();
Mat ExtractBoundary(Mat srcImage, Mat element);
bool ImagesEqual(Mat Image1, Mat Image2);
Mat ImageIntersection(Mat Image1, Mat Image2);
Mat ExtractConnectedComponent(Mat srcImage, Mat initArray, Mat element);
void StatisticsConnectedComponent(Mat srcImage, Mat element);
Mat GetOverlappingPart(Mat srcImage, Mat element, int threshold);