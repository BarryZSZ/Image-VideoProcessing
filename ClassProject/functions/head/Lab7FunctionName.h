#include <opencv2\opencv.hpp>
#include <iostream>

using namespace std;
using namespace cv;

void Lab7Main();
void Lab7Q1();
void Lab7Q2();
void Lab7Q3();
Mat GlobalThreshold(Mat InputArray, float threshold);
Mat RobertsGradient(Mat inputArray);
Mat PrewittGradient(Mat srcImage);
Mat SobelGradient(Mat srcImage);
Mat GetGaussianWindow(float sigma);
Mat LaplacianOfGaussian(Mat srcImage, float sigma);
Mat CannyEdgeDetection(Mat srcImage, float thresholdLower, float thresholdUpper);
void followEdges(int x, int y, Mat &magnitude, int tUpper, int tLower, Mat &edges);
void edgeDetect(Mat &magnitude, int tUpper, int tLower, Mat &edges);