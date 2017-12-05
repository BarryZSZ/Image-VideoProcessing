#include <opencv2\opencv.hpp>
#include <iostream>

using namespace std;
using namespace cv;

void Lab8Main();
void Lab8Q1();
void Lab8Q2();
void Lab8Q3();
void Lab8Q4();
Mat OtsuBinarization(Mat srcImage);
Mat PartitionOtsu(Mat srcImage, int h, int w);
Mat MovingAverageThreshold(Mat srcImage, int n, float b);
