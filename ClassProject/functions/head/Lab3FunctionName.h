#include <opencv2\opencv.hpp>

using namespace cv;

void Lab3Main();
Mat DFT(Mat srcImage);
Mat IDFT(Mat srcImage);
Mat DFTshow(Mat complexImg);
Mat Move2Center(Mat srcImage);
Mat PhaseAngleReconstruct(Mat srcImage);
Mat IDLPF(Mat srcImage, int D0);
Mat IDLP(Mat srcImage, int D0);
Mat BLPF(Mat srcImage, int D0, int n);
Mat BLP(Mat srcImage, int D0, int n);
Mat GLPF(Mat srcImage, int D0);
Mat GLP(Mat srcImage, int D0);