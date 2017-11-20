#include <opencv2\opencv.hpp>

using namespace cv;

void Lab4Main();
Mat IDHPF(Mat srcImage, int D0);
Mat IDHP(Mat srcImage, int D0);
Mat BHPF(Mat srcImage, int D0, int n);
Mat BHP(Mat srcImage, int D0, int n);
Mat GHPF(Mat srcImage, int D0);
Mat GHP(Mat srcImage, int D0);
Mat HomomorphicF(Mat srcImage, float H, float L, float c, int D0);
Mat Homomorphic(Mat srcImage, float H, float L, float c, int D0);
Mat AddSinN(Mat srcImage, int d);
Mat IDBandRejectF(Mat srcImage, float D0, float W);
Mat IDBandReject(Mat srcImage, float D0, float W);
Mat BBandRejectF(Mat srcImage, int D0, int n, float W);
Mat BBandReject(Mat srcImage, int D0, int n, float W);
Mat GBandRejectF(Mat srcImage, int D0, float W);
Mat GBandReject(Mat srcImage, int D0, float W);
Mat CorrelationFrequency(Mat templateI, Mat srcImage);
