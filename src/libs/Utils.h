/*
 * Utils.h
 *
 *  Created on: 18 abr. 2020
 *      Author: claud
 */

#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#include <QtWidgets/QtWidgets>
#include <opencv4/opencv2/opencv.hpp>

// OpenCV Interface
QImage cvMatToQImage(const cv::Mat &inMat);
QString cvMatToString(const cv::Mat &M);
QString cvPoint2fListToString(const std::vector<cv::Point2f> &pointList);
int findBimodalThreshold(const QImage &image);
int findBimodalThreshold(const cv::Mat &gray);
cv::Mat qImageToCvMat(const QImage &inImage);
cv::Mat stringToCvMat(QString str);
std::vector<cv::Point2f> stringToCvPoint2fList(QString str);
QString commonPrefix(QStringList strs);


#endif /* SRC_UTILS_H_ */
