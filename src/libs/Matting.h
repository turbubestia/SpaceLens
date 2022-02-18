/*
 * matting.h
 *
 *  Created on: 2 jun. 2020
 *      Author: claud
 */

#ifndef SRC_MATTING_H_
#define SRC_MATTING_H_

#include <opencv4/opencv2/opencv.hpp>

void solveMultiLevelAlpha32FC3(const cv::Mat &I, const cv::Mat &consts_map, const cv::Mat &consts_vals, cv::Mat &alpha, int level);
void solveAlpha32FC3(const cv::Mat &I, const cv::Mat &consts_map, const cv::Mat &consts_vals, cv::Mat &alpha);

#endif /* SRC_MATTING_H_ */
