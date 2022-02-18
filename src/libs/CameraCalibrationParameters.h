/*
 * CameraCalibrationParameters.h
 *
 *  Created on: 30 jun. 2020
 *      Author: claud
 */

#ifndef SRC_LIBS_CAMERACALIBRATIONPARAMETERS_H_
#define SRC_LIBS_CAMERACALIBRATIONPARAMETERS_H_

#include <opencv4/opencv2/opencv.hpp>

struct CameraCalibrationParameters {

	CameraCalibrationParameters() {
		ccdwidth = 0;
		width = height = 0;
		cx = cy = 0;
		fx = fy = 0;
		k1 = k2 = k3 = 0;
		p1 = p2 = 0;
		isValid = false;
	}

	cv::Mat cameraMatrix() {
		cv::Mat _cameraMatrix = cv::Mat::zeros(cv::Size(3,3), CV_64F);
		_cameraMatrix.at<double>(0,2) = cx;
		_cameraMatrix.at<double>(1,2) = cy;
		_cameraMatrix.at<double>(0,0) = fx;
		_cameraMatrix.at<double>(1,1) = fy;
		_cameraMatrix.at<double>(2,2) = 1;
		return _cameraMatrix;
	}

	cv::Mat distCoeffs() {
		cv::Mat _distCoeffs = cv::Mat::zeros(cv::Size(1,5), CV_64F);
		_distCoeffs.at<double>(0,0) = k1;
		_distCoeffs.at<double>(0,1) = k2;
		_distCoeffs.at<double>(0,2) = p1;
		_distCoeffs.at<double>(0,3) = p2;
		_distCoeffs.at<double>(0,4) = k3;
		return _distCoeffs;
	}

	double ccdwidth;
	int width;
	int height;
	double cx, cy;
	double fx, fy;
	double k1,k2,p1,p2,k3;

	bool isValid;
};



#endif /* SRC_LIBS_CAMERACALIBRATIONPARAMETERS_H_ */
