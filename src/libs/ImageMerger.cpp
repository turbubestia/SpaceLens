/*
 * ImageMerger.cpp
 *
 *  Created on: 26 jun. 2020
 *      Author: claud
 */

#include <QtCore/QtCore>

#include <opencv4/opencv2/opencv.hpp>

#include <Utils.h>

#include "ImageMerger.h"

ImageMerger::ImageMerger(QObject *parent)
	: QObject(parent) {

	int nbWorkers = QThread::idealThreadCount();

	for (int i = 0; i < nbWorkers; i++) {
		ImageMergerWorker *worker = new ImageMergerWorker(i);
		connect(worker, &ImageMergerWorker::merged, this, &ImageMerger::merged);
		connect(worker, &ImageMergerWorker::workerFinished, this, &ImageMerger::workerFinished);
		workerList.append(worker);
	}

	finishedJobs = 0;
	workersSpanwned = 0;
}

ImageMerger::~ImageMerger() {

}

QImage ImageMerger::getMergedImage() {
	return cvMatToQImage(mMerged);
}

void ImageMerger::setImageFilenames(const QStringList &_filenames) {
    filenameList = _filenames;
}

void ImageMerger::setCalibration(const CameraCalibrationParameters &_calibration) {
	calibration = _calibration;
}

void ImageMerger::startMerge() {

    int nbWorkers = workerList.count();
    int groupSize = ceil((double) filenameList.count() / nbWorkers);
    if (groupSize < 1) {
        groupSize = 1;
    }

    finishedJobs = 0;
    workersSpanwned = 0;

    for (int i = 0, k = 0; i < nbWorkers && k < filenameList.count(); i++)
    {
        int currentGroupSize = groupSize;
        int remElements = filenameList.count() - k;
        if (currentGroupSize > remElements) {
            currentGroupSize = remElements;
        }

        QStringList groupFilenames = filenameList.mid(k, currentGroupSize);
        workerList.at(i)->setFilenames(groupFilenames);
        workerList.at(i)->setCalibration(calibration);
        workerList.at(i)->start();

        k += currentGroupSize;

        workersSpanwned++;
    }
}

void ImageMerger::workerFinished() {
    finishedJobs++;

    if (finishedJobs == workersSpanwned) {

    	// Merge local results
    	cv::Mat src = workerList.at(0)->getMerge();

    	int count = workersSpanwned;
    	int rows = src.rows;
    	int cols = src.cols;

    	mMerged = cv::Mat::ones(src.size(), CV_8UC3);
		for (int row = 0; row < rows; row++) {
			cv::Vec3b *mrg_row_ptr = mMerged.ptr<cv::Vec3b>(row);
			for (int col = 0; col < cols; col++) {
				mrg_row_ptr[col] = cv::Vec3b(255,255,255);
			}
		}

    	// traverse through the results of every worker
    	for (int i = 0; i < count; i++) {

    		src = workerList.at(i)->getMerge();

    		// Progressive average
    		for (int row = 0; row < rows; row++) {
				cv::Vec3b *src_row_ptr = src.ptr<cv::Vec3b>(row);
				cv::Vec3b *mrg_row_ptr = mMerged.ptr<cv::Vec3b>(row);
				for (int col = 0; col < cols; col++) {
					if (src_row_ptr[col][0] < mrg_row_ptr[col][0])
						mrg_row_ptr[col][0] = src_row_ptr[col][0];
					if (src_row_ptr[col][1] < mrg_row_ptr[col][1])
						mrg_row_ptr[col][1] = src_row_ptr[col][1];
					if (src_row_ptr[col][2] < mrg_row_ptr[col][2])
						mrg_row_ptr[col][2] = src_row_ptr[col][2];
				}
			}
    	}

        emit mergeFinished();
    }
}

ImageMergerWorker::ImageMergerWorker(int _group, QObject *parent)
	: QThread(parent) {
    group = _group;
}

int ImageMergerWorker::count() {
	return filenames.count();
}

int ImageMergerWorker::getGroup() {
	return group;
}

void ImageMergerWorker::setFilenames(const QStringList &_filenames) {
	filenames = _filenames;
}

void ImageMergerWorker::setCalibration(const CameraCalibrationParameters &_calibration) {
	calibration = _calibration;
}

const QImage &ImageMergerWorker::getMergedBuffer() {
	return qiMergeBuffer;
}

cv::Mat &ImageMergerWorker::getMerge() {
	return mMerged;
}

void ImageMergerWorker::run() {
	cv::Mat src = cv::imread(filenames.at(0).toStdString());

	cv::Mat cameraMatrix = calibration.cameraMatrix();
	cv::Mat distCoeffs = calibration.distCoeffs();

	if (calibration.isValid) {
		undistort(src.clone(), src, cameraMatrix, distCoeffs);
	}

	int len = count();

	if (len == 1) {
	    mMerged = src;
	    emit merged(filenames.at(0), group);
	}
	else {
        int rows = src.rows;
        int cols = src.cols;

        mMerged = cv::Mat::zeros(src.size(), CV_8UC3);
        for (int row = 0; row < rows; row++) {
            cv::Vec3b *mrg_row_ptr = mMerged.ptr<cv::Vec3b>(row);
            for (int col = 0; col < cols; col++) {
                mrg_row_ptr[col] = cv::Vec3b(255,255,255);
            }
        }

        for (int i = 0; i < len; i++) {
            if (i > 0) {
                src = cv::imread(filenames.at(i).toStdString());
                if (calibration.isValid) {
                    undistort(src.clone(), src, cameraMatrix, distCoeffs);
                }
            }

            // Progressive average
            for (int row = 0; row < rows; row++) {
                cv::Vec3b *src_row_ptr = src.ptr<cv::Vec3b>(row);
                cv::Vec3b *mrg_row_ptr = mMerged.ptr<cv::Vec3b>(row);
                for (int col = 0; col < cols; col++) {
                    if (src_row_ptr[col][0] < mrg_row_ptr[col][0])
                        mrg_row_ptr[col][0] = src_row_ptr[col][0];
                    if (src_row_ptr[col][1] < mrg_row_ptr[col][1])
                        mrg_row_ptr[col][1] = src_row_ptr[col][1];
                    if (src_row_ptr[col][2] < mrg_row_ptr[col][2])
                        mrg_row_ptr[col][2] = src_row_ptr[col][2];
                }
            }

            emit merged(filenames.at(i), group);
        }
	}

	emit workerFinished();
}
