/*
 * ImageMerger.cpp
 *
 *  Created on: 26 jun. 2020
 *      Author: claud
 */

#include <QtCore/QtCore>

#include <opencv4/opencv2/opencv.hpp>
#include <exiv2/exiv2.hpp>

#include <Utils.h>

#include "ImageCropper.h"

ImageCropper::ImageCropper(QObject *parent)
	: QObject(parent) {

	int nbWorkers = QThread::idealThreadCount();

	for (int i = 0; i < nbWorkers; i++) {
		ImageCropperWorker *worker = new ImageCropperWorker(i);
		connect(worker, &ImageCropperWorker::cropped, this, &ImageCropper::cropped);
		connect(worker, &ImageCropperWorker::workerFinished, this, &ImageCropper::workerFinished);
		workerList.append(worker);
	}

	finishedJobs = 0;
}

ImageCropper::~ImageCropper() {

}

void ImageCropper::setInputFilenames(const QStringList &filenames) {
    inputFilenameList = filenames;
}

void ImageCropper::setOutputFilenames(const QStringList &filenames) {
	outputFilenameList = filenames;
}

void ImageCropper::setCalibration(const CameraCalibrationParameters &_calibration) {
	calibration = _calibration;
}

void ImageCropper::setRoi(const QRect &_roi) {
	roi.x = _roi.x();
	roi.y = _roi.y();
	roi.width = _roi.width();
	roi.height = _roi.height();
}

void ImageCropper::start() {

    int nbWorkers = workerList.count();
    int groupSize = ceil((double) inputFilenameList.count() / nbWorkers);
    if (groupSize < 1) {
        groupSize = 1;
    }

    finishedJobs = 0;
    workersSpanwned = 0;

    for (int i = 0, k = 0; i < nbWorkers && k < inputFilenameList.count(); i++)
    {
        int currentGroupSize = groupSize;
        int remElements = inputFilenameList.count() - k;
        if (currentGroupSize > remElements) {
            currentGroupSize = remElements;
        }

        QStringList groupFilenames = inputFilenameList.mid(k, currentGroupSize);
        QStringList groupOutputFilenames = outputFilenameList.mid(k, currentGroupSize);
        workerList.at(i)->setInputFilenames(groupFilenames);
        workerList.at(i)->setOutputFilenames(groupOutputFilenames);
        workerList.at(i)->setCalibration(calibration);
        workerList.at(i)->setRoi(roi);
        workerList.at(i)->start();

        k += currentGroupSize;

        workersSpanwned++;
    }
}

void ImageCropper::workerFinished() {
    finishedJobs++;
    if (finishedJobs == workersSpanwned) {
        emit finished();
    }
}

ImageCropperWorker::ImageCropperWorker(int _group, QObject *parent)
	: QThread(parent) {
    group = _group;
}

int ImageCropperWorker::count() {
	return inputFilenames.count();
}

int ImageCropperWorker::getGroup() {
	return group;
}

void ImageCropperWorker::setInputFilenames(const QStringList &_filenames) {
	inputFilenames = _filenames;
}

void ImageCropperWorker::setOutputFilenames(const QStringList &filenames) {
	outputFilenames = filenames;
}

void ImageCropperWorker::setCalibration(const CameraCalibrationParameters &_calibration) {
	calibration = _calibration;
}

void ImageCropperWorker::setRoi(const cv::Rect &_roi) {
	roi = _roi;
}

void ImageCropperWorker::run() {

	cv::Mat cameraMatrix = calibration.cameraMatrix();
	cv::Mat distCoeffs = calibration.distCoeffs();

	int len = count();

	for (int i = 0; i < len; i++) {
		std::string inputFilename = inputFilenames.at(i).toStdString();
		std::string outputFilename = outputFilenames.at(i).toStdString();

		// Read source image
		cv::Mat src = cv::imread(inputFilename);
		if (calibration.isValid) {
			undistort(src.clone(), src, cameraMatrix, distCoeffs);
		}

		// Read source image metadata
		Exiv2::Image::AutoPtr srcMetadata = Exiv2::ImageFactory::open(inputFilename);
		srcMetadata->readMetadata();

		// extract the cropped region
		cv::Mat imageCropped = src(roi);

		// save the image
		std::vector<int> compression_params;
		compression_params.push_back(cv::IMWRITE_JPEG_QUALITY);
		compression_params.push_back(95);
		cv::imwrite(outputFilename,imageCropped,compression_params);

		Exiv2::Image::AutoPtr dstMetadata = Exiv2::ImageFactory::open(outputFilename);

		Exiv2::ExifData &exifData = srcMetadata->exifData();
		exifData["Exif.Photo.PixelYDimension"] = Exiv2::LongValue(imageCropped.rows);
		exifData["Exif.Photo.PixelXDimension"] = Exiv2::LongValue(imageCropped.cols);
		dstMetadata->setExifData(exifData);

		dstMetadata->setIptcData(srcMetadata->iptcData());

		Exiv2::XmpData &xmpData = srcMetadata->xmpData();
		xmpData["Xmp.xmp.CreatorTool"] = Exiv2::StringValue("Ceres-TinyLenses");
		QDate date = QDate::currentDate();
		xmpData["Xmp.xmp.CreatorDate"] = Exiv2::DateValue(date.year(), date.month(), date.day());
		dstMetadata->setXmpData(xmpData);

		dstMetadata->writeMetadata();

		emit cropped(inputFilenames.at(i), group);
	}

	emit workerFinished();
}
