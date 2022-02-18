/*
 * ImageMerger.h
 *
 *  Created on: 26 jun. 2020
 *      Author: claud
 */

#ifndef SRC_LIBS_IMAGECROPPER_H_
#define SRC_LIBS_IMAGECROPPER_H_

#include <QtWidgets/QtWidgets>

#include <opencv4/opencv2/opencv.hpp>

#include <CameraCalibrationParameters.h>

class ImageCropperWorker;

class ImageCropper : public QObject {
	Q_OBJECT

	public:
		ImageCropper(QObject *parent = nullptr);
		virtual ~ImageCropper();

		void setInputFilenames(const QStringList &filenames);
		void setOutputFilenames(const QStringList &filenames);
		void setCalibration(const CameraCalibrationParameters &calibration);
		void setRoi(const QRect &roi);

	public slots:
		void start();

	signals:
		void cropped(const QString filename, int group);
		void finished();

	private slots:
		void workerFinished();

	private:
		QStringList inputFilenameList;
		QStringList outputFilenameList;

		CameraCalibrationParameters calibration;

		QList<ImageCropperWorker *> workerList;
		int finishedJobs;
		int workersSpanwned;

		cv::Rect roi;
};

class ImageCropperWorker : public QThread {
	Q_OBJECT

	public:
		ImageCropperWorker(int group, QObject *parent = nullptr);

		int count();
		int getGroup();
		void setInputFilenames(const QStringList &filenames);
		void setOutputFilenames(const QStringList &filenames);
		void setCalibration(const CameraCalibrationParameters &calibration);
		void setRoi(const cv::Rect &roi);

	signals:
		void cropped(const QString filename, int group);
		void workerFinished();

	protected:
		void run() override;

	private:
		int group;
		CameraCalibrationParameters calibration;
		QStringList inputFilenames;
		QStringList outputFilenames;
		cv::Rect roi;
};


#endif /* SRC_LIBS_IMAGECROPPER_H_ */
