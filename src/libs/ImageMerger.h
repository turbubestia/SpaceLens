/*
 * ImageMerger.h
 *
 *  Created on: 26 jun. 2020
 *      Author: claud
 */

#ifndef SRC_LIBS_IMAGEMERGER_H_
#define SRC_LIBS_IMAGEMERGER_H_

#include <QtWidgets/QtWidgets>

#include <opencv4/opencv2/opencv.hpp>

#include <CameraCalibrationParameters.h>

class ImageMergerWorker;

class ImageMerger : public QObject {
	Q_OBJECT

	public:
		ImageMerger(QObject *parent = nullptr);
		virtual ~ImageMerger();

		QImage getMergedImage();
		void setImageFilenames(const QStringList &filenames);
		void setCalibration(const CameraCalibrationParameters &calibration);

	public slots:
		void startMerge();

	signals:
		void merged(const QString filename, int group);
		void mergeFinished();

	private slots:
		void workerFinished();

	private:
		QStringList filenameList;

		CameraCalibrationParameters calibration;

		QList<ImageMergerWorker *> workerList;
		int finishedJobs;
		int workersSpanwned;

		cv::Mat mMerged;
};

class ImageMergerWorker : public QThread {
	Q_OBJECT

	public:
		ImageMergerWorker(int group, QObject *parent = nullptr);

		int count();
		int getGroup();
		void setFilenames(const QStringList &filenames);
		void setCalibration(const CameraCalibrationParameters &calibration);
		const QImage &getMergedBuffer();
		cv::Mat &getMerge();

	signals:
		void merged(const QString filename, int group);
		void workerFinished();

	protected:
		void run() override;

	private:
		int group;

		CameraCalibrationParameters calibration;

		QImage qiMergeBuffer;
		cv::Mat mMerged;

		QStringList filenames;
};


#endif /* SRC_LIBS_IMAGEMERGER_H_ */
