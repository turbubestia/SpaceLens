/*
 * CameraCalibration.h
 *
 *  Created on: Jun 12, 2020
 *      Author: claud
 */

#ifndef SRC_LIBS_CAMERACALIBRATION_H_
#define SRC_LIBS_CAMERACALIBRATION_H_

#include <QtCore/QtCore>

#include <opencv4/opencv2/opencv.hpp>

class CameraCalibrationWorker;

class CameraCalibration : public QObject {
    Q_OBJECT

    public:
        CameraCalibration(QObject *parent = nullptr);
        virtual ~CameraCalibration();

        void setImageFilenames(const QStringList &filenames);
        void setPattern(QSize dim, QSize corners);
        QImage binaryzeImage(const QImage &image);
        QVector<QPointF> getCorners(int index);
        QImage undistorded(int index);
        void getMatrixCoefficients(double *coeff);
        void getDistortionCoefficients(double *coeff);
        void setMatrixCoefficients(double *coeff);
        void setDistortionCoefficients(double *coeff);

    public slots:
        void startCalibration();
        void setCols(double value);
        void setQuadWidth(double value);
        void setQuadHeight(double value);
        void setRows(double value);
        void setSensorWidth(double value);

    signals:
        void imageCalibrated(const QString filename, bool found);
        void calibrationFinished();

    private slots:
        void workerFinished();

    private:
        double quadWidth;
        double quadHeight;
        double rows;
        double cols;
        double sensorWidth;

        QStringList filenameList;

        QList<CameraCalibrationWorker *> workerList;
        int finishedJobs;

        cv::Mat cameraMatrix;
        cv::Mat distCoeffs;
};

class CameraCalibrationWorker : public QThread {
    Q_OBJECT

    public:
        CameraCalibrationWorker(int group, QObject *parent = nullptr);

        int count();
        int getGroup();
        void setup(int rows, int cols);
        void setFilenames(const QStringList &filenames);
        const std::vector<cv::Point2f> &getCorners(int index);

    signals:
        void imageCalibrated(const QString filename, bool found);
        void workerFinished();

    protected:
        void run() override;

    private:
        int group;
        int rows;
        int cols;
        QStringList filenames;
        QVector<std::vector<cv::Point2f>> corners;
};

#endif /* SRC_LIBS_CAMERACALIBRATION_H_ */
