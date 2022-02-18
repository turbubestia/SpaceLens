/*
 * MainWindow.h
 *
 *  Created on: 22 sep. 2019
 *      Author: claud
 */

#ifndef SRC_MAINWINDOW_H_
#define SRC_MAINWINDOW_H_

#include <QtWidgets/QtWidgets>
#include <opencv4/opencv2/opencv.hpp>

#include <ExifWidget.h>
#include <LayerGraphicWidget.h>
#include <PixelLayer.h>
#include <PointArrayLayer.h>
#include <CameraCalibration.h>
#include <HistogramWidget.h>

#include "Project.h"

class MainWindow : public QMainWindow
{
	Q_OBJECT

	public:
		MainWindow(QWidget *parent = nullptr);
		virtual ~MainWindow();

	public slots:
		void addImages();
		void imageSelected(int index);
	    void newProject();
	    void openProject();
	    void closeProject();
	    void saveProjectAs();
	    void loadProject();
	    void removeImage();

	private:
		// GUI
	    enum AppState {NotReady, SetupNotValid, SetupValid, Calibrated};

		void createActions();
		void createObject();
		void createLayout();
		void createConnections();
		void setState(AppState state);

		QString getLastBrowsedPath();
		void setLastBrowsedPath(QString path);
		bool isValid();
		bool isCalibrated();
		void saveXml(QString filename);
		void updateCalibrationResuts();

	private slots:
		void about();
		void calibrate();
		void imageCalibrated(const QString filename, bool found);
		void calibrationFinished();
		void projectChanged();
		void saveChanges();
		void channelChanged();
		void exportXml();

	private:
		// Menu actions
		QAction *newAction;
		QAction *openAction;
		QAction *closeAction;
		QAction *saveAsAction;

		// Project properties
		Project project;

		// Image list
		QListWidget *imageList;
		QPushButton *imageAdd;
		QPushButton *imageRemove;

		QGroupBox *gbChannels;
		QRadioButton *rbOriginal;
		QRadioButton *rbBinary;
		QRadioButton *rbCorners;
		QRadioButton *rbUndistorded;

		// Tools
		QTabWidget *twTool;

		// Project
		QLineEdit *leName;

		// Calibration Inputs
		QGroupBox *gbInputs;
		QDoubleSpinBox *sbSensorWidth;
		QDoubleSpinBox *sbSensorHeight;
		QSpinBox *sbColCorners;
		QSpinBox *sbRowCorners;
		QDoubleSpinBox *sbQuadWidth;
		QDoubleSpinBox *sbQuadHeight;
		QPushButton *pbCalibrate;

		// Calibration outputs
		QGroupBox *gbOutputs;
		QLabel *lcx;
		QLabel *lcy;
		QLabel *lfx;
		QLabel *lfy;
		QPushButton *pbExportXml;

		// Image viewer
		LayerGraphicWidget *imageViewer;
		PixelLayer *imageLayer;
		PointArrayLayer *patternLayer;

		// Exif metadata viewer
		ExifWidget *ewExif;
		HistogramWidget *histo;

		bool calibrationRunning;

		CameraCalibration *calibrator;

		QTimer triggerSave;
};



#endif /* SRC_MAINWINDOW_H_ */
