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

#include <LayerGraphicWidget.h>
#include <PixelLayer.h>
#include <RoiEditorLayer.h>
#include <HistogramWidget.h>
#include <ExifWidget.h>
#include <ImageMerger.h>
#include <ImageCropper.h>

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
			void mergeImages();
			void removeImage();
			void detectRoi();
			void cropImages();
			void importCalibration();
			void ExportCalibration();
			void clearCalibration();

	private:
		// GUI
		enum AppState {Empty, SetupNotValid, SetupValid, Working};

		void createActions();
		void createObject();
		void createLayout();
		void createConnections();
		void setState(AppState state);

		QString getLastBrowsedPath();
		void setLastBrowsedPath(QString path);
		bool isValid();
		bool isRoiValid();
		bool hasMergedImage();
		bool hasDetectedRoi();

	private slots:
		void saveChanges();
		void about();
		void projectChanged();
		void channelChanged();
		void imageMerged(const QString filename, int group);
		void imageCropped(const QString filename, int group);
		void mergeFinished();
		void cropFinished();
		void updateRoi(QRect roi);
		void updateRoiEditor();

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
		QRadioButton *rbMerged;
		QRadioButton *rbDetection;
		QRadioButton *rbCropped;

		// Tools
		QTabWidget *twTool;

		// Project properties
		QLineEdit *leName;

		// Calibration Inputs
		QGroupBox *gbInputs;
		QPushButton *pbOpenCalibration;
		QPushButton *pbClearCalibration;

		QPushButton *pbMerge;
		QPushButton *pbDetect;
		QPushButton *pbCrop;

		// Calibration outputs
		QGroupBox *gbOutputs;
		QSpinBox *sbRoiX;
		QSpinBox *sbRoiY;
		QSpinBox *sbRoiWidth;
		QSpinBox *sbRoiHeight;
		QPushButton *pbExportXml;

		// Image viewer
		LayerGraphicWidget *imageViewer;
		PixelLayer *imageLayer;
		RoiEditorLayer *roiEditorLayer;

		// Exif metadata viewer
		ExifWidget *ewExif;

		QTimer triggerSave;

		ImageMerger *merger;
		bool isMerging;

		ImageCropper *cropper;
		bool isCropping;
};



#endif /* SRC_MAINWINDOW_H_ */
