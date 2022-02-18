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
#include <PaintLayer.h>

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
		enum AppState {Empty, SetupNotValid, SetupValid, Working};

		void createActions();
		void createObject();
		void createLayout();
		void createConnections();
		QString getCurrentBimapFileName();
		QString getCurrentMaskFileName();
		void setState(AppState state);

		QString getLastBrowsedPath();
		void setLastBrowsedPath(QString path);
		bool isValid();

		void createCvMatrices(cv::Mat &I, cv::Mat &mI, cv::Mat &const_map, cv::Mat &const_val);
		void maskRoi(const cv::Mat &I, cv::Mat &mI, const cv::Mat &const_map, const cv::Mat &const_val, cv::Rect roi);

	private slots:
		void about();
		void imageChanged();
		void projectChanged();
		void saveChanges();
		void saveImage();
		void masking();
		void maskingRoi(QRect roi);
		QImage toMask(const QImage& source);

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

		QGroupBox *gbInputs;

		QPushButton *pbMasking;
		QCheckBox *cbRefine;

		// Tools
		QTabWidget *twTool;

		// Project properties
		QLineEdit *leName;

		// Exif metadata viewer
		ExifWidget *ewExif;

		// Image viewer
        LayerGraphicWidget *imageViewer;
        PixelLayer *imageLayer;
        PaintLayer *maskPainterLayer;

        QTimer triggerSave;
        QTimer triggerImageSave;
};



#endif /* SRC_MAINWINDOW_H_ */
