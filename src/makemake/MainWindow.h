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
#include <ExifWidget.h>
#include <FocusStack.h>

#include "Project.h"

class MainWindow : public QMainWindow
{
	Q_OBJECT

	public:
		MainWindow(QWidget *parent = nullptr);
		virtual ~MainWindow();

	public slots:
		void about();
		void addImages();
		void closeProject();
		void imageSelected(QTreeWidgetItem *current, QTreeWidgetItem *previous);
		void groupCriteriaChanged();
		void loadProject();
		void newProject();
		void openProject();
		void removeImage();
		void saveChanges();
		void saveProjectAs();
		void updateImageList();
		void focus();

	private:
		// GUI
		enum AppState {Empty, SetupNotValid, SetupValid, Working};

		void createActions();
		void createObject();
		void createLayout();
		void createConnections();
		bool isValid();
		void setState(AppState state);
		QString getLastBrowsedPath();
		void setLastBrowsedPath(QString path);

	private slots:
		void projectChanged();
		void imageFocused(int group);
		void imageAligned(int group, QString filename);
		void focusFinish();
		void channelChanged();

	private:
		// Menu actions
		QAction *newAction;
		QAction *openAction;
		QAction *closeAction;
		QAction *saveAsAction;

		// Project properties
		Project project;

		// Image list
		QTreeWidget *imageList;
		QPushButton *imageAdd;
		QPushButton *imageRemove;

		QGroupBox *gbChannels;
		QRadioButton *rbOriginal;
		QRadioButton *rbFocused;

		// Tools
		QTabWidget *twTool;

		// Project properties
		QLineEdit *leName;

		// Simple or Batch sorting
		QGroupBox *gbSorting;
		QRadioButton *rbByName;
		QRadioButton *rbSimple;
		QRadioButton *rbBatch;
		QRadioButton *rbBatchByCount;
		QSpinBox *sbGroupCount;

		// Focus operations
		QGroupBox *gbActions;
		QPushButton *pbStack;

		// Image viewer
		LayerGraphicWidget *imageViewer;
		PixelLayer *imageLayer;

		// Exif metadata viewer
		ExifWidget *ewExif;

		QTimer triggerSave;

		FocusStack *focusStacker;
};



#endif /* SRC_MAINWINDOW_H_ */
