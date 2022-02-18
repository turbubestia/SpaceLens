/*
 * MainWindow.cpp
 *
 *  Created on: 22 sep. 2019
 *      Author: claud
 */

#include <QtWidgets/QtWidgets>

#include <opencv4/opencv2/opencv.hpp>
#include <opencv4/opencv2/xfeatures2d.hpp>
#include <exiv2/exiv2.hpp>

#include <Utils.h>
#include <LayerGraphicWidget.h>
#include <PixelLayer.h>
#include <HistogramWidget.h>
#include <ExifWidget.h>
#include <CameraCalibrationParameters.h>
#include <ImageMerger.h>
#include <ImageCropper.h>

#include "MainWindow.h"
#include "About.h"
#include "Project.h"

#define SETTINGS_FILENAME "ceres.ini"

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	// Resize to a reasonable size
	resize(QGuiApplication::primaryScreen()->availableSize() * 4/5);

	createActions();
	createObject();
	createLayout();
	createConnections();

	setState(Empty);

	triggerSave.setSingleShot(true);
}

MainWindow::~MainWindow() {

}

void MainWindow::about() {
	About dialog(this);
	dialog.show();
	dialog.exec();
}

void MainWindow::addImages() {
	QStringList filenames = QFileDialog::getOpenFileNames(this, "Open Image", getLastBrowsedPath(), "Images (*.jpg)");
	if (filenames.count() == 0) {
	    return;
	}

	setLastBrowsedPath(filenames.at(0));

	for (int i = 0; i < filenames.count(); i++) {
        QString filename = filenames.at(i);
        if (filename.isEmpty()) {
            return;
        }

        // Add image to the project manager
        if (project.addImage(filename)) {
			// Add the entry in the image list
			QFileInfo fileinfo(filename);
			QListWidgetItem *newItem = new QListWidgetItem;
			newItem->setText(QFileInfo(filename).fileName());
			newItem->setIcon(QIcon(QPixmap(":/icons/media/pending.png")));
			imageList->addItem(newItem);
        }
    }

	project.save();

	imageList->setCurrentRow(imageList->count()-1);

	setState(SetupValid);
}

void MainWindow::channelChanged() {
	imageSelected(imageList->currentRow());
}

void MainWindow::clearCalibration() {
	project.calibration = CameraCalibrationParameters();
	imageSelected(imageList->currentRow());
}

void MainWindow::closeProject() {
    project.close();

    imageList->clear();

    rbOriginal->setChecked(true);

    leName->clear();

    sbRoiX->clear();
    sbRoiY->clear();
    sbRoiWidth->clear();
    sbRoiHeight->clear();

    imageLayer->clear();

    setState(Empty);
}

void MainWindow::createActions() {
	QMenu *projectMenu = menuBar()->addMenu("Project");

	newAction = new QAction("New", this);
	openAction = new QAction("Open", this);
	closeAction = new QAction("Close", this);
	saveAsAction = new QAction("Save As", this);

	projectMenu->addAction(newAction);
	projectMenu->addAction(openAction);
	projectMenu->addAction(closeAction);
	projectMenu->addSeparator();
	projectMenu->addAction(saveAsAction);

	connect(newAction, &QAction::triggered, this, &MainWindow::newProject);
	connect(openAction, &QAction::triggered, this, &MainWindow::openProject);
	connect(closeAction, &QAction::triggered, this, &MainWindow::closeProject);
	connect(saveAsAction, &QAction::triggered, this, &MainWindow::saveProjectAs);

	QMenu *helpMenu = menuBar()->addMenu("Help");

	QAction *aboutAction = new QAction("About", this);

	helpMenu->addAction(aboutAction);

	connect(aboutAction, &QAction::triggered, this, &MainWindow::about);
}

void MainWindow::createConnections() {
	connect(&triggerSave, &QTimer::timeout, this, &MainWindow::saveChanges);

	connect(imageAdd, &QPushButton::clicked, this, &MainWindow::addImages);
	connect(imageRemove, &QPushButton::clicked, this, &MainWindow::removeImage);
	connect(imageList, &QListWidget::currentRowChanged, this, &MainWindow::imageSelected);

	connect(roiEditorLayer, &RoiEditorLayer::roiChanged, this, &MainWindow::updateRoi);

	connect(rbOriginal, &QRadioButton::clicked, this, &MainWindow::channelChanged);
	connect(rbMerged, &QRadioButton::clicked, this, &MainWindow::channelChanged);
	connect(rbDetection, &QRadioButton::clicked, this, &MainWindow::channelChanged);
	connect(rbCropped, &QRadioButton::clicked, this, &MainWindow::channelChanged);

	connect(pbMerge, &QPushButton::clicked, this, &MainWindow::mergeImages);
	connect(pbDetect, &QPushButton::clicked, this, &MainWindow::detectRoi);
	connect(pbCrop, &QPushButton::clicked, this, &MainWindow::cropImages);

	connect(leName, &QLineEdit::textEdited, this, &MainWindow::projectChanged);

	connect(pbOpenCalibration, &QPushButton::clicked, this, &MainWindow::importCalibration);
	connect(pbClearCalibration, &QPushButton::clicked, this, &MainWindow::clearCalibration);

	connect(sbRoiX, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::projectChanged);  // @suppress("Invalid arguments")
	connect(sbRoiY, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::projectChanged);  // @suppress("Invalid arguments")
	connect(sbRoiWidth, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::projectChanged);  // @suppress("Invalid arguments")
	connect(sbRoiHeight, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::projectChanged);  // @suppress("Invalid arguments")

	connect(sbRoiX, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::updateRoiEditor);  // @suppress("Invalid arguments")
	connect(sbRoiY, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::updateRoiEditor);  // @suppress("Invalid arguments")
	connect(sbRoiWidth, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::updateRoiEditor);  // @suppress("Invalid arguments")
	connect(sbRoiHeight, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::updateRoiEditor);  // @suppress("Invalid arguments")

	connect(pbExportXml, &QPushButton::clicked, this, &MainWindow::ExportCalibration);

	connect(merger, &ImageMerger::merged, this, &MainWindow::imageMerged);
	connect(merger, &ImageMerger::mergeFinished, this, &MainWindow::mergeFinished);
	connect(cropper, &ImageCropper::cropped, this, &MainWindow::imageCropped);
	connect(cropper, &ImageCropper::finished, this, &MainWindow::cropFinished);
}

void MainWindow::createObject() {
	leName = new QLineEdit;

	imageList = new QListWidget();
	imageAdd = new QPushButton("Add");
	imageRemove = new QPushButton("Remove");

	rbOriginal = new QRadioButton("Original");
	rbMerged = new QRadioButton("Merged");
	rbDetection = new QRadioButton("Detection");
	rbCropped = new QRadioButton("Cropped");

	// ------------------------------------------------------------------------
	// Tools
	twTool = new QTabWidget;

	// ------------------------------------------------------------------------
	// Calibration Inputs
	gbInputs = new QGroupBox("Process Inputs");
	pbOpenCalibration = new QPushButton("Open");
	pbClearCalibration = new QPushButton("Clear");

	// ------------------------------------------------------------------------
	// Action buttons
	pbMerge = new QPushButton("Merge");
	pbDetect = new QPushButton("Detect");
	pbCrop = new QPushButton("Crop");

	pbMerge->setEnabled(false);
	pbDetect->setEnabled(false);
	pbCrop->setEnabled(false);

	// ------------------------------------------------------------------------
	// Calibration Outputs
	gbOutputs = new QGroupBox("Crop Outputs");
	sbRoiX = new QSpinBox;
	sbRoiY = new QSpinBox;
	sbRoiWidth = new QSpinBox;
	sbRoiHeight = new QSpinBox;

	sbRoiX->setMaximum(1000000);
	sbRoiY->setMaximum(1000000);
	sbRoiWidth->setMaximum(1000000);
	sbRoiHeight->setMaximum(1000000);

	pbExportXml = new QPushButton("Export XML");

	rbOriginal->setChecked(true);
	gbOutputs->setEnabled(true);
	pbExportXml->setEnabled(true);

	// ------------------------------------------------------------------------
	// Image viewer
	imageViewer = new LayerGraphicWidget;
	imageLayer = new PixelLayer;
	roiEditorLayer = new RoiEditorLayer;

	imageViewer->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	imageViewer->addLayer(imageLayer);
	imageViewer->addLayer(roiEditorLayer);

	// ------------------------------------------------------------------------
	// Exif viewer
	ewExif = new ExifWidget;

	// ------------------------------------------------------------------------
	// Non-GUI objects
	merger = new ImageMerger(this);
	cropper = new ImageCropper(this);
	isMerging = false;
	isCropping = false;
}

void MainWindow::cropImages() {
	if (!isRoiValid()) {
		return;
	}

	for (int i = 0; i < imageList->count(); i++) {
		QListWidgetItem *item = imageList->item(i);
		item->setIcon(QIcon(QPixmap(":/icons/media/pending.png")));
	}

	isCropping = true;

	int x = sbRoiX->value();
	int y = sbRoiY->value();
	int width = sbRoiWidth->value();
	int height = sbRoiHeight->value();

	project.getProjectDir().mkdir("crop");

	const QStringList &inputFilenames = project.getImagesFilenameList();
	QStringList outputFilenames;
	for (int i = 0; i < inputFilenames.count(); i++) {
		QFileInfo fileinfo(inputFilenames.at(i));
		QString filename = QString("%1/crop/%2_crop.jpg").arg(project.getProjectPath()).arg(fileinfo.baseName());
		outputFilenames.append(filename);
	}

	cropper->setInputFilenames(inputFilenames);
	cropper->setOutputFilenames(outputFilenames);
	cropper->setCalibration(project.calibration);
	cropper->setRoi(QRect(x,y,width,height));
	cropper->start();

	setState(Working);
}

void MainWindow::detectRoi() {
	if (isValid() && hasMergedImage()) {
		setState(Working);

		std::string merged_filename = project.toAbsolutePath("cache/merged.jpg").toStdString();
		cv::Mat merged = cv::imread(merged_filename);

		int cols = merged.cols;
		int rows = merged.rows;

		cv::Mat binary, gray;

		// Use a binary version of the image to detect the contours
		cv::cvtColor(merged, binary, cv::COLOR_BGR2GRAY);
		// Remove random noise
		cv::GaussianBlur(binary, gray, cv::Size(5,5), 0);
		// Bi-mode binarization
		cv::threshold(gray, gray ,0, 255, cv::THRESH_BINARY | cv::THRESH_OTSU);
		cv::bitwise_not(gray, gray);
		// Fill to remove inside noise (white spot due to reflection or bad print quality of the pattern)
		cv::morphologyEx(gray, binary, cv::MORPH_DILATE,
				cv::getStructuringElement( cv::MORPH_ELLIPSE, cv::Size(9,9)), cv::Point(-1,-1), 2);

		// Find contour
		std::vector<std::vector<cv::Point>> contours;
		std::vector<std::vector<cv::Point>> contours_simple;
		std::vector<cv::Vec4i> hierarchy;
		cv::findContours(binary, contours, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE);

		// Contour filtering
		std::vector<cv::Rect> boundRect;

		double cx = cols/2.0;
		double cy = rows/2.0;
		double minDistance = 1e10;

		cv::Rect roi;

		for(int i = 0; i < (int)contours.size(); i++) {
			std::vector<cv::Point> cnt = contours.at(i);
			std::vector<cv::Point> approx_cnt;

			double perimeter = cv::arcLength(cnt, true);
			cv::approxPolyDP(cnt, approx_cnt, 0.005 * perimeter, true);
			cv::Rect recti = cv::boundingRect(approx_cnt);

			contours_simple.push_back(approx_cnt);
			boundRect.push_back(recti);

			// distance of the rectangle to the center of image
			double cxi = recti.x + recti.width/2.0;
			double cyi = recti.y + recti.height/2.0;
			double dx = cxi-cx;
			double dy = cyi-cy;
			double distance = sqrt(dx*dx+dy*dy);

			double areaRoi = recti.width * recti.height;

			// Choose the object closest to the center
			if (distance < minDistance && areaRoi < 0.5*(rows*cols)) {
				minDistance = distance;
				roi = recti;
			}
		}

		// Offset the roi
		roi.x -= 200;
		roi.y -= 200;
		roi.width += 400;
		roi.height += 400;

		cv::Mat detected;
		cv::cvtColor(binary, detected, cv::COLOR_GRAY2BGR);
		for (int i = 0; i < (int) contours_simple.size(); i++) {
			cv::drawContours(detected, contours_simple, i, cv::Scalar(255,0,0), 2, 8);
			cv::rectangle(detected, boundRect[i].tl(), boundRect[i].br(), cv::Scalar(0,255,0), 2, 8, 0);
		}
		cv::rectangle(detected, roi.tl(), roi.br(), cv::Scalar(255,0,255), 2, 8, 0);

		// save the image
		std::vector<int> compression_params;
		compression_params.push_back(cv::IMWRITE_JPEG_QUALITY);
		compression_params.push_back(95);
		std::string detected_filename = project.toAbsolutePath("cache/detection.jpg").toStdString();
		cv::imwrite(detected_filename,detected,compression_params);

		// Set the new ROI in the GUI
		QRect qRoi = QRect(roi.x,roi.y,roi.width,roi.height);
		updateRoi(qRoi);

		setState(SetupValid);
	}
}

void MainWindow::ExportCalibration() {
	QString projectPath = project.getProjectPath();

	QFile xmlCalibration(QString("%1/crop/calibration.xml").arg(projectPath));
	if (!xmlCalibration.open(QIODevice::ReadWrite)) {
		return;
	}

	QString filename;
	filename = project.getImagesFilenameList().at(0);
	filename = QString("%1/crop/%2_crop.jpg").arg(project.getProjectPath()).arg(QFileInfo(filename).baseName());

	QImageReader reader(filename);
	reader.setAutoTransform(true);
	QImage image = reader.read();

	double ccdwidth = project.calibration.ccdwidth * image.width() / (double) project.calibration.width;
	double cx = project.calibration.cx - project.roi.x();
	double cy = project.calibration.cy - project.roi.y();

	QXmlStreamWriter stream(&xmlCalibration);
	stream.setAutoFormatting(true);
	stream.writeStartDocument();
	stream.writeStartElement("calibrations");
	stream.writeStartElement("camera");
	stream.writeAttribute("cameraMaker", ewExif->maker);
	stream.writeAttribute("cameraModel", ewExif->model);
	stream.writeAttribute("lense", ewExif->lense);
	stream.writeAttribute("ccdwidth", QString::number(ccdwidth,'f',6));
	stream.writeAttribute("w", QString::number(image.width()));
	stream.writeAttribute("h", QString::number(image.height()));
	stream.writeAttribute("cx", QString::number(cx,'f',6));
	stream.writeAttribute("cy", QString::number(cy,'f',6));
	stream.writeAttribute("fx", QString::number(project.calibration.fx,'f',6));
	stream.writeAttribute("fy", QString::number(project.calibration.fy,'f',6));
	stream.writeAttribute("k1", "0");
	stream.writeAttribute("k2", "0");
	stream.writeAttribute("p1", "0");
	stream.writeAttribute("p2", "0");
	stream.writeAttribute("k3", "0");
	stream.writeAttribute("skew", "0");
	stream.writeAttribute("name", project.name.replace(' ','_'));
	stream.writeEndElement();
	stream.writeEndDocument();

	xmlCalibration.close();
}

void MainWindow::createLayout() {
	// Image list layout
	QGridLayout *gridlayout1 = new QGridLayout;
	gridlayout1->addWidget(imageList,0,0,1,2);
	gridlayout1->addWidget(imageAdd,1,0);
	gridlayout1->addWidget(imageRemove,1,1);

	gbChannels = new QGroupBox("Channels");
	QGridLayout *gridlayout_channels = new QGridLayout;
	gridlayout_channels->addWidget(rbOriginal, 0, 0);
	gridlayout_channels->addWidget(rbMerged, 0, 1);
	gridlayout_channels->addWidget(rbDetection, 1, 0);
	gridlayout_channels->addWidget(rbCropped, 1, 1);
	gbChannels->setLayout(gridlayout_channels);

	gridlayout1->addWidget(gbChannels, 3,0,1,2);

	// ------------------------------------------------------------------------
	// Calibration Inputs
	QGridLayout *gridlayout2 = new QGridLayout;
	gridlayout2->addWidget(new QLabel("Name"), 0, 0);
	gridlayout2->addWidget(leName, 0, 1, 1, 2);
	gridlayout2->addWidget(new QLabel("Calibration"), 1, 0);
	gridlayout2->addWidget(pbOpenCalibration, 1, 1);
	gridlayout2->addWidget(pbClearCalibration, 1, 2);
	gbInputs->setLayout(gridlayout2);

	// ------------------------------------------------------------------------
	// Action buttons
	QHBoxLayout *boxlayout1 = new QHBoxLayout;
	boxlayout1->addWidget(pbMerge);
	boxlayout1->addWidget(pbDetect);
	boxlayout1->addWidget(pbCrop);

	// ------------------------------------------------------------------------
	// Calibration Outputs
	QGridLayout *gridlayout3 = new QGridLayout;
	gridlayout3->addWidget(new QLabel("ROI Position"), 0,0);
	gridlayout3->addWidget(sbRoiX, 0,1);
	gridlayout3->addWidget(sbRoiY, 0,2);
	gridlayout3->addWidget(new QLabel("ROI Size"), 1,0);
	gridlayout3->addWidget(sbRoiWidth, 1, 1);
	gridlayout3->addWidget(sbRoiHeight, 1, 2);
	gridlayout3->addWidget(pbExportXml, 2,0, 1, 3);
	gbOutputs->setLayout(gridlayout3);

	// ------------------------------------------------------------------------
	// Tools panel widget
	QVBoxLayout *vboxlayout1 = new QVBoxLayout;
	vboxlayout1->addWidget(gbInputs);
	vboxlayout1->addLayout(boxlayout1);
	vboxlayout1->addWidget(gbOutputs);
	vboxlayout1->addItem(new QSpacerItem(0,10, QSizePolicy::Expanding, QSizePolicy::Expanding));
	QWidget *wCalibration = new QWidget;
	wCalibration->setLayout(vboxlayout1);
	twTool->addTab(wCalibration, "Calibration");
	twTool->addTab(ewExif, "Exif");

	// ------------------------------------------------------------------------
	// Left panel widget
	QVBoxLayout *vboxlayout2 = new QVBoxLayout;
	vboxlayout2->addLayout(gridlayout1, 1);
	vboxlayout2->addWidget(twTool, 0);
	QWidget *leftPanelWidget = new QWidget;
	leftPanelWidget->setLayout(vboxlayout2);

	// ------------------------------------------------------------------------
	// Join them with a splitter
	QSplitter *splitter = new QSplitter(this);
	splitter->addWidget(leftPanelWidget);
	splitter->addWidget(imageViewer);
	setCentralWidget(splitter);

	// Allow more space for the images
	double sr = 0.20;
	int w1 = width()*sr;
	int w2 = width()*(1-sr);
	splitter->setSizes(QList<int>({w1,w2}));
}

QString MainWindow::getLastBrowsedPath() {
    QString lastPath; // = project.getProjectPath();
    if (lastPath.isEmpty()) {
        lastPath = QDir::homePath();
    }

    QSettings settings(SETTINGS_FILENAME, QSettings::IniFormat);
    settings.beginGroup("application");
    lastPath = settings.value("lastBrowsedPath", lastPath).toString();
    settings.endGroup();

    QFileInfo fileinfo(lastPath);

    return fileinfo.path();
}

bool MainWindow::hasMergedImage() {
	return QFileInfo::exists(project.toAbsolutePath("cache/merged.jpg"));
}

bool MainWindow::hasDetectedRoi() {
	return QFileInfo::exists(project.toAbsolutePath("cache/detection.jpg"));
}

void MainWindow::imageSelected(int index) {

	if (rbMerged->isChecked() || rbDetection->isChecked()) {
		index = -1;
	}

	QImage image;
	QString filename;

	roiEditorLayer->setVisible(true);

	if (index != -1) {
		// The original image is not stored in memory
		filename = project.getImagesFilenameList().at(index);

		ewExif->setImageFilename(filename);

		// Check if the cropped version of the image exist
		QFileInfo fileinfo(filename);
		QString cropped = QString("%1/crop/%2_crop.jpg").arg(project.getProjectPath()).arg(fileinfo.baseName());
		if (rbCropped->isChecked() && QFileInfo::exists(cropped)) {
			filename = cropped;
			roiEditorLayer->setVisible(false);
		}
	} else {
		// Search the selected processed images only is the project is marked as cropped
		if (rbMerged->isChecked() && hasMergedImage()) {
			filename = project.toAbsolutePath("cache/merged.jpg");
		} else if(rbDetection->isChecked() && hasDetectedRoi()) {
			filename = project.toAbsolutePath("cache/detection.jpg");
		}
	}

	if (rbOriginal->isChecked() && project.calibration.isValid) {
		cv::Mat cameraMatrix = project.calibration.cameraMatrix();
		cv::Mat distCoeffs = project.calibration.distCoeffs();
		cv::Mat src = cv::imread(filename.toStdString()), dst;
		undistort(src, dst, cameraMatrix, distCoeffs);
		image = cvMatToQImage(dst);
	} else {
		QImageReader reader(filename);
		reader.setAutoTransform(true);
		image = reader.read();
	}

	// Set the image in the viewer
	imageLayer->setImage(image);
	roiEditorLayer->resize(image.size());
	// Make the canvas size equal to the image size
	imageViewer->resizeCanvasToLayer(0);
	imageViewer->update();
}

void MainWindow::imageMerged(const QString filename, int group) {
	Q_UNUSED(group)

	QFileInfo fileinfo(filename);
	QList<QListWidgetItem*> items = imageList->findItems(fileinfo.baseName(), Qt::MatchContains);
	if (items.count() > 0) {
		QListWidgetItem* item = items.at(0);
		item->setIcon(QIcon(QPixmap(":/icons/media/valid.png")));
	}
}

void MainWindow::imageCropped(const QString filename, int group) {
	Q_UNUSED(group)

	QFileInfo fileinfo(filename);
	QList<QListWidgetItem*> items = imageList->findItems(fileinfo.baseName(), Qt::MatchContains);
	if (items.count() > 0) {
		QListWidgetItem* item = items.at(0);
		item->setIcon(QIcon(QPixmap(":/icons/media/valid.png")));
	}
}

void MainWindow::importCalibration() {
	QString filename = QFileDialog::getOpenFileName(this, "Open File", getLastBrowsedPath(), "Calibration (*.xml)");
	if (filename.isEmpty()) {
		return;
	}

	setLastBrowsedPath(filename);

	QFile xmlCalibration(filename);
	if (!xmlCalibration.open(QIODevice::ReadWrite)) {
		return;
	}

	QXmlStreamReader xml(&xmlCalibration);

	int state = 0;
	while (!xml.atEnd()) {
		xml.readNextStartElement();

		if (xml.isEndElement()) {
			continue;
		}

		switch(state) {
			case 0:
				// <calibrations>
				if (xml.name() == "calibrations") {
					state = 1;
				}
				break;

			case 1:
				// <camera />
				if (xml.name() == "camera") {
					QXmlStreamAttributes attributes = xml.attributes();
					if (attributes.count() > 0) {
						project.calibration.ccdwidth = attributes.value("ccdwidth").toDouble();
						project.calibration.width = attributes.value("w").toInt();
						project.calibration.height = attributes.value("h").toInt();
						project.calibration.cx = attributes.value("cx").toDouble();
						project.calibration.cy = attributes.value("cy").toDouble();
						project.calibration.fx = attributes.value("fx").toDouble();
						project.calibration.fy = attributes.value("fy").toDouble();
						project.calibration.k1 = attributes.value("k1").toDouble();
						project.calibration.k2 = attributes.value("k2").toDouble();
						project.calibration.p1 = attributes.value("p1").toDouble();
						project.calibration.p2 = attributes.value("p2").toDouble();
						project.calibration.k3 = attributes.value("k3").toDouble();
						project.calibration.isValid = true;

						project.save();
					}
				}
				// anything else
				else {
					xml.skipCurrentElement();
				}
				break;
		}
	}

	imageSelected(imageList->currentRow());
}

bool MainWindow::isValid() {
	return project.imageCount() > 0;
}

bool MainWindow::isRoiValid() {
	int width = sbRoiWidth->value();
	int height = sbRoiHeight->value();
	return width > 0 && height > 0;
}

void MainWindow::loadProject() {
	disconnect(&triggerSave, &QTimer::timeout, this, &MainWindow::saveChanges);

	leName->setText(project.name);

	int x = project.roi.x();
	int y = project.roi.y();
	int width = project.roi.width();
	int height = project.roi.height();
	updateRoi(QRect(x,y,width,height));

	const QStringList &filenames = project.getImagesFilenameList();
	for (int i = 0; i < filenames.count(); i++) {
		QListWidgetItem *newItem = new QListWidgetItem;
		QFileInfo fileinfo(filenames.at(i));
		newItem->setText(fileinfo.fileName());
		newItem->setIcon(QIcon(QPixmap(":/icons/media/pending.png")));
		imageList->addItem(newItem);
	}

	connect(&triggerSave, &QTimer::timeout, this, &MainWindow::saveChanges);
}

void MainWindow::mergeImages() {
	for (int i = 0; i < imageList->count(); i++) {
		QListWidgetItem *item = imageList->item(i);
		item->setIcon(QIcon(QPixmap(":/icons/media/pending.png")));
	}

	isMerging = true;

	merger->setImageFilenames(project.getImagesFilenameList());
	merger->setCalibration(project.calibration);
	merger->startMerge();

	setState(Working);
}

void MainWindow::mergeFinished() {
	isMerging = false;

	QImage image = merger->getMergedImage();

	QString merged_filename = project.toAbsolutePath("cache/merged.jpg");
	QFileInfo merged_fileinfo(merged_filename);
	project.getProjectDir().mkpath(merged_fileinfo.path());

	QImageWriter writer(merged_filename);
	writer.setQuality(95);
	writer.write(image);

	// Set the image in the viewer
	imageLayer->setImage(image);
	// Make the canvas size equal to the image size
	imageViewer->resizeCanvasToLayer(0);
	imageViewer->update();

	setState(SetupValid);
}

void MainWindow::cropFinished() {
	isCropping = false;

	rbCropped->setChecked(true);
	imageSelected(imageList->currentRow());

	if (project.calibration.isValid) {
		ExportCalibration();
	}

	setState(SetupValid);
}

void MainWindow::newProject() {

    QString filename = QFileDialog::getSaveFileName(this, "Open File", getLastBrowsedPath(), "Project (*.ceres)");
    if (filename.isEmpty()) {
        return;
    }

    setLastBrowsedPath(filename);

    project.loadDefault();
    project.setFilename(filename);
    project.save();

    loadProject();

    setState(SetupNotValid);
}

void MainWindow::openProject() {
    QString filename = QFileDialog::getOpenFileName(this, "Open File", getLastBrowsedPath(), "Project (*.ceres)");
    if (filename.isEmpty()) {
        return;
    }

    setLastBrowsedPath(filename);

    closeProject();
    project.open(filename);
    loadProject();

    if (isValid()) {
    	setState(SetupValid);
    } else {
    	setState(SetupNotValid);
    }
}

void MainWindow::projectChanged() {
	triggerSave.start();
}

void MainWindow::removeImage() {
	int index = imageList->currentRow();

	if (index != -1) {
		project.removeImage(index);

		delete imageList->item(index);
		if (index == imageList->count()) {
			index = imageList->count()-1;
		}
		imageList->setCurrentRow(index);

		project.save();
	}
}

void MainWindow::setState(AppState state) {

	switch (state) {
		case Empty:
			// actions
			newAction->setEnabled(true);
			openAction->setEnabled(true);
			closeAction->setEnabled(false);
			saveAsAction->setEnabled(false);
			// Images
			imageList->setEnabled(false);
			imageAdd->setEnabled(false);
			imageRemove->setEnabled(false);
			// channels
			gbChannels->setEnabled(false);
			rbOriginal->setChecked(true);
			// Tools
			twTool->setEnabled(false);
			break;

		case SetupNotValid:
			// actions
			newAction->setEnabled(true);
			openAction->setEnabled(true);
			closeAction->setEnabled(true);
			saveAsAction->setEnabled(true);
			// Images
			imageList->setEnabled(true);
			imageAdd->setEnabled(true);
			imageRemove->setEnabled(true);
			// channels
			gbChannels->setEnabled(true);
			rbOriginal->setEnabled(true);
			rbMerged->setEnabled(false);
			rbDetection->setEnabled(false);
			rbCropped->setEnabled(false);
			// Tools
			twTool->setEnabled(true);
			pbMerge->setEnabled(false);
			pbExportXml->setEnabled(false);
			break;

		case SetupValid:
			// actions
			newAction->setEnabled(true);
			openAction->setEnabled(true);
			closeAction->setEnabled(true);
			saveAsAction->setEnabled(true);
			// Images
			imageList->setEnabled(true);
			imageAdd->setEnabled(true);
			imageRemove->setEnabled(true);
			// channels
			gbChannels->setEnabled(true);
			rbOriginal->setEnabled(true);
			rbMerged->setEnabled(hasMergedImage());
			rbDetection->setEnabled(hasDetectedRoi());
			rbCropped->setEnabled(true);
			// Tools
			twTool->setEnabled(true);
			pbMerge->setEnabled(true);
			pbDetect->setEnabled(hasMergedImage());
			pbCrop->setEnabled(isRoiValid());
			pbExportXml->setEnabled(project.calibration.isValid);
			break;

		case Working:
			// actions
			newAction->setEnabled(false);
			openAction->setEnabled(false);
			closeAction->setEnabled(false);
			saveAsAction->setEnabled(false);
			// Images
			imageList->setEnabled(true);
			imageAdd->setEnabled(false);
			imageRemove->setEnabled(false);
			// channels
			gbChannels->setEnabled(false);
			rbOriginal->setChecked(true);
			// Tools
			twTool->setEnabled(false);
			break;
	}
}

void MainWindow::saveChanges() {
	project.name = leName->text();

	int roiX = sbRoiX->value();
	int roiY = sbRoiY->value();
	int roiWidth = sbRoiWidth->value();
	int roiHeight = sbRoiHeight->value();
	project.roi.setRect(roiX, roiY, roiWidth, roiHeight);

	project.save();
}

void MainWindow::saveProjectAs() {
    QString filename = QFileDialog::getSaveFileName(this, "Open File", getLastBrowsedPath(), "Project (*.ceres)");
    if (filename.isEmpty()) {
        return;
    }

    setLastBrowsedPath(filename);

    project.setFilename(filename);
    project.save();
}

void MainWindow::setLastBrowsedPath(QString filename) {
    QSettings settings(SETTINGS_FILENAME, QSettings::IniFormat);
    settings.beginGroup("application");
    settings.setValue("lastBrowsedPath", QFileInfo(filename).filePath());
    settings.endGroup();
}

void MainWindow::updateRoi(QRect roi) {
	sbRoiX->setValue(roi.x());
	sbRoiY->setValue(roi.y());
	sbRoiWidth->setValue(roi.width());
	sbRoiHeight->setValue(roi.height());

	roiEditorLayer->setRoi(roi);
}

void MainWindow::updateRoiEditor() {
	int x = sbRoiX->value();
	int y = sbRoiY->value();
	int width = sbRoiWidth->value();
	int height = sbRoiHeight->value();

	roiEditorLayer->setRoi(QRect(x,y,width,height));
}
