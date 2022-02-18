/*
 * MainWindow.cpp
 *
 *  Created on: 22 sep. 2019
 *      Author: claud
 */

#include <QtWidgets/QtWidgets>
#include <QtConcurrent/QtConcurrent>

#include <opencv4/opencv2/opencv.hpp>

#include <ExifWidget.h>
#include <LayerGraphicWidget.h>
#include <PixelLayer.h>
#include <PointArrayLayer.h>
#include <CameraCalibration.h>
#include <HistogramWidget.h>

#include "MainWindow.h"
#include "Project.h"
#include "About.h"

#define SETTINGS_FILENAME "eris.ini"

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	// Resize to a reasonable size
	resize(QGuiApplication::primaryScreen()->availableSize() * 4/5);

	calibrationRunning = false;

	createActions();
	createObject();
	createLayout();
	createConnections();

	setState(NotReady);
}

MainWindow::~MainWindow() {

}

void MainWindow::about() {
	About dialog(this);
	dialog.show();
	dialog.exec();
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

	connect(rbOriginal, &QRadioButton::clicked, this, &MainWindow::channelChanged);
	connect(rbBinary, &QRadioButton::clicked, this, &MainWindow::channelChanged);
	connect(rbCorners, &QRadioButton::clicked, this, &MainWindow::channelChanged);
	connect(rbUndistorded, &QRadioButton::clicked, this, &MainWindow::channelChanged);

	connect(leName, &QLineEdit::textEdited, this, &MainWindow::projectChanged);

	connect(sbSensorWidth, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MainWindow::projectChanged);  // @suppress("Invalid arguments")
	connect(sbSensorHeight, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MainWindow::projectChanged); // @suppress("Invalid arguments")
	connect(sbColCorners, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::projectChanged);            // @suppress("Invalid arguments")
	connect(sbRowCorners, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::projectChanged);            // @suppress("Invalid arguments")
	connect(sbQuadWidth, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MainWindow::projectChanged);    // @suppress("Invalid arguments")
	connect(sbQuadHeight, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &MainWindow::projectChanged);   // @suppress("Invalid arguments")

	connect(pbCalibrate, &QPushButton::clicked, this, &MainWindow::calibrate);
	connect(calibrator, &CameraCalibration::imageCalibrated, this, &MainWindow::imageCalibrated);
	connect(calibrator, &CameraCalibration::calibrationFinished, this, &MainWindow::calibrationFinished);

	connect(pbExportXml, &QPushButton::clicked, this, &MainWindow::exportXml);
}

void MainWindow::createObject() {
	leName = new QLineEdit;

	imageList = new QListWidget();
	imageAdd = new QPushButton("Add");
	imageRemove = new QPushButton("Remove");

	histo = new HistogramWidget;

	rbOriginal = new QRadioButton("Original");
	rbBinary = new QRadioButton("Binary");
	rbCorners = new QRadioButton("Corners");
	rbUndistorded = new QRadioButton("Undistorded");

	twTool = new QTabWidget;

	// ------------------------------------------------------------------------
	// Calibration Inputs
	gbInputs = new QGroupBox("Calibration Inputs");
	sbSensorWidth = new QDoubleSpinBox;
	sbSensorHeight = new QDoubleSpinBox;
	sbColCorners = new QSpinBox;
	sbRowCorners = new QSpinBox;
	sbQuadWidth = new QDoubleSpinBox;
	sbQuadHeight = new QDoubleSpinBox;

	sbSensorWidth->setValue(10.0);
	sbSensorHeight->setValue(10.0);
	sbColCorners->setValue(7);
	sbRowCorners->setValue(4);
	sbQuadWidth->setValue(7.5);
	sbQuadHeight->setValue(7.5);

	pbCalibrate = new QPushButton("Calibrate");

	// ------------------------------------------------------------------------
	// Calibration Outputs
	gbOutputs = new QGroupBox("Calibration Outputs");

	pbExportXml = new QPushButton("Export XML");

	rbOriginal->setChecked(true);

	lcx = new QLabel("<b>0.0 px</b>");
	lcy = new QLabel("<b>0.0 px</b>");
	lfx = new QLabel("<b>0.0 px</b>");
	lfy = new QLabel("<b>0.0 px</b>");

	// ------------------------------------------------------------------------
	// Image viewer
	imageViewer = new LayerGraphicWidget;
	imageLayer = new PixelLayer;
	patternLayer = new PointArrayLayer;

	imageViewer->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	imageViewer->addLayer(imageLayer);
	imageViewer->addLayer(patternLayer);

	// ------------------------------------------------------------------------
	// Exif viewer
	ewExif = new ExifWidget;

	triggerSave.setSingleShot(true);
	triggerSave.setInterval(2000);

	calibrator = new CameraCalibration(this);
}

void MainWindow::createLayout() {
	// Image list layout
	QGridLayout *gridlayout1 = new QGridLayout;
	gridlayout1->addWidget(imageList,0,0,1,2);
	gridlayout1->addWidget(imageAdd,1,0);
	gridlayout1->addWidget(imageRemove,1,1);
	gridlayout1->addWidget(histo,2,0,1,2);

	gbChannels = new QGroupBox("Channels");
	QHBoxLayout *hboxlayout0 = new QHBoxLayout;
	hboxlayout0->addWidget(rbOriginal);
	hboxlayout0->addWidget(rbBinary);
	hboxlayout0->addWidget(rbCorners);
	hboxlayout0->addWidget(rbUndistorded);
	gbChannels->setLayout(hboxlayout0);

	gridlayout1->addWidget(gbChannels, 3,0,1,2);

	// ------------------------------------------------------------------------
	// Calibration Inputs
	QGridLayout *gridlayout2 = new QGridLayout;
	gridlayout2->addWidget(new QLabel("Name"), 0, 0);
	gridlayout2->addWidget(leName, 0, 1, 1, 2);
	gridlayout2->addWidget(new QLabel("Horizontal"), 1, 1);
	gridlayout2->addWidget(new QLabel("Vertical"), 1, 2);
	gridlayout2->addWidget(new QLabel("Sensor size [mm]"), 2, 0);
	gridlayout2->addWidget(sbSensorWidth, 2, 1);
	gridlayout2->addWidget(sbSensorHeight, 2, 2);
	gridlayout2->addWidget(new QLabel("Inner corners"), 3, 0);
	gridlayout2->addWidget(sbColCorners, 3, 1);
	gridlayout2->addWidget(sbRowCorners, 3, 2);
	gridlayout2->addWidget(new QLabel("Cell size [mm]"), 4, 0);
	gridlayout2->addWidget(sbQuadWidth, 4, 1);
	gridlayout2->addWidget(sbQuadHeight, 4, 2);
	gbInputs->setLayout(gridlayout2);

	// ------------------------------------------------------------------------
	// Action buttons
	QHBoxLayout *boxlayout1 = new QHBoxLayout;
	boxlayout1->addWidget(pbCalibrate);

	// ------------------------------------------------------------------------
	// Calibration Outputs
	QGridLayout *gridlayout3 = new QGridLayout;

	gridlayout3->addWidget(new QLabel("Principal point x"), 1,0,1,2);
	gridlayout3->addWidget(lcx, 1,2,1,2);
	gridlayout3->addWidget(new QLabel("Principal point y"), 2,0,1,2);
	gridlayout3->addWidget(lcy, 2,2,1,2);
	gridlayout3->addWidget(new QLabel("Focal distance x"), 3,0,1,2);
	gridlayout3->addWidget(lfx, 3,2,1,2);
	gridlayout3->addWidget(new QLabel("Focal distance y"), 4,0,1,2);
	gridlayout3->addWidget(lfy, 4,2,1,2);
	gridlayout3->addWidget(pbExportXml, 5, 2);

	gbOutputs->setLayout(gridlayout3);

	// ------------------------------------------------------------------------
	// Tools panel widget
	QVBoxLayout *vboxlayout1 = new QVBoxLayout;
	vboxlayout1->addWidget(gbInputs);
	vboxlayout1->addLayout(boxlayout1);
	vboxlayout1->addWidget(gbOutputs);
	QWidget *wCalibration = new QWidget;
	wCalibration->setLayout(vboxlayout1);
	twTool->addTab(wCalibration, "Calibration");
	twTool->addTab(ewExif, "Exif");

	// ------------------------------------------------------------------------
	// Left panel widget
	QVBoxLayout *vboxlayout2 = new QVBoxLayout;
	vboxlayout2->addLayout(gridlayout1);
	vboxlayout2->addWidget(twTool);
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
        project.addImage(filename);

        // Add the entry in the image list
        QFileInfo fileinfo(filename);
        QListWidgetItem *newItem = new QListWidgetItem;
        newItem->setText(QFileInfo(filename).fileName());
        newItem->setIcon(QIcon(QPixmap(":/icons/media/pending.png")));
        imageList->addItem(newItem);
    }

	project.save();

	imageList->setCurrentRow(imageList->count()-1);

	if (isCalibrated()) {
		setState(Calibrated);
	} else if (isValid()) {
		setState(SetupValid);
	} else {
		setState(SetupNotValid);
	}
}

void MainWindow::imageSelected(int index) {
	if (calibrationRunning) {
		return;
	}

	if (index != -1) {
		QImage image;

		// The original image is not stored in memory
		QString filename = project.getImagesFilenameList().at(index);

		ewExif->setImageFilename(filename);

		if (rbUndistorded->isChecked() and project.isCalibrated) {
			double matrixCoeff[4] = { project.cx, project.cy, project.fx, project.fy};
			double distortionCoeff[5] = {project.k1, project.k2, project.p1, project.p2, project.k3};
			calibrator->setImageFilenames(project.getImagesFilenameList());
			calibrator->setMatrixCoefficients(matrixCoeff);
			calibrator->setDistortionCoefficients(distortionCoeff);
			image = calibrator->undistorded(index);
		} else {
			QImageReader reader(filename);
			reader.setAutoTransform(true);
			image = reader.read();
		}

		histo->fromImage(image);

		if (rbBinary->isChecked()) {
			image = calibrator->binaryzeImage(image);
		}

		// Draw pattern if selected
		if (rbCorners->isChecked()) {
			patternLayer->setPoints(project.getImagePoints(index));
		} else {
			patternLayer->clear();
		}

		// Set the image in the viewer
		imageLayer->setImage(image);
		// Make the canvas size equal to the image size
		imageViewer->resizeCanvasToLayer(0);

	} else {
		imageLayer->clear();
		patternLayer->clear();
		histo->clear();
	}

	imageViewer->update();
}

QString MainWindow::getLastBrowsedPath() {
    QString lastPath = project.getProjectPath();
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

void MainWindow::setLastBrowsedPath(QString filename) {
    QSettings settings(SETTINGS_FILENAME, QSettings::IniFormat);
    settings.beginGroup("application");
    settings.setValue("lastBrowsedPath", QFileInfo(filename).filePath());
    settings.endGroup();
}

void MainWindow::newProject() {

    QString filename = QFileDialog::getSaveFileName(this, "Open File", getLastBrowsedPath(), "Project (*.eris)");
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
    QString filename = QFileDialog::getOpenFileName(this, "Open File", getLastBrowsedPath(), "Project (*.eris)");
    if (filename.isEmpty()) {
        return;
    }

    setLastBrowsedPath(filename);

    closeProject();
    project.open(filename);
    loadProject();

    if (isCalibrated()) {
    	setState(Calibrated);
    } else if (isValid()) {
    	setState(SetupValid);
    } else {
    	setState(SetupNotValid);
    }
}

void MainWindow::closeProject() {
    project.close();

    imageList->clear();

    rbOriginal->setChecked(true);

    leName->clear();

    sbSensorWidth->clear();
    sbSensorHeight->clear();
    sbColCorners->clear();
    sbRowCorners->clear();
    sbQuadWidth->clear();
    sbQuadHeight->clear();

    lcx->clear();
    lcy->clear();
    lfx->clear();
    lfy->clear();

    imageLayer->clear();
    patternLayer->clear();

    setState(NotReady);
}

void MainWindow::saveProjectAs() {
    QString filename = QFileDialog::getSaveFileName(this, "Open File", getLastBrowsedPath(), "Project (*.eris)");
    if (filename.isEmpty()) {
        return;
    }

    setLastBrowsedPath(filename);

    project.setFilename(filename);
    project.save();
}

void MainWindow::setState(AppState state) {

    switch(state) {
        case NotReady:
        	// Menu
        	closeAction->setEnabled(false);
        	saveAsAction->setEnabled(false);
        	// GUI
            imageAdd->setEnabled(false);
            imageRemove->setEnabled(false);
            gbChannels->setEnabled(false);
            twTool->setEnabled(false);
            imageViewer->setEnabled(false);
            // Default state
            rbOriginal->setChecked(true);
            break;

        case SetupNotValid:
        	// Menu
			closeAction->setEnabled(true);
			saveAsAction->setEnabled(true);
			// GUI
            imageAdd->setEnabled(true);
            imageRemove->setEnabled(true);
            gbChannels->setEnabled(true);
            rbOriginal->setEnabled(true);
            rbBinary->setEnabled(imageList->count() > 0);
            rbCorners->setEnabled(false);
            rbUndistorded->setEnabled(false);
            twTool->setEnabled(true);
            gbInputs->setEnabled(true);
            pbCalibrate->setEnabled(false);
            gbOutputs->setEnabled(false);
            imageViewer->setEnabled(true);
            break;

        case SetupValid:
        	// Menu
			closeAction->setEnabled(true);
			saveAsAction->setEnabled(true);
			// GUI
            imageAdd->setEnabled(true);
            imageRemove->setEnabled(true);
            gbChannels->setEnabled(true);
            rbOriginal->setEnabled(true);
			rbBinary->setEnabled(true);
			rbCorners->setEnabled(false);
			rbUndistorded->setEnabled(false);
            twTool->setEnabled(true);
            gbInputs->setEnabled(true);
            pbCalibrate->setEnabled(true);
            gbOutputs->setEnabled(false);
            imageViewer->setEnabled(true);
            break;

        case Calibrated:
        	// Menu
			closeAction->setEnabled(true);
			saveAsAction->setEnabled(true);
			// GUI
			imageAdd->setEnabled(true);
			imageRemove->setEnabled(true);
			gbChannels->setEnabled(true);
			rbOriginal->setEnabled(true);
			rbBinary->setEnabled(true);
			rbCorners->setEnabled(true);
			rbUndistorded->setEnabled(true);
			twTool->setEnabled(true);
			gbInputs->setEnabled(true);
			pbCalibrate->setEnabled(true);
			gbOutputs->setEnabled(true);
			imageViewer->setEnabled(true);
        	break;
    }
}

void MainWindow::loadProject() {
	disconnect(&triggerSave, &QTimer::timeout, this, &MainWindow::saveChanges);

    leName->setText(project.name);
    sbSensorWidth->setValue(project.sensorWidth);
    sbSensorHeight->setValue(project.sensorHeight);
    sbColCorners->setValue(project.cols);
    sbRowCorners->setValue(project.rows);
    sbQuadWidth->setValue(project.width);
    sbQuadHeight->setValue(project.height);

    const QStringList &filenames = project.getImagesFilenameList();
    for (int i = 0; i < filenames.count(); i++) {
        QListWidgetItem *newItem = new QListWidgetItem;
        QFileInfo fileinfo(filenames.at(i));
        newItem->setText(fileinfo.fileName());

        const QVector<QPointF> &points = project.getImagePoints(i);
        if (points.isEmpty()) {
        	newItem->setIcon(QIcon(QPixmap(":/icons/media/pending.png")));
        } else {
        	newItem->setIcon(QIcon(QPixmap(":/icons/media/valid.png")));
        }
        imageList->addItem(newItem);
    }

    updateCalibrationResuts();

    connect(&triggerSave, &QTimer::timeout, this, &MainWindow::saveChanges);
}

bool MainWindow::isValid() {
    bool valid = true;
    valid &= sbColCorners->value() > 2;
    valid &= sbRowCorners->value() > 2;
    valid &= imageList->count() > 0;
    return valid;
}

bool MainWindow::isCalibrated() {
	return project.isCalibrated;
}

void MainWindow::saveXml(QString filename) {
	qDebug() << filename;

	QFile xmlCalibration(filename);
	if (!xmlCalibration.open(QIODevice::ReadWrite)) {
		return;
	}

	QString imageFilename = project.getImagesFilenameList().at(0);
	cv::Mat image = cv::imread(imageFilename.toLocal8Bit().data());

	QXmlStreamWriter stream(&xmlCalibration);
	stream.setAutoFormatting(true);
	stream.writeStartDocument();
	stream.writeStartElement("calibrations");
	stream.writeStartElement("camera");
	stream.writeAttribute("cameraMaker", ewExif->maker);
	stream.writeAttribute("cameraModel", ewExif->model);
	stream.writeAttribute("lense", ewExif->lense);
	stream.writeAttribute("ccdwidth", QString::number(project.sensorWidth,'f',2));
	stream.writeAttribute("w", QString::number(image.cols));
	stream.writeAttribute("h", QString::number(image.rows));
	stream.writeAttribute("cx", QString::number(project.cx,'f',6));
	stream.writeAttribute("cy", QString::number(project.cy,'f',6));
	stream.writeAttribute("fx", QString::number(project.fx,'f',6));
	stream.writeAttribute("fy", QString::number(project.fy,'f',6));
	stream.writeAttribute("k1", QString::number(project.k1,'f',6));
	stream.writeAttribute("k2", QString::number(project.k2,'f',6));
	stream.writeAttribute("p1", QString::number(project.p1,'f',6));
	stream.writeAttribute("p2", QString::number(project.p2,'f',6));
	stream.writeAttribute("k3", QString::number(project.k3,'f',6));
	stream.writeAttribute("skew", "0");
	stream.writeAttribute("name", project.name.replace(' ','_'));
	stream.writeEndElement();
	stream.writeEndDocument();

	xmlCalibration.close();
}

void MainWindow::calibrate() {
	for (int i = 0; i < imageList->count(); i++) {
		QListWidgetItem *item = imageList->item(i);
		item->setIcon(QIcon(QPixmap(":/icons/media/pending.png")));
	}

	calibrationRunning = true;

	calibrator->setImageFilenames(project.getImagesFilenameList());
	calibrator->setPattern(QSize(project.width, project.height), QSize(project.cols, project.rows));
	calibrator->startCalibration();

	setState(NotReady);
}

void MainWindow::imageCalibrated(const QString filename, bool found) {
	QFileInfo fileinfo(filename);
	QList<QListWidgetItem*> items = imageList->findItems(fileinfo.baseName(), Qt::MatchContains);
	if (items.count() > 0) {
		QListWidgetItem* item = items.at(0);
		if (!found) {
			item->setIcon(QIcon(QPixmap(":/icons/media/bad.png")));
		} else {
			item->setIcon(QIcon(QPixmap(":/icons/media/valid.png")));
		}
	}
}

void MainWindow::calibrationFinished() {
	calibrationRunning = false;
	setState(Calibrated);

	double camera_coeff[4];
	double distortion_coeff[5];

	calibrator->getMatrixCoefficients(camera_coeff);
	calibrator->getDistortionCoefficients(distortion_coeff);

	project.cx = camera_coeff[0];
	project.cy = camera_coeff[1];
	project.fx = camera_coeff[2];
	project.fy = camera_coeff[3];

	project.k1 = distortion_coeff[0];
	project.k2 = distortion_coeff[1];
	project.p1 = distortion_coeff[2];
	project.p2 = distortion_coeff[3];
	project.k3 = distortion_coeff[4];

	project.isCalibrated = true;

	// Update image points
	for (int i = 0; i < project.imageCount(); i++) {
		project.setImagePoints(calibrator->getCorners(i), i);
	}

	project.save();

	saveXml(project.getProjectPath() + "/calibration.xml");

	updateCalibrationResuts();
}

void MainWindow::updateCalibrationResuts() {
	const QStringList &filenames = project.getImagesFilenameList();
	if (filenames.count() > 0) {
		QImageReader reader(filenames.at(0));
		reader.setAutoTransform(true);
		QImage image = reader.read();

		double iw = image.width();
		double ih = image.height();
		double ws = project.sensorWidth / iw;
		double hs = project.sensorHeight / ih;

		lcx->setText(QString("<b>%1 px (%2 mm)</b>").arg(project.cx, 0,'f',1).arg(project.cx*ws, 0,'f',1));
		lcy->setText(QString("<b>%1 px (%2 mm)</b>").arg(project.cy, 0,'f',1).arg(project.cy*hs, 0,'f',1));
		lfx->setText(QString("<b>%1 px (%2 mm)</b>").arg(project.fx, 0,'f',1).arg(project.fx*ws, 0,'f',1));
		lfy->setText(QString("<b>%1 px (%2 mm)</b>").arg(project.fy, 0,'f',1).arg(project.fy*hs, 0,'f',1));
	} else {
		lcx->setText("<b>0 px (0 mm)</b>");
		lcy->setText("<b>0 px (0 mm)</b>");
		lfx->setText("<b>0 px (0 mm)</b>");
		lfy->setText("<b>0 px (0 mm)</b>");
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

void MainWindow::saveChanges() {
	project.name = leName->text();
	project.sensorWidth = sbSensorWidth->value();
	project.sensorHeight = sbSensorHeight->value();
	project.cols = sbColCorners->value();
	project.rows = sbRowCorners->value();
	project.width = sbQuadWidth->value();
	project.height = sbQuadHeight->value();
	project.save();
}

void MainWindow::channelChanged() {
	imageSelected(imageList->currentRow());
}

void MainWindow::exportXml() {
	QString filename = QFileDialog::getSaveFileName(this, "Open File", getLastBrowsedPath(), "Project (*.xml)");
	if (filename.isEmpty()) {
		return;
	}

	setLastBrowsedPath(filename);

	saveXml(filename);
}
