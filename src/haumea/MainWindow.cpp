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

#include <ExifWidget.h>
#include <Utils.h>
#include <Matting.h>

#include "MainWindow.h"
#include "About.h"
#include "Project.h"


#define SETTINGS_FILENAME "haumea.ini"

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	// Resize to a reasonable size
	resize(QGuiApplication::primaryScreen()->availableSize() * 4/5);

	createActions();
	createObject();
	createLayout();
	createConnections();

	triggerSave.setSingleShot(true);
	triggerSave.setInterval(2000);

	triggerImageSave.setSingleShot(true);
	triggerImageSave.setInterval(2500);

	setState(Empty);
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

void MainWindow::closeProject() {
    project.close();

    imageList->clear();
    leName->clear();
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
    connect(&triggerImageSave, &QTimer::timeout, this, &MainWindow::saveImage);

	connect(imageAdd, &QPushButton::clicked, this, &MainWindow::addImages);
	connect(imageRemove, &QPushButton::clicked, this, &MainWindow::removeImage);
	connect(imageList, &QListWidget::currentRowChanged, this, &MainWindow::imageSelected);

	connect(leName, &QLineEdit::textEdited, this, &MainWindow::projectChanged);

	connect(pbMasking, &QPushButton::clicked, this, &MainWindow::masking);
	connect(maskPainterLayer, &PaintLayer::roiChanged, this, &MainWindow::maskingRoi);

	connect(maskPainterLayer, &PaintLayer::changed, this, &MainWindow::imageChanged);
}

void MainWindow::createObject() {
	leName = new QLineEdit;

	imageList = new QListWidget();
	imageAdd = new QPushButton("Add");
	imageRemove = new QPushButton("Remove");

	// ------------------------------------------------------------------------
	// Tools
	twTool = new QTabWidget;

	// ------------------------------------------------------------------------
    // Inputs
    gbInputs = new QGroupBox("Inputs");

    pbMasking = new QPushButton("Masking");
    cbRefine = new QCheckBox("Refine");

	// ------------------------------------------------------------------------
	// Image viewer
	imageViewer = new LayerGraphicWidget;
	imageLayer = new PixelLayer;
	maskPainterLayer = new PaintLayer;

	imageViewer->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	imageViewer->addLayer(imageLayer);
	imageViewer->addLayer(maskPainterLayer);
	imageViewer->setCoreInteration(false);

	// ------------------------------------------------------------------------
	// Exif viewer
	ewExif = new ExifWidget;
}

void MainWindow::createLayout() {
	// Image list layout
	QGridLayout *gridlayout1 = new QGridLayout;
	gridlayout1->addWidget(imageList,0,0,1,2);
	gridlayout1->addWidget(imageAdd,1,0);
	gridlayout1->addWidget(imageRemove,1,1);

	// ------------------------------------------------------------------------
	// Inputs
	QGridLayout *gridlayout2 = new QGridLayout;
	gridlayout2->addWidget(new QLabel("Name"), 0, 0);
	gridlayout2->addWidget(leName, 0, 1, 1, 2);
	gridlayout2->addWidget(pbMasking, 1, 0);
	gridlayout2->addWidget(cbRefine, 1,1);
	gbInputs->setLayout(gridlayout2);

	// ------------------------------------------------------------------------
	// Tools panel widget
	QVBoxLayout *vboxlayout1 = new QVBoxLayout;
	vboxlayout1->addWidget(gbInputs);
	vboxlayout1->addItem(new QSpacerItem(0,10, QSizePolicy::Expanding, QSizePolicy::Expanding));
	QWidget *wControl = new QWidget;
	wControl->setLayout(vboxlayout1);

	twTool->addTab(wControl, "Control");
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

QString MainWindow::getCurrentBimapFileName() {
    int index = imageList->currentRow();
    if (index == -1) {
        return QString();

    }
    QString filename = project.getImagesFilenameList().at(index);
    QFileInfo fileinfo(filename);

    return QString("%1/cache/%2_bimap.png").arg(project.getProjectPath()).arg(fileinfo.baseName());
}

QString MainWindow::getCurrentMaskFileName() {
    int index = imageList->currentRow();
    if (index == -1) {
        return QString();

    }
    QString filename = project.getImagesFilenameList().at(index);
    QFileInfo fileinfo(filename);

    return QString("%1/%2_mask.png").arg(fileinfo.path()).arg(fileinfo.baseName());
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

void MainWindow::imageChanged() {
    qDebug() << "MainWindow::imageChanged()";
    triggerImageSave.start();
}

void MainWindow::imageSelected(int index) {

	QImage image;
	QImage bimap;

	QString filename;

	if (index != -1) {
		// The original image is not stored in memory
		filename = project.getImagesFilenameList().at(index);
		ewExif->setImageFilename(filename);
	}

	if (!filename.isEmpty()) {
		QImageReader reader(filename);
		reader.setAutoTransform(true);
		image = reader.read();

		QString bimap_filename = getCurrentBimapFileName();
		if (QFile(bimap_filename).exists()) {
		    QImageReader reader(bimap_filename);
		    reader.setAutoTransform(true);
		    bimap = reader.read();
		}
	}

	// Set the image in the viewer
	imageLayer->setImage(image);
	imageViewer->resizeCanvasToLayer(0);

	if (bimap.isNull()) {
        maskPainterLayer->clear();
        maskPainterLayer->resizeToCanvas();
	} else {
	    maskPainterLayer->setImage(bimap);
	}

	imageViewer->update();
}

bool MainWindow::isValid() {
	return project.imageCount() > 0;
}

void MainWindow::loadProject() {
	disconnect(&triggerSave, &QTimer::timeout, this, &MainWindow::saveChanges);

	leName->setText(project.name);

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

void MainWindow::createCvMatrices(cv::Mat &I, cv::Mat &mI, cv::Mat &const_map, cv::Mat &const_val) {
    // Convert the image to float in range 0 to 1
    cv::Mat I_d = qImageToCvMat(imageLayer->image());
    I.create(I_d.size(), CV_32FC3);
    I_d.convertTo(I, CV_32FC3, 1.0/255.0);

    // Mask background and foreground
    mI = qImageToCvMat(maskPainterLayer->image());

    // Known pixel values
    const_map = cv::Mat::zeros(mI.size(), CV_32F);
    cv::Mat const_map_FG = cv::Mat::zeros(mI.size(), CV_32F);
    cv::Mat const_map_BG = cv::Mat::zeros(mI.size(), CV_32F);

    // Background (0) or foreground (1) values
    const_val = cv::Mat::zeros(mI.size(), CV_32F);

    for (int i = 0; i < mI.rows; i++) {

        cv::Vec4b *mI_ptr = mI.ptr<cv::Vec4b>(i);

        float *const_val_ptr = const_val.ptr<float>(i);
        float *const_map_FG_ptr = const_map_FG.ptr<float>(i);
        float *const_map_BG_ptr = const_map_BG.ptr<float>(i);

        for (int j = 0; j < mI.cols; j++) {
            const_map_FG_ptr[j] = mI_ptr[j][0] > 0 ? 1.0 : 0.0;
            const_map_BG_ptr[j] = mI_ptr[j][2] > 0 ? 1.0 : 0.0;
            const_val_ptr[j] = mI_ptr[j][0] > 0 ? 1.0 : 0.0;
        }
    }

    cv::Mat A,B,element;

    // For fully computed mask, make a separation between the foreground and background
    element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(5,5));
    cv::erode(const_map_FG,A,element);
    cv::erode(const_map_BG,B,element);
    const_map = A + B;

#if 0
    char title[64];
    sprintf(title, "mI");
    cv::imshow(title, mI);
#endif
}

void MainWindow::maskRoi(const cv::Mat &_I, cv::Mat &_mI, const cv::Mat &_const_map, const cv::Mat &_const_val, cv::Rect roi) {
    const cv::Mat *I_ptr;
    cv::Mat *mI_ptr;
    const cv::Mat *const_map_ptr;
    const cv::Mat *const_val_ptr;

    cv::Mat I_roi;
    cv::Mat mI_roi;
    cv::Mat const_map_roi;
    cv::Mat const_val_roi;

    if (_I.size() == roi.size()) {
        qDebug() << "MainWindow::maskRoi():full image";
        I_ptr = &_I;
        mI_ptr = &_mI;
        const_val_ptr = &_const_val;
        const_map_ptr = &_const_map;
    } else {
        I_roi = _I(roi);
        mI_roi = _mI(roi);
        const_map_roi = _const_map(roi);
        const_val_roi = _const_val(roi);

        I_ptr = &I_roi;
        mI_ptr = &mI_roi;
        const_val_ptr = &const_val_roi;
        const_map_ptr = &const_map_roi;
    }

    cv::Mat A,B,element;

    // Solve for the mask of the entire image
    cv::Mat alpha;
    solveMultiLevelAlpha32FC3(*I_ptr, *const_map_ptr, *const_val_ptr, alpha, 0);
    cv::convertScaleAbs(alpha, A, 255, 0);

    // Close operation
    element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(3,3));
    cv::dilate(A,B,element);
    element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(5,5));
    cv::erode(B,A,element);

    // Alpha mask binarization
    int thr = 255*0.2;
    for (int i = 0; i < mI_ptr->rows; i++) {

        cv::Vec4b *mI_row_ptr = mI_ptr->ptr<cv::Vec4b>(i);
        uchar *B_ptr = B.ptr<uchar>(i);
        uchar *A_ptr = A.ptr<uchar>(i);

        for (int j = 0; j < mI_ptr->cols; j++) {
            mI_row_ptr[j][2] = (A_ptr[j] < thr) ? 255 : 0;
            mI_row_ptr[j][0] = (A_ptr[j] >= thr) ? 255 : 0;
            mI_row_ptr[j][3] = (mI_row_ptr[j][0] || mI_row_ptr[j][2]) ? 100 : 0;

            B_ptr[j] = (A_ptr[j] >= thr) ? 255 : 0;
        }
    }

    if (_I.size() != roi.size()) {
        _mI(roi) = mI_roi;
    }

#if 0
    char title[64];
    sprintf(title, "ROI");
    cv::imshow(title, A);
#endif
}

void MainWindow::masking() {

    triggerImageSave.stop();
    saveImage();

    // Obtain the mask matrices
    cv::Mat I, mI, const_map, const_val;
    createCvMatrices(I, mI, const_map, const_val);

    // Obtain the coarse mask
    maskRoi(I, mI, const_map, const_val, cv::Rect(0,0,I.cols,I.rows));

    // Update the FB/BG mask definition
    QImage mask = cvMatToQImage(mI);
    maskPainterLayer->setImage(mask);
    imageViewer->update();

    saveImage();
}

void MainWindow::maskingRoi(QRect roi) {
    // Only allowed when th refine option is selected
    if (!cbRefine->isChecked()) {
        return;
    }

    triggerImageSave.stop();
    saveImage();

    // Obtain the mask matrices
    cv::Mat I, mI, const_map, const_val;
    createCvMatrices(I, mI, const_map, const_val);

    // Obtain the fine mask
    roi.adjust(-20,-20,20,20);
    if(roi.left() < 0) {
        roi.setLeft(0);
    }
    if (roi.top() < 0) {
        roi.setTop(0);
    }
    if (roi.right() >= I.cols) {
        roi.setRight(I.cols-1);
    }
    if (roi.bottom() >= I.rows) {
        roi.setBottom(I.rows-1);
    }
    maskRoi(I, mI, const_map, const_val, cv::Rect(roi.x(),roi.y(),roi.width(),roi.height()));

    // Update the FB/BG mask definition
    QImage mask = cvMatToQImage(mI);
    maskPainterLayer->setImage(mask);
    imageViewer->update();

    saveImage();
}

void MainWindow::newProject() {

    QString filename = QFileDialog::getSaveFileName(this, "Open File", getLastBrowsedPath(), "Project (*.haum)");
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
    QString filename = QFileDialog::getOpenFileName(this, "Open File", getLastBrowsedPath(), "Project (*.haum)");
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
			// Tools
			twTool->setEnabled(true);
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
			// Tools
			twTool->setEnabled(true);
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
			// Tools
			twTool->setEnabled(false);
			break;
	}
}

void MainWindow::saveChanges() {
	project.name = leName->text();
	project.save();
}

void MainWindow::saveImage() {
    qDebug() << "MainWindow::saveImage()";
    QString bimap_filename = getCurrentBimapFileName();
    const QImage &bimap = maskPainterLayer->image();
    bimap.save(bimap_filename, nullptr, 100);

    QString mask_filename = getCurrentMaskFileName();
    QImage mask = toMask(bimap);
    bool rc = mask.save(mask_filename, nullptr, -1);
    qDebug() << rc;
}

QImage MainWindow::toMask(const QImage& source) {
    int width = source.width();
    int height = source.height();

    qDebug() << width << height;

    QImage mask = QImage(source.size(), QImage::Format_RGB32);

    for (int j=0; j < width; j++) {
        for (int i=0; i < height; i++) {
            QColor color = source.pixel(j,i);
            if (color.red() > 0) {
                mask.setPixelColor(j,i, QColor(255,255,255));
            } else {
                mask.setPixelColor(j,i, QColor(0,0,0));
            }
        }
    }

    return mask;
}

void MainWindow::saveProjectAs() {
    QString filename = QFileDialog::getSaveFileName(this, "Open File", getLastBrowsedPath(), "Project (*.haum)");
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

