/*
 * MainWindow.cpp
 *
 *  Created on: 22 sep. 2019
 *      Author: claud
 */

#include <QtWidgets/QtWidgets>

#include <opencv4/opencv2/opencv.hpp>
#include <exiv2/exiv2.hpp>

#include <Utils.h>
#include <LayerGraphicWidget.h>
#include <PixelLayer.h>
#include <HistogramWidget.h>
#include <ExifWidget.h>

#include "MainWindow.h"
#include "About.h"
#include "Project.h"

#define SETTINGS_FILENAME "makemake.ini"

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
        project.addImage(filenames.at(i));
    }

	project.updateList();

	updateImageList();

	setState(SetupValid);

	triggerSave.start();
}

void MainWindow::channelChanged() {
    imageSelected(imageList->currentItem(), nullptr);
}

void MainWindow::closeProject() {
    project.close();

    imageList->clear();
    imageList->clear();

    rbOriginal->setChecked(true);

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

	connect(imageList, &QTreeWidget::currentItemChanged, this, &MainWindow::imageSelected);

	connect(imageAdd, &QPushButton::clicked, this, &MainWindow::addImages);
	connect(imageRemove, &QPushButton::clicked, this, &MainWindow::removeImage);

	connect(rbByName, &QPushButton::clicked, this, &MainWindow::groupCriteriaChanged);
	connect(rbSimple, &QPushButton::clicked, this, &MainWindow::groupCriteriaChanged);
	connect(rbBatch, &QPushButton::clicked, this, &MainWindow::groupCriteriaChanged);
	connect(rbBatchByCount, &QPushButton::clicked, this, &MainWindow::groupCriteriaChanged);
	connect(sbGroupCount, QOverload<int>::of(&QSpinBox::valueChanged), this, &MainWindow::groupCriteriaChanged);  // @suppress("Invalid arguments")

	connect(rbOriginal, &QRadioButton::clicked, this, &MainWindow::channelChanged);
	connect(rbFocused, &QRadioButton::clicked, this, &MainWindow::channelChanged);

	connect(leName, &QLineEdit::textEdited, this, &MainWindow::projectChanged);

	connect(pbStack, &QPushButton::clicked, this, &MainWindow::focus);

	connect(focusStacker, &FocusStack::imageFocused, this, &MainWindow::imageFocused);
	connect(focusStacker, &FocusStack::imageAligned, this, &MainWindow::imageAligned);
	connect(focusStacker, &FocusStack::finished, this, &MainWindow::focusFinish);
}

void MainWindow::createObject() {
	leName = new QLineEdit;

	imageList = new QTreeWidget();
	imageList->header()->setHidden(true);

	imageAdd = new QPushButton("Add");
	imageRemove = new QPushButton("Remove");

	gbChannels = new QGroupBox("Channels");
	rbOriginal = new QRadioButton("Original");
	rbFocused = new QRadioButton("Focused");

	// ------------------------------------------------------------------------
	// Tools
	twTool = new QTabWidget;

	// ------------------------------------------------------------------------
	// Simple or Batch sorting
	gbSorting = new QGroupBox("Mode");

	rbByName = new QRadioButton("Name");
	QButtonGroup *sortMode = new QButtonGroup;
	sortMode->addButton(rbByName);

	rbSimple = new QRadioButton("Simple");
	rbBatch = new QRadioButton("Batch");
	QButtonGroup *groupMode = new QButtonGroup;
	groupMode->addButton(rbSimple);
	groupMode->addButton(rbBatch);

	rbBatchByCount = new QRadioButton("Count");
	QButtonGroup *batchMode = new QButtonGroup;
	batchMode->addButton(rbBatchByCount);

	sbGroupCount = new QSpinBox;
	sbGroupCount->setMinimum(2);
	sbGroupCount->setMaximum(20);

	// ------------------------------------------------------------------------
	// Action buttons
	gbActions = new QGroupBox("Actions");
	pbStack = new QPushButton("Stack");
	pbStack->setEnabled(false);

	// ------------------------------------------------------------------------
	// Image viewer
	imageViewer = new LayerGraphicWidget;
	imageLayer = new PixelLayer;

	imageViewer->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	imageViewer->addLayer(imageLayer);
	imageViewer->setCoreInteration(true);

	// ------------------------------------------------------------------------
	// Exif viewer
	ewExif = new ExifWidget;

	// ------------------------------------------------------------------------
	// Non-GUI objects

	focusStacker = new FocusStack;
}

void MainWindow::createLayout() {
	// Image list layout
	QGridLayout *gridlayout_imagelist = new QGridLayout;
	gridlayout_imagelist->addWidget(imageList,0,0,1,2);
	gridlayout_imagelist->addWidget(imageAdd,1,0);
	gridlayout_imagelist->addWidget(imageRemove,1,1);

	QGridLayout *gridlayout_channels = new QGridLayout;
	gridlayout_channels->addWidget(rbOriginal, 0, 0);
	gridlayout_channels->addWidget(rbFocused, 0, 1);
	gbChannels->setLayout(gridlayout_channels);

	gridlayout_imagelist->addWidget(gbChannels, 2,0,1,2);

	// ------------------------------------------------------------------------
	// Simple or Batch sorting
	QGridLayout *gridlayout_mode = new QGridLayout;
	gridlayout_mode->addWidget(new QLabel("Sort by"),0,0);
	gridlayout_mode->addWidget(rbByName,0,1);
	//gridlayout_mode->addWidget(rbByDateTime,0,2);
	gridlayout_mode->addWidget(new QLabel("Process"),1,0);
	gridlayout_mode->addWidget(rbSimple,1,1);
	gridlayout_mode->addWidget(rbBatch,1,2);
	gridlayout_mode->addWidget(new QLabel("Batch group by"),2,0);
	gridlayout_mode->addWidget(rbBatchByCount,2,1);
	gridlayout_mode->addWidget(sbGroupCount,2,2);
	gbSorting->setLayout(gridlayout_mode);

	// ------------------------------------------------------------------------
	// Action buttons
	QHBoxLayout *hboxlayout_actions = new QHBoxLayout;
	hboxlayout_actions->addWidget(pbStack);

	// ------------------------------------------------------------------------
	// Tools panel widget
	QVBoxLayout *vboxlayout_tools = new QVBoxLayout;
	vboxlayout_tools->addWidget(gbSorting);
	vboxlayout_tools->addLayout(hboxlayout_actions);
	vboxlayout_tools->addItem(new QSpacerItem(0,10, QSizePolicy::Expanding, QSizePolicy::Expanding));
	QWidget *wTools = new QWidget;
	wTools->setLayout(vboxlayout_tools);
	twTool->addTab(wTools, "Tools");
	twTool->addTab(ewExif, "Exif");

	// ------------------------------------------------------------------------
	// Left panel widget
	QVBoxLayout *vboxlayout_left = new QVBoxLayout;
	vboxlayout_left->addLayout(gridlayout_imagelist, 1);
	vboxlayout_left->addWidget(twTool, 0);
	QWidget *leftPanelWidget = new QWidget;
	leftPanelWidget->setLayout(vboxlayout_left);

	// ------------------------------------------------------------------------
	// Join them with a splitter
	QSplitter *splitter = new QSplitter(this);
	splitter->addWidget(leftPanelWidget);
	splitter->addWidget(imageViewer);
	setCentralWidget(splitter);
}

QString MainWindow::getLastBrowsedPath() {
    QString lastPath;
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

void MainWindow::imageSelected(QTreeWidgetItem *item, QTreeWidgetItem *previous) {
    Q_UNUSED(previous);

    if (item != nullptr) {
        int index = item->data(0, Qt::UserRole).toInt();

        QString filename = project.getImageFilename(index);

        if (rbFocused->isChecked()) {
            QString basename;
            if (project.isBatchEnabled()) {
                if (item->parent() != nullptr) {
                    index = item->parent()->data(0, Qt::UserRole).toInt();
                }
            } else {
                index = 0;
            }

            QString focused = project.getImageFilename(index);
            QFileInfo fileinfo(focused);
            focused = project.getProjectPath() + "/output/" + fileinfo.baseName() + "_focused.jpg";
            fileinfo = QFileInfo(focused);

            if (fileinfo.exists()) {
                filename = focused;
            }
        }

        ewExif->setImageFilename(filename);

        QImageReader reader(filename);
        reader.setAutoTransform(true);
        QImage image = reader.read();

        // Set the image in the viewer
        imageLayer->setImage(image);
        // Make the canvas size equal to the image size
        imageViewer->resizeCanvasToLayer(0);
    } else {
        imageLayer->clear();
    }

    imageViewer->update();
}

void MainWindow::imageFocused(int group) {
    QTreeWidgetItem *top;

    if (project.isBatchEnabled()) {
        top = imageList->invisibleRootItem()->child(group);
    } else {
        top = imageList->invisibleRootItem();
    }

    for(int i = 0; i < top->childCount(); i++) {
        QTreeWidgetItem *item = top->child(i);
        item->setIcon(0,QIcon(QPixmap(":/icons/media/valid.png")));
    }
}

void MainWindow::imageAligned(int group, QString filename) {
    QTreeWidgetItem *item;

    if (project.isBatchEnabled()) {
        int index = project.indexOf(group, filename);
        QTreeWidgetItem *top = imageList->invisibleRootItem()->child(group);
        item = top->child(index);
    } else {
        int index = project.indexOf(-1, filename);
        item = imageList->invisibleRootItem()->child(index);
    }

    item->setIcon(0,QIcon(QPixmap(":/icons/media/pending.png")));
}

void MainWindow::focus() {
    int groups = project.groupCount();

    QTreeWidgetItem *root = imageList->invisibleRootItem();
    for(int i = 0; i < root->childCount(); i++) {
        if (project.isBatchEnabled()) {
            QTreeWidgetItem *top = root->child(i);
            for (int j = 0; j < top->childCount(); j++) {\
                QTreeWidgetItem *item = top->child(j);
                item->setIcon(0,QIcon(QPixmap(":/icons/media/neutral.png")));
            }
        } else {
            QTreeWidgetItem *item = root->child(i);
            item->setIcon(0,QIcon(QPixmap(":/icons/media/neutral.png")));
        }
    }

    focusStacker->clear();

    QString cacheFolder = QString("%1/cache").arg(project.getProjectPath());
    project.getProjectDir().mkdir("cache");

    QString outputFolder = QString("%1/output").arg(project.getProjectPath());
    project.getProjectDir().mkdir("output");

    focusStacker->setCacheFolder(cacheFolder);
    focusStacker->setOutputFolder(outputFolder);

    if (project.isBatchEnabled()) {
        for (int g = 0; g < groups; g++) {
            focusStacker->addStackGroup(project.getImageFilenames(g));
        }
    } else {
        focusStacker->addStackGroup(project.getImageFilenames(-1));
    }

    focusStacker->startFocusStack();

    setState(Working);
}

void MainWindow::focusFinish() {
    setState(SetupValid);
}

void MainWindow::groupCriteriaChanged() {
    project.setBatchEnable(rbBatch->isChecked());

    if (rbByName->isChecked()) {
        project.setSortCriteria(Project::ByName);
    }

    if (rbBatchByCount->isChecked()) {
        project.setGroupByCount(sbGroupCount->value());
    }

    project.updateList();

    updateImageList();

    triggerSave.start();
}

bool MainWindow::isValid() {
	return project.imageCount() > 0;
}

void MainWindow::loadProject() {
	disconnect(&triggerSave, &QTimer::timeout, this, &MainWindow::saveChanges);

	leName->setText(project.name);

	switch (project.getSortCriteria()) {
	    case Project::ByName: rbByName->setChecked(true); break;
	    default: rbByName->setChecked(true); break;
	}

	rbBatch->setChecked(project.isBatchEnabled());

	switch (project.getGroupCriteria()) {
	    case Project::ByCount: rbBatchByCount->setChecked(true); break;
	}

	sbGroupCount->setValue(project.getGroupByCountSize());

	updateImageList();

	connect(&triggerSave, &QTimer::timeout, this, &MainWindow::saveChanges);
}

void MainWindow::newProject() {

    QString filename = QFileDialog::getSaveFileName(this, "Open File", getLastBrowsedPath(), "Project (*.make)");
    if (filename.isEmpty()) {
        return;
    }

    setLastBrowsedPath(filename);

    closeProject();

    project.loadDefault();
    project.setFilename(filename);
    project.save();

    loadProject();

    setState(SetupNotValid);
}

void MainWindow::openProject() {
    QString filename = QFileDialog::getOpenFileName(this, "Open File", getLastBrowsedPath(), "Project (*.make)");
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
    int index = -1;

    QTreeWidgetItem *item = imageList->currentItem();
    if (item != nullptr) {
        index = item->data(0, Qt::UserRole).toInt();
    }

	if (index != -1) {
		project.removeImage(index);
		project.updateList();

		updateImageList();

		triggerSave.start();
	}
}

void MainWindow::setState(AppState state) {

    switch(state) {
        case Empty:
            imageAdd->setEnabled(false);
            imageRemove->setEnabled(false);
            gbChannels->setEnabled(false);
            rbOriginal->setChecked(true);
            leName->setEnabled(false);
            gbSorting->setEnabled(false);
            rbByName->setChecked(true);
            rbSimple->setChecked(true);
            rbBatchByCount->setChecked(true);
            sbGroupCount->setValue(1);
            break;

        case SetupNotValid:
            imageAdd->setEnabled(true);
            imageRemove->setEnabled(true);
            gbChannels->setEnabled(false);
            leName->setEnabled(true);
            gbSorting->setEnabled(false);
            pbStack->setEnabled(false);
            break;

        case SetupValid:
            imageAdd->setEnabled(true);
            imageRemove->setEnabled(true);
            gbChannels->setEnabled(true);
            leName->setEnabled(true);
            gbSorting->setEnabled(true);
            pbStack->setEnabled(true);
            break;

        case Working:
            imageAdd->setEnabled(false);
            imageRemove->setEnabled(false);
            gbChannels->setEnabled(false);
            leName->setEnabled(false);
            gbSorting->setEnabled(false);
            pbStack->setEnabled(false);
            break;
    }
}

void MainWindow::saveChanges() {
	project.name = leName->text();

	project.save();
}

void MainWindow::saveProjectAs() {
    QString filename = QFileDialog::getSaveFileName(this, "Open File", getLastBrowsedPath(), "Project (*.make)");
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

void MainWindow::updateImageList() {

	imageList->clear();

	if (project.isBatchEnabled()) {
		int groupCount = project.groupCount();
		for (int i = 0; i < groupCount; i++) {
			QStringList filenames = project.getImageFilenames(i);

			QString groupName = QString("image_%1").arg(i,3,10,QChar('0'));
			QTreeWidgetItem *newTopItem = new QTreeWidgetItem();
			newTopItem->setText(0,groupName);
			newTopItem->setData(0,Qt::UserRole,project.getImageIndex(0,i));

			for (int j = 0; j < filenames.count(); j++) {
				QTreeWidgetItem *newChildItem = new QTreeWidgetItem;
				QFileInfo fileinfo(filenames.at(j));
				newChildItem->setText(0,fileinfo.fileName());
				newChildItem->setIcon(0,QIcon(QPixmap(":/icons/media/neutral.png")));
				newChildItem->setData(0,Qt::UserRole,project.getImageIndex(j,i));
				newTopItem->addChild(newChildItem);
			}

			imageList->addTopLevelItem(newTopItem);
		}
	} else {
		QStringList filenames = project.getImageFilenames();
		for (int i = 0; i < filenames.count(); i++) {
			QTreeWidgetItem *newItem = new QTreeWidgetItem;
			QFileInfo fileinfo(filenames.at(i));
			newItem->setText(0,fileinfo.fileName());
			newItem->setIcon(0,QIcon(QPixmap(":/icons/media/neutral.png")));
			newItem->setData(0,Qt::UserRole,i);
			imageList->addTopLevelItem(newItem);
		}
	}
}
