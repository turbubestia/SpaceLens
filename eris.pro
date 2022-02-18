CONFIG -= debug release
CONFIG += debug_and_release

QT += core gui widgets

Debug: TARGET = erisd
Release: TARGET = eris

TEMPLATE = app

SOURCES += src/libs/Utils.cpp \
	src/libs/ExifWidget.cpp \
	src/libs/LayerGraphicWidget.cpp \
	src/libs/AbstractLayer.cpp \
	src/libs/PixelLayer.cpp \
	src/libs/PointArrayLayer.cpp \
	src/libs/CameraCalibration.cpp \
	src/libs/HistogramWidget.cpp
	
HEADERS += src/libs/Utils.h \
	src/libs/ExifWidget.h \
	src/libs/LayerGraphicWidget.h \
	src/libs/AbstractLayer.h \
	src/libs/PixelLayer.h \
	src/libs/PointArrayLayer.h \
	src/libs/CameraCalibration.h \
	src/libs/HistogramWidget.h

SOURCES += src/eris/main.cpp
	
SOURCES += src/eris/MainWindow.cpp \
	src/eris/Project.cpp \
	src/eris/About.cpp
	
HEADERS += src/eris/MainWindow.h \
	src/eris/Project.h \
	src/eris/About.h

INCLUDEPATH += src/eris
INCLUDEPATH += src/libs
INCLUDEPATH += $$(MINGW_HOME)/include/opencv4
INCLUDEPATH += $$(MINGW_HOME)/include/opencv4/opencv2

LIBS += -L$$(MINGW_HOME)/lib
LIBS += -llibopencv_core -llibopencv_highgui -llibopencv_imgcodecs
LIBS += -llibopencv_imgproc -llibopencv_features2d -llibopencv_calib3d -llibopencv_photo
LIBS += -llibexiv2

RC_ICONS = media/eris.ico

RESOURCES += qdarkstyle/style.qrc
RESOURCES += resources.qrc