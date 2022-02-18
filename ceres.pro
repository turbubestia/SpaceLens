CONFIG -= debug release
CONFIG += debug_and_release

QT += core gui widgets

Debug: TARGET = ceresd
Release: TARGET = ceres

TEMPLATE = app

SOURCES += src/libs/Utils.cpp \
	src/libs/ExifWidget.cpp \
	src/libs/LayerGraphicWidget.cpp \
	src/libs/AbstractLayer.cpp \
	src/libs/PixelLayer.cpp \
	src/libs/RoiEditorLayer.cpp \
	src/libs/HistogramWidget.cpp \
	src/libs/ImageMerger.cpp \
	src/libs/ImageCropper.cpp
	
HEADERS += src/libs/Utils.h \
	src/libs/ExifWidget.h \
	src/libs/LayerGraphicWidget.h \
	src/libs/AbstractLayer.h \
	src/libs/PixelLayer.h \
	src/libs/RoiEditorLayer.h \
	src/libs/HistogramWidget.h \
	src/libs/CameraCalibrationParameters.h \
	src/libs/ImageMerger.h \
	src/libs/ImageCropper.h

SOURCES += src/ceres/main.cpp
	
SOURCES += src/ceres/MainWindow.cpp \
	src/ceres/About.cpp \
	src/ceres/Project.cpp
	
HEADERS += src/ceres/MainWindow.h \
	src/ceres/About.h \
	src/ceres/Project.h

INCLUDEPATH += src/ceres
INCLUDEPATH += src/libs
INCLUDEPATH += $$(MINGW_HOME)/include/opencv4
INCLUDEPATH += $$(MINGW_HOME)/include/opencv4/opencv2

LIBS += -L$$(MINGW_HOME)/lib
LIBS += -llibopencv_core -llibopencv_highgui -llibopencv_imgcodecs
LIBS += -llibopencv_imgproc -llibopencv_features2d -llibopencv_calib3d -llibopencv_photo
LIBS += -llibexiv2

#RC_ICONS = media/eris.ico

RESOURCES += qdarkstyle/style.qrc
RESOURCES += resources.qrc