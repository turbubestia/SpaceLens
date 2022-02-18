CONFIG -= debug release
CONFIG += debug_and_release

QT += core gui widgets

Debug: TARGET = makemaked
Release: TARGET = makemake

TEMPLATE = app

SOURCES += src/libs/Utils.cpp \
	src/libs/ExifWidget.cpp \
	src/libs/LayerGraphicWidget.cpp \
	src/libs/AbstractLayer.cpp \
	src/libs/PixelLayer.cpp \
	src/libs/FocusStack.cpp
	
HEADERS += src/libs/Utils.h \
	src/libs/ExifWidget.h \
	src/libs/LayerGraphicWidget.h \
	src/libs/AbstractLayer.h \
	src/libs/PixelLayer.h \
	src/libs/FocusStack.h

SOURCES += src/makemake/main.cpp
	
SOURCES += src/makemake/MainWindow.cpp \
	src/makemake/About.cpp \
	src/makemake/Project.cpp
	
HEADERS += src/makemake/MainWindow.h \
	src/makemake/About.h \
	src/makemake/Project.h

INCLUDEPATH += src/ceres
INCLUDEPATH += src/libs
INCLUDEPATH += $$(MINGW_HOME)/include/opencv4
INCLUDEPATH += $$(MINGW_HOME)/include/opencv4/opencv2

LIBS += -L$$(MINGW_HOME)/lib
LIBS += -llibopencv_core -llibopencv_highgui -llibopencv_imgcodecs
LIBS += -llibopencv_imgproc -llibopencv_features2d -llibopencv_xfeatures2d -llibopencv_calib3d -llibopencv_photo -llibopencv_video -llibopencv_ximgproc
LIBS += -llibexiv2

#RC_ICONS = media/eris.ico

RESOURCES += qdarkstyle/style.qrc
RESOURCES += resources.qrc