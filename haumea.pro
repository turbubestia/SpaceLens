CONFIG -= debug release
CONFIG += debug_and_release

QT += core gui widgets

Debug: TARGET = haumead
Release: TARGET = haumea

TEMPLATE = app

SOURCES += src/libs/Utils.cpp \
	src/libs/ExifWidget.cpp \
	src/libs/LayerGraphicWidget.cpp \
	src/libs/AbstractLayer.cpp \
	src/libs/PixelLayer.cpp \
	src/libs/PaintLayer.cpp \
	src/libs/Matting.cpp
	
HEADERS += src/libs/Utils.h \
	src/libs/ExifWidget.h \
	src/libs/LayerGraphicWidget.h \
	src/libs/AbstractLayer.h \
	src/libs/PixelLayer.h \
	src/libs/PaintLayer.h \
	src/libs/Matting.h

SOURCES += src/haumea/main.cpp
	
SOURCES += src/haumea/MainWindow.cpp \
	src/haumea/Project.cpp \
	src/haumea/About.cpp
	
HEADERS += src/haumea/MainWindow.h \
	src/haumea/Project.h \
	src/haumea/About.h

INCLUDEPATH += src/haumea
INCLUDEPATH += src/libs
INCLUDEPATH += $$(MINGW_HOME)/include/opencv4
INCLUDEPATH += $$(MINGW_HOME)/include/opencv4/opencv2

LIBS += -L$$(MINGW_HOME)/lib
LIBS += -llibopencv_core -llibopencv_highgui -llibopencv_imgcodecs
LIBS += -llibopencv_imgproc -llibopencv_features2d -llibopencv_calib3d -llibopencv_photo
LIBS += -llibexiv2

#RC_ICONS = media/haumea.ico

RESOURCES += qdarkstyle/style.qrc
RESOURCES += resources.qrc