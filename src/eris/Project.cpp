/*
 * Project.cpp
 *
 *  Created on: 1 abr. 2020
 *      Author: claud
 */

#include <QtCore/QtCore>
#include <QtWidgets/QtWidgets>

#include <exiv2/exiv2.hpp>
#include <opencv4/opencv2/opencv.hpp>

#include "Project.h"

Project::Project() {
    sensorWidth = 20;
    sensorHeight = 20;
	cols = 7;
	rows = 4;
	width = 10.0;
	height = 10.0;
	cx = cy = 0;
	fx = fy = 0;
	k1 = k2 = k3 = 0;
	p1 = p2 = 0;
	isCalibrated = false;
}

Project::~Project() {

}

bool Project::addImage(QString filename){
	// Internally the image are always stored in absolute mode
	filename = projectpath.absoluteFilePath(filename);
	// Don't add duplicates
	if (!imageFilenames.contains(filename)) {
        imageFilenames.append(filename);
        imagePoints.append(QVector<QPointF>());
	}

	return true;
}

void Project::close() {
	filename.clear();
	imageFilenames.clear();
	imagePoints.clear();
}

const QVector<QPointF> & Project::getImagePoints(int row) {
	return imagePoints.at(row);
}

const QStringList& Project::getImagesFilenameList() {
	return imageFilenames;
}

QString Project::getProjectFilename() {
    QFileInfo fileinfo(projectpath.absoluteFilePath(filename));
    return fileinfo.baseName();
}

QString Project::getProjectPath() {
	return projectpath.path();
}

int Project::imageCount() {
	return imageFilenames.count();
}

bool Project::isSaved() {
	return !filename.isEmpty();
}

void Project::loadDefault() {
    name = "Untitled";
    sensorWidth = 20.0;
    rows = 8;
    cols = 11;
    width = 20;
    height = 20;
    cx = cy = 0;
    fx = fy = 0;
    k1 = k2 = p1 = p2 = k3 = 0;
    isCalibrated = false;
    imageFilenames.clear();
    imagePoints.clear();
}

void Project::open(QString filepath) {
    QFile xmlCalibration(filepath);
    if (!xmlCalibration.open(QIODevice::ReadWrite)) {
        return;
    }

    close();

    setFilename(filepath);
    filename = toRelativePath(filepath);

    QXmlStreamReader xml(&xmlCalibration);
    QXmlStreamAttributes attributes;

    int state = 0;
    int nbImages = 0;
    int nbPoints = 0;

    QVector<QPointF> points;

    while(!xml.atEnd()) {
        xml.readNextStartElement();
        attributes = xml.attributes();

        if (xml.isEndElement()) {
        	continue;
        }

        switch(state) {
            case 0:
                if (xml.name() == u"project") {
                    if (attributes.count() > 0) {
                        name = attributes.value("name").toString();
                    }
                    state = 1;
                }
                break;

            case 1:
                if (xml.name() == u"pattern")
                {
                    rows = attributes.value("rows").toInt();
                    cols = attributes.value("cols").toInt();
                    width = attributes.value("width").toDouble();
                    height = attributes.value("height").toDouble();
                    sensorWidth = attributes.value("sensor_width").toDouble();
                    sensorHeight = attributes.value("sensor_height").toDouble();
                }
                else if(xml.name() == u"calibration")
                {
                    cx = attributes.value("cx").toDouble();
                    cy = attributes.value("cy").toDouble();
                    fx = attributes.value("fx").toDouble();
                    fy = attributes.value("fy").toDouble();

                    k1 = attributes.value("k1").toDouble();
                    k2 = attributes.value("k2").toDouble();
                    p1 = attributes.value("p1").toDouble();
                    p2 = attributes.value("p2").toDouble();
                    k3 = attributes.value("k3").toDouble();

                    isCalibrated = attributes.value("calibrated").toInt();
                }
                else if(xml.name() == u"images")
                {
                    nbImages = attributes.value("count").toInt();
                    state = 2;
                }
                break;

            case 2:
                if (xml.name() == u"image") {
                    addImage(toAbsolutePath(attributes.value("filename").toString()));
                    nbPoints = attributes.value("points").toInt();
                    if (imageFilenames.count() == nbImages && nbPoints == 0) {
                        state = 1;
                    } else if (nbPoints > 0) {
                        state = 3;
                        points.clear();
                    }
                }
                break;

            case 3:
                if (xml.name() == u"point") {
                    double x = attributes.value("x").toDouble();
                    double y = attributes.value("y").toDouble();
                    points.append(QPointF(x,y));

                    if (points.count() == nbPoints) {
                        imagePoints[imageFilenames.count()-1] = points;
                        if (imageFilenames.count() == nbImages) {
                            state = 1;
                        } else {
                            state = 2;
                        }
                    }
                }
                break;
        }
    }

    xmlCalibration.close();
}

void Project::removeImage(int index) {
	imageFilenames.removeAt(index);
	imagePoints.removeAt(index);
}

void Project::save() {

	if (filename.isEmpty()) {
		return;
	}

	// Project absolute filename path
	QString filepath = projectpath.absoluteFilePath(filename);

	QFile xmlCalibration(filepath);
    if (!xmlCalibration.open(QIODevice::ReadWrite)) {
        return;
    }

    QXmlStreamWriter stream(&xmlCalibration);
    stream.setAutoFormatting(true);
    stream.writeStartDocument();

    stream.writeStartElement("project");
    stream.writeAttribute("name", name);

    stream.writeStartElement("pattern");
    stream.writeAttribute("rows", QString::number(rows));
    stream.writeAttribute("cols", QString::number(cols));
    stream.writeAttribute("width", QString::number(width,'f',6));
    stream.writeAttribute("height", QString::number(height,'f',6));
    stream.writeAttribute("sensor_width", QString::number(sensorWidth,'f',6));
    stream.writeAttribute("sensor_height", QString::number(sensorHeight,'f',6));
    stream.writeEndElement();

    stream.writeStartElement("calibration");
    stream.writeAttribute("calibrated", isCalibrated ? "1" : "0");
    stream.writeAttribute("cx", QString::number(cx,'f',6));
    stream.writeAttribute("cy", QString::number(cy,'f',6));
    stream.writeAttribute("fx", QString::number(fx,'f',6));
    stream.writeAttribute("fy", QString::number(fy,'f',6));
    stream.writeAttribute("k1", QString::number(k1,'f',6));
    stream.writeAttribute("k2", QString::number(k2,'f',6));
    stream.writeAttribute("p1", QString::number(p1,'f',6));
    stream.writeAttribute("p2", QString::number(p2,'f',6));
    stream.writeAttribute("k3", QString::number(k3,'f',6));
    stream.writeEndElement();

    stream.writeStartElement("images");
    stream.writeAttribute("count", QString::number(imageFilenames.count()));
    for (int i = 0; i < imageFilenames.count(); i++) {
        stream.writeStartElement("image");
        stream.writeAttribute("filename", toRelativePath(imageFilenames.at(i)));
        int nbPoints = imagePoints.at(i).count();
        stream.writeAttribute("points", QString::number(nbPoints));
        for (int j = 0; j < nbPoints; j++) {
            const QPointF &point = imagePoints.at(i).at(j);
            stream.writeStartElement("point");
            stream.writeAttribute("x", QString::number(point.x(),'f',6));
            stream.writeAttribute("y", QString::number(point.y(),'f',6));
            stream.writeEndElement();
        }
        stream.writeEndElement();
    }
    stream.writeEndElement();

    stream.writeEndElement();
    stream.writeEndDocument();

    xmlCalibration.close();
}

void Project::setFilename(QString filepath) {
	QFileInfo fileInfo(filepath);
	filename = fileInfo.fileName();
	projectpath = fileInfo.absoluteDir();
}

void Project::setImagePoints(const QVector<QPointF> &points, int index) {
	imagePoints[index] = points;
}

QString Project::toAbsolutePath(QString filepath) {
	return projectpath.absoluteFilePath(filepath);
}

QString Project::toRelativePath(QString filepath) {
	return projectpath.relativeFilePath(filepath);
}

