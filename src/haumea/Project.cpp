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
#include "Utils.h"

Project::Project() {

}


Project::~Project() {

}

bool Project::addImage(QString filename){
	// Internally the image are always stored in absolute mode
	filename = projectpath.absoluteFilePath(filename);
	// Don't add duplicates
	if (!imageFilenames.contains(filename)) {
		imageFilenames.append(filename);
		return true;
	}

	return false;
}

void Project::close() {
	filename.clear();
	imageFilenames.clear();
}


const QStringList& Project::getImagesFilenameList() {
	return imageFilenames;
}

QString Project::getProjectPath() {
	return projectpath.path();
}

QDir Project::getProjectDir() {
	return projectpath;
}

QString Project::getProjectFilename() {
	QFileInfo fileinfo(projectpath.absoluteFilePath(filename));
	return fileinfo.baseName();
}

int Project::imageCount() {
	return imageFilenames.count();
}

bool Project::isSaved() {
	return !filename.isEmpty();
}

void Project::loadDefault() {
    name = "Untitled";
    imageFilenames.clear();
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
                if(xml.name() == u"images")
				{
					nbImages = attributes.value("count").toInt();
					state = 2;
				}
				break;

			case 2:
                if (xml.name() == u"image") {
					addImage(toAbsolutePath(attributes.value("filename").toString()));
					if (imageFilenames.count() == nbImages) {
						state = 1;
					}
				}
				break;
		}
	}

	xmlCalibration.close();
}

void Project::removeImage(int index) {
	imageFilenames.removeAt(index);
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

	stream.writeStartElement("images");
	stream.writeAttribute("count", QString::number(imageFilenames.count()));
	for (int i = 0; i < imageFilenames.count(); i++) {
		stream.writeStartElement("image");
		stream.writeAttribute("filename", toRelativePath(imageFilenames.at(i)));
		stream.writeEndElement();
	}
	stream.writeEndElement(); // </images>

	stream.writeEndElement(); // </project>

	stream.writeEndDocument();

	xmlCalibration.close();
}

void Project::setFilename(QString filepath) {
	QFileInfo fileInfo(filepath);
	filename = fileInfo.fileName();
	projectpath = fileInfo.absoluteDir();
}

QString Project::toAbsolutePath(QString filepath) {
	return projectpath.absoluteFilePath(filepath);
}

QString Project::toRelativePath(QString filepath) {
	return projectpath.relativeFilePath(filepath);
}

