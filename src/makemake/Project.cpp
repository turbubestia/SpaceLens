/*
 * Project.cpp
 *
 *  Created on: 1 abr. 2020
 *      Author: claud
 */

#include <QtCore/QtCore>
#include <QtWidgets/QtWidgets>

#include <algorithm>
#include <exiv2/exiv2.hpp>
#include <opencv4/opencv2/opencv.hpp>

#include "Project.h"
#include "Utils.h"

Project::Project() {
	sortCriteria = ByName;

	batch = false;

	groupCriteria = ByCount;
	groupByCountSize = 0;
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
	batch = false;
	imageFilenames.clear();
	imageGroups.clear();
}

void Project::groupByCount() {
	for (int i = 0; i < imageFilenames.count(); ) {
		QVector<int> group;
		for (int j = 0; j < groupByCountSize && i < imageFilenames.count(); j++, i++) {
			group.append(i);
		}
		imageGroups.append(group);
	}
}

int Project::groupCount() {
	return imageGroups.count();
}

int Project::groupOf(QString filename) {
    if (batch) {
        for (int g = 0; g < imageGroups.count(); g++) {
            for (int f = 0; f < imageGroups.at(g).count(); f++) {
                int index = imageGroups.at(g).at(f);
                if (filename == imageFilenames.at(index)) {
                    return g;
                }
            }
        }
    }

    return -1;
}

int Project::indexOf(int group, QString filename) {
    if (group == -1) {
        for (int f = 0; f < imageFilenames.count(); f++) {
            if (filename == imageFilenames.at(f)) {
                return f;
            }
        }
    } else {
        for (int f = 0; f < imageGroups.at(group).count(); f++) {
            int index = imageGroups.at(group).at(f);
            if (filename == imageFilenames.at(index)) {
                return f;
            }
        }
    }

    return -1;
}

QString Project::getImageFilename(int index) {
    return imageFilenames.at(index);
}

QStringList Project::getImageFilenames(int group) {

	if (group == -1) {
		return imageFilenames;
	}
	else if (batch && group < imageGroups.count()) {

	    const QVector<int> &images = imageGroups.at(group);

	    QStringList temp;
	    for (int i = 0; i < images.count(); i++) {
			temp.append(imageFilenames.at(images.at(i)));
		}
	    return temp;
	}

	return QStringList();
}

int Project::getImageIndex(int subgroup, int group) {
    if (group == -1) {
        return subgroup;
    } else {
        return imageGroups.at(group).at(subgroup);
    }
}

Project::GroupCriteria Project::getGroupCriteria() {
    return groupCriteria;
}

int Project::getGroupByCountSize() {
    return groupByCountSize;
}

QString Project::getProjectPath() {
	return projectpath.path();
}

Project::SortCriteria Project::getSortCriteria() {
    return sortCriteria;
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

bool Project::isBatchEnabled() {
	return batch;
}

void Project::loadDefault() {
    name = "Untitled";

    sortCriteria = ByName;

    batch = false;

    groupCriteria = ByCount;
    groupByCountSize = 1;

    imageFilenames.clear();
    imageGroups.clear();
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
            if (xml.name() == u"images" || xml.name() == u"batch") {
				state = 1;
			}
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
                if(xml.name() == u"images") {
					nbImages = attributes.value("count").toInt();
					state = 2;
				}
                else if (xml.name() == u"sortCriteria") {
                    QString sc = attributes.value("value").toString();
                    if (sc == "ByName") {
                        sortCriteria = ByName;
                    } else if (sc == "ByDate") {
                        sortCriteria = ByDate;
                    } else {
                        sortCriteria = ByName;
                    }
                }
                else if(xml.name() == u"batch") {
					batch = attributes.value("enabled").toInt() == 1;
					state = 3;
				}

				break;

			case 2:
                if (xml.name() == u"image") {
					QString filename = attributes.value("filename").toString();
					imageFilenames.append(toAbsolutePath(filename));
				}
				break;

			case 3:
                if (xml.name() == u"groupCriteria") {
					QString gc = attributes.value("value").toString();
					if (gc == "ByCount") {
						groupCriteria = ByCount;
					} else {
						groupCriteria = ByCount;
					}
                } else if (xml.name() == u"groupByCountSize") {
					groupByCountSize = attributes.value("value").toInt();
				}
		}
	}

	xmlCalibration.close();

	updateList();
}

void Project::removeImage(int index) {
	imageFilenames.removeAt(index);

	updateList();
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
		stream.writeAttribute("id", QString::number(i));
		stream.writeAttribute("filename", toRelativePath(imageFilenames.at(i)));
		stream.writeEndElement();
	}
	stream.writeEndElement(); // </images>

	stream.writeStartElement("sortCriteria");
    switch (sortCriteria) {
        case ByName: stream.writeAttribute("value", "ByName"); break;
        case ByDate: stream.writeAttribute("value", "ByDate"); break;
    }
    stream.writeEndElement();

	stream.writeStartElement("batch");
	stream.writeAttribute("enabled", batch ? "1" : "0");
	stream.writeStartElement("groupCriteria");
	switch (groupCriteria) {
		case ByCount: stream.writeAttribute("value", "ByCount"); break;
	}
	stream.writeEndElement();
	stream.writeStartElement("groupByCountSize");
	stream.writeAttribute("value", QString::number(groupByCountSize));
	stream.writeEndElement();
	stream.writeEndElement(); // </batch>

	stream.writeEndElement(); // </project>

	stream.writeEndDocument();

	xmlCalibration.close();
}

void Project::setBatchEnable(bool enable) {
	batch = enable;
}

void Project::setFilename(QString filepath) {
	QFileInfo fileInfo(filepath);
	filename = fileInfo.fileName();
	projectpath = fileInfo.absoluteDir();
}

void Project::setSortCriteria(Project::SortCriteria criteria) {
	sortCriteria = criteria;
}

void Project::setGroupByCount(int count) {
	groupCriteria = ByCount;
	groupByCountSize = count;
}

QString Project::toAbsolutePath(QString filepath) {
	return projectpath.absoluteFilePath(filepath);
}

QString Project::toRelativePath(QString filepath) {
	return projectpath.relativeFilePath(filepath);
}

void Project::updateList() {
    QVector<SortObjectWrapper> sortedFilenames;
    for (int i = 0; i < imageFilenames.count(); i++) {
        sortedFilenames.append(SortObjectWrapper(imageFilenames.at(i), sortCriteria));
    }
    std::sort(sortedFilenames.begin(), sortedFilenames.end());

    imageFilenames.clear();
    imageGroups.clear();

    for (int i = 0; i < sortedFilenames.count(); i++) {
        imageFilenames.append(sortedFilenames.at(i).getFilename());
    }

	if (batch) {
		switch(groupCriteria) {
			case ByCount:
				groupByCount();
				break;
		}
	}
}

