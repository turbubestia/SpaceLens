/*
 * ExifWidget.cpp
 *
 *  Created on: 14 abr. 2020
 *      Author: Claudio Alvarez Barros
 *
 *      This file is part of CameraCalibration.
 *
 *   CameraCalibration is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   CameraCalibration is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with CameraCalibration.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <QtCore/QtCore>
#include <QtWidgets/QtWidgets>
#include <exiv2/exiv2.hpp>
#include <math.h>
#include "ExifWidget.h"

ExifWidget::ExifWidget(QWidget *parent)
	: QWidget(parent) {

	int row = 0;

	lMake = new QLabel("");
	lModel = new QLabel("");
	lLense = new QLabel("");
	lIso = new QLabel("");
	lAperture = new QLabel("");
	lSpeed = new QLabel("");
	lFocalLength = new QLabel("");
	lWidth = new QLabel("");
	lHeight = new QLabel("");

	QGridLayout *layout = new QGridLayout;
	layout->addWidget(new QLabel("<b>Make</b>"), row, 0);
	layout->addWidget(lMake, row++, 1);
	layout->addWidget(new QLabel("<b>Model</b>"), row, 0);
	layout->addWidget(lModel, row++, 1);
	layout->addWidget(new QLabel("<b>Model</b>"), row, 0);
	layout->addWidget(lLense, row++, 1);
	layout->addWidget(new QLabel("<b>ISO</b>"), row, 0);
	layout->addWidget(lIso, row++, 1);
	layout->addWidget(new QLabel("<b>Aperture</b>"), row, 0);
	layout->addWidget(lAperture, row++, 1);
	layout->addWidget(new QLabel("<b>Speed</b>"), row, 0);
	layout->addWidget(lSpeed, row++, 1);
	layout->addWidget(new QLabel("<b>Focal Length</b>"), row, 0);
	layout->addWidget(lFocalLength, row++, 1);
	layout->addWidget(new QLabel("<b>Width</b>"), row, 0);
	layout->addWidget(lWidth, row++, 1);
	layout->addWidget(new QLabel("<b>Height</b>"), row, 0);
	layout->addWidget(lHeight, row++, 1);

	layout->addItem(new QSpacerItem(0,10, QSizePolicy::Expanding, QSizePolicy::Expanding), row, 0);

	setLayout(layout);

	focalLength = 0;
}

ExifWidget::~ExifWidget() {

}

void ExifWidget::setImageFilename(const QString &filename) {
	this->filename = filename;

	Exiv2::Image::AutoPtr image = Exiv2::ImageFactory::open(filename.toStdString());

	image->readMetadata();
	Exiv2::ExifData &exifData = image->exifData();

	Exiv2::ExifData::const_iterator end = exifData.end();
	for (Exiv2::ExifData::const_iterator i = exifData.begin(); i != end; ++i) {
		QString key = QString::fromStdString(i->key());

		if(key == "Exif.Image.Make") {
			maker = QString::fromStdString(i->value().toString());
			lMake->setText(maker);
		}
		else if(key == "Exif.Image.Model") {
			model = QString::fromStdString(i->value().toString());
			lModel->setText(model);
		}
		else if (key == "Exif.Photo.LensModel") {
			lense = QString::fromStdString(i->value().toString());
			lLense->setText(lense);
		}
		else if(key == "Exif.Photo.ISOSpeedRatings") {
			int iso = i->value().toLong();
			lIso->setText(QString("%1").arg(iso));
		}
		else if (key == "Exif.Photo.FNumber") {
			Exiv2::Rational value = i->value().toRational();
			double fstop = (double)value.first / (double)value.second;
			lAperture->setText(QString("1/%1").arg(fstop,0,'f',1));
		}
		else if (key == "Exif.Photo.ExposureTime") {
			Exiv2::Rational value = i->value().toRational();
			lSpeed->setText(QString("%1/%2 s").arg(value.first).arg(value.second));
		}
		else if (key == "Exif.Photo.FocalLength") {
			Exiv2::Rational value = i->value().toRational();
			focalLength = (double)value.first / (double)value.second;
			lFocalLength->setText(QString("%1 mm").arg(focalLength));
		}
		else if (key == "Exif.Photo.PixelXDimension") {
			int value = i->value().toLong();
			lWidth->setText(QString("%1 px").arg(value));
		}
		else if (key == "Exif.Photo.PixelYDimension") {
			int value = i->value().toLong();
			lHeight->setText(QString("%1 px").arg(value));
		}
	}

	Exiv2::XmpData &xmpData = image->xmpData();
	Exiv2::XmpData::const_iterator xend = xmpData.end();
	for (Exiv2::XmpData::const_iterator i = xmpData.begin(); i != xend; ++i) {
		QString key = QString::fromStdString(i->key());
	}
}
