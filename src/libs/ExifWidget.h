/*
 * ExifWidget.h
 *
 *  Created on: 14 abr. 2020
 *      Author: claud
 */

#ifndef SRC_EXIFWIDGET_H_
#define SRC_EXIFWIDGET_H_

#include <QtWidgets/QtWidgets>

class ExifWidget : public QWidget {
		Q_OBJECT

	public:
		ExifWidget(QWidget *parent = nullptr);
		~ExifWidget();

		void setImageFilename(const QString &filename);

	public:
		QString maker;
		QString model;
		QString lense;

		double focalLength;

	private:
		QString filename;

		QLabel *lMake;
		QLabel *lModel;
		QLabel *lLense;
		QLabel *lIso;
		QLabel *lAperture;
		QLabel *lSpeed;
		QLabel *lFocalLength;
		QLabel *lWidth;
		QLabel *lHeight;
};



#endif /* SRC_EXIFWIDGET_H_ */
