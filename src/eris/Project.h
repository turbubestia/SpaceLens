/*
 * Project.h
 *
 *  Created on: 1 abr. 2020
 *      Author: claud
 */

#ifndef SRC_PROJECT_H_
#define SRC_PROJECT_H_

#include <QtCore/QtCore>

#include <opencv4/opencv2/opencv.hpp>

class Project {
	public:
		Project();
		virtual ~Project();

		bool addImage(QString filename);
		void close();
		const QVector<QPointF> & getImagePoints(int row);
		const QStringList& getImagesFilenameList();
		QString getProjectFilename();
		QString getProjectPath();
		int imageCount();
		bool isSaved();
		void loadDefault();
		void open(QString filename);
		void removeImage(int index);
		void save();
		void setFilename(QString filename);
		void setImagePoints(const QVector<QPointF> &points, int index);
		QString toAbsolutePath(QString filepath);
		QString toRelativePath(QString filepath);

	public:
		// Properties to be save on file
		QString name;
		double sensorWidth;
		double sensorHeight;
		int rows;
		int cols;
		double width;
		double height;
		double cx, cy;
		double fx, fy;
		double k1,k2,p1,p2,k3;
		bool isCalibrated;

	private:
		// Properties used to control the state of the project
		QDir projectpath;
		QString filename;

		QStringList imageFilenames;
		QList<QVector<QPointF>> imagePoints;
};


#endif /* SRC_PROJECT_H_ */
