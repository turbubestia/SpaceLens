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
		const QStringList& getImagesFilenameList();
		QDir getProjectDir();
		QString getProjectFilename();
		QString getProjectPath();
		int imageCount();
		bool isSaved();
		void loadDefault();
		void open(QString filename);
		void removeImage(int index);
		void save();
		void setFilename(QString filename);
		QString toAbsolutePath(QString filepath);
		QString toRelativePath(QString filepath);

	public:
		// Properties to be save on file
		QString name;

		QRect roi;

	private:
		// Properties used to control the state of the project
		QDir projectpath;
		QString filename;

		QStringList imageFilenames;
};


#endif /* SRC_PROJECT_H_ */
