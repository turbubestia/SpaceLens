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

#include <CameraCalibrationParameters.h>

class Project {
	public:
		Project();
		virtual ~Project();

		enum SortCriteria {ByName, ByDate};
		enum GroupCriteria {ByCount};

		bool addImage(QString filename);
		void close();
		QString getImageFilename(int index);
		QStringList getImageFilenames(int group = -1);
		int getImageIndex(int subgroup, int group = -1);
		Project::GroupCriteria getGroupCriteria();
		int getGroupByCountSize();
		QDir getProjectDir();
		QString getProjectFilename();
		QString getProjectPath();
		Project::SortCriteria getSortCriteria();
		int groupCount();
		int groupOf(QString filename);
		int imageCount();
		int indexOf(int group, QString filename);
		bool isSaved();
		bool isBatchEnabled();
		void loadDefault();
		void open(QString filename);
		void removeImage(int index);
		void save();
		void setBatchEnable(bool enable);
		void setFilename(QString filename);
		void setSortCriteria(Project::SortCriteria criteria);
		void setGroupByCount(int count);
		QString toAbsolutePath(QString filepath);
		QString toRelativePath(QString filepath);
		void updateList();

	private:
		void groupByCount();

	public:
		// Properties to be save on file
		QString name;

	private:
		// Properties used to control the state of the project
		QDir projectpath;
		QString filename;

		SortCriteria sortCriteria;

		bool batch;

		GroupCriteria groupCriteria;
		int groupByCountSize;

		QStringList imageFilenames;
		QVector<QVector<int>> imageGroups;
};

class SortObjectWrapper {
	public:
		SortObjectWrapper() {
			criteria = Project::ByName;
		}

		SortObjectWrapper(const SortObjectWrapper& other) {
			filename = other.filename;
			criteria = other.criteria;
		}

		SortObjectWrapper(const QString &_filename, Project::SortCriteria _criteria)
			: filename(_filename), criteria(_criteria) {

		}

		const QString &getFilename() const {
			return filename;
		}

		bool operator<(const SortObjectWrapper& other) const {
			QFileInfo a(filename);
			QFileInfo b(other.filename);

			switch (criteria) {
				case Project::ByName:
					return a.baseName() < b.baseName();
				case Project::ByDate:
					return a.lastModified() < b.lastModified();
				default:
					return false;
			}
		}

		SortObjectWrapper &operator=(const SortObjectWrapper& other) {
			filename = other.filename;
			criteria = other.criteria;
			return *this;
		}

	private:
		QString filename;
		Project::SortCriteria criteria;
};


#endif /* SRC_PROJECT_H_ */
