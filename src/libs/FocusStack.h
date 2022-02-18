/*
 * FocusStack.h
 *
 *  Created on: Aug 29, 2020
 *      Author: claud
 */

#ifndef SRC_LIBS_FOCUSSTACK_H_
#define SRC_LIBS_FOCUSSTACK_H_

#include <QtCore/QtCore>

#include <opencv4/opencv2/opencv.hpp>

class FocusStackWorker;

class FocusStack : public QObject {
    Q_OBJECT

    public:
        FocusStack(QObject *parent = nullptr);
        virtual ~FocusStack();

        void addStackGroup(QStringList stack);
        void clear();
        void setCacheFolder(const QString &folder);
        void setOutputFolder(const QString &folder);

    public slots:
        void startFocusStack();

    signals:
        void imageFocused(int group);
        void imageAligned(int group, QString filename);
        void finished();

    private slots:
        void workerFinished();

    private:
        QString cacheFolder;
        QString outputFolder;

        QList<QStringList> groups;
        QList<FocusStackWorker *> workerList;
        int finishedJobs;
        int workersSpanwned;
};

class FocusStackWorker : public QThread {
    Q_OBJECT

    public:
        FocusStackWorker(QObject *parent = nullptr);

        void clear();
        int count();
        void addGroup(int id, const QStringList &filenames);
        void setCacheFolder(const QString &folder);
        void setOutputFolder(const QString &folder);

    signals:
        void imageFocused(int group);
        void imageAligned(int group, QString filename);
        void workerFinished();

    protected:
        void run() override;

    private:
        void align();
        void coarseDepthMap(cv::Mat &dst);
        void downSampleCoarseDepthMap(cv::Mat &src, cv::Mat &dst);
        void fineDepthMap(cv::Mat &src, cv::Mat &dst);
        void upSampleFineDepthMap(cv::Mat &src, cv::Mat &dst, cv::Size size);
        void stack(cv::Mat &src, cv::Mat &dst);

    private:
        QString cacheFolder;
        QStringList cacheFilenames;

        QString outputFolder;

        QVector<int> groupIds;
        QList<QStringList> groups;

        int currentGroup;
        std::string prefix;
};

#endif /* SRC_LIBS_FOCUSSTACK_H_ */
