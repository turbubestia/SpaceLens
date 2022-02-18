/*
 * CameraCalibration.cpp
 *
 *  Created on: Jun 12, 2020
 *      Author: claud
 */

#include <QtCore/QtCore>
#include "CameraCalibration.h"
#include "utils.h"

#include <opencv4/opencv2/opencv.hpp>
#include <opencv4/opencv2/core/types_c.h>

bool CustomFindChessboardCorners(cv::Mat &img, cv::Size pattern_size, std::vector<cv::Point2f> &corners, int flags);
void binarizeCvImage(const cv::Mat &gray, cv::Mat &binary);

void dump_vector_double(const std::vector<double> &src) {
	QFile fid("std_vector.csv");
	fid.open(QFile::WriteOnly);
	if (!fid.isOpen()) {
		return;
	}
	for (int i = 0; i < (int)src.size(); i++) {
		fid.write(QString("%1\n").arg(src[i],0,'f',3).toLocal8Bit());
	}
	fid.close();
}

CameraCalibration::CameraCalibration(QObject *parent)
    : QObject(parent) {
    quadWidth = 10;
    quadHeight = 10;
    rows = 7;
    cols = 11;
    sensorWidth = 20;

    int nbWorkers = QThread::idealThreadCount();

    for (int i = 0; i < nbWorkers; i++) {
        CameraCalibrationWorker *worker = new CameraCalibrationWorker(i);
        connect(worker, &CameraCalibrationWorker::imageCalibrated, this, &CameraCalibration::imageCalibrated);
        connect(worker, &CameraCalibrationWorker::workerFinished, this, &CameraCalibration::workerFinished);
        workerList.append(worker);
    }

    finishedJobs = 0;
}

CameraCalibration::~CameraCalibration() {

}

void CameraCalibration::setCols(double value) {
    cols = value;
}

void CameraCalibration::setImageFilenames(const QStringList &_filenames) {
    filenameList = _filenames;
}

void CameraCalibration::setQuadWidth(double value) {
    quadWidth = value;
}

void CameraCalibration::setPattern(QSize dimensions, QSize corners) {
    quadWidth = dimensions.width();
    quadHeight = dimensions.height();
    cols = corners.width();
    rows = corners.height();
}

QImage CameraCalibration::binaryzeImage(const QImage &image) {
	cv::Mat src = qImageToCvMat(image);

	cv::Mat gray, binary;

	cvtColor(src, gray, cv::COLOR_BGRA2GRAY);
	binarizeCvImage(gray, binary);
	cvtColor(binary, gray, cv::COLOR_GRAY2BGR);

	return cvMatToQImage(gray);
}

QVector<QPointF> CameraCalibration::getCorners(int index) {

	int nbWorkers = workerList.count();
	int groupSize = ceil((double) filenameList.count() / nbWorkers);
	if (groupSize < 1) {
		groupSize = 1;
	}
	int group = index / groupSize;
	int subindex = index % groupSize;

	const std::vector<cv::Point2f> &cv_corners = workerList.at(group)->getCorners(subindex);
	QVector<QPointF> qt_corners(cv_corners.size());

	for (int i = 0; i < (int)cv_corners.size(); i++) {
		qt_corners[i] = QPointF(cv_corners.at(i).x, cv_corners.at(i).y);
	}

	return qt_corners;
}

QImage CameraCalibration::undistorded(int index) {
	cv::Mat src = cv::imread(filenameList.at(index).toStdString());
	cv::Mat dst;

	cv::undistort(src, dst, cameraMatrix, distCoeffs);

	QImage image = cvMatToQImage(dst);

	return image;
}

void CameraCalibration::getMatrixCoefficients(double *coeff) {
	coeff[0] = cameraMatrix.at<double>(0,2);
	coeff[1] = cameraMatrix.at<double>(1,2);
	coeff[2] = cameraMatrix.at<double>(0,0);
	coeff[3] = cameraMatrix.at<double>(1,1);
}

void CameraCalibration::getDistortionCoefficients(double *coeff) {
	coeff[0] = distCoeffs.at<double>(0);
	coeff[1] = distCoeffs.at<double>(1);
	coeff[2] = distCoeffs.at<double>(2);
	coeff[3] = distCoeffs.at<double>(3);
	coeff[4] = distCoeffs.at<double>(4);
}

void CameraCalibration::setMatrixCoefficients(double *coeff) {
	cameraMatrix = cv::Mat::zeros(cv::Size(3,3), CV_64F);

	cameraMatrix.at<double>(0,2) = coeff[0];
	cameraMatrix.at<double>(1,2) = coeff[1];
	cameraMatrix.at<double>(0,0) = coeff[2];
	cameraMatrix.at<double>(1,1) = coeff[3];
	cameraMatrix.at<double>(2,2) = 1;
}

void CameraCalibration::setDistortionCoefficients(double *coeff) {
	distCoeffs = cv::Mat::zeros(cv::Size(1,5), CV_64F);

	distCoeffs.at<double>(0,0) = coeff[0];
	distCoeffs.at<double>(0,1) = coeff[1];
	distCoeffs.at<double>(0,2) = coeff[2];
	distCoeffs.at<double>(0,3) = coeff[3];
	distCoeffs.at<double>(0,4) = coeff[4];
}

void CameraCalibration::setQuadHeight(double value) {
    quadHeight = value;
}

void CameraCalibration::setRows(double value) {
    rows = value;
}

void CameraCalibration::setSensorWidth(double value) {
    sensorWidth = value;
}

void CameraCalibration::startCalibration() {

    int nbWorkers = workerList.count();
    int groupSize = ceil((double) filenameList.count() / nbWorkers);
    if (groupSize < 1) {
        groupSize = 1;
    }

    finishedJobs = 0;

    for (int i = 0, k = 0; i < nbWorkers && k < filenameList.count(); i++)
    {
        int currentGroupSize = groupSize;
        int remElements = filenameList.count() - k;
        if (currentGroupSize > remElements) {
            currentGroupSize = remElements;
        }

        QStringList groupFilenames = filenameList.mid(k, currentGroupSize);
        workerList.at(i)->setFilenames(groupFilenames);
        workerList.at(i)->setup(rows, cols);
        workerList.at(i)->start();

        k += currentGroupSize;
    }
}

void CameraCalibration::workerFinished() {
    finishedJobs++;
    if (finishedJobs == workerList.count()) {

    	// List of point detected in the images
		qDebug() << "List of point detected in the images";
		std::vector<std::vector<cv::Point2f>> image_points_list;
		for (int w = 0; w < workerList.count(); w++) {
			for (int i = 0; i < workerList.at(i)->count(); i++) {
				const std::vector<cv::Point2f> &corners = workerList.at(i)->getCorners(i);
				if (corners.size() > 0) {
					image_points_list.push_back(corners);
				}
			}
		}

		// List of point in the real world
		qDebug() << "List of point in the real world";
		std::vector<cv::Point3f> world_points;
		for (int col = 0; col < cols; col++) {
			for (int row = 0; row < rows ; row++) {
				world_points.push_back(cv::Point3f((float) row * quadHeight, (float) col * quadWidth, 0));
			}
		}

		// Map a world coordinate for each image
		qDebug() << "Map a world coordinate for each image";
		std::vector<std::vector<cv::Point3f>> world_points_list;
		for (int i = 0; i < (int) image_points_list.size(); i++) {
			world_points_list.push_back(world_points);
		}

		// Open the first image to obtain the image size
		cv::Mat image = cv::imread(filenameList.at(0).toStdString());

		// Calibrate
		qDebug() << "Calibrate camera";
		cv::Mat R, T;

		cameraMatrix = cv::Mat::eye(3,3,CV_64F);
		cameraMatrix.at<double>(0,2) = (image.cols-1)/2;
		cameraMatrix.at<double>(1,2) = (image.rows-1)/2;

		int flag = cv::CALIB_FIX_K4 | cv::CALIB_FIX_K5 | cv::CALIB_FIX_PRINCIPAL_POINT | cv::CALIB_FIX_ASPECT_RATIO;
		cv::calibrateCamera(world_points_list, image_points_list, image.size(), cameraMatrix, distCoeffs, R, T, flag);
		//cameraMatrix = cv::getOptimalNewCameraMatrix(cameraMatrix, distCoeffs, image.size(), 1.0, cv::Size(), nullptr, true);

        emit calibrationFinished();
    }
}

CameraCalibrationWorker::CameraCalibrationWorker(int _group, QObject *parent)
    : QThread(parent) {
    group = _group;
    rows = 0;
    cols = 0;
}

int CameraCalibrationWorker::count() {
	return filenames.count();
}

int CameraCalibrationWorker::getGroup() {
    return group;
}

void CameraCalibrationWorker::setup(int _rows, int _cols) {
    rows = _rows;
    cols = _cols;
}

void CameraCalibrationWorker::setFilenames(const QStringList &_filenames) {
    filenames = _filenames;
    corners.resize(filenames.size());
}

const std::vector<cv::Point2f> &CameraCalibrationWorker::getCorners(int index) {
	return corners.at(index);
}

void CameraCalibrationWorker::run() {
    for (int i = 0; i < filenames.count(); i++) {
        cv::Mat src = cv::imread(filenames.at(i).toStdString());

        if (src.data == nullptr) {
            emit imageCalibrated(filenames.at(i), false);
            continue;
        }

        cv::Mat gray, binary;
        cvtColor(src, gray, cv::COLOR_BGRA2GRAY);
        binarizeCvImage(gray, binary);

        std::vector<cv::Point2f> curCorners;

        // Search the coarse position of the chessboard inner curCorners
        bool found = CustomFindChessboardCorners( binary, cv::Size(rows,cols), curCorners, 0);

        if (found) {
            // If the chessboard was found, refine the corner position
            cv::TermCriteria criteria(CV_TERMCRIT_EPS | CV_TERMCRIT_ITER, 30, 0.001);
            cv::cornerSubPix(gray, curCorners, cv::Size(11,11), cv::Size(-1,-1), criteria);
        }

        corners[i] = curCorners;

        emit imageCalibrated(filenames.at(i), found);
    }

    emit workerFinished();
}

void binarizeCvImage(const cv::Mat &gray, cv::Mat &binary) {
	// We need a binary image to find a coarse corner position
	cv::Mat temp;

	// Fill quads to remove inside noise (white spot due to reflection or bad print quality of the pattern)
	cv::morphologyEx(gray, temp, cv::MORPH_OPEN, cv::getStructuringElement( cv::MORPH_ELLIPSE, cv::Size(7,7), cv::Point(3,3)));
	// Remove random noise
	cv::GaussianBlur(temp, binary, cv::Size(5,5), 0);

	int thr = findBimodalThreshold(binary);

	// Bi-mode binarization
	cv::threshold(binary, binary, thr, 255, cv::THRESH_BINARY);
}

//M*//////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                        Intel License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000, Intel Corporation, all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of Intel Corporation may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

/************************************************************************************\
    This is improved variant of chessboard corner detection algorithm that
    uses a graph of connected quads. It is based on the code contributed
    by Vladimir Vezhnevets and Philip Gruebele.
    Here is the copyright notice from the original Vladimir's code:
    ===============================================================

    The algorithms developed and implemented by Vezhnevets Vldimir
    aka Dead Moroz (vvp@graphics.cs.msu.ru)
    See http://graphics.cs.msu.su/en/research/calibration/opencv.html
    for detailed information.

    Reliability additions and modifications made by Philip Gruebele.
    <a href="mailto:pgruebele@cox.net">pgruebele@cox.net</a>

    Some further improvements for detection of partially ocluded boards at non-ideal
    lighting conditions have been made by Alex Bovyrin and Kurt Kolonige

\************************************************************************************/

/************************************************************************************\
  This version adds a new and improved variant of chessboard corner detection
  that works better in poor lighting condition. It is based on work from
  Oliver Schreer and Stefano Masneri. This method works faster than the previous
  one and reverts back to the older method in case no chessboard detection is
  possible. Overall performance improves also because now the method avoids
  performing the same computation multiple times when not necessary.

\************************************************************************************/

/************************************************************************************\
 This version remove the the older method to prevent the algorithm to freeze
 indefinitedly. The image binarization method was also changed.
\************************************************************************************/

#include <stack>
#include <algorithm>
#include <utility>

using namespace cv;

//=====================================================================================
// Implementation for the enhanced calibration object detection
//=====================================================================================

#define MAX_CONTOUR_APPROX  7
#define USE_CV_FINDCONTOURS

struct QuadCountour {
    Point pt[4];
    int parent_contour;

    QuadCountour(const Point pt_[4], int parent_contour_) :
        parent_contour(parent_contour_)
    {
        pt[0] = pt_[0]; pt[1] = pt_[1]; pt[2] = pt_[2]; pt[3] = pt_[3];
    }
};

/** This structure stores information about the chessboard corner.*/
struct ChessBoardCorner
{
    cv::Point2f pt; // Coordinates of the corner
    int row; // Board row index
    int count; // Number of neighbor corners
    struct ChessBoardCorner* neighbors[4]; // Neighbor corners

    ChessBoardCorner(const cv::Point2f& pt_ = cv::Point2f()) :
        pt(pt_), row(0), count(0)
    {
        neighbors[0] = neighbors[1] = neighbors[2] = neighbors[3] = NULL;
    }

    float sumDist(int& n_) const
    {
        float sum = 0;
        int n = 0;
        for (int i = 0; i < 4; ++i)
        {
            if (neighbors[i])
            {
                sum += sqrt(normL2Sqr<float>(neighbors[i]->pt - pt));
                n++;
            }
        }
        n_ = n;
        return sum;
    }
};

/** This structure stores information about the chessboard quadrangle.*/
struct ChessBoardQuad
{
    int count; // Number of quad neighbors
    int group_idx; // quad group ID
    int row, col; // row and column of this quad
    bool ordered; // true if corners/neighbors are ordered counter-clockwise
    float edge_len; // quad edge len, in pix^2
    // neighbors and corners are synced, i.e., neighbor 0 shares corner 0
    ChessBoardCorner *corners[4]; // Coordinates of quad corners
    struct ChessBoardQuad *neighbors[4]; // Pointers of quad neighbors

    ChessBoardQuad(int group_idx_ = -1) :
        count(0),
        group_idx(group_idx_),
        row(0), col(0),
        ordered(0),
        edge_len(0)
    {
        corners[0] = corners[1] = corners[2] = corners[3] = NULL;
        neighbors[0] = neighbors[1] = neighbors[2] = neighbors[3] = NULL;
    }
};

class ChessBoardDetector {
    public:
        cv::Mat binarized_image;
        Size pattern_size;

        cv::AutoBuffer<ChessBoardQuad> all_quads;
        cv::AutoBuffer<ChessBoardCorner> all_corners;

        int all_quads_count;

        ChessBoardDetector(const Size &pattern_size_)
                : pattern_size(pattern_size_), all_quads_count(0) {
        }

        void reset() {
            all_quads.deallocate();
            all_corners.deallocate();
            all_quads_count = 0;
        }

        void generateQuads(const cv::Mat &image_, int flags);
        bool processQuads(std::vector<cv::Point2f> &out_corners, int &prev_sqr_size);
        void findQuadNeighbors();
        void findConnectedQuads(std::vector<ChessBoardQuad*> &out_group, int group_idx);
        int checkQuadGroup(std::vector<ChessBoardQuad*> &quad_group, std::vector<ChessBoardCorner*> &out_corners);
        int cleanFoundConnectedQuads(std::vector<ChessBoardQuad*> &quad_group);
        int orderFoundConnectedQuads(std::vector<ChessBoardQuad*> &quads);
        void orderQuad(ChessBoardQuad &quad, ChessBoardCorner &corner, int common);
        int addOuterQuad(ChessBoardQuad &quad, std::vector<ChessBoardQuad*> &quads);
        void removeQuadFromGroup(std::vector<ChessBoardQuad*> &quads, ChessBoardQuad &q0);
        bool checkBoardMonotony(const std::vector<cv::Point2f> &corners);
};

bool CustomFindChessboardCorners(cv::Mat &img, cv::Size pattern_size, std::vector<cv::Point2f> &corners, int flags)
{

    bool found = false;

    const int min_dilations = 0;
    const int max_dilations = 7;

    if (pattern_size.width <= 2 || pattern_size.height <= 2)
        CV_Error(Error::StsOutOfRange, "Both width and height of the pattern should have bigger than 2");

    int prev_sqr_size = 0;

    ChessBoardDetector detector(pattern_size);

    // Try our standard "1" dilation, but if the pattern is not found, iterate the whole procedure with higher dilations.
    // This is necessary because some squares simply do not separate properly with a single dilation.  However,
    // we want to use the minimum number of dilations possible since dilations cause the squares to become smaller,
    // making it difficult to detect smaller squares.
    for (int dilations = min_dilations; dilations <= max_dilations; dilations++)
    {
        dilate( img, img, Mat(), Point(-1, -1), 1 );

        // So we can find rectangles that go to the edge, we draw a white line around the image edge.
        // Otherwise FindContours will miss those clipped rectangle contours.
        // The border color will be the image mean, because otherwise we risk screwing up filters like cvSmooth()...
        rectangle( img, Point(0,0), Point(img.cols-1, img.rows-1), Scalar(255,255,255), 3, LINE_8);

        Mat binarized_img = img;

        detector.reset();
        detector.generateQuads(binarized_img, flags);

        if (detector.processQuads(corners, prev_sqr_size)) {
            found = true;
            break;
        }
    }

    if (found)
    {
        found = detector.checkBoardMonotony(corners);
    }

    // check that none of the found corners is too close to the image boundary
    if (found)
    {
        const int BORDER = 8;
        for (int k = 0; k < pattern_size.width*pattern_size.height; ++k)
        {
            if( corners[k].x <= BORDER || corners[k].x > img.cols - BORDER ||
                corners[k].y <= BORDER || corners[k].y > img.rows - BORDER )
            {
                found = false;
                break;
            }
        }
    }

    if (found)
    {
        if ((pattern_size.height & 1) == 0 && (pattern_size.width & 1) == 0 )
        {
            int last_row = (pattern_size.height-1)*pattern_size.width;
            double dy0 = corners[last_row].y - corners[0].y;
            if (dy0 < 0)
            {
                int n = pattern_size.width*pattern_size.height;
                for(int i = 0; i < n/2; i++ )
                {
                    cv::Point2f temp = corners[i];
                    corners[i] = corners[n-i-1];
                    corners[n-i-1] = temp;
                }
            }
        }
        cv::cornerSubPix(img, corners, Size(2, 2), Size(-1,-1),
                         cv::TermCriteria(TermCriteria::EPS + TermCriteria::MAX_ITER, 15, 0.1));
    }

    return found;
}

//
// Checks that each board row and column is pretty much monotonous curve:
// It analyzes each row and each column of the chessboard as following:
//    for each corner c lying between end points in the same row/column it checks that
//    the point projection to the line segment (a,b) is lying between projections
//    of the neighbor corners in the same row/column.
//
// This function has been created as temporary workaround for the bug in current implementation
// of cvFindChessboardCornes that produces absolutely unordered sets of corners.
//
bool ChessBoardDetector::checkBoardMonotony(const std::vector<cv::Point2f>& corners)
{
    for (int k = 0; k < 2; ++k)
    {
        int max_i = (k == 0 ? pattern_size.height : pattern_size.width);
        int max_j = (k == 0 ? pattern_size.width: pattern_size.height) - 1;
        for (int i = 0; i < max_i; ++i)
        {
            cv::Point2f a = k == 0 ? corners[i*pattern_size.width] : corners[i];
            cv::Point2f b = k == 0 ? corners[(i+1)*pattern_size.width-1]
                                   : corners[(pattern_size.height-1)*pattern_size.width + i];
            float dx0 = b.x - a.x, dy0 = b.y - a.y;
            if (fabs(dx0) + fabs(dy0) < FLT_EPSILON)
                return false;
            float prevt = 0;
            for (int j = 1; j < max_j; ++j)
            {
                cv::Point2f c = k == 0 ? corners[i*pattern_size.width + j]
                                       : corners[j*pattern_size.width + i];
                float t = ((c.x - a.x)*dx0 + (c.y - a.y)*dy0)/(dx0*dx0 + dy0*dy0);
                if (t < prevt || t > 1)
                    return false;
                prevt = t;
            }
        }
    }
    return true;
}

//
// order a group of connected quads
// order of corners:
//   0 is top left
//   clockwise from there
// note: "top left" is nominal, depends on initial ordering of starting quad
//   but all other quads are ordered consistently
//
// can change the number of quads in the group
// can add quads, so we need to have quad/corner arrays passed in
//
int ChessBoardDetector::orderFoundConnectedQuads(std::vector<ChessBoardQuad*>& quads)
{
    const int max_quad_buf_size = (int)all_quads.size();
    int quad_count = (int)quads.size();

    std::stack<ChessBoardQuad*> stack;

    // first find an interior quad
    ChessBoardQuad *start = NULL;
    for (int i = 0; i < quad_count; i++)
    {
        if (quads[i]->count == 4)
        {
            start = quads[i];
            break;
        }
    }

    if (start == NULL)
        return 0;   // no 4-connected quad

    // start with first one, assign rows/cols
    int row_min = 0, col_min = 0, row_max=0, col_max = 0;

    std::map<int, int> col_hist;
    std::map<int, int> row_hist;

    stack.push(start);
    start->row = 0;
    start->col = 0;
    start->ordered = true;

    // Recursively order the quads so that all position numbers (e.g.,
    // 0,1,2,3) are in the at the same relative corner (e.g., lower right).

    while (!stack.empty())
    {
        ChessBoardQuad* q = stack.top(); stack.pop(); CV_Assert(q);

        int col = q->col;
        int row = q->row;
        col_hist[col]++;
        row_hist[row]++;

        // check min/max
        if (row > row_max) row_max = row;
        if (row < row_min) row_min = row;
        if (col > col_max) col_max = col;
        if (col < col_min) col_min = col;

        for (int i = 0; i < 4; i++)
        {
            ChessBoardQuad *neighbor = q->neighbors[i];
            switch(i)   // adjust col, row for this quad
            {           // start at top left, go clockwise
            case 0:
                row--; col--; break;
            case 1:
                col += 2; break;
            case 2:
                row += 2;   break;
            case 3:
                col -= 2; break;
            }

            // just do inside quads
            if (neighbor && neighbor->ordered == false && neighbor->count == 4)
            {
                //DPRINTF("col: %d  row: %d", col, row);
                CV_Assert(q->corners[i]);
                orderQuad(*neighbor, *(q->corners[i]), (i+2)&3); // set in order
                neighbor->ordered = true;
                neighbor->row = row;
                neighbor->col = col;
                stack.push(neighbor);
            }
        }
    }

    // analyze inner quad structure
    int w = pattern_size.width - 1;
    int h = pattern_size.height - 1;
    int drow = row_max - row_min + 1;
    int dcol = col_max - col_min + 1;

    // normalize pattern and found quad indices
    if ((w > h && dcol < drow) ||
        (w < h && drow < dcol))
    {
        h = pattern_size.width - 1;
        w = pattern_size.height - 1;
    }


    // check if there are enough inner quads
    if (dcol < w || drow < h)   // found enough inner quads?
    {
        return 0;   // no, return
    }

#ifdef ENABLE_TRIM_COL_ROW
    // too many columns, not very common
    if (dcol == w+1)    // too many, trim
    {
        if (col_hist[col_max] > col_hist[col_min]) {
            trimCol(quads, col_min, -1);
        }
        else {
            trimCol(quads, col_max, +1);
        }
    }

    // too many rows, not very common
    if (drow == h+1)    // too many, trim
    {
        if (row_hist[row_max] > row_hist[row_min]) {
            trimRow(quads, row_min, -1);
        }
        else {
            trimRow(quads, row_max, +1);
        }
    }

    quad_count = (int)quads.size(); // update after icvTrimCol/icvTrimRow
#endif

    // check edges of inner quads
    // if there is an outer quad missing, fill it in
    // first order all inner quads
    int found = 0;
    for (int i=0; i < quad_count; ++i)
    {
        ChessBoardQuad& q = *quads[i];
        if (q.count != 4)
            continue;

        {   // ok, look at neighbors
            int col = q.col;
            int row = q.row;
            for (int j = 0; j < 4; j++)
            {
                switch(j)   // adjust col, row for this quad
                {           // start at top left, go clockwise
                    case 0:
                        row--; col--; break;
                    case 1:
                        col += 2; break;
                    case 2:
                        row += 2;   break;
                    case 3:
                        col -= 2; break;
                }
                ChessBoardQuad *neighbor = q.neighbors[j];
                if (neighbor && !neighbor->ordered && // is it an inner quad?
                    col <= col_max && col >= col_min &&
                    row <= row_max && row >= row_min)
                {
                    // if so, set in order
                    found++;
                    CV_Assert(q.corners[j]);
                    orderQuad(*neighbor, *q.corners[j], (j+2)&3);
                    neighbor->ordered = true;
                    neighbor->row = row;
                    neighbor->col = col;
                }
            }
        }
    }

    // if we have found inner quads, add corresponding outer quads,
    //   which are missing
    if (found > 0)
    {
        for (int i = 0; i < quad_count && all_quads_count < max_quad_buf_size; i++)
        {
            ChessBoardQuad& q = *quads[i];
            if (q.count < 4 && q.ordered)
            {
                int added = addOuterQuad(q, quads);
                quad_count += added;
            }
        }

        if (all_quads_count >= max_quad_buf_size)
            return 0;
    }


    // final trimming of outer quads
    if (dcol == w && drow == h) // found correct inner quads
    {
        for (int i = quad_count - 1; i >= 0; i--) // eliminate any quad not connected to an ordered quad
        {
            ChessBoardQuad& q = *quads[i];
            if (q.ordered == false)
            {
                bool outer = false;
                for (int j=0; j<4; j++) // any neighbors that are ordered?
                {
                    if (q.neighbors[j] && q.neighbors[j]->ordered)
                        outer = true;
                }
                if (!outer) // not an outer quad, eliminate
                {
                    removeQuadFromGroup(quads, q);
                }
            }

        }
        return (int)quads.size();
    }

    return 0;
}

// add an outer quad
// looks for the neighbor of <quad> that isn't present, tries to add it in.
// <quad> is ordered
int ChessBoardDetector::addOuterQuad(ChessBoardQuad& quad, std::vector<ChessBoardQuad*>& quads)
{
    int added = 0;
    int max_quad_buf_size = (int)all_quads.size();

    for (int i = 0; i < 4 && all_quads_count < max_quad_buf_size; i++) // find no-neighbor corners
    {
        if (!quad.neighbors[i])    // ok, create and add neighbor
        {
            int j = (i+2)&3;
            int q_index = all_quads_count++;
            ChessBoardQuad& q = all_quads[q_index];
            q = ChessBoardQuad(0);
            added++;
            quads.push_back(&q);

            // set neighbor and group id
            quad.neighbors[i] = &q;
            quad.count += 1;
            q.neighbors[j] = &quad;
            q.group_idx = quad.group_idx;
            q.count = 1;   // number of neighbors
            q.ordered = false;
            q.edge_len = quad.edge_len;

            // make corners of new quad
            // same as neighbor quad, but offset
            const cv::Point2f pt_offset = quad.corners[i]->pt - quad.corners[j]->pt;
            for (int k = 0; k < 4; k++)
            {
                ChessBoardCorner& corner = (ChessBoardCorner&)all_corners[q_index * 4 + k];
                const cv::Point2f& pt = quad.corners[k]->pt;
                corner = ChessBoardCorner(pt);
                q.corners[k] = &corner;
                corner.pt += pt_offset;
            }
            // have to set exact corner
            q.corners[j] = quad.corners[i];

            // now find other neighbor and add it, if possible
            int next_i = (i + 1) & 3;
            int prev_i = (i + 3) & 3; // equal to (j + 1) & 3
            ChessBoardQuad* quad_prev = quad.neighbors[prev_i];
            if (quad_prev &&
                quad_prev->ordered &&
                quad_prev->neighbors[i] &&
                quad_prev->neighbors[i]->ordered )
            {
                ChessBoardQuad* qn = quad_prev->neighbors[i];
                q.count = 2;
                q.neighbors[prev_i] = qn;
                qn->neighbors[next_i] = &q;
                qn->count += 1;
                // have to set exact corner
                q.corners[prev_i] = qn->corners[next_i];
            }
        }
    }
    return added;
}

// trimming routines
#ifdef ENABLE_TRIM_COL_ROW
void ChessBoardDetector::trimCol(std::vector<ChessBoardQuad*>& quads, int col, int dir)
{
    std::vector<ChessBoardQuad*> quads_(quads);
    // find the right quad(s)
    for (size_t i = 0; i < quads_.size(); ++i)
    {
        ChessBoardQuad& q = *quads_[i];

        if (q.ordered && q.col == col)
        {
            if (dir == 1)
            {
                if (q.neighbors[1])
                {
                    removeQuadFromGroup(quads, *q.neighbors[1]);
                }
                if (q.neighbors[2])
                {
                    removeQuadFromGroup(quads, *q.neighbors[2]);
                }
            }
            else
            {
                if (q.neighbors[0])
                {
                    removeQuadFromGroup(quads, *q.neighbors[0]);
                }
                if (q.neighbors[3])
                {
                    removeQuadFromGroup(quads, *q.neighbors[3]);
                }
            }
        }
    }
}

void ChessBoardDetector::trimRow(std::vector<ChessBoardQuad*>& quads, int row, int dir)
{
    std::vector<ChessBoardQuad*> quads_(quads);
    // find the right quad(s)
    for (size_t i = 0; i < quads_.size(); ++i)
    {
        ChessBoardQuad& q = *quads_[i];

        if (q.ordered && q.row == row)
        {
            if (dir == 1)   // remove from bottom
            {
                if (q.neighbors[2])
                {
                    removeQuadFromGroup(quads, *q.neighbors[2]);
                }
                if (q.neighbors[3])
                {
                    removeQuadFromGroup(quads, *q.neighbors[3]);
                }
            }
            else    // remove from top
            {
                if (q.neighbors[0])
                {
                    removeQuadFromGroup(quads, *q.neighbors[0]);
                }
                if (q.neighbors[1])
                {
                    removeQuadFromGroup(quads, *q.neighbors[1]);
                }
            }

        }
    }
}
#endif

//
// remove quad from quad group
//
void ChessBoardDetector::removeQuadFromGroup(std::vector<ChessBoardQuad*>& quads, ChessBoardQuad& q0)
{
    const int count = (int)quads.size();

    int self_idx = -1;

    // remove any references to this quad as a neighbor
    for (int i = 0; i < count; ++i)
    {
        ChessBoardQuad* q = quads[i];
        if (q == &q0)
            self_idx = i;
        for (int j = 0; j < 4; j++)
        {
            if (q->neighbors[j] == &q0)
            {
                q->neighbors[j] = NULL;
                q->count--;
                for (int k = 0; k < 4; ++k)
                {
                    if (q0.neighbors[k] == q)
                    {
                        q0.neighbors[k] = 0;
                        q0.count--;
#ifndef _DEBUG
                        break;
#endif
                    }
                }
                break;
            }
        }
    }
    CV_Assert(self_idx >= 0); // item itself should be found

    // remove the quad
    if (self_idx != count-1)
        quads[self_idx] = quads[count-1];
    quads.resize(count - 1);
}

//
// put quad into correct order, where <corner> has value <common>
//
void ChessBoardDetector::orderQuad(ChessBoardQuad& quad, ChessBoardCorner& corner, int common)
{
    CV_DbgAssert(common >= 0 && common <= 3);

    // find the corner
    int tc = 0;;
    for (; tc < 4; ++tc)
        if (quad.corners[tc]->pt == corner.pt)
            break;

    // set corner order
    // shift
    while (tc != common)
    {
        // shift by one
        ChessBoardCorner *tempc = quad.corners[3];
        ChessBoardQuad *tempq = quad.neighbors[3];
        for (int i = 3; i > 0; --i)
        {
            quad.corners[i] = quad.corners[i-1];
            quad.neighbors[i] = quad.neighbors[i-1];
        }
        quad.corners[0] = tempc;
        quad.neighbors[0] = tempq;
        tc = (tc + 1) & 3;
    }
}

// if we found too many connect quads, remove those which probably do not belong.
int ChessBoardDetector::cleanFoundConnectedQuads(std::vector<ChessBoardQuad*>& quad_group)
{
    // number of quads this pattern should contain
    int count = ((pattern_size.width + 1)*(pattern_size.height + 1) + 1)/2;

    // if we have more quadrangles than we should,
    // try to eliminate duplicates or ones which don't belong to the pattern rectangle...
    int quad_count = (int)quad_group.size();
    if (quad_count <= count)
        return quad_count;
    CV_DbgAssert(quad_count > 0);

    // create an array of quadrangle centers
    cv::AutoBuffer<cv::Point2f> centers(quad_count);

    cv::Point2f center;
    for (int i = 0; i < quad_count; ++i)
    {
        ChessBoardQuad* q = quad_group[i];

        const cv::Point2f ci = (
                q->corners[0]->pt +
                q->corners[1]->pt +
                q->corners[2]->pt +
                q->corners[3]->pt
            ) * 0.25f;

        centers[i] = ci;
        center += ci;
    }
    center.x *= (1.0f / quad_count);

    // If we still have more quadrangles than we should,
    // we try to eliminate bad ones based on minimizing the bounding box.
    // We iteratively remove the point which reduces the size of
    // the bounding box of the blobs the most
    // (since we want the rectangle to be as small as possible)
    // remove the quadrange that causes the biggest reduction
    // in pattern size until we have the correct number
    for (; quad_count > count; quad_count--)
    {
        double min_box_area = DBL_MAX;
        int min_box_area_index = -1;

        // For each point, calculate box area without that point
        for (int skip = 0; skip < quad_count; ++skip)
        {
            // get bounding rectangle
            cv::Point2f temp = centers[skip]; // temporarily make index 'skip' the same as
            centers[skip] = center;            // pattern center (so it is not counted for convex hull)
            std::vector<Point2f> hull;
            Mat points(1, quad_count, CV_32FC2, &centers[0]);
            cv::convexHull(points, hull, true);
            centers[skip] = temp;
            double hull_area = contourArea(hull, true);

            // remember smallest box area
            if (hull_area < min_box_area)
            {
                min_box_area = hull_area;
                min_box_area_index = skip;
            }
        }

        ChessBoardQuad *q0 = quad_group[min_box_area_index];

        // remove any references to this quad as a neighbor
        for (int i = 0; i < quad_count; ++i)
        {
            ChessBoardQuad *q = quad_group[i];
            for (int j = 0; j < 4; ++j)
            {
                if (q->neighbors[j] == q0)
                {
                    q->neighbors[j] = 0;
                    q->count--;
                    for (int k = 0; k < 4; ++k)
                    {
                        if (q0->neighbors[k] == q)
                        {
                            q0->neighbors[k] = 0;
                            q0->count--;
                            break;
                        }
                    }
                    break;
                }
            }
        }

        // remove the quad
        quad_count--;
        quad_group[min_box_area_index] = quad_group[quad_count];
        centers[min_box_area_index] = centers[quad_count];
    }

    return quad_count;
}

void ChessBoardDetector::findConnectedQuads(std::vector<ChessBoardQuad*>& out_group, int group_idx)
{
    out_group.clear();

    std::stack<ChessBoardQuad*> stack;

    int i = 0;
    for (; i < all_quads_count; i++)
    {
        ChessBoardQuad* q = (ChessBoardQuad*)&all_quads[i];

        // Scan the array for a first unlabeled quad
        if (q->count <= 0 || q->group_idx >= 0) continue;

        // Recursively find a group of connected quads starting from the seed all_quads[i]
        stack.push(q);
        out_group.push_back(q);
        q->group_idx = group_idx;
        q->ordered = false;

        while (!stack.empty())
        {
            q = stack.top(); CV_Assert(q);
            stack.pop();
            for (int k = 0; k < 4; k++ )
            {
                ChessBoardQuad *neighbor = q->neighbors[k];
                if (neighbor && neighbor->count > 0 && neighbor->group_idx < 0 )
                {
                    stack.push(neighbor);
                    out_group.push_back(neighbor);
                    neighbor->group_idx = group_idx;
                    neighbor->ordered = false;
                }
            }
        }
        break;
    }
}

int ChessBoardDetector::checkQuadGroup(std::vector<ChessBoardQuad*>& quad_group, std::vector<ChessBoardCorner*>& out_corners)
{
    const int ROW1 = 1000000;
    const int ROW2 = 2000000;
    const int ROW_ = 3000000;

    int quad_count = (int)quad_group.size();

    std::vector<ChessBoardCorner*> corners(quad_count*4);
    int corner_count = 0;
    int result = 0;

    int width = 0, height = 0;
    int hist[5] = {0,0,0,0,0};

    // build dual graph, which vertices are internal quad corners
    // and two vertices are connected iff they lie on the same quad edge
    for (int i = 0; i < quad_count; ++i)
    {
        ChessBoardQuad* q = quad_group[i];

        for (int j = 0; j < 4; ++j)
        {
            if (q->neighbors[j])
            {
                int next_j = (j + 1) & 3;
                ChessBoardCorner *a = q->corners[j], *b = q->corners[next_j];
                // mark internal corners that belong to:
                //   - a quad with a single neighbor - with ROW1,
                //   - a quad with two neighbors     - with ROW2
                // make the rest of internal corners with ROW_
                int row_flag = q->count == 1 ? ROW1 : q->count == 2 ? ROW2 : ROW_;

                if (a->row == 0)
                {
                    corners[corner_count++] = a;
                    a->row = row_flag;
                }
                else if (a->row > row_flag)
                {
                    a->row = row_flag;
                }

                if (q->neighbors[next_j])
                {
                    if (a->count >= 4 || b->count >= 4)
                        goto finalize;
                    for (int k = 0; k < 4; ++k)
                    {
                        if (a->neighbors[k] == b)
                            goto finalize;
                        if (b->neighbors[k] == a)
                            goto finalize;
                    }
                    a->neighbors[a->count++] = b;
                    b->neighbors[b->count++] = a;
                }
            }
        }
    }

    if (corner_count != pattern_size.width*pattern_size.height)
        goto finalize;

{
    ChessBoardCorner* first = NULL, *first2 = NULL;
    for (int i = 0; i < corner_count; ++i)
    {
        int n = corners[i]->count;
        CV_DbgAssert(0 <= n && n <= 4);
        hist[n]++;
        if (!first && n == 2)
        {
            if (corners[i]->row == ROW1)
                first = corners[i];
            else if (!first2 && corners[i]->row == ROW2)
                first2 = corners[i];
        }
    }

    // start with a corner that belongs to a quad with a single neighbor.
    // if we do not have such, start with a corner of a quad with two neighbors.
    if( !first )
        first = first2;

    if( !first || hist[0] != 0 || hist[1] != 0 || hist[2] != 4 ||
        hist[3] != (pattern_size.width + pattern_size.height)*2 - 8 )
        goto finalize;

    ChessBoardCorner* cur = first;
    ChessBoardCorner* right = NULL;
    ChessBoardCorner* below = NULL;
    out_corners.push_back(cur);

    for (int k = 0; k < 4; ++k)
    {
        ChessBoardCorner* c = cur->neighbors[k];
        if (c)
        {
            if (!right)
                right = c;
            else if (!below)
                below = c;
        }
    }

    if( !right || (right->count != 2 && right->count != 3) ||
        !below || (below->count != 2 && below->count != 3) )
        goto finalize;

    cur->row = 0;
    //cvCircle( debug_img, cvPointFrom32f(cur->pt), 3, cvScalar(0,255,0), -1, 8, 0 );

    first = below; // remember the first corner in the next row

    // find and store the first row (or column)
    for (int j = 1; ; ++j)
    {
        right->row = 0;
        out_corners.push_back(right);
        //cvCircle( debug_img, cvPointFrom32f(right->pt), 3, cvScalar(0,255-j*10,0), -1, 8, 0 );
        if( right->count == 2 )
            break;
        if( right->count != 3 || (int)out_corners.size() >= std::max(pattern_size.width,pattern_size.height) )
            goto finalize;
        cur = right;
        for (int k = 0; k < 4; ++k)
        {
            ChessBoardCorner* c = cur->neighbors[k];
            if (c && c->row > 0)
            {
                int kk = 0;
                for (; kk < 4; ++kk)
                {
                    if (c->neighbors[kk] == below)
                        break;
                }
                if (kk < 4)
                    below = c;
                else
                    right = c;
            }
        }
    }

    width = (int)out_corners.size();
    if (width == pattern_size.width)
        height = pattern_size.height;
    else if (width == pattern_size.height)
        height = pattern_size.width;
    else
        goto finalize;

    // find and store all the other rows
    for (int i = 1; ; ++i)
    {
        if( !first )
            break;
        cur = first;
        first = 0;
        int j = 0;
        for (; ; ++j)
        {
            cur->row = i;
            out_corners.push_back(cur);
            //cvCircle( debug_img, cvPointFrom32f(cur->pt), 3, cvScalar(0,0,255-j*10), -1, 8, 0 );
            if (cur->count == 2 + (i < height-1) && j > 0)
                break;

            right = 0;

            // find a neighbor that has not been processed yet
            // and that has a neighbor from the previous row
            for (int k = 0; k < 4; ++k)
            {
                ChessBoardCorner* c = cur->neighbors[k];
                if (c && c->row > i)
                {
                    int kk = 0;
                    for (; kk < 4; ++kk)
                    {
                        if (c->neighbors[kk] && c->neighbors[kk]->row == i-1)
                            break;
                    }
                    if(kk < 4)
                    {
                        right = c;
                        if (j > 0)
                            break;
                    }
                    else if (j == 0)
                        first = c;
                }
            }
            if (!right)
                goto finalize;
            cur = right;
        }

        if (j != width - 1)
            goto finalize;
    }

    if ((int)out_corners.size() != corner_count)
        goto finalize;

    // check if we need to transpose the board
    if (width != pattern_size.width)
    {
        int temp = width;
        width = height;
        height = temp;

        std::vector<ChessBoardCorner*> tmp(out_corners);
        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                out_corners[i*width + j] = tmp[j*height + i];
            }
        }
    }

    // check if we need to revert the order in each row
    {
        cv::Point2f p0 = out_corners[0]->pt,
                    p1 = out_corners[pattern_size.width-1]->pt,
                    p2 = out_corners[pattern_size.width]->pt;
        if( (p1.x - p0.x)*(p2.y - p1.y) - (p1.y - p0.y)*(p2.x - p1.x) < 0 )
        {
            if (width % 2 == 0)
            {
                for (int i = 0; i < height; ++i)
                    for (int j = 0; j < width/2; ++j) {
                        ChessBoardCorner *temp = out_corners[i*width+j];
                        out_corners[i*width+j] = out_corners[i*width+width-j-1];
                        out_corners[i*width+width-j-1] = temp;
                    }
            }
            else
            {
                for (int j = 0; j < width; ++j) {
                    for (int i = 0; i < height/2; ++i) {
                        ChessBoardCorner *temp = out_corners[i*width+j];
                        out_corners[i*width+j] = out_corners[(height - i - 1)*width+j];
                        out_corners[(height - i - 1)*width+j] = temp;
                    }
                }
            }
        }
    }

    result = corner_count;
}

finalize:
    if (result <= 0)
    {
        corner_count = std::min(corner_count, pattern_size.area());
        out_corners.resize(corner_count);
        for (int i = 0; i < corner_count; i++)
            out_corners[i] = corners[i];

        result = -corner_count;

        if (result == -pattern_size.area())
            result = -result;
    }

    return result;
}

void ChessBoardDetector::findQuadNeighbors()
{
    const float thresh_scale = 1.f;
    // find quad neighbors
    for (int idx = 0; idx < all_quads_count; idx++)
    {
        ChessBoardQuad& cur_quad = (ChessBoardQuad&)all_quads[idx];

        // choose the points of the current quadrangle that are close to
        // some points of the other quadrangles
        // (it can happen for split corners (due to dilation) of the
        // checker board). Search only in other quadrangles!

        // for each corner of this quadrangle
        for (int i = 0; i < 4; i++)
        {
            if (cur_quad.neighbors[i])
                continue;

            float min_dist = FLT_MAX;
            int closest_corner_idx = -1;
            ChessBoardQuad *closest_quad = 0;

            cv::Point2f pt = cur_quad.corners[i]->pt;

            // find the closest corner in all other quadrangles
            for (int k = 0; k < all_quads_count; k++)
            {
                if (k == idx)
                    continue;

                ChessBoardQuad& q_k = all_quads[k];

                for (int j = 0; j < 4; j++)
                {
                    if (q_k.neighbors[j])
                        continue;

                    float dist = normL2Sqr<float>(pt - q_k.corners[j]->pt);
                    if (dist < min_dist &&
                        dist <= cur_quad.edge_len*thresh_scale &&
                        dist <= q_k.edge_len*thresh_scale )
                    {
                        // check edge lengths, make sure they're compatible
                        // edges that are different by more than 1:4 are rejected
                        float ediff = cur_quad.edge_len - q_k.edge_len;
                        if (ediff > 32*cur_quad.edge_len ||
                            ediff > 32*q_k.edge_len)
                        {
                            continue;
                        }
                        closest_corner_idx = j;
                        closest_quad = &q_k;
                        min_dist = dist;
                    }
                }
            }

            // we found a matching corner point?
            if (closest_corner_idx >= 0 && min_dist < FLT_MAX)
            {
                CV_Assert(closest_quad);

                if (cur_quad.count >= 4 || closest_quad->count >= 4)
                    continue;

                // If another point from our current quad is closer to the found corner
                // than the current one, then we don't count this one after all.
                // This is necessary to support small squares where otherwise the wrong
                // corner will get matched to closest_quad;
                ChessBoardCorner& closest_corner = *closest_quad->corners[closest_corner_idx];

                int j = 0;
                for (; j < 4; j++)
                {
                    if (cur_quad.neighbors[j] == closest_quad)
                        break;

                    if (normL2Sqr<float>(closest_corner.pt - cur_quad.corners[j]->pt) < min_dist)
                        break;
                }
                if (j < 4)
                    continue;

                // Check that each corner is a neighbor of different quads
                for(j = 0; j < closest_quad->count; j++ )
                {
                    if (closest_quad->neighbors[j] == &cur_quad)
                        break;
                }
                if (j < closest_quad->count)
                    continue;

                // check whether the closest corner to closest_corner
                // is different from cur_quad->corners[i]->pt
                for (j = 0; j < all_quads_count; j++ )
                {
                    ChessBoardQuad* q = &const_cast<ChessBoardQuad&>(all_quads[j]);
                    if (j == idx || q == closest_quad)
                        continue;

                    int k = 0;
                    for (; k < 4; k++ )
                    {
                        if (!q->neighbors[k])
                        {
                            if (normL2Sqr<float>(closest_corner.pt - q->corners[k]->pt) < min_dist)
                                break;
                        }
                    }
                    if (k < 4)
                        break;
                }
                if (j < all_quads_count)
                    continue;

                closest_corner.pt = (pt + closest_corner.pt) * 0.5f;

                // We've found one more corner - remember it
                cur_quad.count++;
                cur_quad.neighbors[i] = closest_quad;
                cur_quad.corners[i] = &closest_corner;

                closest_quad->count++;
                closest_quad->neighbors[closest_corner_idx] = &cur_quad;
            }
        }
    }
}

// returns corners in clockwise order
// corners don't necessarily start at same position on quad (e.g. top left corner)
void ChessBoardDetector::generateQuads(const cv::Mat& image_, int flags)
{
    binarized_image = image_;  // save for debug purposes

    int quad_count = 0;

    all_quads.deallocate();
    all_corners.deallocate();

    // empiric bound for minimal allowed perimeter for squares
    int min_size = 25; //cvRound( image->cols * image->rows * .03 * 0.01 * 0.92 );

    bool filterQuads = (flags & CALIB_CB_FILTER_QUADS) != 0;

    std::vector<std::vector<Point> > contours;
    std::vector<Vec4i> hierarchy;

    cv::findContours(image_, contours, hierarchy, RETR_CCOMP, CHAIN_APPROX_SIMPLE);

    if (contours.empty())
    {
        return;
    }

    std::vector<int> contour_child_counter(contours.size(), 0);
    int boardIdx = -1;

    std::vector<QuadCountour> contour_quads;

    for (int idx = (int)(contours.size() - 1); idx >= 0; --idx)
    {
        int parentIdx = hierarchy[idx][3];
        if (hierarchy[idx][2] != -1 || parentIdx == -1)  // holes only (no child contours and with parent)
            continue;
        const std::vector<Point>& contour = contours[idx];

        Rect contour_rect = boundingRect(contour);
        if (contour_rect.area() < min_size)
            continue;

        std::vector<Point> approx_contour;

        const int min_approx_level = 1, max_approx_level = MAX_CONTOUR_APPROX;
        for (int approx_level = min_approx_level; approx_level <= max_approx_level; approx_level++ )
        {
            approxPolyDP(contour, approx_contour, (float)approx_level, true);
            if (approx_contour.size() == 4)
                break;

            // we call this again on its own output, because sometimes
            // approxPoly() does not simplify as much as it should.
            std::vector<Point> approx_contour_tmp;
            std::swap(approx_contour, approx_contour_tmp);
            approxPolyDP(approx_contour_tmp, approx_contour, (float)approx_level, true);
            if (approx_contour.size() == 4)
                break;
        }

        // reject non-quadrangles
        if (approx_contour.size() != 4)
            continue;
        if (!cv::isContourConvex(approx_contour))
            continue;

        cv::Point pt[4];
        for (int i = 0; i < 4; ++i) {
            pt[i] = approx_contour[i];
        }

        if (filterQuads)
        {
            double p = cv::arcLength(approx_contour, true);
            double area = cv::contourArea(approx_contour, false);

            double d1 = sqrt(normL2Sqr<double>(pt[0] - pt[2]));
            double d2 = sqrt(normL2Sqr<double>(pt[1] - pt[3]));

            // philipg.  Only accept those quadrangles which are more square
            // than rectangular and which are big enough
            double d3 = sqrt(normL2Sqr<double>(pt[0] - pt[1]));
            double d4 = sqrt(normL2Sqr<double>(pt[1] - pt[2]));
            if (!(d3*4 > d4 && d4*4 > d3 && d3*d4 < area*1.5 && area > min_size &&
                d1 >= 0.15 * p && d2 >= 0.15 * p))
                continue;
        }

        contour_child_counter[parentIdx]++;
        if (boardIdx != parentIdx && (boardIdx < 0 || contour_child_counter[boardIdx] < contour_child_counter[parentIdx]))
            boardIdx = parentIdx;

        contour_quads.push_back(QuadCountour(pt, parentIdx));
    }

    size_t total = contour_quads.size();
    size_t max_quad_buf_size = std::max((size_t)2, total * 3);
    all_quads.allocate(max_quad_buf_size);
    all_corners.allocate(max_quad_buf_size * 4);

    // Create array of quads structures
    for (size_t idx = 0; idx < total; ++idx)
    {
        QuadCountour& qc = contour_quads[idx];
        if (filterQuads && qc.parent_contour != boardIdx)
            continue;

        int quad_idx = quad_count++;
        ChessBoardQuad& q = all_quads[quad_idx];

        // reset group ID
        q = ChessBoardQuad();
        for (int i = 0; i < 4; ++i)
        {
            Point2f pt(qc.pt[i]);
            ChessBoardCorner& corner = all_corners[quad_idx * 4 + i];

            corner = ChessBoardCorner(pt);
            q.corners[i] = &corner;
        }
        q.edge_len = FLT_MAX;
        for (int i = 0; i < 4; ++i)
        {
            float d = normL2Sqr<float>(q.corners[i]->pt - q.corners[(i+1)&3]->pt);
            q.edge_len = std::min(q.edge_len, d);
        }
    }

    all_quads_count = quad_count;
}

bool ChessBoardDetector::processQuads(std::vector<cv::Point2f>& out_corners, int &prev_sqr_size)
{
    out_corners.resize(0);
    if (all_quads_count <= 0)
        return false;

    size_t max_quad_buf_size = all_quads.size();

    // Find quad's neighbors
    findQuadNeighbors();

    // allocate extra for adding in orderFoundQuads
    std::vector<ChessBoardQuad*> quad_group;
    std::vector<ChessBoardCorner*> corner_group; corner_group.reserve(max_quad_buf_size * 4);

    for (int group_idx = 0; ; group_idx++)
    {
        findConnectedQuads(quad_group, group_idx);
        if (quad_group.empty())
            break;

        int count = (int)quad_group.size();

        // order the quad corners globally
        // maybe delete or add some
        count = orderFoundConnectedQuads(quad_group);

        if (count == 0)
            continue;       // haven't found inner quads

        // If count is more than it should be, this will remove those quads
        // which cause maximum deviation from a nice square pattern.
        count = cleanFoundConnectedQuads(quad_group);

        count = checkQuadGroup(quad_group, corner_group);

        int n = count > 0 ? pattern_size.width * pattern_size.height : -count;
        n = std::min(n, pattern_size.width * pattern_size.height);
        float sum_dist = 0;
        int total = 0;

        for(int i = 0; i < n; i++ )
        {
            int ni = 0;
            float sum = corner_group[i]->sumDist(ni);
            sum_dist += sum;
            total += ni;
        }
        prev_sqr_size = cvRound(sum_dist/std::max(total, 1));

        if (count > 0 || (-count > (int)out_corners.size()))
        {
            // copy corners to output array
            out_corners.reserve(n);
            for (int i = 0; i < n; ++i)
                out_corners.push_back(corner_group[i]->pt);

            if (count == pattern_size.width*pattern_size.height
                    && checkBoardMonotony(out_corners))
            {
                return true;
            }
        }
    }

    return false;
}
