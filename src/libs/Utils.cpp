/*
 * Utils.cpp
 *
 *  Created on: 18 abr. 2020
 *      Author: claud
 */

#include "Utils.h"

QImage cvMatToQImage(const cv::Mat &inMat) {
	switch(inMat.type()) {
		case CV_8UC4: {
			QImage outImage(inMat.data, inMat.cols, inMat.rows, static_cast<int>(inMat.step), QImage::Format_ARGB32);
			return outImage.rgbSwapped();
		}

		case CV_8UC3: {
			QImage outImage(inMat.data, inMat.cols, inMat.rows, static_cast<int>(inMat.step), QImage::Format_RGB888);
			return outImage.rgbSwapped();
		}

		default:
			return QImage();
	}
}

QString cvMatToString(const cv::Mat &M) {
	QString str;

	int type = M.type();
	str.append(QString("%1;").arg(type, 0, 10));

	for (int row = 0; row < M.rows; row++) {
		for (int col = 0; col < M.cols; col++) {
			QString temp;

			switch(CV_MAT_DEPTH(type)) {
				case CV_8U:
					temp = QString::number(M.at<uint8_t>(row,col), 10);
					break;
				case CV_8S:
					temp = QString::number(M.at<int8_t>(row,col), 10);
					break;
				case CV_16U:
					temp = QString::number(M.at<uint16_t>(row,col), 10);
					break;
				case CV_16S:
					temp = QString::number(M.at<int16_t>(row,col), 10);
					break;
				case CV_32F:
					temp = QString::number((double) M.at<float>(row,col), 'e', 7);
					break;
				case CV_64F:
					temp = QString::number(M.at<double>(row,col), 'e', 16);
					break;
			}

			str.append(temp);
			if (col < M.cols-1) {
				str.append(",");
			}
		}

		if (row < M.rows-1) {
			str.append(";");
		}
	}

	return str;
}

QString cvPoint2fListToString(const std::vector<cv::Point2f> &pointList) {
	QString str;

	for (int i = 0; i < (int) pointList.size(); i++) {
		cv::Point2f p = pointList.at(i);
		QString temp = QString("%1,%2").arg((double) p.x, 0, 'e', 7).arg((double) p.y, 0, 'e', 7);
		str.append(temp);
		if (i < (int) pointList.size() - 1) {
			str.append(";");
		}
	}

	return str;
}

int findBimodalThreshold(const QImage &image) {
	cv::Mat src = qImageToCvMat(image);

	cv::Mat gray;
	cvtColor(src, gray, cv::COLOR_BGRA2GRAY);

	return findBimodalThreshold(gray);
}

void smooth(const std::vector<double> &src, std::vector<double> &dst, int r) {
	int len = src.size();

	std::vector<double> temp(src.size(),0);
	for (int i = 0; i < len; i++) {
		int s = i-r < 0 ? 0 : i-r;
		int e = i+r > len-1 ? len-1 : i+r;

		double m = 0;
		for (int j = s; j <= e; j++ ) {
			m += src[j];
		}

		temp[i] = m / (e-s+1);
	}

	dst = temp;
}

int findBimodalThreshold(const cv::Mat &gray) {

	// compute the sparse histogram
	std::vector<int> histo(256, 0);
	for (int i = 0; i < gray.rows; i++) {
		const uchar *gray_ptr = gray.ptr<uchar>(i);
		for (int j = 0; j < gray.cols; j++) {
			int bin = gray_ptr[j];
			histo[bin]++;
		}
	}

	// compute the accumulated histogram. This provides a smoother distribution
	// that can be refined to obtain a dense histogram
	std::vector<int> acc_histo(256, 0);
	acc_histo[1] = histo[0];
	for (int i = 1; i < 256; i++) {
		acc_histo[i] = acc_histo[i-1] + histo[i];
	}

	// interpolate a fine (dense) histogram and smooth it
	std::vector<double> acc_sub_histo(2550,0);
	for (int i = 0, k = 0; i < 255; i++) {
		double h0 = acc_histo[i];
		double h1 = acc_histo[i+1];
		// use linear interpolation
		for (int j = 0; j < 10; j++, k++) {
			double t = j/10.0;
			acc_sub_histo[k] = h0 * (1 - t) + h1 * t;
		}
	}
	smooth(acc_sub_histo, acc_sub_histo, 19);

	// from the fine histogram reconstruct the histogram and smooth it
	std::vector<double> sub_histo(2550,0);
	for (int i = 1; i < 2549; i++) {
		sub_histo[i] = acc_sub_histo[i+1] - acc_sub_histo[i-1];
	}
	smooth(sub_histo, sub_histo, 19);

	// with the gradient of the histogram we can find the peaks
	std::vector<double> grad_sub_histo(2550,0);
	for (int i = 1; i < 2549; i++) {
		grad_sub_histo[i] = sub_histo[i+1] - sub_histo[i-1];
	}
	smooth(grad_sub_histo, grad_sub_histo, 19);

	// find the peaks and valley

	int peak_index = -1, valley_index = 64;
	for (int i = 5; i < 2549; i++) {
		if (grad_sub_histo[i] > 0 && grad_sub_histo[i+1] <= 0) {
			peak_index = i;
		}
		if (grad_sub_histo[i] <= 0 && grad_sub_histo[i+1] > 0 && peak_index > 0) {
			valley_index = i;
			break;
		}
	}

	return round(valley_index/10.0) + 25;
}

cv::Mat qImageToCvMat(const QImage &inImage) {
	switch (inImage.format()) {
		case QImage::Format_ARGB32:
		case QImage::Format_ARGB32_Premultiplied: {
		    QImage temp = inImage.rgbSwapped();
			cv::Mat outMat(temp.height(), temp.width(), CV_8UC4, const_cast<uchar*>(temp.bits()), static_cast<size_t>(temp.bytesPerLine()));
			return outMat.clone();
		}

		case QImage::Format_RGB32: {
			cv::Mat mat(inImage.height(), inImage.width(), CV_8UC4, const_cast<uchar*>(inImage.bits()),
					static_cast<size_t>(inImage.bytesPerLine()));
			cv::Mat outMat;
			cv::cvtColor(mat, outMat, cv::COLOR_BGRA2BGR);
			return outMat;
		}

		case QImage::Format_RGB888: {
			QImage swapped = inImage.rgbSwapped();
			cv::Mat outMat(swapped.height(), swapped.width(), CV_8UC3,const_cast<uchar*>(swapped.bits()),
					static_cast<size_t>(swapped.bytesPerLine()));
			return outMat.clone();
		}

		default:
			qDebug() << "Unsupported format";
			return cv::Mat();
	}
}

cv::Mat stringToCvMat(QString str) {
	QStringList rowsString = str.split(";");
	if (rowsString.count() == 0) {
		return cv::Mat();
	}

	QStringList colsString = rowsString.at(1).split(',');

	int type = rowsString.at(0).toInt();
	int rows = rowsString.count() - 1;
	int cols = colsString.count();

	cv::Mat M(rows, cols, type);

	for (int row = 0; row < rows; row++) {
		colsString = rowsString.at(1+row).split(',');
		for (int col = 0; col < cols; col++) {
			switch (type) {
				case CV_8U:
					M.at<uint8_t>(row, col) = colsString.at(col).toUInt();
					break;
				case CV_8S:
					M.at<uint8_t>(row, col) = colsString.at(col).toInt();
					break;
				case CV_16U:
					M.at<uint16_t>(row, col) = colsString.at(col).toUInt();
					break;
				case CV_16S:
					M.at<int16_t>(row, col) = colsString.at(col).toInt();
					break;
				case CV_32F:
					M.at<float>(row, col) = colsString.at(col).toFloat();
					break;
				case CV_64F:
					M.at<double>(row, col) = colsString.at(col).toDouble();
					break;
			}
		}
	}

	return M;
}

std::vector<cv::Point2f> stringToCvPoint2fList(QString str) {
	std::vector<cv::Point2f> pointList;

	QStringList points = str.split(";");
	for(int i = 0; i < points.count(); i++) {
		QStringList coord = points.at(i).split(',');
		if (coord.count() == 2) {
			double x = coord.at(0).toDouble();
			double y = coord.at(1).toDouble();
			pointList.push_back(cv::Point2f(x,y));
		} else {
			return std::vector<cv::Point2f>();
		}
	}

	return pointList;
}

QString commonPrefix(QStringList strs) {
    int length = strs.at(0).count();
    for (int i = 1; i < strs.count(); i++) {
        if (strs.at(i).count() < length) {
            length = strs.at(i).count();
        }
    }

    int substr = 0;
    for (int i = 0; i < length; i++) {
        QChar c = strs.at(i).at(0).toLower();
        for (int k = 1; k < strs.count(); k++) {
            if (strs.at(i).at(k).toLower() != c) {
                substr = i;
                i = length;
                break;
            }
        }
    }

    return strs.at(0).mid(0, substr-1);
}


