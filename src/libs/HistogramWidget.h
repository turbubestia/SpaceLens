/*
 * HystogramWidget.h
 *
 *  Created on: 15 jun. 2020
 *      Author: claud
 */

#ifndef SRC_LIBS_HISTOGRAMWIDGET_H_
#define SRC_LIBS_HISTOGRAMWIDGET_H_

#include <QtWidgets/QtWidgets>

class HistogramWidget: public QWidget {
	Q_OBJECT

	public:
		HistogramWidget(QWidget *parent = nullptr);
		virtual ~HistogramWidget();

		void fromImage(const QImage &image);
		void clear();

	protected:
		void paintEvent(QPaintEvent *event) override;

	private:
		int histogram[256];
};

#endif /* SRC_LIBS_HISTOGRAMWIDGET_H_ */
