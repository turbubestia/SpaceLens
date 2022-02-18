/*
 * HystogramWidget.cpp
 *
 *  Created on: 15 jun. 2020
 *      Author: claud
 */

#include <HistogramWidget.h>

HistogramWidget::HistogramWidget(QWidget *parent)
	: QWidget(parent) {

	clear();

	setMinimumHeight(100);
}

HistogramWidget::~HistogramWidget() {

}

void HistogramWidget::fromImage(const QImage &image) {

	clear();

	if (!image.isNull()) {
		for (int i = 0; i < image.height(); i++) {
			for (int j = 0; j < image.width(); j++) {
				int bin = qGray(image.pixel(j,i));
				histogram[bin]++;
			}
		}
	}

	update();
}

void HistogramWidget::clear() {
	memset(histogram, 0, sizeof(int)*256);
}

void HistogramWidget::paintEvent(QPaintEvent *event) {
	QPainter painter(this);
	painter.setRenderHint(QPainter::Antialiasing);

	double w = width();
	double h = height();

	double max = 0;
	for(int i = 0; i < 256; i++) {
		if (histogram[i] > max) {
			max = histogram[i];
		}
	}

    QTransform view;
	view.translate(5,h-5);
	view.scale((w-10)/256.0, (10-h)/max);
    painter.setTransform(view);

	QPainterPath contour;
	contour.moveTo(5,h-5);
	for (int i = 0; i < 256; i++) {
		contour.lineTo(i*(w-10)/255.0 + 5, histogram[i]*(10-h)/max + h-5);
	}
	contour.lineTo(w-5, h-5);

	painter.fillRect(0,0,w,h, QColor(50,50,50));

	painter.setPen(Qt::white);
	painter.setBrush(Qt::white);
	painter.drawPath(contour);

	event->accept();
}

