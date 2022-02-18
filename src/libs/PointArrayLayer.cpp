/*
 * PointArrayLayer.cpp
 *
 *  Created on: 18 jun. 2020
 *      Author: claud
 */

#include "PointArrayLayer.h"

PointArrayLayer::PointArrayLayer() {

}

PointArrayLayer::~PointArrayLayer() {

}

void PointArrayLayer::render(QPainter *painter) {
	if (points.isEmpty()) {
		return;
	}

	painter->save();

	QPen pen;
	pen.setColor(QColor(250,20,20,255));
	pen.setWidthF(2);
	painter->setPen(pen);

	painter->setBrush(QColor(0,0,0,0));

	for (int i = 0; i < points.count(); i++) {
		double x = points.at(i).x();
		double y = points.at(i).y();
		painter->drawEllipse(QPointF(x,y),10.0,10.0);
	}

	for (int i = 0; i < points.count()-1; i++) {
		double x1 = points.at(i).x();
		double y1 = points.at(i).y();
		double x2 = points.at(i+1).x();
		double y2 = points.at(i+1).y();
		painter->drawLine(QPointF(x1,y1),QPointF(x2,y2));
	}

	painter->restore();
}

void PointArrayLayer::setPoints(const QVector<QPointF> &_points) {
	points = _points;
}

void PointArrayLayer::clear() {
	points.clear();
}
