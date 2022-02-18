/*
 * RoiEditor.cpp
 *
 *  Created on: 28 jun. 2020
 *      Author: claud
 */

#include "LayerGraphicWidget.h"
#include "RoiEditorLayer.h"

RoiEditorLayer::RoiEditorLayer() {
	x = 0;
	y = 0;
	width = 0;
	height = 0;
	drag = false;
	corner = 0;
	r = 7;
}

RoiEditorLayer::~RoiEditorLayer() {

}

QRect RoiEditorLayer::getRoi() {
	return QRect(x,y,width,height);
}

bool RoiEditorLayer::mousePressEvent(QMouseEvent *event) {
	lastDragPoint = parent()->mapWidgetToCanvas(event->pos());

	int cx = lastDragPoint.x();
	int cy = lastDragPoint.y();

	int vx[4] = {x,x+width,x+width,x};
	int vy[4] = {y,y,y+height,y+height};

	double rr = r/parent()->getScale();

	corner = -1;
	for(int i = 0; i < 4; i++) {
		double dx = cx - vx[i];
		double dy = cy - vy[i];
		double d = sqrt(dx*dx+dy*dy);
		if (d < rr) {
			corner = i;
			break;
		}
	}

	drag = corner != -1;

	return drag;
}

bool RoiEditorLayer::mouseReleaseEvent(QMouseEvent *event) {
	Q_UNUSED(event);

	if (drag) {
		drag = false;
		return true;
	} else {
		return false;
	}
}

bool RoiEditorLayer::mouseMoveEvent(QMouseEvent *event) {
	if (drag) {
		QPoint pos = parent()->mapWidgetToCanvas(event->pos());
		if (corner == 0) {
			x = pos.x();
			y = pos.y();
		} else if(corner == 1) {
			width = abs(x - pos.x());
		} else if(corner == 2) {
			width = abs(x - pos.x());
			height = abs(y - pos.y());
		} else {
			height = abs(y - pos.y());
		}

		emit roiChanged(QRect(x,y,width,height));

		return true;
	} else {
		return false;
	}
}

void RoiEditorLayer::render(QPainter *painter) {

	if (!isVisible() || width <= 0 || height <= 0) {
		return;
	}

	painter->save();

	QPen pen;
	pen.setColor(QColor(250,20,20,255));
	pen.setWidthF(1.5/parent()->getScale());
	painter->setPen(pen);

	painter->setBrush(QColor(0,0,0,0));

	double rr = r/parent()->getScale();

	painter->drawEllipse(QPointF(x,y),rr,rr);
	painter->drawEllipse(QPointF(x+width,y),rr,rr);
	painter->drawEllipse(QPointF(x+width,y+height),rr,rr);
	painter->drawEllipse(QPointF(x,y+height),rr,rr);

	painter->drawLine(QPointF(x+rr,y),QPointF(x+width-rr,y));
	painter->drawLine(QPointF(x+width,y+rr),QPointF(x+width,y+height-rr));
	painter->drawLine(QPointF(x+width-rr,y+height),QPointF(x+rr,y+height));
	painter->drawLine(QPointF(x,y+height-rr),QPointF(x,y+rr));

	painter->restore();
}

void RoiEditorLayer::setRoi(QRect roi) {
	x = roi.x();
	y = roi.y();
	width = roi.width();
	height = roi.height();

	parent()->update();
}
