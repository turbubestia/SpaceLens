/*
 * PaintLayer.cpp
 *
 *  Created on: Jan 7, 2021
 *      Author: claud
 */

#include "LayerGraphicWidget.h"
#include "PaintLayer.h"

PaintLayer::PaintLayer() {
    drag = false;
    button = 0;
    penSize = 10;
}

PaintLayer::~PaintLayer() {

}

void PaintLayer::clear() {
    _image = QImage(100,100,QImage::Format_ARGB32);
    _image.fill(QColor(0,0,0,0));
    _width = 100;
    _height = 100;
}

const QImage & PaintLayer::image() {
    return _image;
}

bool PaintLayer::mousePressEvent(QMouseEvent *event) {
    lastDragPoint = parent()->mapWidgetToCanvas(event->pos());
    drag = true;
    button = event->button();

    QPen pen;
    switch (button) {
        case Qt::LeftButton: pen.setColor(QColor(255,0,0,100)); break;
        case Qt::RightButton: pen.setColor(QColor(0,0,255,100)); break;
        case Qt::MiddleButton: pen.setColor(QColor(0,0,0,0)); break;
        default: pen.setColor(QColor(127,127,127,0)); break;
    }
    pen.setWidth(penSize);
    pen.setCapStyle(Qt::RoundCap);
    pen.setJoinStyle(Qt::RoundJoin);

    QPainter painter(&_image);
    painter.setCompositionMode(QPainter::CompositionMode_Source);
    painter.setPen(pen);
    painter.drawPoint(lastDragPoint);

    int x = lastDragPoint.x() - penSize/2;
    int y = lastDragPoint.y() - penSize/2;
    int width = penSize;
    int height = penSize;

    roi = QRect(x,y,width,height);

    return true;
}

bool PaintLayer::mouseReleaseEvent(QMouseEvent *event) {
    Q_UNUSED(event);

    if (drag) {
        drag = false;
        qDebug() << "PaintLayer::mouseReleaseEvent()" << roi;
        emit changed();
        emit roiChanged(roi);
    }

    return true;
}

bool PaintLayer::mouseMoveEvent(QMouseEvent *event) {
    if (drag) {
        QPoint currentDragPoint = parent()->mapWidgetToCanvas(event->pos());

        QPen pen;
        switch (button) {
            case Qt::LeftButton: pen.setColor(QColor(255,0,0,100)); break;
            case Qt::RightButton: pen.setColor(QColor(0,0,255,100)); break;
            case Qt::MiddleButton: pen.setColor(QColor(0,0,0,0)); break;
            default: pen.setColor(QColor(127,127,127,0)); break;
        }
        pen.setWidth(penSize);
        pen.setCapStyle(Qt::RoundCap);
        pen.setJoinStyle(Qt::RoundJoin);

        QPainter painter(&_image);
        painter.setCompositionMode(QPainter::CompositionMode_Source);
        painter.setPen(pen);
        painter.drawLine(lastDragPoint, currentDragPoint);

        lastDragPoint = currentDragPoint;

        int x = lastDragPoint.x() - penSize/2;
        int y = lastDragPoint.y() - penSize/2;
        int width = penSize;
        int height = penSize;

        if (x < roi.left()) {
            roi.setLeft(x);
        }
        if (y < roi.top()) {
            roi.setTop(y);
        }
        if (x+width > roi.right()) {
            roi.setRight(x+width);
        }
        if (y+height > roi.bottom()) {
            roi.setBottom(y+height);
        }
    }

    parent()->update();

    return true;
}

void PaintLayer::render(QPainter *painter) {

    if (!isVisible()) {
        return;
    }

    QPoint pos = QCursor::pos();
    pos = parent()->mapFromGlobal(pos);
    pos = parent()->mapWidgetToCanvas(pos);

    painter->drawImage(0,0,_image);
    painter->setPen(QColor(0,0,0,255));
    painter->drawEllipse(QPointF(pos),penSize/2.0,penSize/2.0);
}

void PaintLayer::resize(QSize size) {
    if (!_image.isNull() && _image.size() != size) {
        _image = _image.scaled(size);
        _width = size.width();
        _height = size.height();
    }
}

void PaintLayer::resizeToCanvas() {
    if (!_image.isNull()) {
        _image = _image.scaled(parent()->canvasSize());
    }
}

void PaintLayer::setImage(const QImage &image) {
    QSize size = image.size();
    _image = image;
    if (!image.isNull()) {
        _width = size.width();
        _height = size.height();
    } else {
        _width = _height = 100;
    }
}

bool PaintLayer::wheelEvent(QWheelEvent *event) {
    if (event->angleDelta().y() < 0) {
        penSize -= 2;
        if (penSize < 2) {
            penSize = 2;
        }
    } else {
        penSize += 2;
        if (penSize > 100) {
            penSize = 100;
        }
    }

    parent()->update();

    return true;
}
