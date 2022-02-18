/*
 * PaintLayer.h
 *
 *  Created on: Jan 7, 2021
 *      Author: claud
 */

#ifndef SRC_LIBS_PAINTLAYER_H_
#define SRC_LIBS_PAINTLAYER_H_

#include <QtWidgets/QtWidgets>

#include "AbstractLayer.h"

class PaintLayer: public AbstractLayer {
    Q_OBJECT

    public:
        PaintLayer();
        virtual ~PaintLayer();

        void clear();
        bool mousePressEvent(QMouseEvent *event) override;
        bool mouseReleaseEvent(QMouseEvent *event) override;
        bool mouseMoveEvent(QMouseEvent *event) override;
        void render(QPainter *painter) override;
        void resize(QSize size) override;
        void resizeToCanvas() override;
        void setImage(const QImage &image);
        const QImage &image();
        bool wheelEvent(QWheelEvent *event) override;

    signals:
        void changed();
        void roiChanged(QRect roi);

    private:
        QImage _image;

        bool drag;
        QPoint lastDragPoint;
        int button;
        int penSize;
        QRect roi;
};

#endif /* SRC_LIBS_PAINTLAYER_H_ */
