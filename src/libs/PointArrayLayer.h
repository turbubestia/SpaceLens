/*
 * PointArrayLayer.h
 *
 *  Created on: 18 jun. 2020
 *      Author: claud
 */

#ifndef SRC_LIBS_POINTARRAYLAYER_H_
#define SRC_LIBS_POINTARRAYLAYER_H_

#include <QtWidgets/QtWidgets>

#include "AbstractLayer.h"

class PointArrayLayer: public AbstractLayer {
	public:
		PointArrayLayer();
		virtual ~PointArrayLayer();

		void render(QPainter *painter) override;
		void setPoints(const QVector<QPointF> &points);
		void clear();

	private:
		QVector<QPointF> points;
};

#endif /* SRC_LIBS_POINTARRAYLAYER_H_ */
