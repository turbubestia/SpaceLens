/*
 * RoiEditor.h
 *
 *  Created on: 28 jun. 2020
 *      Author: claud
 */

#ifndef SRC_LIBS_ROIEDITORLAYER_H_
#define SRC_LIBS_ROIEDITORLAYER_H_

#include <QtWidgets/QtWidgets>

#include "AbstractLayer.h"

class RoiEditorLayer : public AbstractLayer {
	Q_OBJECT

	public:
		RoiEditorLayer();
		virtual ~RoiEditorLayer();

		QRect getRoi();
		bool mousePressEvent(QMouseEvent *event) override;
		bool mouseReleaseEvent(QMouseEvent *event) override;
		bool mouseMoveEvent(QMouseEvent *event) override;
		void render(QPainter *painter) override;
		void setRoi(QRect roi);

	signals:
		void roiChanged(QRect roi);

	private:
		int x;
		int y;
		int width;
		int height;

		bool drag;
		QPoint lastDragPoint;

		int corner;

		double r;
};

#endif /* SRC_LIBS_ROIEDITORLAYER_H_ */
