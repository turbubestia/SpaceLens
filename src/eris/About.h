/*
 * About.h
 *
 *  Created on: 21 jun. 2020
 *      Author: claud
 */

#ifndef SRC_ERIS_ABOUT_H_
#define SRC_ERIS_ABOUT_H_

#include <QtWidgets/QtWidgets>

class About: public QDialog {
	Q_OBJECT

	public:
		About(QWidget *parent = nullptr, Qt::WindowFlags f = Qt::WindowFlags());
		virtual ~About();

	private:
		QTabWidget *sections;
		QLabel *about;
};

#endif /* SRC_ERIS_ABOUT_H_ */
