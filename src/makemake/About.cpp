/*
 * About.cpp
 *
 *  Created on: 21 jun. 2020
 *      Author: claud
 */

#include "About.h"

About::About(QWidget *parent, Qt::WindowFlags f)
	: QDialog(parent, f) {

	about = new QLabel;
	about->setText(QString(
	"<h3>MakeMake - Focus Stacker</h3> <br>" \
	"<b>MakeMake</b> is part of the <b>TinyLenses</b> software kit<br><br>" \

	"Developed by Claudio Alvarez - <a href=\"https://www.lantern.cl\">www.lantern.cl</a><br><br>" \

	"<pre>This program is free software: you can redistribute it and/or modify\n" \
    "it under the terms of the GNU General Public License as published by\n" \
    "the Free Software Foundation, either version 3 of the License, or\n" \
    "(at your option) any later version.\n\n" \
    "This program is distributed in the hope that it will be useful,\n" \
    "but WITHOUT ANY WARRANTY; without even the implied warranty of\n" \
    "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n" \
    "GNU General Public License for more details.\n\n" \
    "You should have received a copy of the GNU General Public License\n" \
    "along with this program.  If not, see <a href=\"https://www.gnu.org/licenses/\">www.gnu.org/licenses/</a>.</pre>"));

	sections = new QTabWidget;
	sections->addTab(about, "About Makemake");

	QLabel *qt5 = new QLabel;
	qt5->setText(QString(
	"This software uses Qt5 using the GNU General Public License 3. <br><br>"
	"<pre>GNU General Public License Usage \n"
	"Alternatively, this file may be used under the terms of the GNU\n"
	"General Public License version 2.0 or (at your option) the GNU General\n"
	"Public license version 3 or any later version approved by the KDE Free\n"
	"Qt Foundation. The licenses are as published by the Free Software\n"
	"Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3\n"
	"included in the packaging of this file. Please review the following\n"
	"information to ensure the GNU General Public License requirements will\n"
	"be met: <a href=\"https://www.gnu.org/licenses/gpl-2.0.html\">www.gnu.org/licenses/gpl-2.0.html</a> and\n"
	"<a href=\"https://www.gnu.org/licenses/gpl-3.0.html\">www.gnu.org/licenses/gpl-3.0.html.</a>"
	));

	sections->addTab(qt5, "Qt 5");

	QLabel *openCV = new QLabel;
	openCV->setText(QString(
	"This software uses OpenCV 4. Get a copy of the license at <a href=\"https://opencv.org/license/\">opencv.org/license/</a>"
	));
	sections->addTab(openCV, "OpenCV 4");

	QVBoxLayout *vboxlayout = new QVBoxLayout;
	vboxlayout->addWidget(sections, 0);

	setLayout(vboxlayout);
}

About::~About() {

}

