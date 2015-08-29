////// Using FLTK 
#include "GUIFrame.h"

//int main(int argc, char **argv) 
//{  
//	Fl_Window *window = new Fl_Window(340,180);  
//	Fl_Box *box = new Fl_Box(20,40,300,100,"Hello, World!");  
//	box->box(FL_UP_BOX);  
//	box->labelfont(FL_BOLD+FL_ITALIC);  
//	box->labelsize(36);  
//	box->labeltype(FL_SHADOW_LABEL);  
//	GLCanvas *can = new GLCanvas(0, 0, 800, 800);
//	can->show();
//	window->end();  
//	window->resize(0, 0, 1000, 1000);
//	window->show(argc, argv);  
//	return Fl::run();
//}

int main(int argc, char **argv)
{
	//GLCanvas *can = new GLCanvas(0, 0, 800, 800);
	//can->show();
	//return Fl::run();

	GUIFrame app("DGP Studio (Subdivision Surface Fitting)", 800, 640);
	app.add_mesh_model(new DecimationModel);
	return app.run();

}

//
////
////// Using QT4.5
//#include <QTextCodec>
//#include <QApplication>
//#pragma  comment (lib, "QtOpenGL4.lib")
//#pragma  comment (lib, "QtGui4.lib")
//#pragma  comment (lib, "QtCore4.lib")
//#pragma  comment (lib, "qtmain.lib")
//#include "mainwindow.h"
//
//#include <iostream>
//
//int main(int argc, char *argv[])
//{
//	//QTextCodec::setCodecForTr(QTextCodec::codecForName("GB2312"));
//	//QTextCodec::setCodecForLocale(QTextCodec::codecForName("GB2312"));
//	//QTextCodec::setCodecForCStrings(QTextCodec::codecForName("GB2312"));
//
// //	QTextCodec::setCodecForTr(QTextCodec::codecForName("UTF8"));
// //	QTextCodec::setCodecForCStrings(QTextCodec::codecForName("UTF8"));
//
//	QTextCodec::setCodecForTr(QTextCodec::codecForLocale());
//	QTextCodec::setCodecForCStrings(QTextCodec::codecForLocale());
//
//	QApplication app(argc, argv);
//	MainWindow mainWin;
//	mainWin.show(); 
//	return app.exec(); 
//}
