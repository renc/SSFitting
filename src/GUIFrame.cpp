// ***************************************************************
//  
//  GUIFrame 
//  Copyright (C) 2007 - by rencanjiang All Rights Reserved
//  rencanjiang@163.com	
//	-------------------------------------------------------------
//  
//  version: 1.0
//  data: 12/18/2007
// ***************************************************************
//  
//  File description:
//  Imitate M.Garland's libgfx and qslim to build up a framework.
// ***************************************************************

#include <cstdio>
#include <cstdarg>
#include <fstream>

#include <OpenMesh/Tools/Utils/Timer.hh>
//#include "../../../src/CourseExamples/01-MeshViewer/gl.hh" // renc 20150829 //

#include "GUIFrame.h"

GUIFrame *GUIFrame::current = NULL;
// ------------------------
// class GLCanvas
// ------------------------
void GLCanvas::resize(int x, int y, int w, int h)
{	//std::cout << "GLCanvas.resize().\n";
    Fl_Gl_Window::resize(x, y, w, h);

    if( shown() ) {
        make_current();// api from fltk
        //glViewport(0, 0, w, h);//�ӿڱ任 //
		if (app) { 
			app->reshape(w, h); 
		}
        invalidate();// api from fltk 
    }
}

void GLCanvas::draw()
{	 
    if( !valid() ) {//��Ч������ "valid() is turned on by FLTK after draw() returns"
        valid(1);//ʹ����Ч
        if(app) app->setup_for_drawing();//��ʼ��������ͶӰ�任
		// ��������ʹ��ÿ��resize function֮�����µ���set_for_drawing(), ����invoke init_display() function, 
		// ������glut��, ���init_display()Ӧ��ֻ�ǿ�ʼ������һ�ε�.֮��even if resize, there's no need to init_display() again.
    }

    if(app) app->draw_contents();
}

int GLCanvas::handle(int event)
{
    bool need_redraw = false;

    int where[2];  where[0] = Fl::event_x();  where[1] = Fl::event_y();//��굱ǰ�����¼���λ��

    // NOTE: Normally, we examine event_state() rather than
    // event_button() because it is valid no matter what the generating
    // event whereas event_button() is only valid for PUSH and RELEASE
    // events.  However, since event_state() only tells us what buttons
    // are *pressed*, we need to revert to event_button() on RELEASE
    // events.
    //
    int which = Fl::event_button();//which��ȡֵ��1,2,3

    if( event != FL_RELEASE ) { //�м�������ʱ�ж����ĸ���
		if( Fl::event_state(FL_BUTTON1) ) {
			// emulate middle button by combination of left & right buttons
			if( Fl::event_state(FL_BUTTON3) )  which = 2; //�������ֻ���������Ļ�, ͬʱ���¾ͱ�ʾ�м�
			else                               which = 1;
		} else if( Fl::event_state(FL_BUTTON2) ) {
			which = 2;
		} else if( Fl::event_state(FL_BUTTON3) ) {
			which = 3;
		}
    }

    switch( event ) {
		case FL_FOCUS:
		case FL_UNFOCUS:
		// Do nothing special �����¼���Ӧ��û�¸ɵ�
		break;
		
		case FL_PUSH:
		need_redraw = app && app->mouse_down(where, which);
		last_click[0]=where[0];  last_click[1]=where[1];
		break;

		case FL_RELEASE:
		need_redraw = app && app->mouse_up(where, which);
		break;

		case FL_DRAG:               //��굱ǰλ��, �ϴ�λ��,   �ĸ�����
		need_redraw = app && app->mouse_drag(where, last_click, which);
		last_click[0]=where[0];  last_click[1]=where[1];
		break;

	    case FL_KEYBOARD:
		if( !app || !app->key_press(Fl::event_key()) )
			return 0;
		break;

		default:
		return Fl_Gl_Window::handle(event);
    }

    if( need_redraw ) redraw();

    return 1;
}

// --------------------------------
// class GUIFrame
// --------------------------------
GUIFrame::GUIFrame(const char* title, const int width, const int height)
:canvas(NULL), status_line(NULL) {
	Fl::scheme("plastic");		// optional
	// ��������㴰�� ---------------
	// ������width, height��ʵ�������������������canvas��͸�, ��������ڵĿ�͸�Ҫ�����
	// ����pad��canvas�͵��µ�״̬���뿪��������߿�ľ���, Ŀ���Ǻÿ�����
	int pad=5;
    //Fl::visual(FL_RGB8);

	int yfill = 0;//�봰�ڶ���y�߶�, �����ʼΪ��Ϳ��Ե���
	const unsigned int iMenuBarHeight = 30;
	const unsigned int iCanvasHeight = height; 
	const unsigned int iStatusBarHeight = 30; 

	Fl_Window *w = new Fl_Window(width+2*pad, iMenuBarHeight + iCanvasHeight + iStatusBarHeight + pad * 2);
    toplevel = w;//�ǵ�Ч�������������е���toplevel = create_window(,,);
	toplevel->label(title);
	w->box(FL_UP_BOX);
		 
		// �˵� -------------�������������в˵���Ϊ��ͳһ�ù���.
		menu_bar = new Fl_Menu_Bar(0, yfill, w->w(), iMenuBarHeight);
		Fl_Menu_Item menu_layout[] = {
			{ "&File", FL_ALT + 'f', 0, 0, FL_SUBMENU },
				{ "&Open File...",		FL_CTRL + 'o',	cb_open},
				{ "&Save as...",		0,				cb_file_save_as},
				{ "E&xit",				FL_CTRL + 'x',	cb_exit, 0 }, 
			{ 0 },

			{ "&View", 0, 0, 0, FL_SUBMENU },
				{ "&DrawModle", 0, 0, 0, FL_SUBMENU },
					{ "Wireframe",		0,				cb_view_wireframe},
					{ "HiddenLine",		0,				cb_view_hiddenline},
					{ "SolidFlat",		0,				cb_view_solidflat},
					{ "SolidSmooth",	0,				cb_view_solidsmooth},
					{ "Line&SolidFlat",	0,				cb_view_lineandsolidflat},					
				{ 0 },
				{ "&DrawModle", 0, 0, 0, FL_SUBMENU },
					//{ "��ʾ������",		0,				cb_draw_coordinate_system_true},
					//{ "����������",		0,				cb_draw_coordinate_system_false},
					{ "Turn on Coordinates",		0,				cb_draw_coordinate_system_true},
					{ "Turn off Coordiante",		0,				cb_draw_coordinate_system_false},
				{ 0 },

				//
				//{ "Curvature",		0,				cb_view_curvature},
				//{ "Mean Saliency",	0, 0, 0, FL_SUBMENU },
				//	{ "Saliency 1",		0,			cb_view_saliency_1},
				//	{ "Saliency 2",		0,			cb_view_saliency_2},
				//	{ "Saliency 3",		0,			cb_view_saliency_3},
				//	{ "Saliency 4",		0,			cb_view_saliency_4},
				//	{ "Saliency 5",		0,			cb_view_saliency_5},
				//	{ "MultiScales",	0,			cb_view_saliency_multiscales},
				//{ 0 },
			{ 0 },
			//
			//{ "Smoothing", 0, 0, 0, FL_SUBMENU },
			//	{ "Laplacian Smoothing",0,				cb_smoothing_laplaciansmoothing },
			//{ 0 }, 

			{ "&Edit", 0, 0, 0, FL_SUBMENU },
				{ "Decimation QEM", 0, 0, 0, FL_SUBMENU },
					//{ "To Percentage",	0,				cb_fitting_decimation_percentage },
					{ "To Vertices Num",0,				cb_fitting_decimation_vertices_num },
					{ "draw on/off mesh_collapsed_", 0, cb_draw_mesh_collapsed_or_not },
				{ 0 },
				//{ "Refinement", 0, 0, 0, FL_SUBMENU },
				//	{ "Sequence Refine",0,				cb_fitting_sequence_refine },
				//	{ "SelectiveRefine",0,				cb_fitting_selective_refine },
				//{ 0 },
				//{ "�Բ������ֲ�ƽ��",	0,				cb_local_smooth_parameterization },
				{ "Resample midpoint on simplified_mesh_", 0, cb_resample_edge_midpoint }, //{ "��ģ�͵��е����", 0,				cb_resample_edge_midpoint }, 
				{ "Midpoint Subdivision one time",	0,					cb_midpointsubdivision_1 },//{ "�е�ϸ��һ��",	0,					cb_midpointsubdivision_1 },
				{ "Get initial_control_mesh_", 0,				cb_create_initialcontrolmesh },//{ "�������Ŀ�������", 0,				cb_create_initialcontrolmesh },
				{ "Feature perseving Loop subd one time", 0,			cb_loop_subdivision },//{ "��������Loopϸ��һ��", 0,			cb_loop_subdivision },
				{ "Create limit surface",	0,				cb_create_limit_surface },//{ "����limit surface",	0,				cb_create_limit_surface },
				//{ "һ������fitting",	0,				cb_fitting },
				//{ "Error Driven refinement", 0,			cb_error_driven_selective_refine }, 
				{ "Fitting Error", 0, 0, 0, FL_SUBMENU },
					{ "Evaluation&&Calculate error", 0,			cb_evaluation },
					{ "View fitting error of mesh2_", 0,		cb_view_fitting_error },
					{ "View fitting error of simplified_mesh_", 0,		cb_view_fitting_error2 },
					{ "Refine", 0, cb_refine },
				{ 0 },
			{ "====================", 0, NULL },
				{ "try bounding",		0,				cb_try_boudnding },
				{"draw on/off mesh2_",	0,				cb_draw_mesh2},
			{ 0 }, 

			{ "&Help", 0, 0, 0, FL_SUBMENU },
				{ "&About...",		FL_CTRL + 'a',	cb_help_about},
			{ 0 },

			{ 0 }
		};
		menu_bar->copy(menu_layout, this);//�˵�λ����x=0, y=yfill=0, width=w->w(), height=30;
		//����this�����Ǻ���Ҫ�ķ���callback�����оͲ�����ʹ��void *dataת��ΪGUIFrame������, 
		//�����Ļ���ֻ��ʹ��GUIFrame�е�currentָ��ָ��ǰ������.

		add_menu("&File/&try", 0, (Fl_Callback *)NULL);
		yfill += menu_bar->h();

		yfill += pad;
		 		
		// canvas���� -------------------------------------
		canvas = new GLCanvas(pad, yfill, width, iCanvasHeight);//
		canvas->box(FL_DOWN_FRAME);
		canvas->attach_app(this);

		int glmode = 0;//��Ⱦģʽ, ���ظ�ʽpixel format
		if(canvas->can_do(FL_RGB8))    glmode|=FL_RGB8;
		if(canvas->can_do(FL_DOUBLE))  glmode|=FL_DOUBLE;
		if(canvas->can_do(FL_DEPTH))   glmode|=FL_DEPTH;
		if(canvas->can_do(FL_ALPHA))   glmode|=FL_ALPHA;// mode(int )��class Fl_Gl_Window�ĺ���
		if(glmode)
			canvas->mode(glmode);//��ʵ����Fl_Gl_Window��draw��������ǰ���õ�

		yfill += canvas->h();

		yfill += pad;
		
		// ���µ�״̬��Ϣ�� --------------------------
		status_line = new Fl_Output(pad, yfill, width, iStatusBarHeight);//
		status_line->color(48);
		status_line->labeltype(FL_NO_LABEL);
		yfill += status_line->h();
		
	w->end();// end() ��ʾ�Ѿ����ؼ�����������

    w->size(w->w(), yfill+pad);	// adjust window height ���ں�����ʼʱ��߶Ȼ�û�ж�
    w->resizable(*canvas);		// resize canvas with window, cavas�Ͷ��㴰��һ�������С
	
    // These are used by resize_canvas() to resize the window based on
    // the target size of the canvas.
    w_offset = w->w() - width;
    h_offset = w->h() - height;

	// init mouse buttons
	for (int i=0; i<10; ++i)
		button_down_[i] = false;
	
	GUIFrame::current = this;

	has_init_display_ = false;
	background_color_ = white_color_glframe;
} 
////////////////////////////////////////////////////////////////////////
int GUIFrame::add_menu(const std::string& s, int key, Fl_Callback *cb, int flags)
{
    return menu_bar->add(s.c_str(), key, cb, this, flags);
	//�����this, ʹ���������е���Ӧ����void GUIFrame::cb_function(Fl_Widget *w, void *data)
	//��void *data������ת��ΪGUIFrame *gui = static_cast<GUIFrame*>(data)����ʹ����GUIFrame�еĳ�Ա����.
}

int GUIFrame::add_toggle_menu(const std::string& s, int key, bool& val, int flags)
{	// ����һ���ر��ѡ��ť�˵�(��ѡ�ɲ�ѡ), ����SMFView��Ŀ��֪����
	bool will_draw_surface;
	add_toggle_menu("&View/Will draw/Surface", FL_CTRL+'s', will_draw_surface);//����ǲ�����bool����
    return menu_bar->add(s.c_str(), key, (Fl_Callback *)cb_toggle, &val, FL_MENU_TOGGLE|(val?FL_MENU_VALUE:0)|flags);
}
void GUIFrame::cb_toggle(Fl_Menu_ *m, bool *flag)
{   //��add_toggle_menu����ʹ��, ��������ʹ�õ�
    *flag = m->mvalue()->value()!=0;
    current->canvas->redraw();
}

///////////////////////////////////////////////////////////////////////
//�������������״̬��Ϣ��������
int GUIFrame::status(const char *fmt, ...)
{	
	static char strbuf[1000];

    va_list args;
    va_start(args, fmt);
    int n = vsprintf(strbuf, fmt, args);

    status_line->value(strbuf);
    return n;
}

void GUIFrame::lock_size()
{
    toplevel->size_range(toplevel->w(), toplevel->h(), toplevel->w(), toplevel->h());
    toplevel->resizable(NULL);
}
void GUIFrame::unlock_size()
{
    toplevel->resizable(*canvas);
    toplevel->size_range(100, 100, 0, 0);
}
void GUIFrame::resize_canvas(int width, int height)
{	//��ʵ��һ��Ϊ�˷���ĺ���, �������ö��㴰�ڵĴ�С, �����canvas�Ĵ�СҲ���Ÿı�
    toplevel->size(width+w_offset, height+h_offset);
    toplevel->redraw();
}

////////////////////////////////////////////////////////////////////////
// ��Ӧ�˵�����
// Default menu system.  Although the application can completely
// override this, it's expected that most programs will just add their
// own additional entries.
// �������漸���������Բ���override��, ����cb_open()���ǿ��ԸĽ�������������ģ���ļ�
//	menu File -------------
void GUIFrame::cb_open(Fl_Widget *w, void *data) { 
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->open();
}
void GUIFrame::open() {
	//fl_alert("This application has not defined this function.");	
	std::cout << std::endl;
	std::cout << "========�û�ѡ������ģ�͵��ļ�==========" << std::endl;//for test
	std::cout << "open file..." << std::endl;//for test
	Fl_File_Chooser		*fc;
	fc = new Fl_File_Chooser(".", "*", Fl_File_Chooser::SINGLE, "Choose the model file.");
	fc->callback(NULL);//��filechoose�е��ĳ���ļ�/�ļ���,���߰�oK/cancle��������õ�

	fc->show();
	while (fc->visible())
		Fl::wait();

	const char *filename = NULL;
	filename = fc->value();
	if (filename != NULL) { ////�û�ѡ��ģ���ļ�ʱ����OK����
		std::cout << "file name choosed: " << filename << std::endl;//for test
		this->status(filename);
		if (mesh_model_->open_mesh(filename)) {
			std::cout << "file open succeed.\n";
		} else {
			// open file fail,������ѡ��Ĳ��Ǻ��ʵ�ģ���ļ�(off, obj)
			std::cout << "file open failed.\n";
		}
	} else {//�û�ѡ��ģ���ļ�ʱ����Cancel����
		std::cout << "User didn't choose file." << std::endl;
	}
	this->canvas->redraw();
	std::cout << "========�û�ѡ���ļ��������============" << std::endl;//for test

}
void GUIFrame::cb_file_save_as(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->save_as();
	
}
void GUIFrame::save_as() {
	char msg[80], pat[8], name[16];
	sprintf(msg, "Save this model to file: ");
	sprintf(pat, "*.%s", "off");//��׺ 
	sprintf(name, "new_model.%s", "off");//��ʾ�Ŀ�ѡ���� //����Ӧ�ü��뱣֤�û�����off������objΪ�����ʽ��
	char *file_output_name = fl_file_chooser(msg, pat, name);//"zzzzzzzzz.off";
	if (file_output_name) {
		mesh_model_->save_as(file_output_name);
		std::cout << "GUIFrame::model save as: " << file_output_name << ".\n";//
	}
}
void GUIFrame::cb_exit(Fl_Widget *, void *data)
{
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->exitt();
}
void GUIFrame::exitt() {
	if (!fl_choice("Are you sure you want to quit?", "Cancel", "Quit", 0L)) return;
	exit(0); 
}
void GUIFrame::cb_help_about(Fl_Widget *w, void *data) {
	//fl_alert("This application has not defined.");//!����
	fl_message("DGP Application!\n����:�βӽ�. ����:2007.10.");//i����
}

//-----------------------------------------------------------------------------
void GUIFrame::cb_view_wireframe(Fl_Widget *w, void *data) { 
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->set_view_drawmodel(MeshModel::WIREFRAME);
	//gui->canvas->redraw();
}
void GUIFrame::cb_view_hiddenline(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->set_view_drawmodel(MeshModel::HIDDENLINE);
}
void GUIFrame::cb_view_solidflat(Fl_Widget *w, void *data) { 
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->set_view_drawmodel(MeshModel::SOLIDFLAT);
	//gui->canvas->redraw();
}
void GUIFrame::cb_view_solidsmooth(Fl_Widget *w, void *data) { 
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->set_view_drawmodel(MeshModel::SOLIDSMOOTH);
	//gui->canvas->redraw();
}
void GUIFrame::cb_view_lineandsolidflat(Fl_Widget *w, void *data) { 
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->set_view_drawmodel(MeshModel::LINEANDSOLIDFLAT);
	//gui->canvas->redraw();
}
void GUIFrame::cb_draw_coordinate_system_true(Fl_Widget *w, void *data) { 
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->set_draw_coordinate_system(true);
	//gui->canvas->redraw();
}
void GUIFrame::cb_draw_coordinate_system_false(Fl_Widget *w, void *data) { 
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->set_draw_coordinate_system(false);
	//gui->canvas->redraw();
}


//-----------------------------------------------------------------------------
/*void GUIFrame::cb_view_curvature(Fl_Widget *w, void *data) { 
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->view_curvature();
	//gui->canvas->redraw();
}
void GUIFrame::cb_view_saliency_1(Fl_Widget *w, void *data) { //ֻ����ʾscale 1Ҳ����sigma1 = 2 * epsilon�����saliencyֵ.
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->view_saliency(1);
}
void GUIFrame::cb_view_saliency_2(Fl_Widget *w, void *data) { //ֻ����ʾscale 2Ҳ����sigma2 = 3 * epsilon�����saliencyֵ.
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->view_saliency(2);
}
void GUIFrame::cb_view_saliency_3(Fl_Widget *w, void *data) { //ֻ����ʾscale 3Ҳ����sigma3 = 4 * epsilon�����saliencyֵ.
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->view_saliency(3);
}
void GUIFrame::cb_view_saliency_4(Fl_Widget *w, void *data) { //ֻ����ʾscale 4Ҳ����sigma4 = 5 * epsilon�����saliencyֵ.
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->view_saliency(4);
}
void GUIFrame::cb_view_saliency_5(Fl_Widget *w, void *data) { //ֻ����ʾscale 5Ҳ����sigma5 = 6 * epsilon�����saliencyֵ.
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->view_saliency(5);
}
void GUIFrame::cb_view_saliency_multiscales(Fl_Widget *w, void *data) {// ��5��scales����non-linear suppression, ����adding.
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->view_saliency(6);// ע��������6��ʾ��5��scales���ۺ�
}
void GUIFrame::cb_smoothing_laplaciansmoothing(Fl_Widget *w, void *data) { 	
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->smoothing_laplaciansmoothing();
	//gui->canvas->redraw();
}*/
//-----------------------------------------------------------------------------
void GUIFrame::cb_fitting_decimation_percentage(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;

	const char *input = fl_input("��������������ϣ���򻯵�ʲô�ٷֱ�. \n����0��100��Ч:", "");// �ڶ�����������ʾ������������
	double per = 0.0;
	if (input) { //����Ӧ�����ж�input�ǲ��������ַ���,��ת����0��100�İٷֱ�
		per = atof(input) / 100.0;//std::cout << input << ", " << per << "\n";//for test

		gui->decimation_process((unsigned int)(gui->mesh_model_->n_vertices() * per));
	}	
}

void GUIFrame::cb_fitting_decimation_vertices_num(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;

	std::string label("��������������ϣ���򻯺�Ķ�����.\n����0��");
	char n[100];
	itoa(gui->mesh_model_->n_vertices(), n, 10);
	label += n;label +="��Ч.";

	const char *input = fl_input(label.c_str(), "");// �ڶ�����������ʾ������������
	unsigned int num = 0;
	if (input) { //����Ӧ�����ж�input�ǲ��������ַ���,��ת�����������
		num = atoi(input);//std::cout << input << ", " << num << ", " << n << "\n";//for test

		gui->decimation_process(num);	
	}	
}
void GUIFrame::cb_draw_mesh_collapsed_or_not(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->draw_mesh_collapsed_or_not_process();
}

void GUIFrame::cb_local_smooth_parameterization(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame *>(data);
	gui->local_smooth_parameterization_process();
}
void GUIFrame::cb_fitting_sequence_refine(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame *>(data);
	gui->sequence_refine();
}

void GUIFrame::cb_fitting_selective_refine(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->selective_refine();
}
void GUIFrame::cb_resample_edge_midpoint(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->resample_edge_midpoint_process();
}

void GUIFrame::cb_midpointsubdivision_1(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->midpointsubdivision_1_process();
}
void GUIFrame::cb_create_initialcontrolmesh(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->create_initialcontrolmesh_process();
}// ��initial control mesh
void GUIFrame::cb_loop_subdivision(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->loop_subdivision_process();
}
void GUIFrame::cb_create_limit_surface(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->create_limit_surface_process();
}
void GUIFrame::cb_fitting(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame *> (data);
	gui->fitting_process();
}
void GUIFrame::cb_error_driven_selective_refine (Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->error_driven_selective_refine_process();
}
void GUIFrame::cb_evaluation(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->evaluation_process();
}
void GUIFrame::cb_view_fitting_error(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->view_fitting_error_mesh2_process();
}
void GUIFrame::cb_view_fitting_error2(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->view_fitting_error_sm_process();
}
void GUIFrame::cb_refine(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->refine_process();
}
void GUIFrame::cb_try_boudnding(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->try_bounding_process(); 
}
void GUIFrame::cb_draw_mesh2(Fl_Widget *w, void *data) {
	GUIFrame *gui = static_cast<GUIFrame*>(data);//(GUIFrame *)data;
	gui->draw_mesh2_process(); 
}
////////////////////////////////////////////////////////////////////////
// Override these methods to control the contents of the GL canvas
void GUIFrame::setup_for_drawing() {
	if (has_init_display_ == false)
	init_display();

	has_init_display_ = true;
}
void GUIFrame::draw_contents() {
	display();
}
// Override these methods to receive events from the GL canvas
bool GUIFrame::mouse_down(int *where, int which) {
	// ����������ǽ�fltk�µ���Ӧת����Ҵ�GlutViewer�ù�����glut�汾����Ӧ
	mouse(which - 1, 0, where[0], where[1]);//whichȡֵ��1,2,3; �ڶ�������0��ʾGLUT_DOWN
	return true;//which��1����Ϊfltk�汾��whichȡֵ1,2,3�ֱ��ʾ3������,������glut�汾����0,1,2�ֱ��ʾ3������
}
bool GUIFrame::mouse_up(int *where, int which) {
	mouse(which - 1, 1, where[0], where[1]);//whichȡֵ��1,2,3; �ڶ�������0��ʾGLUT_UP
	return true;
}
bool GUIFrame::mouse_drag(int *where, int *last, int which) {
	motion(where[0], where[1]);//ǰ���mouse_down() invoke mouse() has remenber the "which"
	return true;
}
bool GUIFrame::key_press(int key) {
	switch(key) {
		// �����Լ�����Ӧ����,����������keyboard����,��ʵҲû�б�Ҫ,ֻ��Ϊ���ö���
		case FL_Escape: exitt(); break;
		default: keyboard(key,0,0); break;
	}
	return true;
}
/////////////////////////////////////////////////////////////////////////////
//=============================================================================
//-----------------------------------------------------------------------------
// protected, virtual, invoked by constructor, 
// set the background color and lighting and material properies
void GUIFrame::init_display()
{	//std::cout << "GUIFrame.init_display().\n";
	// OpenGL state
	glClearColor(background_color_[0], background_color_[1], background_color_[2], 0.0);//background is black color
	glEnable( GL_DEPTH_TEST );

	// some performance settings
	//   glEnable( GL_CULL_FACE );
	glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE );
	glLightModeli( GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE );

	// material
	GLfloat mat_a[] = {0.2f, 0.2f, 0.2f, 1.0};
	GLfloat mat_d[] = {0.4f, 0.4f, 0.4f, 1.0};
	GLfloat mat_s[] = {0.8f, 0.8f, 0.8f, 1.0};
	GLfloat shine[] = {128.0};//0~10, һ������, 10~128����
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_a);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_d);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_s);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);

	// lighting
	glLoadIdentity();

	GLfloat pos1[] = { 0.1f, 0.1f, -0.02f, 0.0};//��������3��lighting
	GLfloat pos2[] = {-0.1f, 0.1f, -0.02f, 0.0};
	GLfloat pos3[] = { 0.0f, 0.0f, 0.1f, 0.0};
	GLfloat col1[] = {.05f, .05f, .6f, 1.0};
	GLfloat col2[] = {.6f, .05f, .05f, 1.0};
	GLfloat col3[] = {1.0, 1.0, 1.0, 1.0};

	glEnable(GL_LIGHT0);    
	glLightfv(GL_LIGHT0,GL_POSITION, pos1);
	glLightfv(GL_LIGHT0,GL_DIFFUSE,  col1);
	glLightfv(GL_LIGHT0,GL_SPECULAR, col1);

	glEnable(GL_LIGHT1);  
	glLightfv(GL_LIGHT1,GL_POSITION, pos2);
	glLightfv(GL_LIGHT1,GL_DIFFUSE,  col2);
	glLightfv(GL_LIGHT1,GL_SPECULAR, col2);

	glEnable(GL_LIGHT2);  
	glLightfv(GL_LIGHT2,GL_POSITION, pos3);
	glLightfv(GL_LIGHT2,GL_DIFFUSE,  col3);
	glLightfv(GL_LIGHT2,GL_SPECULAR, col3);

	// scene pos and size
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix_);
	set_scene(OpenMesh::Vec3f(0.0, 0.0, 0.0), 1.0);//center, radius
}


//-----------------------------------------------------------------------------
void GUIFrame::display(void) 
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// ������һ�䱾���ǲ���Ҫ��, ��������ɫ����ı�Ļ�
  	glClearColor(background_color_[0], background_color_[1], background_color_[2], 0.0);

	//draw();
	mesh_model_->draw(); 
	this->canvas->redraw();
}

void GUIFrame::draw() {	
	//return;
}

// protected, virtual, invoked when the window's size is changed 
void GUIFrame::reshape(int _w, int _h)
{	//�⺯������glut����GlutViewer��ʹ�õ�,��fltk���Ȳ�����
	//std::cout << "GUIFrame.reshape().\n";
	width_  = _w; 
	height_ = _h;
	glViewport(0, 0, _w, _h);
	update_projection_matrix();
	//glutPostRedisplay();
}

//----------------------------------------------------------------------------
// protected, ֻ��reshape������set_scene��������
void GUIFrame::update_projection_matrix()
{
	width_ = canvas->w();
	height_ = canvas->h();//��������fltk�¸����
	
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	gluPerspective(fovy_, (GLfloat)width_/(GLfloat)height_, near_, far_);//����glut�º�reshape����һ�µ�
	glGetDoublev( GL_PROJECTION_MATRIX, projection_matrix_);
	glMatrixMode( GL_MODELVIEW );
}


//----------------------------------------------------------------------------
// public, ��init����������, and invoked right after load a new model
void GUIFrame::set_scene( const OpenMesh::Vec3f& _cog, float _radius )
{
	center_ = _cog;//(0,0,0)
	radius_ = _radius;//1

	near_  = 0.01f * radius_;
	far_   = 100.0f * radius_;
	fovy_ = 45.0;
	update_projection_matrix();

	view_all();
}

//----------------------------------------------------------------------------
// ֻ�Ǳ�set_scene��������
void GUIFrame::view_all()
{  // modelview_matrix_ * center_
	translate( OpenMesh::Vec3f( -(modelview_matrix_[0]*center_[0] + 
						modelview_matrix_[4]*center_[1] +
						modelview_matrix_[8]*center_[2] + 
						modelview_matrix_[12] + 0*1),
					  -(modelview_matrix_[1]*center_[0] + 
						modelview_matrix_[5]*center_[1] +
						modelview_matrix_[9]*center_[2] + 
						modelview_matrix_[13] + 0*1),
					  -(modelview_matrix_[2]*center_[0] + 
						modelview_matrix_[6]*center_[1] +
						modelview_matrix_[10]*center_[2] + 
						modelview_matrix_[14] + 3.0f*radius_*1) 
						) ); //z�������ƶ��ĸ���, ϵ��3�ĳɸ���Ļ�(����6.0)����ȥ�����С,��Ϊ���ӵ��Զ��
}

//----------------------------------------------------------------------------
// protected, used in virtual traceball,
bool GUIFrame::map_to_sphere( const OpenMesh::Vec2i& _v2D, OpenMesh::Vec3f& _v3D )
{
	width_ = canvas->w();//����֮ǰGlutViewer�²��õ�,��Ϊ��width_������reshape�����»������ֵ,
	height_ = canvas->h();//������fltk�汾��,��ʹ��reshape������.����width_��ֵ��Ҫ��ʹ��ǰ�ֶ�����

	if ((_v2D[0] >= 0) && (_v2D[0] <= width_) &&
		(_v2D[1] >= 0) && (_v2D[1] <= height_)) 
	{
		double x  = (double)(_v2D[0] - 0.5*width_)  / (double)width_;
		double y  = (double)(0.5*height_ - _v2D[1]) / (double)height_;
		double sinx         = sin(M_PI * x * 0.5);
		double siny         = sin(M_PI * y * 0.5);
		double sinx2siny2   = sinx * sinx + siny * siny;

		_v3D[0] = sinx;
		_v3D[1] = siny;
		_v3D[2] = sinx2siny2 < 1.0f ? sqrt(1.0 - sinx2siny2) : 0.0f;

		return true;
	}
	else return false;
}

//-----------------------------------------------------------------------------
// protected, virtual, invoked when the mouse is pressed or released
void GUIFrame::mouse(int button, int state, int x, int y)
{
	// mouse press
	if (state == GLUT_DOWN) // #define GLUT_DOWN 0
	{
		last_point_2D_ = OpenMesh::Vec2i(x,y);
		last_point_ok_ = map_to_sphere( last_point_2D_, last_point_3D_ );
		button_down_[button] = true;
	}

	// mouse release
	else // #define GLUT_UP 1
	{
		last_point_ok_ = false;
		button_down_[button] = false;

		// GLUT: button 3 or 4 -> mouse wheel clicked
		if (button == 3) {
			width_ = canvas->w();//����֮ǰGlutViewer�²��õ�,��Ϊ��width_������reshape�����»������ֵ,
			zoom(0, (int)(last_point_2D_[1] - 0.05*width_));			
		} else if (button == 4) {
			width_ = canvas->w();//����֮ǰGlutViewer�²��õ�,��Ϊ��width_������reshape�����»������ֵ,
			zoom(0, (int)(last_point_2D_[1] + 0.05*width_));
		}
	}

	//glutPostRedisplay();// fltk�汾���Լ������Fl_Gl_Window��redraw������,���Բ���Glut�汾�µ�glutPostRedisplay������
}
//-----------------------------------------------------------------------------
void GUIFrame::motion(int x, int y)
{
	// zoom // �����Ҽ����ڰ��˲˵�������һ������ģ������� #define GLUT_RIGHT_BUTTON 2
	// if (button_down_[0] && button_down_[1])
	if(button_down_[2])//��������fltk�汾�¼������Ӧ�Ҽ�
	{
		zoom(x, y);
	}
	// rotation
	else if (button_down_[0]) // ��� #define GLUT_LEFT_BUTTON 0
	{
		rotation(x, y);
	}
	// translation
	else if (button_down_[1]) // �м�,#define GLUT_MIDDLE_BUTTON 1
	{
		translation(x, y);
	}

	// remeber points
	last_point_2D_ = OpenMesh::Vec2i(x, y);
	last_point_ok_ = map_to_sphere(last_point_2D_, last_point_3D_);

	//glutPostRedisplay();// fltk�汾���Լ������Fl_Gl_Window��redraw������,���Բ���Glut�汾�µ�glutPostRedisplay������
}

//-----------------------------------------------------------------------------
void GUIFrame::rotation(int x, int y)
{
	if (last_point_ok_) 
	{
		OpenMesh::Vec3f  new_point_3D;
		
		OpenMesh::Vec2i  new_point_2D = OpenMesh::Vec2i(x, y);
		bool   new_point_ok = map_to_sphere(new_point_2D, new_point_3D);

		if (new_point_ok)
		{
			OpenMesh::Vec3f axis      = (last_point_3D_ % new_point_3D);
			float cos_angle = (last_point_3D_ | new_point_3D);

			if (fabs(cos_angle) < 1.0) 
			{
				float angle = 2.0f*acos(cos_angle) * 180.0f / M_PI;
				rotate(axis, angle);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void GUIFrame::translation(int x, int y)
{ // x, y are the current mouse location.
	float dx = x - last_point_2D_[0];
	float dy = y - last_point_2D_[1];

	float z = - ((modelview_matrix_[ 2]*center_[0] + 
		modelview_matrix_[ 6]*center_[1] + 
		modelview_matrix_[10]*center_[2] + 
		modelview_matrix_[14]) /
		(modelview_matrix_[ 3]*center_[0] + 
		modelview_matrix_[ 7]*center_[1] + 
		modelview_matrix_[11]*center_[2] + 
		modelview_matrix_[15]));

	float aspect = (float)width_ / (float)height_;
	float up     = tan(fovy_/2.0f*M_PI/180.f) * near_;
	float right  = aspect*up;

	translate(OpenMesh::Vec3f(2.0f*dx/width_*right/near_*z, 
		-2.0f*dy/height_*up/near_*z, 
		0.0f));
}


//-----------------------------------------------------------------------------
void GUIFrame::zoom(int x, int y)
{
	float dy = y - last_point_2D_[1];
	float h  = height_;
	translate(OpenMesh::Vec3f(0.0f, 0.0f, radius_ * dy * 3.0f / h));
}

//----------------------------------------------------------------------------
// support the above rotation and translation function.
void GUIFrame::translate( const OpenMesh::Vec3f& _trans )
{
	glLoadIdentity();
	glTranslated( _trans[0], _trans[1], _trans[2] );
	glMultMatrixd( modelview_matrix_ );
	glGetDoublev( GL_MODELVIEW_MATRIX, modelview_matrix_);
}

//----------------------------------------------------------------------------
void GUIFrame::rotate( const OpenMesh::Vec3f& _axis, float _angle )
{
	OpenMesh::Vec3d t( modelview_matrix_[0]*center_[0] + 
		modelview_matrix_[4]*center_[1] +
		modelview_matrix_[8]*center_[2] + 
		modelview_matrix_[12],
		modelview_matrix_[1]*center_[0] + 
		modelview_matrix_[5]*center_[1] +
		modelview_matrix_[9]*center_[2] + 
		modelview_matrix_[13],
		modelview_matrix_[2]*center_[0] + 
		modelview_matrix_[6]*center_[1] +
		modelview_matrix_[10]*center_[2] + 
		modelview_matrix_[14] );

	glLoadIdentity();
	glTranslated(t[0], t[1], t[2]);
	glRotated( _angle, _axis[0], _axis[1], _axis[2]);
	glTranslated(-t[0], -t[1], -t[2]); 
	glMultMatrixd(modelview_matrix_);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix_);
}


//-----------------------------------------------------------------------------
void GUIFrame::keyboard(int key, int x, int y) 
{
	switch (key)
	{
	case '0'://����0.
		{
			std::cerr << "Performance test: ";
			double fps = measure_fps();
			std::cerr << fps << " FPS\n";
			break;
		}

	case 'q': case 'Q': case 27: 
		exit(0); 
		break;

	case '1':
		fitting_process(); break;
	case '2':
		error_driven_selective_refine_process(); break;
	case '3':
		evaluation_process(); break;
	case '9':
		mesh_model_->test9();
	default:
		{
			break;
		}
	}
	this->canvas->redraw();
}
//-----------------------------------------------------------------------------
double GUIFrame::measure_fps()
{
	double fps(0.0);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	unsigned int  frames = 90;
	const float   angle  = 360.0f/(float)frames;
	unsigned int  i;
	OpenMesh::Vec3f         axis;

	OpenMesh::Utils::Timer timer; timer.start();

	for (i=0, axis=OpenMesh::Vec3f(1,0,0); i<frames; ++i)
	{ rotate(axis, angle); display(); }
	for (i=0, axis=OpenMesh::Vec3f(0,1,0); i<frames; ++i)
	{ rotate(axis, angle); display(); }
	for (i=0, axis=OpenMesh::Vec3f(0,0,1); i<frames; ++i)
	{ rotate(axis, angle); display(); }

	glFinish();

	timer.stop();
	fps = (1000.0 / timer.mseconds() * (3.0 * frames));

	glPopMatrix();
	//glutPostRedisplay();// fltk�汾���Լ������Fl_Gl_Window��redraw������,���Բ���Glut�汾�µ�glutPostRedisplay������

	return fps;
}


//=============================================================================
