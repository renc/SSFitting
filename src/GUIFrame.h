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
#ifndef dgpstudio_guiframe_h
#define dgpstudio_guiframe_h

#include <windows.h>
// fltk 1.3 contain gl/glu/glut, using some freeglut code. check the doc. 
#include <FL/gl.h>//#include <GL/gl.h>
#include <FL/glu.h>
#include <fl/glut.h> // // for glutSolidSphere glutWireTeapot 
//#pragma comment(lib, "glut32.lib") // renc 20150829
#pragma comment (lib, "opengl32.lib")
#pragma comment (lib, "glu32.lib")
#include <FL/Fl.H>
#include <FL/fl_file_chooser.H>
#include <FL/filename.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/Fl_Output.H>
#include <fl/fl_ask.h> // for fl_alert(...)

#include <string>
#include <vector>

//#include <OpenMesh/Core/Math/VectorT.hh> // OpenMesh 1.x
#include <OpenMesh/Core/Geometry/VectorT.hh> // OpenMesh 2.0 RC3


// -----------------------------------------
// GLCanvas:
// ������Ѿ�����, �û��������಻�ÿ��ǹ����.
// ������ʵ�־���Ĺ���, ���Ĺ��ܶ��ǵ������ϲ��GUIFrame *app��ʵ�ֵ�.
// GUIFrame: basic gui used fltk, contain the menu, the canvas and the status bar.
class GUIFrame;

static OpenMesh::Vec3f black_color_glframe(0.0f, 0.0f, 0.0f);
static OpenMesh::Vec3f white_color_glframe(1.0f, 1.0f, 1.0f);

class GLCanvas : public Fl_Gl_Window
{
private:
    GUIFrame *app;//���canvas�������ϲ�, ����ͨ����������ʵ�ʹ�����
	int last_click[2];//����ϴε�λ��    

public:
    // Override selected FLTK window methods
    //
    virtual void draw(); //in Fl_Window class, draw() is proteted.
    virtual int handle(int event);//in Fl_Window class, handle() is also public.
    virtual void resize(int x, int y, int w, int h);//in Fl_Window class, resize() is alse public.

public:
    GLCanvas(int x, int y, int w, int h, const char *label=NULL)
		: Fl_Gl_Window(x, y, w, h, label), app(NULL){ last_click[0] = last_click[1] = -1; }
    void attach_app(GUIFrame *a){ if(!app)	app = a; }
};

// --------------------------------------------------
//#include "MeshModel.h"
#include "DecimationModel.h"

class GUIFrame
{
private:
    int w_offset, h_offset;//�������ڵĳߴ��������Ⱦ��canvas���ڵĳߴ�Ŀ�͸߶Ȳ��
    
    Fl_Window *toplevel;
    Fl_Menu_Bar *menu_bar;
    //Fl_Menu_Item *menu_layout;// ��Ч��Fl_Menu_Item menu_layout[];Ҳ����һ��Fl_Menu_Item���͵�����
    GLCanvas *canvas;
    Fl_Output *status_line;//���µ�״̬��Ϣ�� a one line text output field
    
	DecimationModel *mesh_model_;

public:
    // This is the public interface of GUIFrame available to the application.
    //
    static GUIFrame *current;	// There should only be one.

	GUIFrame(const char *title = "Basic GUI Frame", const int width = 800, const int height = 640);//640 * 480, 800 * 640
    virtual ~GUIFrame() {}//��ʾ���м̳е� 

	void add_mesh_model(DecimationModel *_mesh_model) { mesh_model_ = _mesh_model; mesh_model_->add_guiframe(this); };
	int run() {
		toplevel->show();
		return Fl::run();
	}

public:
    void title(const char *l) { toplevel->label(l); }//���ö��㴰�ڱ���

	// Menu construction and standard callbacks
    int add_menu(const std::string&, int key, Fl_Callback *cb, int flags=0);//һ��ֻ���һ�͵�����Ϳ���, �ֱ��ǲ˵����ֺͻص�����.
    int add_toggle_menu(const std::string&, int key, bool& val, int flags=0);
    static void cb_toggle(Fl_Menu_ *m, bool *flag);
 
    int status(const char *fmt, ...);//�������������״̬��Ϣ��������
    
	void lock_size();
    void unlock_size();
	void resize_canvas(int width, int height);

public:
    // ��Ӧ�˵�, Callback functions that get executed in response to menu commands.
    static void cb_open(Fl_Widget *w, void *data);
	void open();
	static void cb_file_save_as(Fl_Widget *w, void *data); //���Ϊoff������obj�ļ�
	void save_as();
    static void cb_exit(Fl_Widget *w, void *data);//����Ͳ��ø���
	void exitt();

	static void cb_view_wireframe(Fl_Widget *w, void *data);
	static void cb_view_hiddenline(Fl_Widget *w, void *data);
	static void cb_view_solidflat(Fl_Widget *w, void *data);
	static void cb_view_solidsmooth(Fl_Widget *w, void *data);
	static void cb_view_lineandsolidflat(Fl_Widget *w, void *data);
	void set_view_drawmodel(MeshModel::DrawMode _dm) { mesh_model_->set_drawmodel(_dm); }
	static void cb_draw_coordinate_system_true(Fl_Widget *w, void *data);
	static void cb_draw_coordinate_system_false(Fl_Widget *w, void *data);
	void set_draw_coordinate_system(bool _draw) { mesh_model_->set_draw_coordinate_system_process(_draw); }

/*
	static void cb_view_curvature(Fl_Widget *w, void *data); 
	void view_curvature() { mesh_model_->view_curvature(); }//
	static void cb_view_saliency_1(Fl_Widget *w, void *data);
	static void cb_view_saliency_2(Fl_Widget *w, void *data);
	static void cb_view_saliency_3(Fl_Widget *w, void *data);
	static void cb_view_saliency_4(Fl_Widget *w, void *data);
	static void cb_view_saliency_5(Fl_Widget *w, void *data);
	static void cb_view_saliency_multiscales(Fl_Widget *w, void *data);// ��5��scales����non-linear suppression, ����adding.
	void view_saliency(int _scale_type) { mesh_model_->view_saliency(_scale_type); }//

	static void cb_smoothing_laplaciansmoothing(Fl_Widget *w, void *data); 
	void smoothing_laplaciansmoothing() { mesh_model_->lalplacian_smoothing();} //
*/
	static void cb_fitting_decimation_percentage(Fl_Widget *w, void *data); //
	static void cb_fitting_decimation_vertices_num(Fl_Widget *w, void *data); //
	static void cb_draw_mesh_collapsed_or_not(Fl_Widget *w, void *data); //
	static void cb_fitting_sequence_refine(Fl_Widget *w, void *data); //
	static void cb_fitting_selective_refine(Fl_Widget *w, void *data);
	void decimation_process(unsigned int _num) { mesh_model_->decimation_process(_num); }
	void draw_mesh_collapsed_or_not_process() { mesh_model_->draw_mesh_collapsed__process(); }
	void sequence_refine() { mesh_model_->sequence_refine(); }
	void selective_refine() { mesh_model_->selective_refine(); }

	static void cb_local_smooth_parameterization(Fl_Widget *w, void *data); // �Լ򻯵Ľ����һ���е����.
	void local_smooth_parameterization_process() { mesh_model_->local_smooth_parameterization_process(); }
	static void cb_resample_edge_midpoint(Fl_Widget *w, void *data); // �Լ򻯵Ľ����һ���е����.
	void resample_edge_midpoint_process() { mesh_model_->resample_edge_midpoint_process(); }
	static void cb_midpointsubdivision_1(Fl_Widget *w, void *data); // �Լ򻯵Ľ����һ���е�ϸ��,�е����.
	void midpointsubdivision_1_process() { mesh_model_->midpointsubdivision_1_process(); } //
	static void cb_create_initialcontrolmesh(Fl_Widget *w, void *data); // ��initial control mesh
	void create_initialcontrolmesh_process() { mesh_model_->create_initialcontrolmesh_process(); } //
	static void cb_loop_subdivision(Fl_Widget *w, void *data); // 
	void loop_subdivision_process() { mesh_model_->loop_subdivision_process(); }
	static void cb_create_limit_surface(Fl_Widget *w, void *data);
	void create_limit_surface_process() { mesh_model_->create_limit_surface_process(); }
	static void cb_fitting(Fl_Widget *w, void *data);
	void fitting_process() { mesh_model_->fitting_process(); }
	static void cb_error_driven_selective_refine(Fl_Widget *w, void *data);
	void error_driven_selective_refine_process() {mesh_model_->error_driven_selective_refine_process(); }
	static void cb_evaluation(Fl_Widget *w, void *data);
	void evaluation_process() {mesh_model_->evaluation_process(); }
	static void cb_view_fitting_error(Fl_Widget *w, void *data);
	void view_fitting_error_mesh2_process() {mesh_model_->view_fitting_error_mesh2_process(); }
	static void cb_view_fitting_error2(Fl_Widget *w, void *data);
	void view_fitting_error_sm_process() {mesh_model_->view_fitting_error_sm_process(); }
	static void cb_refine(Fl_Widget *w, void *data);
	void refine_process() {mesh_model_->refine_process(); }

	static void cb_try_boudnding(Fl_Widget *w, void *data);
	void try_bounding_process() { mesh_model_->try_bounding_process(); }
	static void cb_draw_mesh2(Fl_Widget *w, void *data);
	void draw_mesh2_process() {mesh_model_->draw_mesh2_process(); }


	static void cb_help_about(Fl_Widget *w, void *data);

public:
	///////////////////////////////////////////////////////////////////////////////
	// ============================================================================
    // Applications are customized by overriding the following methods.
	// Override these methods to control the contents of the GL canvas.
	// ��Щ�����Ǳ�GLCanvas�����õ�, ���ڽ�event���ݸ�Model.
    virtual void setup_for_drawing();
    virtual void draw_contents();
    virtual void update_animation() {};

    // Override these methods to receive events from the GL canvas
    virtual bool mouse_down(int *where, int which);
    virtual bool mouse_up(int *where, int which);
    virtual bool mouse_drag(int *where, int *last, int which);
    virtual bool key_press(int key);

    // Override this method to free memory, close files, etc.
    virtual void cleanup_for_exit() {}

	bool has_init_display_;//����fltk��glut��ͬ����Ҫ�����.

public: 
	////////////////////////////////////////////////////////////////////////////////
	// =============================================================================
	// Ϊ��ʵ�������Ӧ�Ͱ�����Ӧ,��ֻ�ǽ�fltk�汾����Ҫ���Լ�ʵ�ֵ������Ӧ����ת�ӵ������glut�汾��.
	// ������Щʵ�����(trackball����)�Ͱ�����Ӧ�ĺ���,������GlutViewer and GlutExaminer���ù�����glut�汾.

	int  width_, height_;//���ڱ���cavas�Ŀ�͸߶�

	void   set_scene(const OpenMesh::Vec3f& _center, float _radius);
	void   view_all();
	double measure_fps();

	// virtual����,��̬����Ҫ��
	virtual void init_display();// ��ʼ������init()��,����������̫������,�����ױ������init�������������,���л�һ�����������
	virtual void draw();

	// overloaded glut functions
	virtual void display(void);//ÿ����Ҫ��Ⱦ�Ĳ�ͬ����override��Ϳ�����
	virtual void reshape(int w, int h); 
	virtual void motion(int x, int y);                       //����ƶ���Ӧ
	virtual void mouse(int button, int state, int x, int y); //��갴����Ӧ
	virtual void keyboard(int key, int x, int y);


protected:
	// ��Щʵ����������ƽ��ת�������ǲ���Ҫ�ĵ���, ������Ҫ����ѡ��ѡ��Ĳ���
	// updates projection matrix
	void update_projection_matrix();
	// translate the scene and update modelview matrix
	void translate(const OpenMesh::Vec3f& _trans);
	// rotate the scene (around its center) and update modelview matrix
	void rotate(const OpenMesh::Vec3f& _axis, float _angle);

	// virtual trackball: map 2D screen point to unit sphere
	bool map_to_sphere(const OpenMesh::Vec2i& _point, OpenMesh::Vec3f& _result);

	// mouse processing functions
	void rotation(int x, int y);
	void translation(int x, int y);
	void zoom(int x, int y);

public:
	OpenMesh::Vec3f center() { return center_; }
	float radius() { return radius_; }

protected://�Ᵽ�����Ա�ʾ�������ǿ���ʹ����Щ��Ϣ��

	// scene position and dimension
	OpenMesh::Vec3f    center_;
	float    radius_;

	// projection parameters
	float    near_, far_, fovy_;

	// OpenGL matrices// 0,1,2,3 is the first column, 0,4,8,12 is the first row
	double   projection_matrix_[16], modelview_matrix_[16];

	// trackball helpers
	OpenMesh::Vec2i    last_point_2D_;
	OpenMesh::Vec3f    last_point_3D_;
	bool     last_point_ok_;

	bool     button_down_[10];//��ʵֻ�õ���ǰ3��
	//������Щʵ�����(trackball����)�Ͱ�����Ӧ�ĺ���,������GlutViewer and GlutExaminer���ù�����glut�汾.
	// =============================================================================

private:
	OpenMesh::Vec3f background_color_;
};

////////////////////////////////////////////////////////////////////////
//
// This template makes it easier to create FLTK-compliant callbacks.
// In particular, its purpose is to construct static thunks for
// calling member functions of MxGUI-derived classes.
// ������fltk��Ϊ�˵������callback function
//template<class Gui>
//struct MxBinder
//{
//	typedef void (GUIFrame::*GuiCommand)();
//	typedef void (GUIFrame::*GuiCommand1)(int);
//	typedef void (GUIFrame::*GuiCommand2)(Fl_Menu_ *);
//
//	template<GuiCommand cmd>
//	static void to(Fl_Widget *, void *data)
//	{
//		GUIFrame *gui = static_cast<GUIFrame*>(data);
//		(gui->*cmd)();//û�в���
//		gui->canvas->redraw();
//	}
//
//	template<GuiCommand1 cmd, int i>
//	static void to_arg(Fl_Widget *, void *data)
//	{
//		GUIFrame *gui = static_cast<GUIFrame*>(data);
//		(gui->*cmd)(i);//һ������
//		gui->canvas->redraw();
//	}
//
//	template<GuiCommand2 cmd>
//	static void to_menu(Fl_Widget *w, void *data)
//	{
//		GUIFrame *gui = static_cast<GUIFrame*>(data);
//		(gui->*cmd)(static_cast<Fl_Menu_ *>(w));
//		gui->canvas->redraw();
//	}
//};

#endif //dgpstudio_guiframe_h