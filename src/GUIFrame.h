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
// 这个类已经完整, 用户和其它类不用考虑过多的.
// 它并不实现具体的功能, 它的功能都是调用其上层的GUIFrame *app来实现的.
// GUIFrame: basic gui used fltk, contain the menu, the canvas and the status bar.
class GUIFrame;

static OpenMesh::Vec3f black_color_glframe(0.0f, 0.0f, 0.0f);
static OpenMesh::Vec3f white_color_glframe(1.0f, 1.0f, 1.0f);

class GLCanvas : public Fl_Gl_Window
{
private:
    GUIFrame *app;//这个canvas所属的上层, 都是通过它来进行实际工作的
	int last_click[2];//鼠标上次的位置    

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
    int w_offset, h_offset;//整个窗口的尺寸和用来渲染的canvas窗口的尺寸的宽和高度差别
    
    Fl_Window *toplevel;
    Fl_Menu_Bar *menu_bar;
    //Fl_Menu_Item *menu_layout;// 等效于Fl_Menu_Item menu_layout[];也就是一个Fl_Menu_Item类型的数组
    GLCanvas *canvas;
    Fl_Output *status_line;//底下的状态信息栏 a one line text output field
    
	DecimationModel *mesh_model_;

public:
    // This is the public interface of GUIFrame available to the application.
    //
    static GUIFrame *current;	// There should only be one.

	GUIFrame(const char *title = "Basic GUI Frame", const int width = 800, const int height = 640);//640 * 480, 800 * 640
    virtual ~GUIFrame() {}//表示是有继承的 

	void add_mesh_model(DecimationModel *_mesh_model) { mesh_model_ = _mesh_model; mesh_model_->add_guiframe(this); };
	int run() {
		toplevel->show();
		return Fl::run();
	}

public:
    void title(const char *l) { toplevel->label(l); }//设置顶层窗口标题

	// Menu construction and standard callbacks
    int add_menu(const std::string&, int key, Fl_Callback *cb, int flags=0);//一般只填第一和第三项就可以, 分别是菜单名字和回调函数.
    int add_toggle_menu(const std::string&, int key, bool& val, int flags=0);
    static void cb_toggle(Fl_Menu_ *m, bool *flag);
 
    int status(const char *fmt, ...);//设置输出到底下状态信息栏的内容
    
	void lock_size();
    void unlock_size();
	void resize_canvas(int width, int height);

public:
    // 响应菜单, Callback functions that get executed in response to menu commands.
    static void cb_open(Fl_Widget *w, void *data);
	void open();
	static void cb_file_save_as(Fl_Widget *w, void *data); //另存为off或者是obj文件
	void save_as();
    static void cb_exit(Fl_Widget *w, void *data);//这个就不用改了
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
	static void cb_view_saliency_multiscales(Fl_Widget *w, void *data);// 将5个scales先做non-linear suppression, 再做adding.
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

	static void cb_local_smooth_parameterization(Fl_Widget *w, void *data); // 对简化的结果做一次中点采样.
	void local_smooth_parameterization_process() { mesh_model_->local_smooth_parameterization_process(); }
	static void cb_resample_edge_midpoint(Fl_Widget *w, void *data); // 对简化的结果做一次中点采样.
	void resample_edge_midpoint_process() { mesh_model_->resample_edge_midpoint_process(); }
	static void cb_midpointsubdivision_1(Fl_Widget *w, void *data); // 对简化的结果做一次中点细分,中点采样.
	void midpointsubdivision_1_process() { mesh_model_->midpointsubdivision_1_process(); } //
	static void cb_create_initialcontrolmesh(Fl_Widget *w, void *data); // 求initial control mesh
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
	// 这些函数是被GLCanvas所调用的, 用于将event传递给Model.
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

	bool has_init_display_;//这是fltk和glut不同所需要加入的.

public: 
	////////////////////////////////////////////////////////////////////////////////
	// =============================================================================
	// 为了实现鼠标响应和按键响应,我只是将fltk版本下需要我自己实现的鼠标响应函数转接到下面的glut版本下.
	// 下面这些实现鼠标(trackball功能)和按键响应的函数,都是在GlutViewer and GlutExaminer中拿过来的glut版本.

	int  width_, height_;//用于保存cavas的宽和高度

	void   set_scene(const OpenMesh::Vec3f& _center, float _radius);
	void   view_all();
	double measure_fps();

	// virtual函数,多态所需要的
	virtual void init_display();// 初始名字是init()的,但是这名字太常用了,很容易被子类的init函数错误的重载,所有换一个具体的名字
	virtual void draw();

	// overloaded glut functions
	virtual void display(void);//每次需要渲染的不同就在override这就可以了
	virtual void reshape(int w, int h); 
	virtual void motion(int x, int y);                       //鼠标移动响应
	virtual void mouse(int button, int state, int x, int y); //鼠标按键响应
	virtual void keyboard(int key, int x, int y);


protected:
	// 这些实现鼠标的缩放平移转动操作是不需要改的了, 除非想要增添选面选点的操作
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

protected://这保护属性表示其子类是可以使用这些信息的

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

	bool     button_down_[10];//其实只用到了前3项
	//上面这些实现鼠标(trackball功能)和按键响应的函数,都是在GlutViewer and GlutExaminer中拿过来的glut版本.
	// =============================================================================

private:
	OpenMesh::Vec3f background_color_;
};

////////////////////////////////////////////////////////////////////////
//
// This template makes it easier to create FLTK-compliant callbacks.
// In particular, its purpose is to construct static thunks for
// calling member functions of MxGUI-derived classes.
// 方便在fltk中为菜单项加入callback function
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
//		(gui->*cmd)();//没有参数
//		gui->canvas->redraw();
//	}
//
//	template<GuiCommand1 cmd, int i>
//	static void to_arg(Fl_Widget *, void *data)
//	{
//		GUIFrame *gui = static_cast<GUIFrame*>(data);
//		(gui->*cmd)(i);//一个参数
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