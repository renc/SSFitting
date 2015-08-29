#ifndef dgp_studio_common_trackball_h
#define dgp_studio_common_trackball_h

#include "OpenMeshAll.h"
#include "gl.h"


namespace DGP {

//Usage:
//	(1). init();
//	(2). set_scene(...); if you want another center and radius;
//	(3). update_reshape(...); when you change the window's width and height;
//	(4). mouse input and manipulation.
//  (5). 区别于在opengl draw函数中还需要将trackball的结果(变换矩阵)应用到opengl的modelview中,
//       这里是直接在trackball的函数中就应用了, 所有opengl draw函数中只需要画自己的东西就ok了.
//   具体的用法示例在class GlutWindow.

class TrackballManipulator {
	// ----------下面是来自GlutExaminer的trackball功能实现
public:
	TrackballManipulator() {
		// init mouse buttons
		for (int i=0; i<10; ++i)
			button_down_[i] = false;
	}

	// init scene pos and size
	void init(int _w, int _h);

	// After read in new mesh model. 
	void set_scene(const OpenMesh::Vec3f& _center, float _radius);

	// updates window's width and height, for example in void reshape(int _w, int _h);
	void update_reshape(int _w, int _h);
	
	// button 0:left, 1. middle, 2:right. 
	// state: 0: press, 1:release
	//  /* Mouse buttons. */
	//#define GLUT_LEFT_BUTTON		0
	//#define GLUT_MIDDLE_BUTTON	1
	//#define GLUT_RIGHT_BUTTON		2 
	//	/* Mouse button  state. */
	//#define GLUT_DOWN				0
	//#define GLUT_UP				1
	void mouse(int button, int state, int x, int y);  
	void motion(int x, int y);

	// 我加入的.
	// To transfrom one position of the mesh to position according to the world coordinate system.
	// 1. Here, the camera/viewpoint is always at origin.
	// 2. From translation and rotate functions, we can see how to computer the pos in world using the modelview_matrix_.
	OpenMesh::Vec3d object_to_eye_coordiante(const OpenMesh::Vec3d & _p);  
public:
	OpenMesh::Vec3f center() { return center_; }
	float radius() { return radius_; }

protected:	
	void view_all(); //只是被set_scene函数调用

	// updates projection matrix, 只被reshape函数和set_scene函数调用
	void update_projection_matrix();

	// virtual trackball: map 2D screen point to unit sphere, invoked by rotation
	bool map_to_sphere(const OpenMesh::Vec2i& _point, OpenMesh::Vec3f& _result);

	// 这些实现鼠标的缩放平移转动操作是不需要改的了, 除非想要增添选面选点的操作
	// mouse processing functions
	void rotation(int x, int y);
	void translation(int x, int y);
	void zoom(int x, int y);
	// translate the scene and update modelview matrix, invoked by rotation(x,y).
	void translate(const OpenMesh::Vec3f& _trans);
	// rotate the scene (around its center) and update modelview matrix, invoked by translation/zoom(x,y).
	void rotate(const OpenMesh::Vec3f& _axis, float _angle);

protected://这保护属性表示其子类是可以使用这些信息的

	// scene position and dimension
	OpenMesh::Vec3f    center_;
	float    radius_;

	// projection parameters
	float    near_, far_, fovy_;

	// OpenGL matrices// 0,1,2,3 is the first column, 0,4,8,12 is the first row
	double   projection_matrix_[16], modelview_matrix_[16];

	int width_, height_;
	// trackball helpers
	OpenMesh::Vec2i    last_point_2D_;
	OpenMesh::Vec3f    last_point_3D_;
	bool     last_point_ok_;

	bool     button_down_[10];//其实只用到了前3项
	//上面这些实现鼠标(trackball功能)和按键响应的函数,都是在GlutViewer and GlutExaminer中拿过来的glut版本.
	// =============================================================================

}; // end of class TrackballManipulator

} // end of namespace DGP
#endif //dgp_studio_common_trackball_h