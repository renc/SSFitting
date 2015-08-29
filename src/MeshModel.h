#ifndef dgpstudio_meshmodel_h
#define dgpstudio_meshmodel_h

//---------------------
// MeshModel
// 目的是实现类似于CourseExamples下的MeshViewer.
// 相对于父类GUIFrame, 添加的功能: 使用OpenMesh载入模型

//#include "../common/DGPcommon.h"
#include "DGPcommon.h"

//class MainWindow; // renc: qt version 
class GUIFrame; 

class MeshModel
{
public:
	enum DrawMode {
		NONE = 0, // 什么都不显示, 这是很特殊的只是方便调试而已.
		WIREFRAME = 1, HIDDENLINE = 2, SOLIDFLAT = 3, SOLIDSMOOTH = 4,	LINEANDSOLIDFLAT = 5,
		COLOR = 6
	};	
	
	MeshModel();
	virtual ~MeshModel() {}
	//void add_guiframe(MainWindow *_guiframe);renc qt version 
	void add_guiframe(GUIFrame *_guiframe);

public:
	//////////////////////////////////////////////////////////////////
	// 下面都是从CourseExamples中的MeshViewer中拷贝过来的,
	// virtual函数是多态所必需的, 暗示着这些函数都会被子类overrided
	/// open mesh
	virtual bool open_mesh(const char* _filename);
	/// update buffer with face indices, 为了加速mesh_toberendered_所指向的那个网格的渲染.
	void update_renderface_indices(); // invoked by open_mesh(const char* _filename);in orther to fasten the render
	/// draw the scene
	virtual void draw();  //override父类中的函数,实现渲染模型文件
	///////////////////////////////////////////////////////////////////

	void save_as(std::string _file_name){ OpenMesh::IO::write_mesh(*mesh_toberendered_, _file_name); }//mesh_

	void set_drawmodel(DrawMode _dm) { draw_mode_ = _dm; } ;
	
	// when the mesh is modified, like simplified and rifined, update the properties
	virtual void update();
	int n_vertices() { return mesh_.n_vertices(); }//以便View/Controller的询问.

	void set_draw_coordinate_system_process(bool _draw) { draw_coordiante_system_ = _draw; } //外部接口

protected:
	std::string filename_;
	//typedef OpenMesh::TriMesh_ArrayKernelT<>  TriMesh;
	TriMesh                       mesh_;
	TriMesh *mesh_toberendered_; // mesh_是原始网格, 而这个是要渲染的网格.
	std::vector<unsigned int>  indices_;//用于方便快速渲染网格面的int array

	bool is_modified;//用于标记这个mesh_是不是被修改过了,例如被简化过,例如被增加过点等等.
	////////////////////////////////////////////////////////////////////////
	DGP::BoundingBox bb;

	//MainWindow *guiframe_;//作为view and controller. renc: qt version
	GUIFrame *guiframe_; 

	DrawMode draw_mode_; 
	bool draw_coordiante_system_;
	OpenMesh::Vec3f line_color_;
};

#endif // dgpstudio_meshmodel_h