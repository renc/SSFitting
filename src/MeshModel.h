#ifndef dgpstudio_meshmodel_h
#define dgpstudio_meshmodel_h

//---------------------
// MeshModel
// Ŀ����ʵ��������CourseExamples�µ�MeshViewer.
// ����ڸ���GUIFrame, ��ӵĹ���: ʹ��OpenMesh����ģ��

//#include "../common/DGPcommon.h"
#include "DGPcommon.h"

//class MainWindow; // renc: qt version 
class GUIFrame; 

class MeshModel
{
public:
	enum DrawMode {
		NONE = 0, // ʲô������ʾ, ���Ǻ������ֻ�Ƿ�����Զ���.
		WIREFRAME = 1, HIDDENLINE = 2, SOLIDFLAT = 3, SOLIDSMOOTH = 4,	LINEANDSOLIDFLAT = 5,
		COLOR = 6
	};	
	
	MeshModel();
	virtual ~MeshModel() {}
	//void add_guiframe(MainWindow *_guiframe);renc qt version 
	void add_guiframe(GUIFrame *_guiframe);

public:
	//////////////////////////////////////////////////////////////////
	// ���涼�Ǵ�CourseExamples�е�MeshViewer�п���������,
	// virtual�����Ƕ�̬�������, ��ʾ����Щ�������ᱻ����overrided
	/// open mesh
	virtual bool open_mesh(const char* _filename);
	/// update buffer with face indices, Ϊ�˼���mesh_toberendered_��ָ����Ǹ��������Ⱦ.
	void update_renderface_indices(); // invoked by open_mesh(const char* _filename);in orther to fasten the render
	/// draw the scene
	virtual void draw();  //override�����еĺ���,ʵ����Ⱦģ���ļ�
	///////////////////////////////////////////////////////////////////

	void save_as(std::string _file_name){ OpenMesh::IO::write_mesh(*mesh_toberendered_, _file_name); }//mesh_

	void set_drawmodel(DrawMode _dm) { draw_mode_ = _dm; } ;
	
	// when the mesh is modified, like simplified and rifined, update the properties
	virtual void update();
	int n_vertices() { return mesh_.n_vertices(); }//�Ա�View/Controller��ѯ��.

	void set_draw_coordinate_system_process(bool _draw) { draw_coordiante_system_ = _draw; } //�ⲿ�ӿ�

protected:
	std::string filename_;
	//typedef OpenMesh::TriMesh_ArrayKernelT<>  TriMesh;
	TriMesh                       mesh_;
	TriMesh *mesh_toberendered_; // mesh_��ԭʼ����, �������Ҫ��Ⱦ������.
	std::vector<unsigned int>  indices_;//���ڷ��������Ⱦ�������int array

	bool is_modified;//���ڱ�����mesh_�ǲ��Ǳ��޸Ĺ���,���类�򻯹�,���类���ӹ���ȵ�.
	////////////////////////////////////////////////////////////////////////
	DGP::BoundingBox bb;

	//MainWindow *guiframe_;//��Ϊview and controller. renc: qt version
	GUIFrame *guiframe_; 

	DrawMode draw_mode_; 
	bool draw_coordiante_system_;
	OpenMesh::Vec3f line_color_;
};

#endif // dgpstudio_meshmodel_h