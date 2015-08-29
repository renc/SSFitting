
// ***************************************************************
//  
//  DecimationModel 
//  Copyright (C) 2007 - by rencanjiang All Rights Reserved
//  rencanjiang@163.com	
//	-------------------------------------------------------------
//  
//  version: 1.0
//  data: 12/10/2007
// ***************************************************************
//  
//  File description:
// ***************************************************************

#ifndef dgpstudio_decimationviewer_h
#define dgpstudio_decimationviewer_h

// 2007-11-06
// 这个类主要是做简化模型网格的,初始来自于CourseExamples中的06-Decimation例子

//== INCLUDES =================================================================

//#include "../common/DGPcommon.h"
#include "DGPcommon.h"
#include <iostream>
#include <fstream>
#include <set>
#include <float.h> // for FLT_MAX
#include <cstring>
#include <string>
#include <map>
#include <algorithm> // for max_element algorithm
#include <functional> // for greater function object in std::set
#include "MeshModel.h" 

//#include "../QuarticBezierTriangle/QuadraticBezierTriangles.h"// renc, 20150829
#include "QuadraticBezierTriangles.h" // renc, 20150829

class DecimationModel: public MeshModel
{
	// constructor and deconstructor
public:
	DecimationModel();
	
	// override functions
public:		
	virtual void update(); 
	virtual void draw(); // 继承

public:
	
	static DecimationModel *currentDV;//static变量来.cpp文件中需要声明.
	
	// 存储一些和某个VertexHandle相关的信息
	OpenMesh::VPropHandleT<DGP::Quadricd>                vp_quadric_;
	OpenMesh::VPropHandleT<double>                   vp_cost_;//存放收缩的代价, 原来名字是vp_prio_.
	OpenMesh::VPropHandleT<TriMesh::HalfedgeHandle>    vp_target_;//存放和这点相关的最优点对(万一收缩时代价最小)

	// access priority of vertex _vh
	double& priority(TriMesh::VertexHandle _vh) { return mesh_.property(vp_cost_, _vh); }

	struct VertexCmp { //VertexCmp function object
		bool operator()(TriMesh::VertexHandle _v0, TriMesh::VertexHandle _v1) const
		{
			// std::set needs UNIQUE keys -> handle equal priorities
			return ((DecimationModel::currentDV->priority(_v0) == DecimationModel::currentDV->priority(_v1)) ? 
				(_v0.idx() < _v1.idx()) : // 优先权相等怎么办, 则以在数组中的索引判断根据
				(DecimationModel::currentDV->priority(_v0) < DecimationModel::currentDV->priority(_v1)));//_v0比_v1优先返回
		} // 代价越小优先级越高. 
	};
	std::set<TriMesh::VertexHandle, VertexCmp>  queue; //less(返回真)的结果是递增序列
	//队列, 以代价递增, 里面放的只是顶点handle, 而顶点的代价等属性放在其property中 

	// 默认的halfedge collapse是(from, to)->to, 和edge collapse是有一点区别的,
	bool vertex_optimal_placement;//是否使用在collapse后计算最优的点的位置
	bool area_weighted_error;//是否使用area weighted error matric

	OpenMesh::VPropHandleT<DGP::v_feature_type> vp_type_; // enum这个关键字在这里可以不用了, 这与C有区别

	TriMesh simplified_mesh_; //
	OpenMesh::VPropHandleT<DGP::v_feature_type> vp_type_sm_;//复杂一个模型的点的类型,
	OpenMesh::VPropHandleT<int > vp_node_handle_sm_;//复制与vertex hierarchy的联系.
	void add_required_properties();
	void remove_required_properties();
	void identify_all_sharp_features();  //准备去掉的.
	double priority_qem(TriMesh::HalfedgeHandle _heh);
	void enqueue_vertex(TriMesh::VertexHandle _vh);
	void decimate(unsigned int _n_vertices);

	OpenMesh::FPropHandleT<TriMesh::FHandle> fp_fh_; //mesh_.property(), mesh_的一个面对应于simplified_mesh_上的哪一个面.
	OpenMesh::FPropHandleT<TriMesh::FHandle> fp_fh_sm_; //simplified_mesh_.property(), simplified_mesh_的一个面对应于mesh_上的哪一个面.
	
public:
	// 响应菜单
	void decimation_process(unsigned int num, bool _feature_or_not = false);
	
public:
	// 响应菜单
	void sequence_refine(); //这两个refinement函数暂时只是用来说明有这个功能而已,还没有实际用途.
	void selective_refine();

private:
	DGP::VHierarchy vhierarchy_;
	OpenMesh::HPropHandleT<TriMesh::VertexHandle> heovh_; // 用于保存每一个halfedge's opposite vertex handle
	OpenMesh::VPropHandleT<int> vp_node_handle_; //mesh_中顶点找vhierarchy上的node.
	OpenMesh::VPropHandleT<int> vp_leaf_node_handle_; // point to the leaf node, never change
	std::vector<TriMesh::Point> vpoints_;
	std::vector<TriMesh::VHandle> voldhandles_; //和vpoints_对应的old vertex handle.记录的是garbage collection前的索引
	std::vector<int> ecol_sequence_list_;

public:
	void draw_mesh2_process() { draw_mesh2_ = !draw_mesh2_; }//开关显示原始网格, 用于参考.

private:
	TriMesh mesh2_; //复制original mesh, 不能对mesh2_有任何的修改. 只为了在deleted一些vertices之后还需要查询original的connectivity.
	OpenMesh::VPropHandleT<double> varea_;
	OpenMesh::HPropHandleT<double> hmvc_;// 每一个半边对应的mean value coordinate.
	TriMesh::VertexHandle test_vertex_, test_vertex1_, test_vertex2_;
	bool draw_mesh2_;
	std::vector<unsigned int>  mesh2_indices_for_render_;//用于方便快速渲染网格面的int array.
	void calc_backupmesh_vertex_area();
	void calc_backupmesh_halfedge_meanvaluecoor();// 计算每一个半边对应的mean value coordinate.
	void copy_backup_mesh();

public:
	// mesh_toberedered_在MeshModel中定义, 当读入一个新网格时就指向mesh_
	enum MeshTypeToBeRendered { ORIGINAL =0 , SIMPLIFIED, RESAMPLE_MIDPOINT, REFINED_SIMPLIFIED, 
		INITIAL_CONTROL, LOOP_ONE, LIMIT };
	void set_mesh_toberendered(TriMesh *_mesh_tobe_rendered, const MeshTypeToBeRendered _meshtype_tobe_rendered) { 
		_mesh_tobe_rendered->update_normals();
		mesh_toberendered_ = _mesh_tobe_rendered; 
		meshtype_toberendered_ = _meshtype_tobe_rendered; 
		update_renderface_indices();
	}
	
private:
	MeshTypeToBeRendered meshtype_toberendered_; //当面需要render的是哪一个mesh

	// -------------------------------------------------
private:
	bool do_parameterization_;//简化的同时是否进行参数化
	bool do_smooth_;                                   // property of mesh_
	OpenMesh::FPropHandleT<std::vector<TriMesh::VertexHandle> > fvset_;//保存了参数化到这个面上面的顶点集合
	OpenMesh::VPropHandleT<TriMesh::FaceHandle> vf_; // vf_(vf0_, vf1_, vf2_)
	OpenMesh::VPropHandleT<TriMesh::VertexHandle> vf0_;// the deleted vertex is parameterized onto this face handle(vf0_, vf1_, vf2_)
	OpenMesh::VPropHandleT<TriMesh::VertexHandle> vf1_;// the deleted vertex is parameterized onto this face handle(vf0_, vf1_, vf2_)
	OpenMesh::VPropHandleT<TriMesh::VertexHandle> vf2_;// the deleted vertex is parameterized onto this face handle(vf0_, vf1_, vf2_)
	OpenMesh::VPropHandleT<OpenMesh::Vec3d> vbc_;// the deleted vertex's barycentric coordinates of the according face
	OpenMesh::VPropHandleT<OpenMesh::Vec2d> vuv_;//the deleted vertex's parameter value uv, 这是在参数化过程中使用的.

	OpenMesh::HPropHandleT<std::vector<TriMesh::HalfedgeHandle> > hep_heset_;//保存参数化到这个*特征*半边上的特征原始的特征半边.
	OpenMesh::VPropHandleT<TriMesh::EdgeHandle> vp_feature_edge_;//特征顶点参数化到哪一个特征边上. 这一对信息应该是一起修改的.

	void initial_parameterization(TriMesh::HalfedgeHandle _hh);
	size_t n_local_sp_; // 在initial parameterize之后进行多少次的local smooth parameterize
	void local_smooth_parameterization_triangledomain();//每次以四个等边三角形作为平滑区域
	void local_smooth_parameterization_triangledomain(TriMesh::FaceHandle _fh, int &, int &);//每次以四个等边三角形作为平滑区域
	void local_smooth_parameterization_circledomain();//以一个圆(circle domain)做平滑区域
	void local_smooth_parameterization_circledomain(TriMesh::VertexHandle _vh);
	void resample_original_feature_edge_midpoint();//对原始边和特征边采样
	void resample_boundary_edge_midpoint(TriMesh &_mesh, TriMesh::HHandle _hh, TriMesh::VHandle _vstart, TriMesh::VHandle _vend);//对边界边采样
	void flatten_vertex_one_ring_resample_edges_midpoint();
	void flatten_vertex_one_ring_resample_edges_midpoint(TriMesh::VertexHandle vh);
	void flatten_face_resample_3edge_midpoint();
	void flatten_face_resample_3edge_midpoint(TriMesh::FaceHandle _fh);
	void project_midpoint_to_triangle();
	void use_closepoints(TriMesh::EHandle _eh);
  

	OpenMesh::EPropHandleT<bool > empl_;//  the mid-point of edge is located? 这个边的中点确定了没有?
	OpenMesh::EPropHandleT<TriMesh::Point> emp_;// the mid-point of edges, which not be deleted.
	OpenMesh::EPropHandleT<TriMesh::VertexHandle> emp_closest_vh_;//这个边中点靠近那个原始顶点.
	bool midpoint_sample_succeed_;//对于中点细分增加的新点是否在原始网格中找到了采样点. 
	TriMesh mesh_collapsed_;//主要是用于显示参数化效果, 和原始网格的连接关系一样, 只是点坐标不一样.
	bool mesh_collapsed_render_or_not_;	
	
private:
	TriMesh refined_simplified_mesh_;//对base mesh做一次midpoint subdivision获得的refined mesh
	OpenMesh::VPropHandleT<DGP::v_feature_type> vp_type_rsm_; 
	void split_edge_rsm(TriMesh::EHandle& );
	void split_face_rsm(TriMesh::FHandle& );
	TriMesh::FHandle corner_cutting_rsm(TriMesh::HHandle& );//返回值是我修改的.
	
	OpenMesh::VPropHandleT<TriMesh::FaceHandle> vf_of_rsm_; // mesh_上的点的参数化被映射到refined_simplified_mesh_上
	OpenMesh::VPropHandleT<TriMesh::VertexHandle> vf0_of_rsm_;//  
	OpenMesh::VPropHandleT<TriMesh::VertexHandle> vf1_of_rsm_;//  
	OpenMesh::VPropHandleT<TriMesh::VertexHandle> vf2_of_rsm_;//  
	OpenMesh::VPropHandleT<OpenMesh::Vec3d> vbc_of_rsm_;//  
	OpenMesh::FPropHandleT<std::vector<TriMesh::VHandle> > fvset_rsm_; //refined_simplified_mesh_上每一个面有什么原始点被参数化到这个面上.

	std::vector<TriMesh::VHandle> Xm_vertex_;//前n个元素对应simplified mesh的顶点, 后面的元素是refined simplified mesh中新增的顶点
	std::vector<TriMesh::EHandle> es_vec_;//存放simplified mesh中的所有边, 其实只是为了方便求Xm_vertex_和Bmn
	std::map<TriMesh::EHandle, OpenMesh::Vec3d> crease_edge_2_mid_pos_; //simplified mesh上特征边的采样中点坐标. 
	int n_corv_of_sm_; //simplified mesh(因为角度不发生变化, 因此也同时是refined simplified mesh)中的corner vertex的数目

	TriMesh initial_control_mesh_;
	OpenMesh::VPropHandleT<DGP::v_feature_type> vp_type_icm_; 
	void calc_AmmBmn(std::vector<std::map<int, double> > &Amm_data, std::vector<std::map<int, double> > &Bmn_data);
	void calc_crease_first(int m, int n);
	
public:
	bool do_parameterization() { return do_parameterization_; }
	void set_do_parameterization(bool _p) { do_parameterization_ = _p; }
	void set_n_local_sp(size_t _n) { n_local_sp_ = _n; }
	// 响应菜单函数
	void draw_mesh_collapsed__process() { mesh_collapsed_render_or_not_ = !mesh_collapsed_render_or_not_; }
	void local_smooth_parameterization_process();
	void resample_edge_midpoint_process();
	void midpointsubdivision_1_process();//对简化的模型进行一次中点细分, 中点采样.
	void create_initialcontrolmesh_process();//setup up the fitting equation and 求解出initial control mesh
	void create_initialcontrolmesh_process2();//setup up the fitting equation and 求解出initial control mesh
	void create_initialcontrolmesh_process3();
	void loop_subdivision_process();//对initial control mesh做一次带feature的loop subdivision
	void create_limit_surface_process();
	bool fitting_process();
    
	bool test9() {
		open_mesh(std::string("D:/cjren/Libraries/MeshCourse07-Code/MSVC/release/fan_t_6475.off").c_str());
		decimation_process(55);
		resample_edge_midpoint_process();
		midpointsubdivision_1_process();
		create_initialcontrolmesh_process2();
		return true;
	}
	//------------------------------------------------------
public:
	void evaluation_process();
	void view_fitting_error_mesh2_process();
	void view_fitting_error_sm_process();
	void refine_process();

	static double calc_step_side1(double _h, void *params); // stam(v,w) step side.
	static OpenMesh::Vec4d calc_step_side1_parames_; static OpenMesh::Vec3d calc_step_side1_sk_;
	static std::vector<OpenMesh::Vec3d> calc_step_side1_c_;
	Stam98ClassicLoopEval & get_cloop_evaluation() { return cloop_evaluation; }

	const double mesh_v_fe(TriMesh::VHandle _v) const { return mesh_.property(vp_fitting_error_, _v); }
	struct FErr_v 
	{
		bool operator() (TriMesh::VHandle v0, TriMesh::VHandle v1) {
			return (DecimationModel::currentDV->mesh_v_fe(v0) < DecimationModel::currentDV->mesh_v_fe(v1));
		}
	};
	std::vector<OpenMesh::Vec3d > test_array_;
	std::vector<TriMesh::VHandle> vsplit_array_;
	const double simplified_mesh_v_fe(TriMesh::VHandle _v) const { return simplified_mesh_.property(vp_fitting_error_sm_, _v); }
	struct FErr_v_sm 
	{
		bool operator() (TriMesh::VHandle v0, TriMesh::VHandle v1) {
			return (DecimationModel::currentDV->simplified_mesh_v_fe(v0) < DecimationModel::currentDV->simplified_mesh_v_fe(v1));
		}
	};
private:
	bool evaluation_get_all_footpoints();//将evaluation的结果作为下一步找最近点的初始值.
	bool evaluation_get_footpoint(const TriMesh& _mesh, TriMesh::FaceHandle _fh_of_initial_control_mesh_, 
		OpenMesh::Vec3d Sk, OpenMesh::Vec3d p, OpenMesh::Vec3d bc,
									OpenMesh::Vec3d &pp, OpenMesh::Vec3d &clp, const int cross_patch_depth);
	void evaluation_get_face_type(TriMesh::FaceHandle fh, bool &_is_feature, bool &_is_regular, int &feature_type, TriMesh::HalfedgeHandle &hh);
	Stam98ClassicLoopEval cloop_evaluation;
	Hoppe94LoopEval h94loop_evaluation;
	const static int Max_Cross_Patch_Depth = 5;
	OpenMesh::VPropHandleT<bool> vp_eval_;//是否求出来了// mesh_.property
	OpenMesh::VPropHandleT<OpenMesh::Vec3d> vp_eval_pos_;//mesh_上原始点在limit surface上的evaluation position
	OpenMesh::VPropHandleT<bool> vp_closest_;//是否求出来了 mesh_.property
	OpenMesh::VPropHandleT<OpenMesh::Vec3d> vp_closest_pos_;//mesh_上原始点在limit surface上的最近点position
	OpenMesh::VPropHandleT<bool> vp_project_;// whether get to Sk's projection point on the limit surfaces.
	OpenMesh::VPropHandleT<OpenMesh::Vec3d> vp_project_pos_;// Sk's projection point on the limit surfaces.
	void solve_func_get_h(const OpenMesh::Vec3d dStam_dv, const OpenMesh::Vec3d dStam_dw, const OpenMesh::Vec3d Sk, const OpenMesh::Vec3d clp, double & delta_v, double & delta_w);
	bool cloop_stam98_get_closest_points(const std::vector<OpenMesh::Vec3d> & cN6, 
										const double v, const double w, OpenMesh::Vec3d Sk, OpenMesh::Vec3d & clp, 
										const TriMesh& _mesh, TriMesh::HHandle h12, int cross_patch_depth); //h12是用来跨界时候找邻面的.
	bool h94loop_get_closest_points(Hoppe94LoopEval &h94loop_evaluation, bool _is_regular, int feature_type, 
									const std::vector<OpenMesh::Vec3d> & _c, const std::vector<bool> &_h,
									const double v, const double w, OpenMesh::Vec3d Sk, OpenMesh::Vec3d & clp, 
									const TriMesh& _mesh, TriMesh::HHandle h84orh12, int cross_patch_depth);
	void update_uvw_when_v_new_is_less_than_zero(const double v, const double delta_v, const double w, const double delta_w, double &v_new, double &w_new);
	void update_uvw_when_w_new_is_less_than_zero(const double v, const double delta_v, const double w, const double delta_w, double &v_new, double &w_new);
	void update_uvw_when_v_new_plus_w_new_is_greaterr_than_one(const double v, const double delta_v, const double w, const double delta_w, double &v_new, double &w_new);
	OpenMesh::VPropHandleT<double> vp_fitting_error_;//mesh_.property
	double standard_error_;
	bool is_view_fitting_error_of_mesh2_;
	OpenMesh::VPropHandleT<TriMesh::Color> vp_color_fitting_error_mesh2_;//mesh2_.property 
	void traverse_binary_tree_count_leaves(int _node_handle, int &_n_leaves);
	void traverse_binary_tree_collect_leaves(int _node_handle, std::vector<int> &leaf_node_handle_array);
	double traverse_binary_tree_add_leaves_fitting_error(int _node_handle);
	void calc_fitting_error();  
	bool is_view_fitting_error_of_sm_;//
	OpenMesh::VPropHandleT<double> vp_fitting_error_sm_;
	OpenMesh::VPropHandleT<TriMesh::Color> vp_color_fitting_error_sm_;//initial_control_mesh_.property 
	double standard_error_sm_;
	std::vector<TriMesh::VHandle> vsplit_array_sm_; // the vertices of simplified_mesh_ are going to be split.
	int count_split;
	
	// 先用QAS曲面逼近拟合曲面再求原始模型上采样点Sk到QAS曲面上的最近点(closest point or footpoint).
	bool evaluation_get_footpoint_byQAS(const TriMesh& _mesh, TriMesh::FaceHandle _fh_of_initial_control_mesh_, 
		OpenMesh::Vec3d Sk,  OpenMesh::Vec3d bc,
		OpenMesh::Vec3d &clp, int cross_patch_depth);
	OpenMesh::FPropHandleT<TriMesh::HalfedgeHandle> fp_hh_icm_; // icm上每一个面记录一个边.
	DGP::QuadraticBezierTriangles *qas_;

	std::set<double, std::greater<double> > appro_error_queue_;// queue of approximation error, 从大到小递减排序
	//std::set<TriMesh::VertexHandle, Appro_Error_Comp> appro_error_queue_;
	bool test_;
	
	
public:
	void error_driven_selective_refine_process();		

private:
	OpenMesh::VPropHandleT<TriMesh::VertexHandle> vp_closestp_icm_;//sample points at limit surface <---> foot points at original surface
	OpenMesh::VPropHandleT<double> vp_distances_;//distances = sample point - foot point
	void find_closestp_icm();
	bool test_is_leaf_of_vertex_hierarchy(TriMesh::VHandle _vh_in_simplified_mesh);
	bool test_is_triangle_flipping_after_vsplit(TriMesh::VHandle _vh_in_simplified_mesh);//在error driven selective refinement中预防生成的三角形出现翻转.
	bool vsplit_refine(TriMesh::VertexHandle v_to);
	bool error_driven_selective_refine(); //找一个可以分裂的顶点, 分裂.


	std::vector<OpenMesh::Vec3d> laplacian_;
	std::vector<OpenMesh::Vec3d> laplacian_rcm_;

	bool fitting_terminate_; //不断迭代, 一直到什么时候终止呢. 2009-07-07
	bool feature_considered_;//
public:
	void pick(float _p0, float _p1, float _p2);
private: // --------------------
	TriMesh controlmesh_forbounding_;//对initial control mesh经过几次loop subd后获得的.
	TriMesh upperbounding_mesh_, lowerbounding_mesh_;// 
public:
	void try_bounding_process();
};
/**/
#endif // dgpstudio_decimationviewer_h