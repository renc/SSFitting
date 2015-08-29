// ***************************************************************
//  
//  MeshType 
//  Copyright (C) 2007 - by rencanjiang All Rights Reserved
//  rencanjiang@163.com	
//	-------------------------------------------------------------
//  
//  version: 1.0
//  data: 12/16/2007
// ***************************************************************
//  
//  File description:
//  ���������ʹ�õ�Mesh������
// ***************************************************************

#ifndef dgp_meshtype_h
#define dgp_meshtype_h
#include "OpenMeshAll.h"

// specify a mesh //������������������
struct MyTraits : public OpenMesh::DefaultTraits
{
	typedef double Scalar; //����������ʹ�õ�����, ������float��double//

	// ---------------
	typedef OpenMesh::Vec3d  Point;//��Point�����ͺ�Normal�������ǿ��Զ��ĳ�Vec3d��

	/// The default normal type is OpenMesh::Vec3f. Normal��Point������Ӧ�ö�����ͬ��
	typedef OpenMesh::Vec3d  Normal;

	/// The default 1D texture coordinate type is float.
	typedef double  TexCoord1D;
	/// The default 2D texture coordinate type is OpenMesh::Vec2f.
	typedef OpenMesh::Vec2d  TexCoord2D;
	/// The default 3D texture coordinate type is OpenMesh::Vec3f.
	typedef OpenMesh::Vec3d  TexCoord3D;

	// ----------------
	// �ⶫ��Ӧ��ȥ����, ͳһ��runtime property��
	// �������runtimeʱ��mesh_.add_property(...)������standard properties,
	// ��ʱֻʣ��Decimation project�е�����.
	// ����ʹ�õ�����compile time����custom property�ķ�ʽ
	VertexTraits //��compile time�ı�mesh items�Ľṹ
	{
	public:
		int node_handle()
		{ return node_handle_; }

		void set_node_handle(int _node_handle)
		{ node_handle_ = _node_handle; }

	private:
		//�����Ӧ�Ľڵ�Node��handle,��ʵ�����Ǹ��ڵ���vertex hierarchy��nodes�����е��±�.
		int  node_handle_;//���ڱ���ÿһ����������Ӧ����vertex hierarchy�еĵ�ǰactive node. 

	};//ʹ��mesh_.vertex(vh).node_handle()������ʹ��.

	// ----------------
	// �������request_..._...()��runtimeʱ������standard properties, ����ʹ����compile time����ķ�ʽ
	// ��Ϊ��Щstatus̫������.
	VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color); // ��compile time �ﵽmesh.request_vertex_status();��Ч��
	HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);//������Ĭ��ʹ�õ���
	EdgeAttributes(OpenMesh::Attributes::Status);
	FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
	// ��compile time �ﵽmesh.request_vertex_status();��Ч��	
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> Mesh;//��������������MeshViewer�е�,����Ϊ�˺����������Ǻ���������ʹ�õ�,���и���������
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> TriMesh;
//����ʹ��Mesh::Scalar, Mesh::Point ��Ϊһ��������, ��Ϊfloat or double���͵Ļ���

// ����������������
struct MyTraits1 : public OpenMesh::DefaultTraits //�������MyTraitsһ����, ֻ�Ǽ򵥰汾.
{
	typedef double Scalar; //����������ʹ�õ�����, ������float��double//

	typedef OpenMesh::Vec3d  Point;//��Point�����ͺ�Normal�������ǿ��Զ��ĳ�Vec3d��
	typedef OpenMesh::Vec3d  Normal;

	typedef double  TexCoord1D;
	typedef OpenMesh::Vec2d  TexCoord2D;
	typedef OpenMesh::Vec3d  TexCoord3D;

	VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color); // ��compile time �ﵽmesh.request_vertex_status();��Ч��
	HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);//������Ĭ��ʹ�õ���
	EdgeAttributes(OpenMesh::Attributes::Status);
	FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);	
};
typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits1> PolyMesh;

namespace DGP {
struct BoundingBox {
	OpenMesh::Vec3d bbMin;
	OpenMesh::Vec3d bbMax;
	OpenMesh::Vec3d center; //= (bb.bbMin + bb.bbMax)*0.5;
	double diagonal_lenght;// �Խ��߳���= (bb.bbMin - bb.bbMax).norm();
};

// ��������
template <typename MeshT>
inline void mesh_copy(MeshT & _dest_mesh, const MeshT & _sour_mesh) {
	//if (_dest_mesh.empty() == false) {
		_dest_mesh.clear();
	//}
	_dest_mesh = _sour_mesh;
}

template <typename TriMeshT> // ��valence = 6������·�����ʰ��.
inline typename TriMeshT::HalfedgeHandle forward3_halfedge_handle(const TriMeshT & _m, const typename TriMeshT::HHandle& _h) {
	return _m.next_halfedge_handle(_m.opposite_halfedge_handle(
		_m.next_halfedge_handle(_m.opposite_halfedge_handle(_m.next_halfedge_handle(_h)))));
}
// 
template <typename TriMeshT>
inline typename TriMeshT::HalfedgeHandle backward3_halfedge_handle(const TriMeshT & _m, const typename TriMeshT::HHandle& _h) {
	return _m.prev_halfedge_handle(_m.opposite_halfedge_handle(
		_m.prev_halfedge_handle(_m.opposite_halfedge_handle(_m.prev_halfedge_handle(_h)))));
}
// ÿһ���ߵ��е�
template <class MeshT> 
inline OpenMesh::Vec3d edge_midpos(const MeshT &_m, const typename MeshT::EHandle &_eh) {
	typename MeshT::HHandle h0(_m.halfedge_handle(_eh, 0)); 
	typename MeshT::HHandle h1(_m.halfedge_handle(_eh, 1));
	return (_m.point(_m.to_vertex_handle(h0)) + _m.point(_m.to_vertex_handle(h1)) )*0.5;
}
// ÿһ�����ж��ٸ��ǵ�,(������ڻ��İ�ߵĳ�����). 
template <class MeshT>
inline int n_face_corners(const MeshT& _m, const typename MeshT::FHandle& _fh) {
	int count = 0;
	for (typename MeshT::ConstFaceVertexIter cfv_it(_m, _fh); cfv_it; ++cfv_it) ++count;
	return count;
}
// ��һ������������.
template <class MeshT> 
inline typename MeshT::Point face_center(const MeshT &_m, const typename MeshT::FHandle &_fh) {
	int count = 0;
	typename MeshT::Point pos_sum(0, 0, 0);
	for (typename MeshT::ConstFaceVertexIter cfv_it(_m, _fh); cfv_it; ++cfv_it) {
		++count;
		pos_sum += _m.point(cfv_it.handle());
	}
	return pos_sum / count; 
}
// �ж�ģ�������Ƿ���б߽�
template <typename MeshT>
inline bool has_boundary(const MeshT &_m) {
	for (typename MeshT::ConstHalfedgeIter ch_it(_m.halfedges_begin()), ch_end(_m.halfedges_end());
		ch_it != ch_end; ++ch_it) {
			if (_m.is_boundary(ch_it.handle())) return true;
	}
	return false;
}

// ����ģ�����������, ֻ�������λ��ı��Σ����ı��λ�ϣ������(������������ϵ�).
enum MeshTypes { MPY_POLY, MPY_QUAD, MPY_TRI, MPY_QUAD_TRI }; 
template <typename MeshT>
inline MeshTypes identify_mesh_poly_type(const MeshT & _mesh, bool _debug = false) {
	int n_count_tris = 0, n_count_quads = 0, n_count_polys = 0;
	for (typename MeshT::ConstFaceIter cf_it(_mesh.faces_begin()), cf_end(_mesh.faces_end()); cf_it != cf_end; ++cf_it) {
		int count_vertices = 0;
		for (typename MeshT::CFVIter cfv_it(_mesh, cf_it); cfv_it; ++cfv_it) { // same as: cfv_it = _mesh.cfv_iter(cf_it.handle());
			 ++count_vertices;
		} 
		if (count_vertices == 3) ++n_count_tris;
		else if (count_vertices == 4) ++n_count_quads;
		else if (count_vertices > 4) ++n_count_polys;
		else std::cerr << "Error: ԭʼģ��������.\n";
	}
	if (n_count_polys != 0) {
		if (_debug ) std::cout << "Info: MeshTypes polygons.\n";
		return MPY_POLY;
	} else {
		if (n_count_quads == 0 && n_count_tris != 0) { 
			if (_debug) std::cout << "Info: MeshTypes triangles.\n"; return MPY_TRI; 
		} else if (n_count_quads != 0 && n_count_tris == 0) { 
			if (_debug) std::cout << "Info: MeshTypes quads.\n"; return MPY_QUAD; 
		} else if (n_count_quads != 0 && n_count_tris != 0) { 
			if (_debug) std::cout << "Info: MeshTypes tri+quads.\n"; return MPY_QUAD_TRI; 
		}
	}
}

// ����OpenMesh 1.9.6�е�PolyConnectivity���е�remove_edge(EdgeHandle _eh):
/** Removes the edge _eh. Its adjacent faces are merged. _eh and one of the 
adjacent faces are set deleted. The handle of the remaining face is 
returned (InvalidFaceHandle is returned if _eh is a boundary edge).

\precondition is_simple_link(_eh). This ensures that there are no hole faces
or isolated vertices appearing in result of the operation. 
����汾��û��is_simple_link�������.

\attention Needs the Attributes::Status attribute for edges and faces.

\note This function does not perform a garbage collection. It
only marks items as deleted. ����ǧ��ǵ���֮��Ҫ_mesh.garbage_collection();
*/
template <class PolyMeshT>
inline typename PolyMeshT::FaceHandle
PolyConnectivity_remove_edge(PolyMeshT &_m, typename PolyMeshT::EdgeHandle _eh) {
	//don't allow "dangling" vertices and edges
	assert(!_m.status(_eh).deleted());

	typename PolyMeshT::HalfedgeHandle heh0 = _m.halfedge_handle(_eh, 0);
	typename PolyMeshT::HalfedgeHandle heh1 = _m.halfedge_handle(_eh, 1);

	//deal with the faces //      remain                         delete
	typename PolyMeshT::FaceHandle rem_fh = _m.face_handle(heh0), del_fh = _m.face_handle(heh1);
	if (!del_fh.is_valid())
	{//boundary case - we must delete the rem_fh
		std::swap(del_fh, rem_fh);
	}
	assert(del_fh.is_valid());
	/*  for (FaceHalfedgeIter fh_it = fh_iter(del_fh); fh_it; ++fh_it)
	{//set the face handle of the halfedges of del_fh to point to rem_fh
	set_face_handle(fh_it, rem_fh);  
	} */
	//fix the halfedge relations
	typename PolyMeshT::HalfedgeHandle prev_heh0 = _m.prev_halfedge_handle(heh0);
	typename PolyMeshT::HalfedgeHandle prev_heh1 = _m.prev_halfedge_handle(heh1);

	typename PolyMeshT::HalfedgeHandle next_heh0 = _m.next_halfedge_handle(heh0);
	typename PolyMeshT::HalfedgeHandle next_heh1 = _m.next_halfedge_handle(heh1);

	_m.set_next_halfedge_handle(prev_heh0, next_heh1);
	_m.set_next_halfedge_handle(prev_heh1, next_heh0); 
	//correct outgoing vertex handles for the _eh vertices (if needed)
	typename PolyMeshT::VertexHandle vh0 = _m.to_vertex_handle(heh0);
	typename PolyMeshT::VertexHandle vh1 = _m.to_vertex_handle(heh1);

	if (_m.halfedge_handle(vh0) == heh1)
	{
		_m.set_halfedge_handle(vh0, next_heh0);
	}
	if (_m.halfedge_handle(vh1) == heh0)
	{
		_m.set_halfedge_handle(vh1, next_heh1);
	} 

	//correct the hafledge handle of rem_fh if needed and preserve its first vertex
	if (_m.halfedge_handle(rem_fh) == heh0)
	{//rem_fh is the face at heh0
		_m.set_halfedge_handle(rem_fh, prev_heh1);
	}
	else if (_m.halfedge_handle(rem_fh) == heh1)
	{//rem_fh is the face at heh1
		_m.set_halfedge_handle(rem_fh, prev_heh0);
	}  
	for (typename PolyMeshT::FaceHalfedgeIter fh_it = _m.fh_iter(rem_fh); fh_it; ++fh_it) //cjren 2009-01 ���rem_fh�Ҿ��ÿ�����Ҫ�ĳ�del_fh,��δ��֤��.
	{//set the face handle of the halfedges of del_fh to point to rem_fh
		_m.set_face_handle(fh_it, rem_fh);   
	}  

	_m.status(_eh).set_deleted(true);  
	_m.status(del_fh).set_deleted(true);  
	return rem_fh;//returns the remaining face handle
	// cjren 2009-01-05 Ϊʲô����_m.adjust_outgoing_halfedge(vh)����Ϊ���ﲻ�����޸��Ǹ��߽�ߡ�
} // end of function PolyConnectivity_remove_edge

// ��һ�������߽�������dual mesh
template <typename MeshT>
inline MeshT dual_mesh(MeshT& _m) {
	MeshT m;
	// step one, find the center of each face of _m
	OpenMesh::FPropHandleT<MeshT::VertexHandle> fp_new_vertex;
	_m.add_property(fp_new_vertex);
	for (typename MeshT::FaceIter f_it(_m.faces_begin()), f_end(_m.faces_end());
		f_it != f_end; ++f_it) {			 
			_m.property(fp_new_vertex, f_it) = m.add_vertex(DGP::face_center(_m, f_it.handle()));
	}
	// step two, add each face of m, according to every vertex of _m.
	for (typename MeshT::CVIter cv_it(_m.vertices_begin()), cv_end(_m.vertices_end());
		cv_it != cv_end; ++cv_it) {
			std::vector<typename MeshT::VHandle> vec; vec.clear();
			typename MeshT::VOHIter voh_it(_m, cv_it);
			typename MeshT::HHandle next = _m.ccw_rotated_halfedge_handle(voh_it);
			while (next != voh_it.handle()) {
				vec.push_back(_m.property(fp_new_vertex, _m.face_handle(next)));
				next = _m.ccw_rotated_halfedge_handle(next);
			}
			vec.push_back(_m.property(fp_new_vertex, _m.face_handle(voh_it)));
			m.add_face(vec); 
	}
	_m.remove_property(fp_new_vertex);
	return m;
}

} // end of namespace DGP

#endif //dgp_meshtype_h