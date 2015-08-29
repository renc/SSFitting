
#include "CLoop.h"

namespace DGP {

//
//void cloop_subdi_1(TriMesh &_m) {
//	cloop_subdi(_m, 1);
//}
void cloop_subdi(TriMesh &_m, int _n) {
	if (_n <= 0) return;
	typedef OpenMesh::Subdivider::Uniform::LoopT<TriMesh, double> Loop;
	Loop subdivide_obj;

	//std::cout << "Subdivide " << n  << " times with '" << subdivide_obj.name() << "'\n";

	// use interface 1
	subdivide_obj(_m, _n);

	// compute face & vertex normals
	_m.update_normals(); 
}
OpenMesh::Vec3d cloop_limit_pos(const TriMesh &_m, const TriMesh::VertexHandle &_vh) {
	int k = _m.valence(_vh);
	double inv_k = 1.0/double(k); double t = (3.0 + 2.0 * cos( 2.0 * M_PI * inv_k) ); // double beta = inv_k * (40.0 - t * t)/64.0;
	double k_alpha = 1.0 / ( 24.0/(40.0-t*t) + 1); double alpha = inv_k * k_alpha;

	TriMesh::Point limit_pos = OpenMesh::Vec3d(0, 0, 0);
	for (TriMesh::ConstVertexVertexIter vv_it(_m, _vh); vv_it; ++vv_it) {
		limit_pos += _m.point(vv_it.handle()) * alpha; 
	}
	limit_pos += _m.point(_vh) * (1.0 - k_alpha); 
	return limit_pos;
}
OpenMesh::Vec3d cloop_tangent_vector_along_edge(const TriMesh &_m, const TriMesh::HalfedgeHandle & _hh) {
	int k = _m.valence(_m.from_vertex_handle(_hh)); // valence of the vertex
	OpenMesh::Vec3d t(0, 0, 0);
	TriMesh::HHandle h = _hh;

	double PI_2_inv_k = 2 * M_PI / k;
	for (int i = 0; i < k; ++i) {
		t += (2.0 / k) * cos(PI_2_inv_k * i) * _m.point(_m.to_vertex_handle(h));
		h = _m.opposite_halfedge_handle(_m.prev_halfedge_handle(h));
	}
	if (h != _hh) std::cout << "Error: h shoud be = _hh.\n";
	return t;
}
OpenMesh::Vec3d cloop_limit_nor(const TriMesh &_m, const TriMesh::VertexHandle &_vh) {
	OpenMesh::Vec3d t1(0, 0, 0); // tangent vector
	OpenMesh::Vec3d t2(0, 0, 0);
	OpenMesh::Vec3d t3(0, 0, 0);
	int i = 0;
	int k = _m.valence(_vh);
	double PI_2_inv_k = 2 * M_PI / k;
	for (TriMesh::ConstVertexOHalfedgeIter voh_it(_m, _vh); voh_it; ++voh_it, ++i) {
		t1 += cos(PI_2_inv_k * i) * _m.point(_m.to_vertex_handle(voh_it)); 
		t2 += sin(PI_2_inv_k * i) * _m.point(_m.to_vertex_handle(voh_it)); 
		t3 += cos(PI_2_inv_k * i) * _m.point(_m.to_vertex_handle(_m.opposite_halfedge_handle(_m.prev_halfedge_handle(voh_it)))); 
	} 
	// cross(t1, t2).normalize() == (-1) * cross(t1, t3).normalize() 奇怪. 
	// return cross(t1, t3).normalize(); //正方向. //cross(t1, t2).normalize();//反向.
	return cross(t2, t1).normalize();
}
void cloop_limit_surf(TriMesh &_m) {
	OpenMesh::VPropHandleT<OpenMesh::Vec3d> vp_limit_pos, vp_limit_nor;
	_m.add_property(vp_limit_pos); _m.add_property(vp_limit_nor);

	for (TriMesh::VIter v_it(_m.vertices_begin()), v_end(_m.vertices_end()); v_it != v_end; ++v_it) {
		_m.property(vp_limit_pos, v_it) = cloop_limit_pos(_m, v_it.handle());		
	}
	for (TriMesh::VIter v_it(_m.vertices_begin()), v_end(_m.vertices_end()); v_it != v_end; ++v_it) {
		_m.property(vp_limit_nor, v_it) = cloop_limit_nor(_m, v_it.handle());		
	} 

	for (TriMesh::VIter v_it(_m.vertices_begin()), v_end(_m.vertices_end()); v_it != v_end; ++v_it) {
		_m.set_point(v_it, _m.property(vp_limit_pos, v_it));		
	}
	for (TriMesh::VIter v_it(_m.vertices_begin()), v_end(_m.vertices_end()); v_it != v_end; ++v_it) {
		_m.set_normal(v_it, _m.property(vp_limit_nor, v_it));		
	} 
	_m.remove_property(vp_limit_pos); _m.remove_property(vp_limit_nor);

	//for (TriMesh::VIter v_it(_m.vertices_begin()), v_end(_m.vertices_end()); v_it != v_end; ++v_it) {
	//	_m.set_point(v_it, cloop_limit_pos(_m, v_it.handle()));		
	//} //这是错误的 因为前面已经改了顶点坐标了 所以后面的法线就求得不对了
	//for (TriMesh::VIter v_it(_m.vertices_begin()), v_end(_m.vertices_end()); v_it != v_end; ++v_it) {
	//	_m.set_normal(v_it, cloop_limit_nor(_m, v_it.handle()));		
	//} 

} // end of cloop_limit_surf(...)
} // end of namespace DGP