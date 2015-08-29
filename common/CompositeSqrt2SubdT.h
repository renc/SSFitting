//*************************************************************
// ��ϸ�Ĺ���������Guiqing Li and Weiyin Ma, 2006
// Composite sqrt2 subdivision surfaces. 
// 1. composite subdi as generalization of 4-direction box splines surfaces.
// 2. composite subdi as extension of B-spline surfaces(2-direction box splines surfaces).
// 
// rencanjiang@163.com 2008-12
//*************************************************************

#ifndef dgp_sqrt2subd_h
#define dgp_sqrt2subd_h

#include "OpenMeshAll.h"
 
namespace DGP{

template <typename QuadMeshT, typename RealType = float>
class CompositeSqrt2SubdT : private OpenMesh::Utils::Noncopyable
{
public:
	typedef RealType                                real_t;
	typedef QuadMeshT                               mesh_t;	

public:
	CompositeSqrt2SubdT(void): attached_(NULL) //�������Կ�����һ��operator, ���԰󶨲�ͬ��mesh.
	{ }

	~CompositeSqrt2SubdT()   { detach(); }
public:
	const char *name() const { return "Composite Sqrt2 Subdi for Quads"; }

	bool attach(mesh_t &_m) {
		if ( attached_ == &_m )
			return true;

		detach();
		if (prepare( _m )) {
			attached_ = &_m;
			return true;
		}
		return false; 
	}
	// Primal unified sqrt2 subdivision. 
	bool primal_unified_sqrt2_subdivide(mesh_t& _mesh, size_t _m) {
		// Ԥ����, Sf(0)
		for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) {
			_mesh.property(fp_pos_, f_it) = mesh_t::Point(0,0,0);
		}
		// Sqar2 splitting operator,  
		for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end()); e_it != e_end; ++e_it) {
			_mesh.status(e_it).set_tagged(true);
		}
		for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) {
			mesh_t::VHandle new_fv = _mesh.new_vertex(_mesh.property(fp_pos_, f_it.handle())); 
			_mesh.split(f_it, new_fv);  
		}
		for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end()); e_it != e_end; ++e_it) {
			if (_mesh.status(e_it).tagged() == true) {  
				_mesh.status(e_it).set_tagged(false);
				DGP::PolyConnectivity_remove_edge(_mesh, e_it.handle()); 
			} 
		} 
		_mesh.garbage_collection();	 
		// Sv(2)
		for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end()); v_it != v_end; ++v_it) {
			_mesh.point(v_it) *= 2.0;
		}   
		// (Afv*Avf)^m
		for (int i = 0; i < _m; ++i) {
			// Avf
			for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) {
				mesh_t::Point pos(0,0,0);
				//for (typename mesh_t::FVIter fv_it(_mesh, f_it); fv_it; ++fv_it) {
				//	pos += _mesh.point(fv_it); 
				//}
				typename mesh_t::FVIter fv_it(_mesh, f_it);
				pos += _mesh.point(fv_it); 
				pos += _mesh.point(++fv_it); 
				pos += _mesh.point(++fv_it); 
				pos += _mesh.point(++fv_it); 
				_mesh.property(fp_pos_, f_it) = pos * 0.25;//����Sqrt2�õ�������4���ǵ�.
			}  
			// Afv
			for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end()); v_it != v_end; ++v_it) {
				mesh_t::Point pos(0, 0, 0);
				int count = 0;
				for (typename mesh_t::VFIter vf_it(_mesh, v_it); vf_it; ++vf_it) {
					pos += _mesh.property(fp_pos_, vf_it);  
					++count;
				} 
				_mesh.property(vp_pos_, v_it) = pos / count; 
			} 	 
		}  
		for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end()); v_it != v_end; ++v_it) {
			_mesh.set_point(v_it, _mesh.property(vp_pos_, v_it));
		}  
		_mesh.update_normals();

		return true;
	}
	// Primal unified 1-4 subdivision based on the sqrt2 splitting.
	bool primal_unified_1to4_subdivide(mesh_t &_mesh, size_t _m) {

		// Ԥ����, ���嶥������ͺ���ĵ�Ϊ������.
		for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end()); v_it != v_end; ++v_it) {
			_mesh.property(vp_flag_, v_it) = VTF_ZERO;
		} //Sf(0)
		for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) {
			_mesh.property(fp_pos_, f_it) = mesh_t::Point(0,0,0);
		}
		// Sqar2 splitting operator,  
		for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end()); e_it != e_end; ++e_it) {
			_mesh.status(e_it).set_tagged(true);
		}
		for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) {
			mesh_t::VHandle new_fv = _mesh.new_vertex(_mesh.property(fp_pos_, f_it.handle())); 
			_mesh.split(f_it, new_fv); 
			_mesh.property(vp_flag_, new_fv) = VTF_ONE;
		}
		for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end()); e_it != e_end; ++e_it) {
			if (_mesh.status(e_it).tagged() == true) {  
				_mesh.status(e_it).set_tagged(false);
				DGP::PolyConnectivity_remove_edge(_mesh, e_it.handle()); 
			} 
		}
		_mesh.garbage_collection();	
		// Avf
		for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) {
			mesh_t::Point pos(0,0,0);
			for (typename mesh_t::FVIter fv_it(_mesh, f_it); fv_it; ++fv_it) {
				pos += _mesh.point(fv_it);
			}
			_mesh.property(fp_pos_, f_it) = pos * 0.25;//����Sqrt2�õ�������4���ǵ�.
		}
		// Sf(2)
		for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) {
			_mesh.property(fp_pos_, f_it) *= 2;
		}
		// A'fv(1)
		for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end()); v_it != v_end; ++v_it) {
			if (_mesh.property(vp_flag_, v_it) == VTF_ONE) {
				mesh_t::Point pos(0, 0, 0);
				int count = 0;
				for (typename mesh_t::VVIter vv_it(_mesh, v_it); vv_it; ++vv_it) {
					pos += _mesh.point(vv_it);  
					++count;
				} 
				if (count != 4) std::cout << "Error: ����Ӧ�������ֻ��4����������.\n";
				_mesh.property(vp_pos_, v_it) = pos * 0.25; 
			} else {
				_mesh.property(vp_pos_, v_it) = _mesh.point(v_it);  
			}
		} 
		for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end()); v_it != v_end; ++v_it) {
			_mesh.set_point(v_it, _mesh.property(vp_pos_, v_it));
		}
		// Sqar2 splitting operator,  
		for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end()); e_it != e_end; ++e_it) {
			_mesh.status(e_it).set_tagged(true);
		}
		for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) {
			mesh_t::VHandle new_fv = _mesh.new_vertex(_mesh.property(fp_pos_, f_it.handle())); 
			_mesh.split(f_it, new_fv); 
			_mesh.property(vp_flag_, new_fv) = VTF_TWO;
		}
		for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end()); e_it != e_end; ++e_it) {
			if (_mesh.status(e_it).tagged() == true) {  
				_mesh.status(e_it).set_tagged(false);
				DGP::PolyConnectivity_remove_edge(_mesh, e_it.handle()); 
			} 
		}
		_mesh.garbage_collection();	
		// ����, �����һ��linear 1-4 splitting.
		// (Afv * Avf)^m
		for (int i = 0; i < _m; ++i) {
			// Avf
			for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) {
				mesh_t::Point pos(0,0,0);
				for (typename mesh_t::FVIter fv_it(_mesh, f_it); fv_it; ++fv_it) {
					pos += _mesh.point(fv_it);
				}
				_mesh.property(fp_pos_, f_it) = pos * 0.25;//����Sqrt2�õ�������4���ǵ�.
			}
			// Afv
			for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end()); v_it != v_end; ++v_it) {
				mesh_t::Point pos(0, 0, 0);
				int count = 0;
				for (typename mesh_t::VFIter vf_it(_mesh, v_it); vf_it; ++vf_it) {
					pos += _mesh.property(fp_pos_, vf_it);  
					++count;
				} 
				_mesh.property(vp_pos_, v_it) = pos / count; 
			} 			
		}
		for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end()); v_it != v_end; ++v_it) {
			_mesh.set_point(v_it, _mesh.property(vp_pos_, v_it));
		}
		_mesh.update_normals();

		return true;
	} // end of primal_unified_1to4_subdivide().

	void detach(void)
	{
		if ( attached_ ) {
			cleanup( *attached_ );
			attached_ = NULL;
		}
	}
protected:
	bool prepare( mesh_t& _m )
	{
		_m.add_property( vp_pos_ );// for quad.
		_m.add_property( fp_pos_ );// for quad.
		_m.add_property(vp_flag_);
		return true;
	}
	bool cleanup( mesh_t& _m )
	{
		_m.remove_property( vp_pos_ );
		_m.remove_property( fp_pos_ );// for quad.
		_m.remove_property( vp_flag_);
		return true;
	} 
	// �����������ֻ����Ϊ�Է���pure virtual func, 
	// ��Ϊ�������û��ʵ��������,����Ҳ����bool operator()������ϸ����.
	bool subdivide(const mesh_t& _m, size_t _n ) {
		std::cout << "Error: Nothing is done.\n"; 
	}
private: // topological modifiers

private: // geometry helper
	
private: // data
	mesh_t *attached_;
	 
	OpenMesh::VPropHandleT< typename mesh_t::Point > vp_pos_;//ÿһ����ĸ���֮���λ��.
	OpenMesh::EPropHandleT< typename mesh_t::Point > ep_pos_;//ÿһ�����������ӵĵ��λ��.
	OpenMesh::FPropHandleT< typename mesh_t::Point > fp_pos_;//�������ı���֮��ÿһ���涼��һ�����ĵ�.

	enum VTypeFlag { VTF_ZERO, VTF_ONE, VTF_TWO}; //Ϊ����������sqrt2ϸ�ָ��Լ���Ķ���.
	OpenMesh::VPropHandleT< unsigned int > vp_flag_;//0: old control pos; 1: ; 2: 
}; // end of class CompositeSqrt2SubdT

// usage function, ��_m����_n��ϸ��.
template <class QuadMeshT>
inline void sqrt2_subdi(QuadMeshT &_m, int _n) { 
	for (int i = 0; i < _n; ++i) {
		// a VF-type averaging operator, �²����ÿ�����.
		OpenMesh::FPropHandleT<OpenMesh::Vec3d> fp_cen;
		_m.add_property(fp_cen);
		for (typename QuadMeshT::FIter f_it(_m.faces_begin()), f_end(_m.faces_end()); f_it != f_end; ++f_it) {
			OpenMesh::Vec3d cen(0, 0, 0);
			int count = 0;
			for (typename QuadMeshT::FVIter fv_it(_m, f_it); fv_it; ++fv_it) {
				++ count;
				cen += _m.point(fv_it);
			}
			_m.property(fp_cen, f_it.handle()) = cen / count;
		}
		// a FV-type averaging operator, ÿ���ɵ㱻FV���͸���֮���������.
		OpenMesh::VPropHandleT<OpenMesh::Vec3d> vp_smoothpos;
		_m.add_property(vp_smoothpos);
		for (typename QuadMeshT::VIter v_it(_m.vertices_begin()), v_end(_m.vertices_end()); v_it != v_end; ++v_it) {
			int n_count_fv = 0;
			OpenMesh::Vec3d sum_fv(0, 0, 0);
			for (typename QuadMeshT::VFIter vf_it(_m, v_it); vf_it; ++vf_it) {
				++ n_count_fv;
				sum_fv += _m.property(fp_cen, vf_it);
			} 
			_m.property(vp_smoothpos, v_it) = sum_fv / n_count_fv; 
			_m.point(v_it) = _m.property(vp_smoothpos, v_it);
		}
		// Sqar2 splitting operator,  
		for (typename QuadMeshT::EIter e_it(_m.edges_begin()), e_end(_m.edges_end()); e_it != e_end; ++e_it) {
			_m.status(e_it).set_tagged(true);
		}
		for (typename QuadMeshT::FIter f_it(_m.faces_begin()), f_end(_m.faces_end()); f_it != f_end; ++f_it) {
			//typename QuadMeshT::VertexHandle new_fv = _m.new_vertex(_m.property(fp_cen, f_it.handle())); 
			//_m.split(f_it, new_fv); 
			_m.split(f_it, _m.property(fp_cen, f_it.handle())); 
		}
		for (typename QuadMeshT::EIter e_it(_m.edges_begin()), e_end(_m.edges_end()); e_it != e_end; ++e_it) {
			if (_m.status(e_it).tagged() == true) {  
				_m.status(e_it).set_tagged(false);
				DGP::PolyConnectivity_remove_edge(_m, e_it.handle()); 
			} 
		}//
		_m.remove_property(fp_cen);
		_m.remove_property(vp_smoothpos);

		_m.garbage_collection();		
	}
 
	// compute face & vertex normals
	_m.update_normals();
} 
} // end of namespace DGP
#endif //dgp_sqrt2subd_h