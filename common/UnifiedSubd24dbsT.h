// 详细的工作内容在Guiqing Li and Weiyin Ma, 2007
// Unified Subdivision Generalizing 2- and 4-Direction Box Spline. 
//  
// rencanjiang@163.com 2008-12
//*************************************************************


#ifndef dgp_unifiedSubd24dbs
#define dgp_unifiedSubd24dbs
 
#include "OpenMeshAll.h"

namespace DGP{

	template <typename QuadMeshT, typename RealType = float>
	class UnifiedSubd24dbsT : private OpenMesh::Utils::Noncopyable
	{
	public:
		typedef RealType                                real_t;
		typedef QuadMeshT                               mesh_t;	

	public:
		UnifiedSubd24dbsT(void)  //这个类可以看做是一个operator, 可以绑定不同的mesh.
		{ }

		~UnifiedSubd24dbsT()   {  }
	public:
		const char *name() const 
		{ return "Unified Subdivision Generalizing 2- and 4-Direction Box Spline for Quads"; }
 
		// Unified subdivision generalize 2 and 4 direction box spline
		bool unified_subdivide(mesh_t &_mesh, size_t _m, size_t _n) {
			//// 预处理, 为下面是sqrt2 split做准备的.
			//for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) {
			//	_mesh.property(fp_pos, f_it) = mesh_t::Point(0,0,0);
			//}
			 
			for (int time = 0; time < 2; ++time) {
				int m_n = 0; 
				if (time == 0) m_n = _n;
				else m_n = _m;

					// Sqrt2 splitting operator,  
					for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end()); e_it != e_end; ++e_it) {
						_mesh.status(e_it).set_tagged(true);
					}
					for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) {
						_mesh.split(f_it, typename typename mesh_t::Point(0, 0, 0));
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
					//增加了四边形之后每一个面都有一个中心点.
					OpenMesh::FPropHandleT< typename mesh_t::Point > fp_pos;
					_mesh.add_property(fp_pos);
					for (int i = 0; i < int(m_n / 2); ++i) {
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
							_mesh.property(fp_pos, f_it) = pos * 0.25;//猜想Sqrt2得到的面是4个角的.
						}  
						// Afv
						for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end()); 
							v_it != v_end; ++v_it) {
								mesh_t::Point pos(0, 0, 0);
								int count = 0;
								for (typename mesh_t::VFIter vf_it(_mesh, v_it); vf_it; ++vf_it) {
									pos += _mesh.property(fp_pos, vf_it);  
									++count;
								} 
								_mesh.set_point(v_it, pos / count);
						} 	  
					}
					_mesh.remove_property(fp_pos);
					for (int i = 0; i < m_n % 2; ++i) {
						_mesh.add_property(fp_pos);
						// Avf
						for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) {
							mesh_t::Point pos(0,0,0);
							//可以猜想Sqrt2得到的面是4个角的, 但是经过Td之后就不是了.
							int count_v = 0;
							for (typename mesh_t::FVIter fv_it(_mesh, f_it); fv_it; ++fv_it) {
								pos += _mesh.point(fv_it);
								++ count_v;
							}
							_mesh.property(fp_pos, f_it) = pos / count_v;
						}  
						// Td
						mesh_t dual_mesh;
						OpenMesh::FPropHandleT<mesh_t::VHandle > fp_new_vertex;
						_mesh.add_property(fp_new_vertex);
						for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end());
							f_it != f_end; ++f_it) {
								_mesh.property(fp_new_vertex, f_it) 
									= dual_mesh.add_vertex(_mesh.property(fp_pos, f_it));
						}
						for (typename mesh_t::CVIter cv_it(_mesh.vertices_begin()), cv_end(_mesh.vertices_end());
							cv_it != cv_end; ++cv_it) {
								std::vector<typename mesh_t::VHandle> vec; vec.clear();
								typename mesh_t::VOHIter voh_it(_mesh, cv_it);
								typename mesh_t::HHandle next = _mesh.ccw_rotated_halfedge_handle(voh_it);
								while (next != voh_it.handle()) {
									vec.push_back(_mesh.property(fp_new_vertex, _mesh.face_handle(next)));
									next = _mesh.ccw_rotated_halfedge_handle(next);
								}
								vec.push_back(_mesh.property(fp_new_vertex, _mesh.face_handle(voh_it)));
								dual_mesh.add_face(vec); 
						}
						_mesh.remove_property(fp_pos);
						_mesh = dual_mesh;
					}
			} // end of for. 2times
			_mesh.update_normals();

			return true;
		}  // end of unified_subdivide(...). 
		bool unified_subdivide_features(mesh_t &_mesh, 
			std::vector<int > & _crease_edges_arr,
			std::vector<int > & _crease_vertices_arr, 
			size_t _m, size_t _n) 
		{
			// 设置特征信息
			
			if (_crease_vertices_arr.size() != _crease_edges_arr.size()) {
				std::cout << "Error: 暂时的实现是只支持crease vertex & crease/boundary edge, not corner or dart.";
				std::cout << " v " << _crease_vertices_arr.size() << ", e " << _crease_edges_arr.size() << ". \n";
			}
			for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end());
				e_it != e_end; ++e_it) {
					_mesh.status(e_it).set_feature(false);
			}
			for (int i = 0; i < _crease_edges_arr.size(); ++i) {
				_mesh.status(_mesh.edge_handle(_crease_edges_arr[i])).set_feature(true);
			}
			OpenMesh::VPropHandleT<DGP::v_feature_type> vp_type;
			_mesh.add_property(vp_type);
			for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end());
				v_it != v_end; ++v_it) {
					_mesh.property(vp_type, v_it) = DGP::SMOOTH_VFT;
			}
			for (int i = 0; i < _crease_vertices_arr.size(); ++i) {
				_mesh.property(vp_type, _mesh.vertex_handle(_crease_vertices_arr[i])) = DGP::CREASE_VFT;
			}
			_crease_edges_arr.clear();
			_crease_vertices_arr.clear();
			 //当考虑到crease modeling时候需要分odd steps and even steps of refinement.
			// 具体请参考文章.
			// 1. odd steps of refinement
			// Sqrt2 splitting operator,  
			for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end()); e_it != e_end; ++e_it) {
				if (_mesh.status(e_it).feature() == false) // 一般的边, 非特征边
					_mesh.status(e_it).set_tagged(true); // will be removed.
			}
			for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) {
				typename mesh_t::VHandle v_new = _mesh.add_vertex(typename typename mesh_t::Point(0, 0, 0));
				_mesh.split(f_it, v_new);

				_mesh.property(vp_type, v_new) = DGP::SMOOTH_VFT;
				for (typename mesh_t::VOHIter voh_it(_mesh, v_new); voh_it; ++voh_it) {
					_mesh.status(_mesh.edge_handle(voh_it)).set_tagged(false);
				}
			}
			for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end()); e_it != e_end; ++e_it) {
				if (_mesh.status(e_it).tagged() == true) {  
					_mesh.status(e_it).set_tagged(false);
					DGP::PolyConnectivity_remove_edge(_mesh, e_it.handle()); 
				} 
			} 
			_mesh.garbage_collection();	 
			OpenMesh::FPropHandleT<bool> fp_special;// special or not.
			_mesh.add_property(fp_special); 
			for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end());
				f_it != f_end; ++f_it) {
					_mesh.property(fp_special, f_it) = false;
			}  
			//// 在上面的基础上特征边一分为二, 参考Loop, CC, 等细分实现中的split_edge()函数.
			for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end()); e_it != e_end; ++e_it) {
				if (_mesh.status(e_it).feature() ) {//  特征边
					split_edge_add_feature(_mesh, vp_type, fp_special, e_it.handle());  
				} 
			} // end of operator Sqrt2.
            // -----这里_mesh 含有的property包括了 edge feature, vertex type, face special.
			// 还要注意face special flag是要在_mesh.garbage_collection()之后才能设置的, 否则后者会修改.
			std::cout << "Odd Sqrt2, ";
			 
 
			// Sv(2), 这个效果在后面的Afv*Avf后会被消除，但是涉及特征边之后就不是了，所以考虑特征时就别尝试R^(1,1).
			// 因为当_n=1时候就不会进入Afv*Avf，Sv(2)的效果没有消去，对后面有影响.
			for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end()); v_it != v_end; ++v_it) {
				_mesh.point(v_it) *= 2.0;
			}   
			// Afv * Avf
			OpenMesh::FPropHandleT< typename mesh_t::Point > fp_pos;//增加了四边形之后每一个面都有一个中心点.
			_mesh.add_property(fp_pos);
			for (int i = 0; i < int(_n / 2); ++i) {
				// Avf
				for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); 
					f_it != f_end; ++f_it) 
				{
					mesh_t::Point pos(0,0,0);
					if (_mesh.property(fp_special, f_it) == false) {
						typename mesh_t::FVIter fv_it(_mesh, f_it);//猜想Sqrt2得到的面是4个角的.	 
						pos += _mesh.point(fv_it); 
						pos += _mesh.point(++fv_it); 
						pos += _mesh.point(++fv_it); 
						pos += _mesh.point(++fv_it); 
						_mesh.property(fp_pos, f_it) = pos * 0.25;
					} 
				}  
				//if (i == 0) { //取消Sv(2)对特征点的放大效果.
				//	for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end());
				//		v_it != v_end; ++v_it) {
				//			if (_mesh.property(vp_type, v_it) == DGP::CREASE_VFT) {
				//				_mesh.set_point(v_it, _mesh.point(v_it)*0.5);
				//			}
				//	}
				//}
				//for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); 
				//	f_it != f_end; ++f_it) 
				//{
				//	if (_mesh.property(fp_special, f_it)) {
				//		typename mesh_t::HHandle h0 = _mesh.halfedge_handle(f_it);
				//		typename mesh_t::HHandle h1 = _mesh.next_halfedge_handle(h0);
				//		while (!_mesh.status(_mesh.edge_handle(h0)).feature()
				//			|| !_mesh.status(_mesh.edge_handle(h1)).feature()) {
				//				h0 = h1; 
				//				h1 = _mesh.next_halfedge_handle(h0);
				//		}
				//		if (!_mesh.status(_mesh.edge_handle(h0)).feature())
				//			std::cerr << "Error: the while above, h0.\n";
				//		if (!_mesh.status(_mesh.edge_handle(h1)).feature())
				//			std::cerr << "Error: the while above, h1.\n";
				//		typename mesh_t::VHandle v0 = _mesh.to_vertex_handle(h0); //就是那个valence=2的特征点.
				//		if (_mesh.valence(v0) != 2) std::cout << "Error: v0 should be 2.\n";

				//		_mesh.set_point(v0, (_mesh.point(_mesh.from_vertex_handle(h0))
				//			+ _mesh.point(_mesh.to_vertex_handle(h1)))*0.5); //两边中点
				//		_mesh.property(fp_pos, f_it) = _mesh.point(v0); 
				//	}
				//} end of Avf
				OpenMesh::EPropHandleT<OpenMesh::Vec3d> ep_mid;
				_mesh.add_property(ep_mid); //每一段特征边的中点
				for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end());
					e_it != e_end; ++e_it) {
						if (_mesh.status(e_it).feature()) {
							_mesh.property(ep_mid, e_it) = DGP::edge_midpos(_mesh, e_it.handle());
						}
				}
				// Afv
				// 首先更新特征边上的特征点
				for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end()); 
					v_it != v_end; ++v_it) {
						if (_mesh.property(vp_type, v_it) == DGP::CREASE_VFT) { 
							std::vector<typename mesh_t::EHandle> e;
							for (typename mesh_t::VEIter ve_it(_mesh, v_it); ve_it; ++ve_it) {
								if (_mesh.status(ve_it.handle()).feature()) e.push_back(ve_it.handle());
							}
							if (e.size() != 2) 
								std::cout << "Error: crease vertex only has 2 crease edges, not " << e.size() << "\n";
							_mesh.set_point(v_it, 
								(_mesh.property(ep_mid, e[0]) + _mesh.property(ep_mid, e[1]))*0.5);
							//然后这个点作为两边的两个面的中点
							if (_mesh.valence(v_it) == 2) {
								_mesh.property(fp_pos, _mesh.face_handle(_mesh.halfedge_handle(e[0], 0))) = _mesh.point(v_it);
								_mesh.property(fp_pos, _mesh.face_handle(_mesh.halfedge_handle(e[0], 1))) = _mesh.point(v_it);
							}

						}
				} 
				_mesh.remove_property(ep_mid);
				// 更新一般的normal vertices
				for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end()); 
					v_it != v_end; ++v_it) {
						if (_mesh.property(vp_type, v_it) == DGP::SMOOTH_VFT) {
							mesh_t::Point pos(0, 0, 0);
							int count = 0;
							for (typename mesh_t::VFIter vf_it(_mesh, v_it); vf_it; ++vf_it) {
								pos += _mesh.property(fp_pos, vf_it);  
								++count;
							} 
							_mesh.set_point(v_it, pos / count);
						}
				} 	
				// end of Afv
			} 
			_mesh.remove_property(fp_pos);
			std::cout << "Odd Afv*Avf end, ";
			// (Td*Avf)mod(n, 2)
			for (int i = 0; i < _n % 2; ++i) {
				OpenMesh::FPropHandleT<OpenMesh::Vec3d> fp_pos; //这里是自给自足
				_mesh.add_property(fp_pos);
				// Avf, 求每一个面的中心点.
				for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); 
					f_it != f_end; ++f_it) 
				{
					if (_mesh.property(fp_special, f_it) == false) { //非特殊面的中心点为所有角点的平均.
						//可以猜想Sqrt2得到的面是4个角的, 但是经过Td之后就不是了.
						_mesh.property(fp_pos, f_it) = DGP::face_center(_mesh, f_it);
					} 		
				}  
				if (_n == 1 && i == 0) { //_n==1时前面的(Afv*Avf)^int(_n/2), 因此需要取消Sv(2)对特征点的放大效果.
					for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end());
						v_it != v_end; ++v_it) {
							if (_mesh.property(vp_type, v_it) == DGP::CREASE_VFT) {
								_mesh.set_point(v_it, _mesh.point(v_it)*0.5);
							}
					}
					for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); 
						f_it != f_end; ++f_it) 
					{
						if (_mesh.property(fp_special, f_it)) {
							typename mesh_t::HHandle h0 = _mesh.halfedge_handle(f_it);
							typename mesh_t::HHandle h1 = _mesh.next_halfedge_handle(h0);
							while (!_mesh.status(_mesh.edge_handle(h0)).feature()
								|| !_mesh.status(_mesh.edge_handle(h1)).feature()) {
									h0 = h1; 
									h1 = _mesh.next_halfedge_handle(h0);
							}
							if (!_mesh.status(_mesh.edge_handle(h0)).feature())
								std::cerr << "Error: the while above, h0.\n";
							if (!_mesh.status(_mesh.edge_handle(h1)).feature())
								std::cerr << "Error: the while above, h1.\n";
							typename mesh_t::VHandle v0 = _mesh.to_vertex_handle(h0); //就是那个valence=2的特征点.
							if (_mesh.valence(v0) != 2) std::cout << "Error: should be 2.\n";

							_mesh.set_point(v0, (_mesh.point(_mesh.from_vertex_handle(h0))
								+ _mesh.point(_mesh.to_vertex_handle(h1)))*0.5); //两边中点
							_mesh.property(fp_pos, f_it) = _mesh.point(v0); 
						}					
					}  
				}
				
				// end of Avf  
				//_mesh.update_normals();return true;// for test only.

				// Td
				mesh_t dual_mesh; //下面谨记dual_mesh也需要设置3个property,
				// 分别是crease edge, crease vertex, special face.
				// 于是每当加入一个新点时候都要设置其特征性质, 每当加入一个新面时候也设置其特征和边的特征性质.
				OpenMesh::VPropHandleT<DGP::v_feature_type> vp_type_dual;//crease edge feature直接设.
				OpenMesh::FPropHandleT<bool > fp_special_dual;
				dual_mesh.add_property(vp_type_dual);
				dual_mesh.add_property(fp_special_dual);
				// 首先添加顶点
				OpenMesh::FPropHandleT<mesh_t::VHandle > fp_dual_vertex; //这个面对于dual_mesh上的哪个顶点.
				_mesh.add_property(fp_dual_vertex);
				for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end());
					f_it != f_end; ++f_it) {
						if (_mesh.property(fp_special, f_it) == false) {
							typename mesh_t::VHandle v = dual_mesh.add_vertex(_mesh.property(fp_pos, f_it));
							dual_mesh.property(vp_type_dual, v) = DGP::SMOOTH_VFT;//dual_mesh每加一点都要设置.
							_mesh.property(fp_dual_vertex, f_it) = v;							
						}
				} 
				//这个特征边对应的valence=2的点是要被加入了dual_mesh.
				OpenMesh::EPropHandleT<bool> ep_check;
				_mesh.add_property(ep_check);
				for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end());
					e_it != e_end; ++ e_it) {
						_mesh.property(ep_check, e_it) = false;
				}
				for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end());
					e_it != e_end; ++ e_it) {
						if (_mesh.status(e_it).feature() && _mesh.property(ep_check, e_it) == false) {
							// 找valence = 2的顶点.
							// /   f0             \   上三角面
							// ===h0===>v0===h1===>,  特征边e_it
							// \            f1    /   下三角面
							// v0 is the mid vertex of the crease edge, which valence equals 2.
							// v0 只加入一次，但是被上下两个三角面共享
							typename mesh_t::HHandle h0 = _mesh.halfedge_handle(e_it, 0), h1;
							if (_mesh.valence(_mesh.to_vertex_handle(h0)) == 2)
								h1 = _mesh.next_halfedge_handle(h0);
							else if (_mesh.valence(_mesh.from_vertex_handle(h0)) == 2) {
								h1 = h0;
								h0 = _mesh.prev_halfedge_handle(h0);
							} else std::cerr << "Error: 定位两个特征半边时候错了.\n";
							
							typename mesh_t::VHandle new_v = dual_mesh.add_vertex(_mesh.point(_mesh.to_vertex_handle(h0)));
							dual_mesh.property(vp_type_dual, new_v) = DGP::CREASE_VFT;//dual_mesh每加一点都要设置.
							typename mesh_t::FHandle f0 = _mesh.face_handle(h0);
							typename mesh_t::FHandle f1 = _mesh.face_handle(_mesh.opposite_halfedge_handle(h1));
							if (f0 == f1) std::cout << "Error: 上面的h0和h1的错了.\n";
							if (f0.is_valid()) _mesh.property(fp_dual_vertex, f0) = new_v;
							if (f1.is_valid()) _mesh.property(fp_dual_vertex, f1) = new_v;
							_mesh.property(ep_check, _mesh.edge_handle(h0)) = true;
							_mesh.property(ep_check, _mesh.edge_handle(h1)) = true;
						}
				} 
				_mesh.remove_property(ep_check);
				//至此, 每一个面都有一个新的顶点用于之后生成dual mesh.
				
				//// 然后就是添加面, 先添加normal vertex对应的面， 再添加special vertex对应的面
				for (typename mesh_t::CVIter cv_it(_mesh.vertices_begin()), cv_end(_mesh.vertices_end());
					cv_it != cv_end; ++cv_it) {
						if (_mesh.property(vp_type, cv_it.handle()) == DGP::SMOOTH_VFT) {
							std::vector<typename mesh_t::VHandle> vec; vec.clear();
							typename mesh_t::VOHIter voh_it(_mesh, cv_it);
							typename mesh_t::HHandle next = _mesh.ccw_rotated_halfedge_handle(voh_it);
							while (next != voh_it.handle()) {
								vec.push_back(_mesh.property(fp_dual_vertex, _mesh.face_handle(next)));
								next = _mesh.ccw_rotated_halfedge_handle(next);
							}
							vec.push_back(_mesh.property(fp_dual_vertex, _mesh.face_handle(voh_it)));  
							typename mesh_t::FHandle new_f = dual_mesh.add_face(vec);  
							dual_mesh.property(fp_special_dual, new_f) = false; //每加一个面都要设置一下.
							for (typename mesh_t::FEIter fe_it(dual_mesh, new_f); fe_it; ++fe_it) {
								dual_mesh.status(fe_it.handle()).set_feature(false);
							}
						}						
				}
				for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end());
					v_it != v_end; ++v_it) {
						if (_mesh.property(vp_type, v_it) == DGP::CREASE_VFT && _mesh.valence(v_it) > 2) {
								// 找两个特征半边
								//                  e_arr[0]
								//<-----------v_it---------->
								//	e_arr[1]
								std::vector<typename mesh_t::HHandle> e_arr;
								for (typename mesh_t::VOHIter vhe_it(_mesh, v_it); vhe_it; ++vhe_it) {
									if (_mesh.status(_mesh.edge_handle(vhe_it.handle())).feature()) 
										e_arr.push_back(vhe_it.handle());
								}
								if (e_arr.size() != 2) std::cout << "Error: e_arr has 2 crease edges.\n";
								
								// 自己也要加入到dual mesh中成为其一个点.
								typename mesh_t::VHandle new_v = dual_mesh.add_vertex(_mesh.point(v_it));
								dual_mesh.property(vp_type_dual, new_v) = DGP::CREASE_VFT;//dual_mesh每加一点都要设置.
								
								if (_mesh.face_handle(e_arr[0]).is_valid()) { // 先加入一个半面, this test is for boundary.
									std::vector<typename mesh_t::VHandle> vec; vec.clear();
									vec.push_back(new_v);
									typename mesh_t::HHandle  h = e_arr[0]; 
									while (h != e_arr[1]) {
										vec.push_back(_mesh.property(fp_dual_vertex, _mesh.face_handle(h)));
										h = _mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(h));
									}
									//std::cout << "@";std::cout << vec.size();
									//if (vec.size)
									typename mesh_t::FHandle new_f = dual_mesh.add_face(vec);  //std::cout << "@";
									dual_mesh.property(fp_special_dual, new_f) = true; //每加一个面都要设置一下.
									for (typename mesh_t::FEIter fe_it(dual_mesh, new_f); fe_it; ++fe_it) {
										dual_mesh.status(fe_it.handle()).set_feature(false);
									} 
								}
								if (_mesh.face_handle(e_arr[1]).is_valid()) { // 再加入另一个半面, this test is for boundary.
									std::vector<typename mesh_t::VHandle> vec; vec.clear();
									vec.push_back(new_v);
									typename mesh_t::HHandle  h = e_arr[1]; 
									while (h != e_arr[0]) {
										vec.push_back(_mesh.property(fp_dual_vertex, _mesh.face_handle(h)));
										h = _mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(h));
									}//std::cout << "#";std::cout<<vec.size();
									typename mesh_t::FHandle new_f = dual_mesh.add_face(vec); //std::cout << "#"; 
									dual_mesh.property(fp_special_dual, new_f) = true; //每加一个面都要设置一下.
									for (typename mesh_t::FEIter fe_it(dual_mesh, new_f); fe_it; ++fe_it) {
										dual_mesh.status(fe_it.handle()).set_feature(false);
									} 
								}
								if (_mesh.face_handle(e_arr[0]).is_valid())
									dual_mesh.status(dual_mesh.edge_handle(dual_mesh.find_halfedge(new_v, _mesh.property(fp_dual_vertex, _mesh.face_handle(e_arr[0]))))).set_feature(true);
								else 
									dual_mesh.status(dual_mesh.edge_handle(dual_mesh.find_halfedge(new_v, _mesh.property(fp_dual_vertex, _mesh.face_handle(_mesh.opposite_halfedge_handle(e_arr[0])))))).set_feature(true);
								if (_mesh.face_handle(e_arr[1]).is_valid())
									dual_mesh.status(dual_mesh.edge_handle(dual_mesh.find_halfedge(new_v, _mesh.property(fp_dual_vertex, _mesh.face_handle(e_arr[1]))))).set_feature(true);
								else 
									dual_mesh.status(dual_mesh.edge_handle(dual_mesh.find_halfedge(new_v, _mesh.property(fp_dual_vertex, _mesh.face_handle(_mesh.opposite_halfedge_handle(e_arr[1])))))).set_feature(true);
								
						}
				}
				_mesh.remove_property(fp_pos);
				_mesh.remove_property(fp_dual_vertex);
				_mesh.remove_property(vp_type);
				_mesh.remove_property(fp_special);
				_mesh = dual_mesh;
				_mesh.add_property(vp_type);
				_mesh.add_property(fp_special);
				for (typename mesh_t::EIter e_it(dual_mesh.edges_begin()), e_end(dual_mesh.edges_end());
					e_it != e_end; ++e_it) {
						_mesh.status(e_it).set_feature(dual_mesh.status(e_it).feature());
				}
				for (typename mesh_t::VIter v_it(dual_mesh.vertices_begin()), v_end(dual_mesh.vertices_end());
					v_it != v_end; ++v_it) {
						_mesh.property(vp_type, v_it) = dual_mesh.property(vp_type_dual, v_it);
				}
				for (typename mesh_t::FIter f_it(dual_mesh.faces_begin()), f_end(dual_mesh.faces_end());
					f_it != f_end; ++f_it) {
						_mesh.property(fp_special, f_it) = dual_mesh.property(fp_special_dual, f_it); 
				}
				// end of Td
			}// end fo (Td*Avf)mod(n, 2)
			
			// check , 只考虑非边界crease edges（带环无dart无corner）的情况下应该满足 v=e=f。
			int n_crease_v=0, n_crease_e=0, n_crease_f=0;
			for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end());
				e_it != e_end; ++e_it) {
					if (_mesh.status(e_it).feature()) ++ n_crease_e;
			}
			for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end());
				v_it != v_end; ++v_it) {
					if (_mesh.property(vp_type, v_it) == DGP::CREASE_VFT) 
						++ n_crease_v;
			}
			for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end());
				f_it != f_end; ++f_it) {
					if (_mesh.property(fp_special, f_it))
						++n_crease_f;
			}
			std::cout << "Odd Td*Avf, ";
			std::cout << "Test: after odd step, "<< n_crease_v << "v, " << n_crease_e << "e, " << n_crease_f << "f.\n";
			

			// 2. Even steps of refinement
			// Sqrt2 splitting operator,  
			for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end()); e_it != e_end; ++e_it) {
				if (_mesh.status(e_it).feature() == false) // 一般的边, 非特征边
					_mesh.status(e_it).set_tagged(true);   // 所有非特征边都是准备将要被去掉的.
			}  
			for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); f_it != f_end; ++f_it) {
				if (_mesh.property(fp_special, f_it) == false)
					_mesh.split(f_it, typename typename mesh_t::Point(0, 0, 0));
					//_mesh.split(f_it, DGP::face_center(_mesh, f_it));//这只是测试用的, 算法中加入的点都是坐标为0的.
			}  
			OpenMesh::EPropHandleT<bool> ep_check;
			_mesh.add_property(ep_check);
			for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end());
				e_it != e_end; ++ e_it) {
					_mesh.property(ep_check, e_it) = false;
			} //std::cout << "for test, begin Even Sqrt2.\n"; 
			for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end());
				e_it != e_end; ++ e_it) {
					if (_mesh.status(e_it).feature() && _mesh.property(ep_check, e_it) == false) {
						// 找valence = 2的顶点.
						// /   f0             \   上三角特征面
						// ===h0===>v0===h1===>,  特征边e_it
						// \            f1    /   下三角特征面
						// v0 is the mid vertex of the crease edge, which valence equals 2.
						// v0 只加入一次，但是被上下两个三角面共享
						typename mesh_t::HHandle h0 = _mesh.halfedge_handle(e_it, 0), h1;
						if (_mesh.valence(_mesh.to_vertex_handle(h0)) == 2)
							h1 = _mesh.next_halfedge_handle(h0);
						else if (_mesh.valence(_mesh.from_vertex_handle(h0)) == 2) {
							h1 = h0;
							h0 = _mesh.prev_halfedge_handle(h0);
						} else std::cerr << "Error: 定位两个特征半边时候错了.\n";
						if (_mesh.face_handle(h0).is_valid())
							if (_mesh.property(fp_special, _mesh.face_handle(h0)) == false) 
								std::cout << "Error: should be special 1.\n";
						if (_mesh.face_handle(_mesh.opposite_halfedge_handle(h0)).is_valid())
							if (_mesh.property(fp_special, _mesh.face_handle(_mesh.opposite_halfedge_handle(h0))) == false) 
								std::cout << "Error: should be special 2.\n";

						//std::cout << "for test, split_face_1to2 a.\n"; 
						if (_mesh.face_handle(h0).is_valid())
							split_face_1to2(_mesh, h0, h1);
						//std::cout << "for test, split_face_1to2 b.\n"; 
						if (_mesh.face_handle(_mesh.opposite_halfedge_handle(h0)).is_valid())
							split_face_1to2(_mesh, _mesh.opposite_halfedge_handle(h1), _mesh.opposite_halfedge_handle(h0));
						//std::cout << "for test, split_face_1to2 c.\n"; 
						_mesh.property(ep_check, _mesh.edge_handle(h0)) = true;
						_mesh.property(ep_check, _mesh.edge_handle(h1)) = true;
						
					}
			}  
			_mesh.remove_property(ep_check);
			//std::cout << "for test, after split_face_1to2.\n";
			for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end()); e_it != e_end; ++e_it) {
				if (_mesh.status(e_it).tagged() == true) {  
					_mesh.status(e_it).set_tagged(false); //std::cout << "a";
					DGP::PolyConnectivity_remove_edge(_mesh, e_it.handle()); //std::cout << "b";
				} 
			}   //std::cout << "for test, before garbage_collection().\n";
			_mesh.garbage_collection();	 
			// end of Sqrt2. 往下好像不再需要special face flag了,但是crease edge and crease vertices还是需要记录的.
			std::cout << "Even Sqrt2.\n";
			// Sv(2) opertor
			for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end()); v_it != v_end; ++v_it) {
				_mesh.point(v_it) *= 2.0;
			}  //std::cout << "Even Sv2.\n";

			// (Afv*Avf)^int(_m/2)
			if (fp_pos.is_valid()) _mesh.remove_property(fp_pos);// fp's center pos for the Avf in even steps.
			_mesh.add_property(fp_pos);
			OpenMesh::EPropHandleT<OpenMesh::Vec3d> ep_pos; // edge's mid pos.
			_mesh.add_property(ep_pos);
			for (int i = 0; i < int(_m/2); ++i) {
				std::cout << "(Afv*Avf)^int(_m/2), loop i=" << i << ".\n";
				// Avf
				for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end());
					f_it != f_end; ++f_it) { //一般面的中点.
						_mesh.property(fp_pos, f_it) = DGP::face_center(_mesh, f_it.handle());
				} 
				for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end());
					e_it != e_end; ++e_it) {
						if (_mesh.status(e_it).feature()) { //一般边的中点.
							_mesh.property(ep_pos, e_it) = DGP::edge_midpos(_mesh, e_it.handle());
							 
							if (i == 0 && _n != 0) { // 由于之前的Sv(2), 第一次Avf要特殊处理.
								_mesh.property(ep_pos, e_it) *= 0.5;//首先边中点要再缩小一倍.
								
								// 特征边两边的面的中点第一次要特殊处理, 因为有一个点坐标为0.
								typename mesh_t::HHandle h0 = _mesh.halfedge_handle(e_it, 0);
								typename mesh_t::HHandle h1 = _mesh.opposite_halfedge_handle(h0);
									
								if (_mesh.face_handle(h0).is_valid()) { //证明h0不是边界.
									typename mesh_t::HHandle h0_next = _mesh.next_halfedge_handle(h0);
									typename mesh_t::HHandle h0_prev = _mesh.prev_halfedge_handle(h0); 
									if (DGP::is_zero(_mesh.point(_mesh.to_vertex_handle(h0_next)))) {
										_mesh.property(fp_pos, _mesh.face_handle(h0)) //前提是上面的is_zero正常,否则死循环.
											= (_mesh.point(_mesh.from_vertex_handle(h0_prev)) 
											+ _mesh.point(_mesh.to_vertex_handle(h0))) * 0.25; 
									} else if (DGP::is_zero(_mesh.point(_mesh.from_vertex_handle(h0_prev)))) {
										_mesh.property(fp_pos, _mesh.face_handle(h0)) 
											= (_mesh.point(_mesh.from_vertex_handle(h0)) 
											+ _mesh.point(_mesh.to_vertex_handle(h0_next))) * 0.25; 
									} else { // 没有零坐标的点, 当_n为偶数并且在一个面的所有边都是crease edges时候就可以出现.

										_mesh.property(fp_pos, _mesh.face_handle(h0)) 
											= DGP::face_center(_mesh, _mesh.face_handle(h0)) * 0.5;
									} 
								}  
								if (_mesh.face_handle(h1).is_valid()) { //证明h1不是边界.
									typename mesh_t::HHandle h1_next = _mesh.next_halfedge_handle(h1);
									typename mesh_t::HHandle h1_prev = _mesh.prev_halfedge_handle(h1); 
									if (DGP::is_zero(_mesh.point(_mesh.to_vertex_handle(h1_next)))) {
										_mesh.property(fp_pos, _mesh.face_handle(h1)) //前提是上面的is_zero正常,否则死循环.
											= (_mesh.point(_mesh.from_vertex_handle(h1_prev)) 
											+ _mesh.point(_mesh.to_vertex_handle(h1))) * 0.25; 
									} else if (DGP::is_zero(_mesh.point(_mesh.from_vertex_handle(h1_prev)))) {
										_mesh.property(fp_pos, _mesh.face_handle(h1)) 
											= (_mesh.point(_mesh.from_vertex_handle(h1)) 
											+ _mesh.point(_mesh.to_vertex_handle(h1_next))) * 0.25; 
									} else { // 没有零坐标的点, 当_n为偶数时候就可以出现.
										_mesh.property(fp_pos, _mesh.face_handle(h1)) 
											= DGP::face_center(_mesh, _mesh.face_handle(h1)) * 0.5;
									}   
								}   
							}
							else if (i == 0 && _n == 0) {
								_mesh.property(ep_pos, e_it) *= 0.5;//首先边中点要再缩小一倍.
							}
						}
				} // end of Avf
				 
				// Afv
				for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end()); 
					v_it != v_end; ++v_it) {
						if (_mesh.property(vp_type, v_it) == DGP::SMOOTH_VFT) {
							mesh_t::Point pos(0, 0, 0);
							int count = 0;
							for (typename mesh_t::VFIter vf_it(_mesh, v_it); vf_it; ++vf_it) {
								pos += _mesh.property(fp_pos, vf_it);  
								++count;
							} 
							_mesh.set_point(v_it, pos / count);
						}
				} 
				for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end()); 
					v_it != v_end; ++v_it) {
						if (_mesh.property(vp_type, v_it) == DGP::CREASE_VFT) { 
							std::vector<typename mesh_t::EHandle> e;
							for (typename mesh_t::VEIter ve_it(_mesh, v_it); ve_it; ++ve_it) {
								if (_mesh.status(ve_it.handle()).feature()) e.push_back(ve_it.handle());
							}
							if (e.size() != 2) std::cout << "Error: even Afv, crease vertex only has 2 crease edges.\n";
							_mesh.set_point(v_it, 
								(_mesh.property(ep_pos, e[0]) + _mesh.property(ep_pos, e[1]))*0.5);
						}
				} // end of Afv 
			} // end of  (Afv*Avf)^int(_m/2)
			_mesh.remove_property(ep_pos);
			_mesh.remove_property(fp_pos);
			std::cout << "Even Afv*Avf, ";



			// (Td*Avf)mod(m, 2)
			for (int i = 0; i < _m % 2; ++i) {
				std::cout << "In (Td*Avf)mod(m, 2) loop, i=" << i << ".\n";
				OpenMesh::FPropHandleT<OpenMesh::Vec3d> fp_pos; //这里是自给自足, center of each face.
				_mesh.add_property(fp_pos);
				OpenMesh::EPropHandleT<OpenMesh::Vec3d> ep_pos;//center of each crease edge.
				_mesh.add_property(ep_pos);
				// Avf, 求每一个面的中心点.
				for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end()); 
					f_it != f_end; ++f_it) 
				{
					if (_mesh.property(fp_special, f_it) == false) { //非特殊面的中心点为所有角点的平均.
						//可以猜想Sqrt2得到的面是4个角的, 但是经过Td之后就不是了.
						_mesh.property(fp_pos, f_it) = DGP::face_center(_mesh, f_it);
					} 		
				}  
				for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end());
					e_it != e_end; ++e_it) {
						if (_mesh.status(e_it).feature()) { //一般边的中点.
							_mesh.property(ep_pos, e_it) = DGP::edge_midpos(_mesh, e_it.handle());

							if (_m == 1 && i == 0 && _n != 1) { // _m=1表示之前的(Afv*Avf)无进行, 由于之前的Sv(2), 第一次Avf要特殊处理.
								_mesh.property(ep_pos, e_it) *= 0.5;//首先边中点要再缩小一倍.
								 
								// 特征边两边的面的中点要特殊处理.
								typename mesh_t::HHandle h0 = _mesh.halfedge_handle(e_it, 0);
								typename mesh_t::HHandle h1 = _mesh.opposite_halfedge_handle(h0);

								if (_mesh.face_handle(h0).is_valid()) { //证明h0不是边界.
									typename mesh_t::HHandle h0_next = _mesh.next_halfedge_handle(h0);
									typename mesh_t::HHandle h0_prev = _mesh.prev_halfedge_handle(h0); 
									if (DGP::is_zero(_mesh.point(_mesh.to_vertex_handle(h0_next)))) {
										_mesh.property(fp_pos, _mesh.face_handle(h0)) //前提是上面的is_zero正常,否则死循环.
											= (_mesh.point(_mesh.from_vertex_handle(h0_prev)) 
											+ _mesh.point(_mesh.to_vertex_handle(h0))) * 0.25; 
									} else if (DGP::is_zero(_mesh.point(_mesh.from_vertex_handle(h0_prev)))) {
										_mesh.property(fp_pos, _mesh.face_handle(h0)) 
											= (_mesh.point(_mesh.from_vertex_handle(h0)) 
											+ _mesh.point(_mesh.to_vertex_handle(h0_next))) * 0.25; 
									} else { // 没有零坐标的点, 当_n为偶数时候就可以出现.
										_mesh.property(fp_pos, _mesh.face_handle(h0)) 
											= DGP::face_center(_mesh, _mesh.face_handle(h0)) * 0.5;
									} 
								}  
								if (_mesh.face_handle(h1).is_valid()) { //证明h1不是边界.
									typename mesh_t::HHandle h1_next = _mesh.next_halfedge_handle(h1);
									typename mesh_t::HHandle h1_prev = _mesh.prev_halfedge_handle(h1); 
									if (DGP::is_zero(_mesh.point(_mesh.to_vertex_handle(h1_next)))) {
										_mesh.property(fp_pos, _mesh.face_handle(h1)) //前提是上面的is_zero正常,否则死循环.
											= (_mesh.point(_mesh.from_vertex_handle(h1_prev)) 
											+ _mesh.point(_mesh.to_vertex_handle(h1))) * 0.25; 
									} else if (DGP::is_zero(_mesh.point(_mesh.from_vertex_handle(h1_prev)))) {
										_mesh.property(fp_pos, _mesh.face_handle(h1)) 
											= (_mesh.point(_mesh.from_vertex_handle(h1)) 
											+ _mesh.point(_mesh.to_vertex_handle(h1_next))) * 0.25; 
									} else { // 没有零坐标的点, 当_n为偶数时候就可以出现.
										_mesh.property(fp_pos, _mesh.face_handle(h1)) 
											= DGP::face_center(_mesh, _mesh.face_handle(h1)) * 0.5;
									}   
								}  
							} else if (_m == 1 && i == 0) {
								_mesh.property(ep_pos, e_it) *= 0.5;//首先边中点要再缩小一倍. 
							} 
						}  
				} // end of Avf  
				
				// Td
				mesh_t dual_mesh; //下面谨记dual_mesh也需要设置2个property,
				// 分别是crease edge, crease vertex .
				// 于是每当加入一个新点时候都要设置其特征性质, 每当加入一个新面时候也设置其边的特征性质.
				OpenMesh::VPropHandleT<DGP::v_feature_type> vp_type_dual;//crease edge feature直接设. 
				dual_mesh.add_property(vp_type_dual); 
				// 首先添加顶点, 并且设置顶点的类型.
				OpenMesh::FPropHandleT<mesh_t::VHandle > fp_dual_vertex; //这个面对于dual_mesh上的哪个顶点.
				_mesh.add_property(fp_dual_vertex);
				OpenMesh::EPropHandleT<mesh_t::VHandle > ep_dual_vertex; //这个边对于dual_mesh上的哪个顶点.
				_mesh.add_property(ep_dual_vertex);
				for (typename mesh_t::FIter f_it(_mesh.faces_begin()), f_end(_mesh.faces_end());
					f_it != f_end; ++f_it) {
						typename mesh_t::VHandle v = dual_mesh.add_vertex(_mesh.property(fp_pos, f_it));
						dual_mesh.property(vp_type_dual, v) = DGP::SMOOTH_VFT;//dual_mesh每加一点都要设置.
						_mesh.property(fp_dual_vertex, f_it) = v;
				}
				for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end());
					e_it != e_end; ++e_it) {
						if (_mesh.status(e_it).feature()) {
							typename mesh_t::VHandle v = dual_mesh.add_vertex(_mesh.property(ep_pos, e_it));
							dual_mesh.property(vp_type_dual, v) = DGP::CREASE_VFT;//dual_mesh每加一点都要设置.
							_mesh.property(ep_dual_vertex, e_it) = v; 
						}
				}
				_mesh.remove_property(fp_pos); //在往dual_mesh中加完顶点之后这两个位置信息就不再需要了.
				_mesh.remove_property(ep_pos);
				// 再添加面
				for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end());
					v_it != v_end; ++v_it) {
						if (_mesh.property(vp_type, v_it) == DGP::SMOOTH_VFT) {
							std::vector<typename mesh_t::VHandle> vec; vec.clear();
							typename mesh_t::VOHIter voh_it(_mesh, v_it);
							typename mesh_t::HHandle next = _mesh.ccw_rotated_halfedge_handle(voh_it);
							while (next != voh_it.handle()) {
								vec.push_back(_mesh.property(fp_dual_vertex, _mesh.face_handle(next)));
								next = _mesh.ccw_rotated_halfedge_handle(next);
							}
							vec.push_back(_mesh.property(fp_dual_vertex, _mesh.face_handle(voh_it)));  
							typename mesh_t::FHandle new_f = dual_mesh.add_face(vec);  
							for (typename mesh_t::FEIter fe_it(dual_mesh, new_f); fe_it; ++fe_it) {//每加一个面都要设置一下.
								dual_mesh.status(fe_it.handle()).set_feature(false);
							}    
						}
				}  
				for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end());
					v_it != v_end; ++v_it) {
						if (_mesh.property(vp_type, v_it) == DGP::CREASE_VFT ) {
							// 找两个特征半边
							//                  e_arr[0]
							//<-----------v_it----------> 如果是边界的话确保e_arr[0]是内部的半边.
							//	e_arr[1]                  每一个边界顶点有且只有一个出半边是在外部的.
							std::vector<typename mesh_t::HHandle> e_arr;
							for (typename mesh_t::VOHIter vhe_it(_mesh, v_it); vhe_it; ++vhe_it) {
								if (_mesh.status(_mesh.edge_handle(vhe_it.handle())).feature()) 
									e_arr.push_back(vhe_it.handle());
							}
							if (e_arr.size() != 2) std::cout << "Error: e_arr should have 2 crease edges.\n";

							if (_mesh.face_handle(e_arr[0]).is_valid() == false && _mesh.face_handle(e_arr[1]).is_valid()) {
								typename mesh_t::HHandle temp = e_arr[0];
								e_arr[0] = e_arr[1]; e_arr[1] = temp; //保证e_arr[0]是内部的, 而e_arr[1]可能是外部的.
							}
							// 这个特征点v_it有两种情况，1是像普通点一样需要去掉的，不加入到dual_mesh中。
							// 2.很特殊不能删去，并加入到dual_mesh中成为一角，这出现在v_it的两个特征边直接交一个角.
							// 可能是小角点
							//if (_mesh.face_handle(e_arr[0]) == _mesh.face_handle(_mesh.opposite_halfedge_handle(e_arr[1]))
							//	|| _mesh.face_handle(e_arr[1]) == _mesh.face_handle(_mesh.opposite_halfedge_handle(e_arr[0]))) 
							//{
							//	//  /|
							//	// e1|   内部，特殊的当v_it是边界角同时也只有两个边,都是特征边.
							//	//  v_it---->e0
							//	// 自己也要加入到dual mesh中成为其一个点.
							//	typename mesh_t::VHandle v_dual; 
							//	if (_m == 1) { //因为当_m==1时候之前的(Afv*Avf)^int(_m,2)并没有起效,这样要人为的消去Sv(2)遗留的效果.
							//		v_dual = dual_mesh.add_vertex(_mesh.point(v_it) * 0.5);
							//	} else v_dual = dual_mesh.add_vertex(_mesh.point(v_it));
							//	dual_mesh.property(vp_type_dual, v_dual) = DGP::CREASE_VFT;//dual_mesh每加一点都要设置.
							//	 
							//	std::vector<typename mesh_t::VHandle> vec; vec.clear();
							//	vec.push_back(v_dual);
							//	vec.push_back(_mesh.property(ep_dual_vertex, _mesh.edge_handle(e_arr[0])));
							//	typename mesh_t::HHandle h = e_arr[0];
							//	while (h != e_arr[1]) {
							//		vec.push_back(_mesh.property(fp_dual_vertex, _mesh.face_handle(h)));  
							//		h = _mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(h));
							//	}
							//	vec.push_back(_mesh.property(ep_dual_vertex, _mesh.edge_handle(e_arr[1])));
							//	typename mesh_t::FHandle new_f = dual_mesh.add_face(vec);   
							//	for (typename mesh_t::FEIter fe_it(dual_mesh, new_f); fe_it; ++fe_it) {//每加一个面都要设置一下.
							//		dual_mesh.status(fe_it.handle()).set_feature(false);
							//	} 
							//	dual_mesh.status(dual_mesh.edge_handle(dual_mesh.find_halfedge(vec[0], vec[1]))).set_feature(true);
							//	dual_mesh.status(dual_mesh.edge_handle(dual_mesh.find_halfedge(vec[0], vec[vec.size()-1]))).set_feature(true);

							//	if (_mesh.face_handle(e_arr[1]).is_valid()) { //加入背面，从e_arr[1] ...往e_arr[0]这一面.

							//		std::vector<typename mesh_t::VHandle> vec; vec.clear();
							//		vec.push_back(v_dual);
							//		vec.push_back(_mesh.property(ep_dual_vertex, _mesh.edge_handle(e_arr[1])));
							//		typename mesh_t::HHandle h = e_arr[1];
							//		while (h != e_arr[0]) {
							//			vec.push_back(_mesh.property(fp_dual_vertex, _mesh.face_handle(h)));  
							//			h = _mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(h));
							//		}
							//		vec.push_back(_mesh.property(ep_dual_vertex, _mesh.edge_handle(e_arr[0])));
							//		typename mesh_t::FHandle new_f = dual_mesh.add_face(vec);  
							//		for (typename mesh_t::FEIter fe_it(dual_mesh, new_f); fe_it; ++fe_it) {//每加一个面都要设置一下.
							//			dual_mesh.status(fe_it.handle()).set_feature(false);
							//		} 
							//		dual_mesh.status(dual_mesh.edge_handle(dual_mesh.find_halfedge(vec[0], vec[1]))).set_feature(true);
							//		dual_mesh.status(dual_mesh.edge_handle(dual_mesh.find_halfedge(vec[0], vec[vec.size()-1]))).set_feature(true);
							//	}
							//}  
							//else { //这是第一种情况，v_it需要被去掉的.
								if (_mesh.face_handle(e_arr[0]).is_valid()) { // 先加入一个半面, this test is for boundary.
									typename mesh_t::HHandle h = e_arr[0]; //其实上面已经保证了e_arr[0]不会是边界的半边了.
									std::vector<typename mesh_t::VHandle> vec; vec.clear();
									vec.push_back(_mesh.property(ep_dual_vertex, _mesh.edge_handle(e_arr[0])));
									while (_mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(h)) != e_arr[1]) 
									{
										vec.push_back(_mesh.property(fp_dual_vertex, _mesh.face_handle(h)));
										h = _mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(h));
									}
									vec.push_back(_mesh.property(fp_dual_vertex, _mesh.face_handle(h))); 
									vec.push_back(_mesh.property(ep_dual_vertex, _mesh.edge_handle(e_arr[1])));   
									typename mesh_t::FHandle new_f = dual_mesh.add_face(vec);   
									for (typename mesh_t::FEIter fe_it(dual_mesh, new_f); fe_it; ++fe_it) { //每加一个面都要设置一下.
										dual_mesh.status(fe_it.handle()).set_feature(false);
									} 
									dual_mesh.status(dual_mesh.edge_handle(dual_mesh.find_halfedge(vec[0], vec[vec.size()-1]))).set_feature(true);
								}
								if (_mesh.face_handle(e_arr[1]).is_valid()) { // 再加入另一个半面, this test is for boundary.
									typename mesh_t::HHandle h = e_arr[1]; 
									std::vector<typename mesh_t::VHandle> vec; vec.clear();
									vec.push_back(_mesh.property(ep_dual_vertex, _mesh.edge_handle(e_arr[1])));
									while (_mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(h)) != e_arr[0]) 
									{
										vec.push_back(_mesh.property(fp_dual_vertex, _mesh.face_handle(h)));
										h = _mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(h));
									}
									vec.push_back(_mesh.property(fp_dual_vertex, _mesh.face_handle(h))); 
									vec.push_back(_mesh.property(ep_dual_vertex, _mesh.edge_handle(e_arr[0])));   
									typename mesh_t::FHandle new_f = dual_mesh.add_face(vec);   
									for (typename mesh_t::FEIter fe_it(dual_mesh, new_f); fe_it; ++fe_it) { //每加一个面都要设置一下.
										dual_mesh.status(fe_it.handle()).set_feature(false);
									} //虽然这里会冲掉上面一个半面的设置，但是下面那句会保证(e_arr[0],  e_arr[1])这个边还是特征边.
									dual_mesh.status(dual_mesh.edge_handle(dual_mesh.find_halfedge(vec[0], vec[vec.size()-1]))).set_feature(true);
								}
							//}
							 
						}
				} // end of 增加面.
				_mesh.remove_property(fp_dual_vertex);
				_mesh.remove_property(ep_dual_vertex);//在给dual_mesh添加完面之后这个两个记录信息就不再需要了.

				_mesh.remove_property(vp_type); 
				_mesh = dual_mesh;
				_mesh.add_property(vp_type); 
				for (typename mesh_t::EIter e_it(dual_mesh.edges_begin()), e_end(dual_mesh.edges_end());
					e_it != e_end; ++e_it) {
						_mesh.status(e_it).set_feature(dual_mesh.status(e_it).feature());
				}
				for (typename mesh_t::VIter v_it(dual_mesh.vertices_begin()), v_end(dual_mesh.vertices_end());
					v_it != v_end; ++v_it) {
						_mesh.property(vp_type, v_it) = dual_mesh.property(vp_type_dual, v_it);
				} // end of Td

			} //end of (Td*Avf)^mod(_m, 2).
			
			// 
			std::cout << "Test: after even step, " << _mesh.n_vertices() << "v, " << _mesh.n_edges() << "e, " << _mesh.n_faces() << "f.\n";
			// 重新设置需要输出的特征信息.
			_crease_edges_arr.clear();
			_crease_vertices_arr.clear();
			for (typename mesh_t::EIter e_it(_mesh.edges_begin()), e_end(_mesh.edges_end());
				e_it != e_end; ++e_it) {
					if (_mesh.status(e_it).feature())
						_crease_edges_arr.push_back(e_it.handle().idx());
			}  
			for (typename mesh_t::VIter v_it(_mesh.vertices_begin()), v_end(_mesh.vertices_end());
				v_it != v_end; ++v_it) {
					if (_mesh.property(vp_type, v_it) == DGP::CREASE_VFT)
						_crease_vertices_arr.push_back(v_it.handle().idx());
			} 
			std::cout << "Test: end R^(m,n), crease vertex and crease edge 's size: "
				<< _crease_vertices_arr.size() << ", " << _crease_edges_arr.size() << ".\n";
			_mesh.remove_property(vp_type);
			// after R^(_m, _n)带特征的处理.
			_mesh.update_normals();
			return true;
		} // end of function unified_subdivide_features().
	protected: 
	private: // topological modifiers
		//将一个特征边分裂为两个特征边, 加入的点是crease vertice, 加入的新边也是crease edge, 两边的面是special face.
		void split_edge_add_feature(mesh_t& _m, 
			OpenMesh::VPropHandleT<DGP::v_feature_type> & _vp_type, // 分裂出来的顶点是crease vertex
			OpenMesh::FPropHandleT<bool> & _fp_special, //分裂出来的面也是special face 
			const typename mesh_t::EdgeHandle& _eh) //这个特征边分裂
		{
			// 这个函数是Loop, CC, 等细分实现中的split_edge()函数的修改版.
			typename mesh_t::HalfedgeHandle heh     = _m.halfedge_handle(_eh, 0);
			typename mesh_t::HHandle opp_heh = _m.halfedge_handle(_eh, 1);

			if (_m.face_handle(heh).is_valid() == false
				&& _m.face_handle(opp_heh).is_valid()) { // confirm the heh's face is valid.
					typename mesh_t::HHandle temp = heh;
					heh = opp_heh;
					opp_heh = temp;
			}

			typename mesh_t::HalfedgeHandle new_heh, opp_new_heh, t_heh;
			typename mesh_t::VertexHandle   vh;
			typename mesh_t::VertexHandle   vh1(_m.to_vertex_handle(heh));
			typename mesh_t::Point          zero(0,0,0);

			// new vertex //这里显示了怎么加入一个新点并修改其拓扑连接。
			vh = _m.new_vertex( zero );

			// Re-link mesh entities
			if (_m.is_boundary(_eh))
			{
				for (t_heh = heh;
					_m.next_halfedge_handle(t_heh) != opp_heh;
					t_heh = _m.opposite_halfedge_handle(_m.next_halfedge_handle(t_heh)))
				{}
			}
			else
			{
				for (t_heh = _m.next_halfedge_handle(opp_heh);
					_m.next_halfedge_handle(t_heh) != opp_heh;
					t_heh = _m.next_halfedge_handle(t_heh) ) 
				{} // 此循环对tri/quad都有效.
			}

			new_heh     = _m.new_edge(vh, vh1);
			opp_new_heh = _m.opposite_halfedge_handle(new_heh);
			_m.set_vertex_handle( heh, vh );

			_m.set_next_halfedge_handle(t_heh, opp_new_heh);  // 接下面的半边
			_m.set_next_halfedge_handle(opp_new_heh, opp_heh);// 接下面的半边
			_m.set_next_halfedge_handle(new_heh, _m.next_halfedge_handle(heh)); // 接上面的半边
			_m.set_next_halfedge_handle(heh, new_heh);							// 接上面的半边

			if (_m.face_handle(opp_heh).is_valid())
			{
				_m.set_face_handle(opp_new_heh, _m.face_handle(opp_heh));
				_m.set_halfedge_handle(_m.face_handle(opp_new_heh), opp_new_heh);
			}

			_m.set_face_handle( new_heh, _m.face_handle(heh) );
			_m.set_halfedge_handle( vh, new_heh);
			_m.set_halfedge_handle( _m.face_handle(heh), heh );
			_m.set_halfedge_handle( vh1, opp_new_heh );

			// Never forget this, when playing with the topology
			_m.adjust_outgoing_halfedge( vh );
			_m.adjust_outgoing_halfedge( vh1 );

			_m.status(_m.edge_handle(new_heh)).set_feature(true);
			if (_vp_type.is_valid() == false) std::cout << "Error: _vp_type";
			if (_fp_special.is_valid() == false ) std::cout << "Error: _fp_special";
			_m.property(_vp_type, vh) = DGP::CREASE_VFT;
			_m.property(_fp_special, _m.face_handle(heh)) = true;  
			if (_m.face_handle(opp_heh).is_valid()) {
				_m.property(_fp_special, _m.face_handle(opp_heh)) = true; 
			}			 
		} // end of split_edge_add_feature

		// 分别将一个面分裂为两个面, 顶点都在了，主要加入相应的连接.
		void split_face_1to2(mesh_t &_mesh, typename mesh_t::HHandle &_h0, typename mesh_t::HHandle &_h1) {
			//
			//            * 2
			//           /|\
			//          / | \
			//      h3 /  |  \ h2
			//        / h4|h5 \
			//       /_   |/   \
			//    3 *---->*---->* 1 
			//        h0  0  h1
			typename mesh_t::HHandle h0 = _h0, h1 = _h1;
			typename mesh_t::HHandle h2 = _mesh.next_halfedge_handle(h1);
			typename mesh_t::HHandle h3 = _mesh.prev_halfedge_handle(h0);
			if (h2 == h3) {
				if (DGP::n_face_corners(_mesh, _mesh.face_handle(h0)) == 3) return;
				else std::cerr << "Error: 这应该是三角形的情况, 不需要分裂.\n";
			}
			typename mesh_t::VHandle v0 = _mesh.to_vertex_handle(h0);
			typename mesh_t::VHandle v2 = _mesh.to_vertex_handle(h2);
			typename mesh_t::HHandle h4 = _mesh.new_edge(v0, v2);
			typename mesh_t::HHandle h5 = _mesh.opposite_halfedge_handle(h4);
			typename mesh_t::FHandle f_old = _mesh.face_handle(h0);
			typename mesh_t::FHandle f_new = _mesh.new_face();

			// Re-Set Handles around old Face
			_mesh.set_next_halfedge_handle(h0, h4);
			_mesh.set_next_halfedge_handle(h4, h3);
			_mesh.set_face_handle(h4, f_old);
			_mesh.set_halfedge_handle(f_old, h4);
			// Re-Set Handles around new Face
			_mesh.set_next_halfedge_handle(h2, h5);
			_mesh.set_next_halfedge_handle(h5, h1);
			_mesh.set_face_handle(h1, f_new);
			_mesh.set_face_handle(h2, f_new);
			_mesh.set_face_handle(h5, f_new);
			_mesh.set_halfedge_handle(f_new, h5);

			_mesh.status(_mesh.edge_handle(h4)).set_feature(false);
			_mesh.status(_mesh.edge_handle(h4)).set_tagged(false);
		}
	private: // geometry helper

	private: // data 
		
	}; // end of class UnifiedSubd24dbsT
 
} // end of namespace DGP
#endif //dgp_unifiedSubd24dbs