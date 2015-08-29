#ifndef dgp_studio_common_qem_decimation_h
#define dgp_studio_common_qem_decimation_h

#include "QuadricT.h"
#include "VertexFeatureType.h"
#include "OpenMeshAll.h"

namespace DGP {
	 

template <typename TriMeshT>
class QEMDecimationT {

public:
	//计算每一个顶点的quadric 
	void init_calc_quadric(TriMeshT &_mesh, OpenMesh::VPropHandleT<Quadricd> & _vp_quadric , bool _area_weighted_error) 
	{	
		if (_vp_quadric.is_valid() == false) {
			std::cerr << "Error: DGP::QEMDecimation::init_calc_quadric(), _vp_quadric invalid.\n";
			return; //这是因为在调用端还没有_mesh.add_property(_vp_quadric);所造成的.
		}
		// 利用面的法向量, 求出计算每一个顶点的Quadric
		_mesh.update_face_normals();// compute face normals

		if (_area_weighted_error == false) {
			// 用顶点来做循环
			TriMesh::VertexIter  v_it, v_end(_mesh.vertices_end());//v_end = _mesh.vertices_end();

			for (v_it=_mesh.vertices_begin(); v_it != v_end; ++v_it) // 用顶点来做循环求每一个顶点的Quadric
			{
				TriMesh::Point       n;
				double				a, b, c, d;

				TriMesh::VertexHandle _vh(v_it.handle());

				_mesh.property(_vp_quadric, _vh).clear();

				if (_mesh.is_boundary(v_it)) { //是边界点
					// calc vertex quadrics from incident triangles 用这个顶点四周的面的信息,所以用到circulator
					//for (TriMesh::VFIter vf_it(_mesh, _vh); vf_it; ++vf_it)
					for (TriMesh::VFIter vf_it = _mesh.vf_iter(_vh); vf_it; ++vf_it)
					{
						// plane equation // {x-x0, y-y0, z-z0}dot{n0, n1, n2} = 0. 
						n = _mesh.normal(vf_it);//前提是已经_mesh.request_face_normals()和_mesh.update_face_normals();//这样获得的法向量是已经单位化的了.

						a = n[0];
						b = n[1];
						c = n[2];
						d = -(n | _mesh.point(_vh));//dot product

						// plane -> quadric, 这是原来的做法, 好像没有使用area-weighted quadric
						_mesh.property(_vp_quadric, _vh) += DGP::Quadricd(a, b, c, d);//将这个面对应的fundamental quadric加起来	

						// ////////////////////////////////////////////////////////////////////////////////////
						// 假如这个面是边界面, 生成一个垂直的面, 这是在QEM原论文的第6节中的Preserving Boundaries.
						TriMesh::EdgeHandle boundary_edge_handle;
						for (TriMesh::ConstFaceEdgeIter cfeit = _mesh.cfe_iter( vf_it ); cfeit; ++cfeit) {
							//  //这个面只要有一个边是边界就将这个面当作是在边界上.而且一个面可能不止只有一个边是边界边(对三角形来说可能是最多两个边是边界边)
							if (_mesh.is_boundary( cfeit.handle() ) ) {
								boundary_edge_handle = cfeit.handle();
								TriMesh::HalfedgeHandle heh = _mesh.halfedge_handle(boundary_edge_handle, 0);
								if (_mesh.is_boundary(heh)) {
									heh = _mesh.opposite_halfedge_handle(heh);
								}
								TriMesh::VHandle v0 = _mesh.from_vertex_handle(heh);
								TriMesh::VHandle v1 = _mesh.to_vertex_handle(heh);
								TriMesh::VHandle v2 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh));
								TriMesh::Point p0 = _mesh.point(v0), p1 = _mesh.point(v1), p2 = _mesh.point(v2);

								TriMesh::Point p5 = p0 + (p1 - p0).normalize() * ((p2 - p0)|(p1 - p0).normalize());
								TriMesh::Point p2_mirror = p5 + (p2 - p5).norm() * n;

								OpenMesh::Vec3d n_mirror = (p2_mirror -p0) % (p1 -p0);
								double length = n_mirror.norm(); //长度
								if (length > FLT_MIN)
								{
									n_mirror /= length;//单位化								
								}

								a = n_mirror[0];
								b = n_mirror[1];
								c = n_mirror[2];
								d = -(OpenMesh::vector_cast<OpenMesh::Vec3d>(p0)|n_mirror);

								_mesh.property(_vp_quadric, _vh) += DGP::Quadricd(a, b, c, d);
								// break;
							}/**/
						} // end of for each edge
						/////////////////////////////////////// end of 在QEM原论文的第6节中的Preserving Boundaries.
					} // end of for each face 
				} else { //这个点不是边界点
					// calc vertex quadrics from incident triangles 用这个顶点四周的面的信息,所以用到circulator
					//for (TriMesh::VFIter vf_it(_mesh, _vh); vf_it; ++vf_it)
					for (TriMesh::VFIter vf_it = _mesh.vf_iter(_vh); vf_it; ++vf_it)
					{
						// plane equation // {x-x0, y-y0, z-z0}dot{n0, n1, n2} = 0. 
						n = _mesh.normal(vf_it);//前提是已经_mesh.request_face_normals()和_mesh.update_face_normals();//这样获得的法向量是已经单位化的了.

						a = n[0];
						b = n[1];
						c = n[2];
						d = -(n | _mesh.point(_vh));//dot product

						// plane -> quadric, 这是原来的做法, 好像没有使用area-weighted quadric
						_mesh.property(_vp_quadric, _vh) += DGP::Quadricd(a, b, c, d);//将这个面对应的fundamental quadric加起来	
					}
				} // end of if-else vertex is at boundary
			} // end of for each vertex

		} else { // _area_weighted_error == true
			// 用面来做循环  参考OpenMesh/Tools/Decimater/ModQuadricT.hh
			// calc (normal weighted) quadric
			TriMesh::FaceIter f_it  = _mesh.faces_begin(), f_end = _mesh.faces_end();

			TriMesh::FaceVertexIter fv_it;
			TriMesh::VertexHandle vh0, vh1, vh2;

			double a,b,c,d, area;

			// use the faces to iterator, to calualate all the vertex's fundermental quadric
			for (; f_it != f_end; ++f_it)
			{	f_it.handle().idx();
			fv_it = _mesh.fv_iter(f_it.handle());
			vh0 = fv_it.handle();  ++fv_it;
			vh1 = fv_it.handle();  ++fv_it;
			vh2 = fv_it.handle();

			OpenMesh::Vec3d v0, v1, v2;
			{
				using namespace OpenMesh;

				v0 = vector_cast<Vec3d>(_mesh.point(vh0));
				v1 = vector_cast<Vec3d>(_mesh.point(vh1));
				v2 = vector_cast<Vec3d>(_mesh.point(vh2));
			}

			OpenMesh::Vec3d n = (v1-v0) % (v2-v0);
			area = n.norm(); //长度
			if (area > FLT_MIN)
			{
				n /= area;//单位化
				area *= 0.5;
			}

			a = n[0];
			b = n[1];
			c = n[2];
			d = -(OpenMesh::vector_cast<OpenMesh::Vec3d>(_mesh.point(vh0))|n); //d表示原点到vh0的距离.

			DGP::Quadricd q(a, b, c, d);
			q *= area;//here uses the area-weighted error metric

			_mesh.property(_vp_quadric, vh0) += q;
			_mesh.property(_vp_quadric, vh1) += q;
			_mesh.property(_vp_quadric, vh2) += q;

			// ////////////////////////////////////////////////////////////////////////////////////
			// 假如这个面是边界面, 生成一个垂直的面, 这是在QEM原论文的第6节中的Preserving Boundaries.
			TriMesh::EdgeHandle boundary_edge_handle;
			for (TriMesh::ConstFaceEdgeIter cfeit = _mesh.cfe_iter( f_it ); cfeit; ++cfeit) {
				//  //这个面只要有一个边是边界就将这个面当作是在边界上.而且一个面可能不止只有一个边是边界边(对三角形来说可能是最多两个边是边界边)
				if (_mesh.is_boundary( cfeit.handle() ) ) {
					boundary_edge_handle = cfeit.handle();
					TriMesh::HalfedgeHandle heh = _mesh.halfedge_handle(boundary_edge_handle, 0);
					if (_mesh.is_boundary(heh)) {
						heh = _mesh.opposite_halfedge_handle(heh);
					}
					TriMesh::VHandle v0 = _mesh.from_vertex_handle(heh);
					TriMesh::VHandle v1 = _mesh.to_vertex_handle(heh);
					TriMesh::VHandle v2 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh));
					TriMesh::Point p0 = _mesh.point(v0), p1 = _mesh.point(v1), p2 = _mesh.point(v2);

					TriMesh::Point p5 = p0 + (p1 - p0).normalize() * ((p2 - p0)|(p1 - p0).normalize());
					TriMesh::Point p2_mirror = p5 + (p2 - p5).norm() * (OpenMesh::Vec3d)n;

					OpenMesh::Vec3d n_mirror = (p2_mirror -p0) % (p1 -p0);
					double length = n_mirror.norm(); //长度
					if (length > FLT_MIN)
					{
						n_mirror /= length;//单位化								
					}

					a = n_mirror[0];
					b = n_mirror[1];
					c = n_mirror[2];
					d = -(OpenMesh::vector_cast<OpenMesh::Vec3d>(p0)|n_mirror);

					_mesh.property(_vp_quadric, v0) += DGP::Quadricd(a, b, c, d);
					_mesh.property(_vp_quadric, v1) += DGP::Quadricd(a, b, c, d);
				} // end of if /**/
			} // end of for each edge 
			////////////////////////////////////// end of 在QEM原论文的第6节中的Preserving Boundaries.
			}
		} // end of if-else area_weighted_error == false/true

	} // end of init_calc_quadric


	// This function can do two job. (1)preserving boundaries by judge whether one 
	// of the halfedge which is supposed to be collapsed, although the original paper
	// do this by weighted by a large penalty factor. This means that halfedge collapses
	// can deal with manifold models with boundaries. (2)preventing mesh inversion,
	// by consider the normal of the faces flips whether or not, before and after edge
	// collapses. More information can be found in the part 6 of original paper Garland97.
	static bool is_collapse_legal(TriMeshT &_mesh, typename TriMeshT::HalfedgeHandle _hh, const OpenMesh::VPropHandleT<DGP::v_feature_type> &_vp_type) 
	{	// (from, to) -> to,也就是(v0, v1)->v1,
		if (_vp_type.is_valid() == false) {
			std::cerr << "Error: DGP::QEMDecimation::is_collapse_legal(), _vp_type invalid.\n";
			return false; //这是因为在调用端还没有_mesh.add_property(_vp_type);所造成的.
		}

		// collect vertices	半边的两个顶点, 若能收缩的话情况将是(from, to)->to
		TriMesh::VertexHandle v0 = _mesh.from_vertex_handle(_hh);
		TriMesh::VertexHandle v1 = _mesh.to_vertex_handle(_hh);

		// collect faces 半边的对应的面
		TriMesh::FaceHandle fl = _mesh.face_handle(_hh);
		TriMesh::FaceHandle fr = _mesh.face_handle(_mesh.opposite_halfedge_handle(_hh));

		// backup point positions 这两个变量的值不被修改
		TriMesh::Point p0 = _mesh.point(v0);
		TriMesh::Point p1 = _mesh.point(v1);

		// (1) topological test // 这半边的两个顶点都没有被deleted,而且都不是在边界上而它的边又不是边界
		if (!_mesh.is_collapse_ok(_hh)) 
			return false;

		// (2) from vertex is at boundary,
		if (_mesh.is_boundary(v0)) {
			if (_mesh.valence(v0) == 2) return false;//(2.3) boundary corner

			if (_mesh.is_boundary(v1) == false) { //(2.1) v0是在边界上, 而v1不是在边界上, 
				// 这是无论如何不可以的, 避免了将boundary往内部收缩(这样的后果之一是洞变大).
				return false; // 也就是说边界边只能沿着边界边来收缩.
			} else { // (2.2) v0是在边界上, v1也是在边界上, 也就是说e(v0, v1)是边界边, 这就讨论v1点的类型
				if (_mesh.property(_vp_type, v0) == DGP::DART_VFT || 
					_mesh.property(_vp_type, v0) == DGP::CORNER_VFT || 
					_mesh.property(_vp_type, v0) == DGP::CREASE_VFT) { 
					//实际上v0只能是smooth点, 是dart/corner当然不能动, 
					//是creast呢?也不行,原因是边界边不是crease边,所以如果v0是crease的话,v1只可能内部点.
					return false;													
				} 
			}
		} else { // (3) from vertex is not at boundary, 现在from点是inner vertex内部点了.

			if (_mesh.property(_vp_type, v0) == DGP::DART_VFT || _mesh.property(_vp_type, v0) == DGP::CORNER_VFT) {// (3.1) from vertex is dart和corner的顶点的是不予被收缩掉的
				return false;
			} else if (_mesh.property(_vp_type, v0) == DGP::CREASE_VFT) {//(3.2) crease vertex只能沿着crease edge收缩
				if (_mesh.status(_mesh.edge_handle(_hh)).feature() == false) { //_hh不是特征边
					return false;
				} else { //_hh是特征边了.
					if (_mesh.status(_mesh.edge_handle(_mesh.prev_halfedge_handle(_hh))).feature() == true
						|| _mesh.status(_mesh.edge_handle(_mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(_hh)))).feature() == true)
					{
						//std::cout << "Crease Angle.\n";
						return false;					
					}
				}

			} // (3.3) v0是smooth vertex
		}

		// test normal flipping: 判断法向量会否被翻转, 从而避免三角形的overlapping
		//   if normal vector of a (non-degenerate) triangle changes by 
		//   more than pi degrees, return false.
		//std::cout << v0 << ", " << v1 << std::endl;//for test
		// simulate collapse
		_mesh.set_point(v0, p1);

		// check for flipping normals
		//double max_deviation = (double)(90.0 / 180.0 * M_PI); //默认情况下是_f = 90, 那么max_deviation = 0.5*M_PI //以90度为界限
		double min_cos = 1.0e-6;//0.0;//= cos(max_deviation); //为了计算简单,min_cos简单可以用0来表示吧 //const float min_cos = (float) cos(0.25*M_PI);
		double c(1.0);
		bool zero_area = false;
		bool mini_angle = false;
		for (TriMesh::ConstVertexFaceIter vf_it(_mesh.cvf_iter(v0)); vf_it; ++vf_it) 
		{
			if (vf_it.handle() != fl && vf_it.handle() != fr)
			{
				// (1) 判断是否有三角形翻转.
				TriMesh::Point n0 = _mesh.normal(vf_it); //这是没有修改v0点坐标前的面的法向, 下面从新计算后就变了.
				TriMesh::Point n1 = _mesh.calc_face_normal(vf_it); //注意这方法求出来的normal都是经过归一化处理的, 长度都是1.
				//std::cout << "n0: " << n0 << ", n1: " << n1 << std::endl;

				c = dot(n0, n1); //std::cout << "c: " << c << std::endl;// c = n0 | n1;
				if (c < min_cos)
					break; //只要from的一邻域中有一个面反转了也是不和要求的.

				// (2) 这里其实可以加入判断,要是n=(0, 0, 0)的话,就否决.也就是不许三点共线.
				// -----自己求一次法向,用来表示面积的 
				TriMesh::ConstFaceVertexIter fv_it(_mesh.cfv_iter(vf_it.handle()));

				const TriMesh::Point p0(_mesh.point(fv_it));  ++fv_it;
				const TriMesh::Point p1(_mesh.point(fv_it));  ++fv_it;
				const TriMesh::Point p2(_mesh.point(fv_it));
				TriMesh::Normal p1p0(p0);  p1p0 -= OpenMesh::vector_cast<TriMesh::Normal>(p1);
				TriMesh::Normal p1p2(p2);  p1p2 -= OpenMesh::vector_cast<TriMesh::Normal>(p1);

				TriMesh::Normal n    = cross(p1p2, p1p0);			
				////if ((fabs(n1[0]) < 1.e-0) || (fabs(n1[1]) < 1.e-0) || (fabs(n1[2]) < 1.e-0)) { zero_area = true; break; }			
				double len = n.norm();
				if (fabs(len) < 1.e-6) { 
					zero_area = true; //std::cout << "zero area ";
					break;
				}
				// (3)阻止生成太尖的triangle
				double min_angle = 6*M_PI/180, max_angle = 174*M_PI/180;
				double tmp = acos(dot(p1p0, p1p2)/(p1p0.norm() * p1p2.norm()));
				if ( tmp < min_angle || tmp > max_angle) { mini_angle = true; break; }
				tmp = acos(dot(p1p0, p0-p2)/(p1p0.norm() * (p0-p2).norm()));
				if ( tmp < min_angle || tmp > max_angle) { mini_angle = true; break; }
				tmp = acos(dot(p1p2, p2-p0)/(p1p2.norm() * (p2-p0).norm()));
				if ( tmp < min_angle || tmp > max_angle) { mini_angle = true; break; }
			}
		}

		// undo simulation changes
		_mesh.set_point(v0, p0);

		if (c < min_cos || zero_area == true || mini_angle == true) { //std::cout << "false.\n";//for test
			return false; //if (c < 0.0) return false;
		}

		/*
		int valences = (_mesh.valence(v0) + _mesh.valence(v1));
		if (_mesh.is_boundary(_hh) || _mesh.is_boundary(_mesh.opposite_halfedge_handle(_hh))) {
		valences -= 3;
		} else valences -= 4;
		if (valences > 10) return false;*/

		// collapse passed all tests -> ok
		return true;
	}


}; // class QEMDecimationT


} // namespace DGP
#endif //dgp_studio_common_qem_decimation_h