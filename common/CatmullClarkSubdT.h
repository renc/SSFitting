// ***************************************************************
//  
//  CatmullClarkSubdT 
//  Copyright (C) 2008 - by rencanjiang All Rights Reserved
//  rencanjiang@163.com	
//	-------------------------------------------------------------
//  
//  data: 2009-01-07
// ***************************************************************
//  
//  File description: 1. ģ�� LoopT.hh �� MidpointLinearSubdT.h��ʵ�ֵ�. 
//  2. ע����������в���ֹ�����ı��λ�����������(������arbitrary polygonal)��
//	��Ҫ���ø�general��CC rules����һ��ϸ�ֵõ����ж���quadrilateral��mesh.
//  3. ���ԶԺ��б߽��ģ�ͽ���ϸ��, ���Ǳ߽��ϵ�ϸ��ϵ��mark�����м�����ͬ�ģ�Ҫ�����������. 
//            *                   *
//           / \                 / \
//          /   \               /   \
//         /     \   ====>     *     *
//        /       \           /  \*/  \
//       /         \         /    |    \
//      *-----------*       *-----*-----* 
//      |           |       |     |     | 
//      |           |       |     |     | 
//      |           |       *-----*-----* 
//      |           |       |     |     | 
//      |           |       |     |     | 
//      |___________|       *-----*-----* 
// ***************************************************************
#ifndef dgp_catmullclarksubd_h
#define dgp_catmullclarksubd_h

#include "OpenMeshAll.h"

namespace DGP {

	template <typename MeshType, typename RealType = float>
	class CatmullClarkSubdT : public OpenMesh::Subdivider::Uniform::SubdividerT<MeshType, RealType>
	{
	public:

		typedef RealType                                real_t;
		typedef MeshType								mesh_t;
		typedef SubdividerT< mesh_t, real_t >           parent_t;

	public:

		CatmullClarkSubdT(void) : parent_t()
		{ }
		~CatmullClarkSubdT() {}
	public:
		const char *name() const { return "Catmull Clark subdi for Quad"; }

		// ����������¼���ģ�����siggraph2000 course notes 'subdivision for modeling and animation':
		// CC schemes are defined for quad mesh. Arbitrary polygonal meshes can be 
		// reduced to a quadrilateral mesh using a more general form of Catmull-Clark rules.
		// Notes: this function should not be invoked before the attached(..) function.
		bool reduce_to_quadmesh( mesh_t &_m) {
			// 1. a face control point for an n-gon is computed as the average of the corners.
			for (typename mesh_t::FIter f_it(_m.faces_begin()), f_end(_m.faces_end());
				f_it != f_end; ++f_it) {
					compute_poly_center(_m, f_it.handle());
			}
			// 2. an edge control point as the average of the endpoints of the edge
			// and newly computed face control points of adjacent faces.
			for (typename mesh_t::EIter e_it(_m.edges_begin()), e_end(_m.edges_end()); 
				e_it != e_end; ++e_it) {
					compute_edge_control_point(_m, e_it.handle());
			}
			// 3. even(old) control points update.
			for (typename mesh_t::VIter v_it(_m.vertices_begin()), v_end(_m.vertices_end());
				v_it != v_end; ++v_it) {
					compute_vertex_control_point(_m, v_it.handle());
			}

			// change the topology.
			for (typename mesh_t::EIter e_it(_m.edges_begin()), e_end(_m.edges_end()); 
				e_it != e_end; ++e_it)
				split_edge(_m, e_it.handle() );
			for (typename mesh_t::FIter f_it(_m.faces_begin()), f_end(_m.faces_end());
				f_it != f_end; ++f_it)
				split_face_quad(_m, f_it.handle());
			// update the geometry.
			for (typename mesh_t::VIter v_it(_m.vertices_begin()), v_end(_m.vertices_end());
				v_it != v_end; ++v_it)
				_m.set_point(v_it, _m.property( vp_pos_, v_it));

			return true; 
		} // end of function reduce_to_quadmesh();
	protected:
		bool prepare( mesh_t& _m )
		{
			_m.add_property( vp_pos_ );// for tri/quad.
			_m.add_property( ep_pos_ );// for tri/quad.
			_m.add_property( fp_pos_ );// for quad.
			return true;
		}
		bool cleanup( mesh_t& _m )
		{
			_m.remove_property( vp_pos_ );
			_m.remove_property( ep_pos_ );
			_m.remove_property( fp_pos_ );// for quad.
			return true;
		}

		bool subdivide( mesh_t& _m, size_t _n)
		{
			typename mesh_t::FaceIter   fit, f_end;
			typename mesh_t::EdgeIter   eit, e_end;
			typename mesh_t::VertexIter vit;

			// Do _n subdivisions
			for (size_t i=0; i < _n; ++i)
			{
				// compute new positions for old vertices, ��tri/quad ����Ч.
				for ( vit  = _m.vertices_begin(); vit != _m.vertices_end(); ++vit)
					smooth( _m, vit.handle() );//���ÿһ���ɵĶ��㱻���º��������.


				// Compute position for new vertices and store them in the edge property, ��tri/quad����Ч.
				for (eit=_m.edges_begin(), e_end = _m.edges_end(); eit != e_end; ++eit)
					compute_midpoint( _m, eit.handle() );//���ÿһ�����������ӵĵ������

				// �����ı�����˵������Ҫ����ÿһ�����¼�����е�. ֻ��quad��Ч.
				for (fit = _m.faces_begin(), f_end = _m.faces_end(); fit != f_end; ++fit) {
					//if (n_face_corners(_m, fit.handle()) == 4) 
						compute_quad_center(_m, fit.handle());
				}

				// Split each edge at midpoint and store precomputed positions (stored in
				// edge property ep_pos_) in the vertex property vp_pos_;

				// Attention! Creating new edges, hence make sure the loop ends correctly. ��tri/quad����Ч.
				for (eit=_m.edges_begin(), e_end = _m.edges_end(); eit != e_end; ++eit)
					split_edge(_m, eit.handle() );
				
				// Commit changes in topology and reconsitute consistency

				// Attention! Creating new faces, hence make sure the loop ends correctly.
				for (fit = _m.faces_begin(), f_end   = _m.faces_end(); fit != f_end; ++fit)
					if (n_face_corners(_m, fit.handle()) == 8) split_face_quad(_m, fit.handle());// for quad.
					//else if (n_face_corners(_m, fit.handle()) == 6) split_face(_m, fit.handle() ); // for tri.
					else std::cout << "Error: CatmullClarkSubdT, split edge error.\n";


				// Commit changes in geometry
				for ( vit  = _m.vertices_begin(); vit != _m.vertices_end(); ++vit)
					_m.set_point(vit, _m.property( vp_pos_, vit ) );

			} // end of for

			return true;
		} // end of function subdivide

	private: // topological modifiers

		void split_face(mesh_t& _m, const typename mesh_t::FaceHandle& _fh)
		{
			typename mesh_t::HalfedgeHandle
				heh1(_m.halfedge_handle(_fh)),
				heh2(_m.next_halfedge_handle(_m.next_halfedge_handle(heh1))),
				heh3(_m.next_halfedge_handle(_m.next_halfedge_handle(heh2)));

			// Cutting off every corner of the 6_gon
			corner_cutting( _m, heh1 );
			corner_cutting( _m, heh2 );
			corner_cutting( _m, heh3 );
		}

		void corner_cutting(mesh_t& _m, const typename mesh_t::HalfedgeHandle& _he)
		{
			// Define Halfedge Handles
			typename mesh_t::HalfedgeHandle
				heh1(_he),
				heh5(heh1),
				heh6(_m.next_halfedge_handle(heh1));

			// Cycle around the polygon to find correct Halfedge
			for (; _m.next_halfedge_handle(_m.next_halfedge_handle(heh5)) != heh1;
				heh5 = _m.next_halfedge_handle(heh5))
			{}

			typename mesh_t::VertexHandle
				vh1 = _m.to_vertex_handle(heh1),
				vh2 = _m.to_vertex_handle(heh5);

			typename mesh_t::HalfedgeHandle
				heh2(_m.next_halfedge_handle(heh5)),
				heh3(_m.new_edge( vh1, vh2)),
				heh4(_m.opposite_halfedge_handle(heh3));

			// Intermediate result
			//
			//            *
			//         5 /|\
			//          /_  \
			//    vh2> *     *
			//        /|\3   |\
			//       /_  \|4   \
			//      *----\*----\*
			//          1 ^   6
			//            vh1 (adjust_outgoing halfedge!)
			//

			// Old and new Face             //����ܺõ���ʾ����ô���һ�����棬�Ծ�����ô�����˵��޸�.
			typename mesh_t::FaceHandle     fh_old(_m.face_handle(heh6));
			typename mesh_t::FaceHandle     fh_new(_m.new_face());


			// Re-Set Handles around old Face
			_m.set_next_halfedge_handle(heh4, heh6);
			_m.set_next_halfedge_handle(heh5, heh4);

			_m.set_face_handle(heh4, fh_old);
			_m.set_face_handle(heh5, fh_old);
			_m.set_face_handle(heh6, fh_old);
			_m.set_halfedge_handle(fh_old, heh4);

			// Re-Set Handles around new Face
			_m.set_next_halfedge_handle(heh1, heh3);
			_m.set_next_halfedge_handle(heh3, heh2);

			_m.set_face_handle(heh1, fh_new);
			_m.set_face_handle(heh2, fh_new);
			_m.set_face_handle(heh3, fh_new);

			_m.set_halfedge_handle(fh_new, heh1);
		}

		void split_face_quad(mesh_t& _m, typename mesh_t::FaceHandle& _fh)
		{ 
		typename mesh_t::HalfedgeHandle
			h1(_m.halfedge_handle(_fh)),
			h2(_m.next_halfedge_handle(_m.next_halfedge_handle(h1))),
			h3(_m.next_halfedge_handle(_m.next_halfedge_handle(h2))),
			h4(_m.next_halfedge_handle(_m.next_halfedge_handle(h3))),
			h5(_m.next_halfedge_handle(h1)),
			h6(_m.next_halfedge_handle(h2)),
			h7(_m.next_halfedge_handle(h3)),
			h8(_m.next_halfedge_handle(h4));
		if (h1 != _m.next_halfedge_handle(_m.next_halfedge_handle(h4))) std::cerr << "Error: CatmullClarkSubdT split_face_quad, halfedges.\n";
		typename mesh_t::VertexHandle
			v1(_m.to_vertex_handle(h1)), v2(_m.to_vertex_handle(h2)), 
			v3(_m.to_vertex_handle(h3)), v4(_m.to_vertex_handle(h4));

		// add new vertex (face vertex)
		typename mesh_t::VHandle vh = _m.new_vertex(typename mesh_t::Point(0, 0, 0));
		_m.property(vp_pos_, vh) = _m.property(fp_pos_, _fh);  

		// �ȼ���4���ߣ������ư��֮������ӹ�ϵ.
		typename mesh_t::HHandle 
			v1h = _m.new_edge(vh, v1), v2h = _m.new_edge(vh, v2),
			v3h = _m.new_edge(vh, v3), v4h = _m.new_edge(vh, v4);
		_m.set_next_halfedge_handle(v1h, h5); // ��ʱ�����޸�ǰh1����һ����߻���ԭ����.
		_m.set_next_halfedge_handle(h1, _m.opposite_halfedge_handle(v1h));//��ʱh1����һ����߾ͱ��޸���.
		_m.set_next_halfedge_handle(v2h,h6);  
		_m.set_next_halfedge_handle(h2, _m.opposite_halfedge_handle(v2h)); 
		_m.set_next_halfedge_handle(v3h, h7);  
		_m.set_next_halfedge_handle(h3, _m.opposite_halfedge_handle(v3h)); 
		_m.set_next_halfedge_handle(v4h, h8);  
		_m.set_next_halfedge_handle(h4, _m.opposite_halfedge_handle(v4h));

		_m.set_next_halfedge_handle(_m.opposite_halfedge_handle(v1h), v4h);
		_m.set_next_halfedge_handle(_m.opposite_halfedge_handle(v2h), v1h);
		_m.set_next_halfedge_handle(_m.opposite_halfedge_handle(v3h), v2h);
		_m.set_next_halfedge_handle(_m.opposite_halfedge_handle(v4h), v3h);

		// �Ǹ��ɵ����������Ϣ��Ҫ������
		_m.set_halfedge_handle(_fh, h1);  
		_m.set_face_handle(h1, _fh); _m.set_face_handle(h8, _fh); 
		_m.set_face_handle(_m.opposite_halfedge_handle(v1h), _fh); _m.set_face_handle(v4h, _fh); 
		// ��һ���¼�������������Ϣ��Ҫ������
		typename mesh_t::FaceHandle f2(_m.new_face());
		_m.set_halfedge_handle(f2, h2);  
		_m.set_face_handle(h2, f2); _m.set_face_handle(h5, f2); 
		_m.set_face_handle(_m.opposite_halfedge_handle(v2h), f2); _m.set_face_handle(v1h, f2); 
		// �ڶ����¼�������������Ϣ��Ҫ������
		typename mesh_t::FaceHandle f3(_m.new_face());
		_m.set_halfedge_handle(f3, h3);  
		_m.set_face_handle(h3, f3); _m.set_face_handle(h6, f3); 
		_m.set_face_handle(_m.opposite_halfedge_handle(v3h), f3); _m.set_face_handle(v2h, f3); 
		// �������¼�������������Ϣ��Ҫ������
		typename mesh_t::FaceHandle f4(_m.new_face());
		_m.set_halfedge_handle(f4, h4);  
		_m.set_face_handle(h4, f4); _m.set_face_handle(h7, f4); 
		_m.set_face_handle(_m.opposite_halfedge_handle(v4h), f4); _m.set_face_handle(v3h, f4); 

		_m.set_halfedge_handle(vh, v1h);
		_m.adjust_outgoing_halfedge(vh);//�����ʵ���Բ�Ҫ, ��Ϊ����vh�����е㣬�������ڱ߽��ϡ�
		}
		void split_edge(mesh_t& _m, const typename mesh_t::EdgeHandle& _eh)
		{
			typename mesh_t::HalfedgeHandle
				heh     = _m.halfedge_handle(_eh, 0),
				opp_heh = _m.halfedge_handle(_eh, 1);

			typename mesh_t::HalfedgeHandle new_heh, opp_new_heh, t_heh;
			typename mesh_t::VertexHandle   vh;
			typename mesh_t::VertexHandle   vh1(_m.to_vertex_handle(heh));
			typename mesh_t::Point          zero(0,0,0);

			// new vertex //������ʾ����ô����һ���µ㲢�޸����������ӡ�
			vh = _m.new_vertex( zero );

			// memorize position, will be set later
			_m.property( vp_pos_, vh ) = _m.property( ep_pos_, _eh );


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
				{} // ��ѭ����tri/quad����Ч.
			}

			new_heh     = _m.new_edge(vh, vh1);
			opp_new_heh = _m.opposite_halfedge_handle(new_heh);
			_m.set_vertex_handle( heh, vh );

			_m.set_next_halfedge_handle(t_heh, opp_new_heh);  // ������İ��
			_m.set_next_halfedge_handle(opp_new_heh, opp_heh);// ������İ��
			_m.set_next_halfedge_handle(new_heh, _m.next_halfedge_handle(heh)); // ������İ��
			_m.set_next_halfedge_handle(heh, new_heh);							// ������İ��

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
		}

	private: // geometry helper
		
		// for arbitrary polygonal mesh to reduce to quad mesh.
		void compute_poly_center(mesh_t& _m, const typename mesh_t::FaceHandle &_fh) {
			int count = 0;
			OpenMesh::Vec3d pos(0,0,0);
			for (typename mesh_t::ConstFaceVertexIter cfv_it(_m, _fh); cfv_it; ++cfv_it) { 
				++count;
				pos += _m.point(cfv_it.handle());
			}
			_m.property(fp_pos_, _fh) = pos / count;
		} 
		void compute_edge_control_point(mesh_t &_m, const typename mesh_t::EdgeHandle &_eh) {
			typename mesh_t::HalfedgeHandle h(_m.halfedge_handle(_eh, 0)); 
			typename mesh_t::HalfedgeHandle oh(_m.halfedge_handle(_eh, 1)); 
			typename mesh_t::FaceHandle f(_m.face_handle(h));
			typename mesh_t::FaceHandle of(_m.face_handle(oh));
			typename mesh_t::Point pos(0, 0, 0);
			pos += _m.point(_m.to_vertex_handle(h)) + _m.point(_m.to_vertex_handle(oh)) + 
				_m.property(fp_pos_, f) + _m.property(fp_pos_, of);
			pos *= 0.25;
			_m.property(ep_pos_, _eh) = pos;
		}
		void compute_vertex_control_point(mesh_t &_m, const typename mesh_t::VertexHandle &_vh) {
			int k = _m.valence(_vh); 
			double invk2 = 1.0 / (k*k);
			typename mesh_t::Point pos(0, 0, 0);
			for (typename mesh_t::VOHIter voh_it(_m, _vh); voh_it; ++voh_it) {
				pos += (_m.point(_m.to_vertex_handle(voh_it))
					+ _m.point(_m.to_vertex_handle(_m.next_halfedge_handle(voh_it))) //������ܴ���, ������ʦ֮���ٸĻ���. P76 sig2000.
					) * invk2; 
			}
			pos += _m.point(_vh) * (k-2.0)/k;
			_m.property(vp_pos_, _vh) = pos;
		}
		// for general quad mesh CC subdi.
		void compute_quad_center(mesh_t& _m, const typename mesh_t::FaceHandle &_fh) {
			int count = 0;
			OpenMesh::Vec3d pos(0,0,0);
			for (typename mesh_t::ConstFaceVertexIter cfv_it(_m, _fh); cfv_it; ++cfv_it) { 
				++count;
				pos += _m.point(cfv_it.handle());
			}
			if (count != 4) std::cerr << "Error: �ı��β��б�Ҫ������е�.\n";
			_m.property(fp_pos_, _fh) = pos / count;
		} 
		void compute_midpoint(mesh_t& _m, const typename mesh_t::EdgeHandle& _eh)
		{	//����ÿһ�����¼�����е�����(odd vertices = new vertices).

			typename mesh_t::HalfedgeHandle heh, opp_heh;

			heh      = _m.halfedge_handle( _eh, 0);
			opp_heh  = _m.halfedge_handle( _eh, 1);

			typename mesh_t::Point pos(0, 0, 0); 
			if (_m.is_boundary(_eh)) { //����ֻ����򵥵ģ���siggraph2000 course notes���и����ӵ�.
				pos = (_m.point(_m.to_vertex_handle(heh)) + _m.point(_m.to_vertex_handle(opp_heh)) ) * 0.5;
			} 
			else {
				pos = (_m.point(_m.to_vertex_handle(heh)) + _m.point(_m.to_vertex_handle(opp_heh)) ) * 0.375; // 3/8.
				typename mesh_t::HHandle
					h1 = _m.next_halfedge_handle(_m.next_halfedge_handle(heh)),
					h2 = _m.next_halfedge_handle(_m.next_halfedge_handle(opp_heh));
				pos += (_m.point(_m.to_vertex_handle(h1)) + _m.point(_m.from_vertex_handle(h1))
					+ _m.point(_m.to_vertex_handle(h2)) + _m.point(_m.from_vertex_handle(h2))) * 0.0625;// 1/16.
			} 
			_m.property( ep_pos_, _eh ) = pos;

		}

		void smooth(mesh_t& _m, const typename mesh_t::VertexHandle& _vh)
		{	//���������㱻����֮�����λ��.
			typename mesh_t::Point            pos(0.0,0.0,0.0);
			typename mesh_t::HHandle h(_m.halfedge_handle(_vh));
			if (_m.is_boundary(h)) {
				typename mesh_t::HHandle prev_h(_m.prev_halfedge_handle(h));
				pos += _m.point(_vh) * 0.75 + 
					(_m.point(_m.to_vertex_handle(h)) + _m.point(_m.from_vertex_handle(prev_h))) * 0.125;
			}
			else {
				int k = _m.valence(_vh);
				double beta = 3.0/(2.0 * k), gamma = 1.0/(4.0 * k);// 
				double beta_invk = beta / k, gamma_invk = gamma / k; 
				for (typename mesh_t::VertexOHalfedgeIter voh_it(_m, _vh); voh_it; ++voh_it) {
					pos += _m.point(_m.to_vertex_handle(voh_it.handle())) * beta_invk
						+ _m.point(_m.to_vertex_handle(_m.next_halfedge_handle(voh_it.handle()))) * gamma_invk; 
				}
				pos += _m.point(_vh) * (1 - beta - gamma);
			}

			_m.property( vp_pos_, _vh ) = pos;
		} // end of function smooth.

	private: // data

		OpenMesh::VPropHandleT< typename mesh_t::Point > vp_pos_;//ÿһ����ĸ���֮���λ��.
		OpenMesh::EPropHandleT< typename mesh_t::Point > ep_pos_;//ÿһ�����������ӵĵ��λ��.
		OpenMesh::FPropHandleT< typename mesh_t::Point > fp_pos_;//�������ı���֮��ÿһ���涼��һ�����ĵ�.
	}; // end of class CatmullClarkSubdT
	 
	// usage function, ��_m����_n��ϸ��.
	template <class MeshT>
	inline void catmull_clark_subdi(MeshT &_m, int _n) {
		typedef DGP::CatmullClarkSubdT<MeshT, double>  CcS;
		CcS subdivide_obj;
		// use interface 2
		subdivide_obj.attach(_m); 
		// 1. please check that all the faces are quads,
		// if it has tri or arbitary poly, then do one more general subdi firstly. 
		MeshTypes mpt = identify_mesh_poly_type(_m);
		if (mpt != MeshTypes::MPY_QUAD) {
			std::cout << "Info: Not only quads in the mesh, 1 general cc done firstly.\n";
			// do it here.
			if (subdivide_obj.reduce_to_quadmesh(_m) == false) {
				std::cout << "Error: Arbitrary polygonal meshes reduce to quadmesh failed in CC.\n";
			} else {
				std::cout << "Arbitrary polygonal mesh reduce to quadrilateral mesh succeeded befor CC.\n";
			}
		}
		subdivide_obj(_n);  
		subdivide_obj.detach();

		// compute face & vertex normals
		_m.update_normals();
	} 
} // end of namespace DGP
#endif //dgp_catmullclarksubd_h