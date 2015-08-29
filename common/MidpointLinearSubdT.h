// ***************************************************************
//  
//  MidpointLinearSubdT 
//  Copyright (C) 2008 - by rencanjiang All Rights Reserved
//  rencanjiang@163.com	
//	-------------------------------------------------------------
//  
//  version: 1.0
//  data: 2008-04-10 
// ***************************************************************
//  
//  File description: 1. 只是三角形线性中点细分, 每个边加入一个中点而已, 
//  模拟 <OpenMesh/Tools/Subdivider/Uniform/LoopT.hh>来实现的. 
//  2. 加入了对四边形的中点线性细分.2009-01-05.
//  3. 注意这里只是支持三角形和四边形.
//            *                   *
//           / \                 / \
//          /   \               /   \
//         /     \   ====>     *-----*
//        /       \           /\     /\
//       /         \         /  \   /  \
//      *-----------*       *----\*/----* 
//      |           |       |     |     | 
//      |           |       |     |     | 
//      |           |       *-----*-----* 
//      |           |       |     |     | 
//      |           |       |     |     | 
//      |___________|       *-----*-----* 
// ***************************************************************
#ifndef dgp_midpointlinearsub_h
#define dgp_midpointlinearsub_h

#include "OpenMeshAll.h"

namespace DGP {
	
template <typename MeshType, typename RealType = float>
class MidpointLinearSubdT : public OpenMesh::Subdivider::Uniform::SubdividerT<MeshType, RealType>
{
public:

	typedef RealType                                real_t;
	typedef MeshType								mesh_t;
	typedef SubdividerT< mesh_t, real_t >           parent_t;

public:

	MidpointLinearSubdT(void) : parent_t()
	{ }
	~MidpointLinearSubdT() {}
public:
	const char *name() const { return "Minpoint linear Subdivision for Tri/Quad"; }

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
			// compute new positions for old vertices, 对tri/quad 都有效.
			for ( vit  = _m.vertices_begin(); vit != _m.vertices_end(); ++vit)
				smooth( _m, vit.handle() );//求出每一个旧的顶点被更新后的新坐标.


			// Compute position for new vertices and store them in the edge property, 对tri/quad都有效.
			for (eit=_m.edges_begin(), e_end = _m.edges_end(); eit != e_end; ++eit)
				compute_midpoint( _m, eit.handle() );//求出每一个边上先增加的点的坐标

			// 对于四边形来说，还需要计算每一个面新加入的中点. 只对quad有效.
			for (fit = _m.faces_begin(), f_end = _m.faces_end(); fit != f_end; ++fit) {
				if (n_face_corners(_m, fit.handle()) == 4) compute_quad_center(_m, fit.handle());
			}
			 
			// Split each edge at midpoint and store precomputed positions (stored in
			// edge property ep_pos_) in the vertex property vp_pos_;

			// Attention! Creating new edges, hence make sure the loop ends correctly. 对tri/quad都有效.
			for (eit=_m.edges_begin(), e_end = _m.edges_end(); eit != e_end; ++eit)
				split_edge(_m, eit.handle() );
			
			// Commit changes in topology and reconsitute consistency

			// Attention! Creating new faces, hence make sure the loop ends correctly.
			for (fit = _m.faces_begin(), f_end   = _m.faces_end(); fit != f_end; ++fit) 
				if (n_face_corners(_m, fit.handle()) == 8) split_face_quad(_m, fit.handle());// for quad. 因为上面split edge后一个边变了2个.
				else if (n_face_corners(_m, fit.handle()) == 6) split_face(_m, fit.handle() ); // for tri.
				else std::cout << "Error: MidpointLinearSubdT, split edge error.\n";

			
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
		//      *---->*---->*
		//          1 ^   6
		//            vh1 (adjust_outgoing halfedge!)
		//

		// Old and new Face             //下面很好地显示了怎么添加一个新面，对旧面怎么做拓扑的修改.
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
		if (h1 != _m.next_halfedge_handle(_m.next_halfedge_handle(h4))) std::cerr << "Error: MidpointLinearSubdT split_face_quad, halfedges.\n";
		typename mesh_t::VertexHandle
			v1(_m.to_vertex_handle(h1)), v2(_m.to_vertex_handle(h2)), 
			v3(_m.to_vertex_handle(h3)), v4(_m.to_vertex_handle(h4));
		
		// add new vertex (face vertex)
		typename mesh_t::VHandle vh = _m.new_vertex(typename mesh_t::Point(0, 0, 0));
		_m.property(vp_pos_, vh) = _m.property(fp_pos_, _fh);  
		
		// 先加入4个边，并完善半边之间的连接关系.
		typename mesh_t::HHandle 
			v1h = _m.new_edge(vh, v1), v2h = _m.new_edge(vh, v2),
			v3h = _m.new_edge(vh, v3), v4h = _m.new_edge(vh, v4);
		_m.set_next_halfedge_handle(v1h, h5); // 此时被下修改前h1的下一个半边还是原来的.
		_m.set_next_halfedge_handle(h1, _m.opposite_halfedge_handle(v1h));//此时h1的下一个半边就被修改了.
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

		// 那个旧的面的连接信息需要被更新
		_m.set_halfedge_handle(_fh, h1);  
		_m.set_face_handle(h1, _fh); _m.set_face_handle(h8, _fh); 
		_m.set_face_handle(_m.opposite_halfedge_handle(v1h), _fh); _m.set_face_handle(v4h, _fh); 
		// 第一个新加入的面的连接信息需要被加入
		typename mesh_t::FaceHandle f2(_m.new_face());
		_m.set_halfedge_handle(f2, h2);  
		_m.set_face_handle(h2, f2); _m.set_face_handle(h5, f2); 
		_m.set_face_handle(_m.opposite_halfedge_handle(v2h), f2); _m.set_face_handle(v1h, f2); 
		// 第二个新加入的面的连接信息需要被加入
		typename mesh_t::FaceHandle f3(_m.new_face());
		_m.set_halfedge_handle(f3, h3);  
		_m.set_face_handle(h3, f3); _m.set_face_handle(h6, f3); 
		_m.set_face_handle(_m.opposite_halfedge_handle(v3h), f3); _m.set_face_handle(v2h, f3); 
		// 第三个新加入的面的连接信息需要被加入
		typename mesh_t::FaceHandle f4(_m.new_face());
		_m.set_halfedge_handle(f4, h4);  
		_m.set_face_handle(h4, f4); _m.set_face_handle(h7, f4); 
		_m.set_face_handle(_m.opposite_halfedge_handle(v4h), f4); _m.set_face_handle(v3h, f4); 

		_m.set_halfedge_handle(vh, v1h);
		_m.adjust_outgoing_halfedge(vh);//这句其实可以不要, 因为这里vh是面中点，不会是在边界上。
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

		// new vertex //这里显示了怎么加入一个新点并修改其拓扑连接。
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
	}

private: // geometry helper
	void compute_quad_center(mesh_t& _m, const typename mesh_t::FaceHandle &_fh) {
		int count = 0;
		OpenMesh::Vec3d pos(0,0,0);
		for (typename mesh_t::ConstFaceVertexIter cfv_it(_m, _fh); cfv_it; ++cfv_it) { 
			++count;
			pos += _m.point(cfv_it.handle());
		}
		if (count != 4) std::cerr << "Error: 四边形才有必要求面的中点.\n";
		_m.property(fp_pos_, _fh) = pos / count;
	} 
	void compute_midpoint(mesh_t& _m, const typename mesh_t::EdgeHandle& _eh)
	{	//计算每一个边新加入的中点坐标

		typename mesh_t::HalfedgeHandle heh, opp_heh;

		heh      = _m.halfedge_handle( _eh, 0);
		opp_heh  = _m.halfedge_handle( _eh, 1);

		typename mesh_t::Point pos(_m.point(_m.to_vertex_handle(heh)));

		pos += ( _m.point(_m.to_vertex_handle(opp_heh)) );

		pos *= 0.5;
		_m.property( ep_pos_, _eh ) = pos;

	}

	void smooth(mesh_t& _m, const typename mesh_t::VertexHandle& _vh)
	{	//求出这个顶点被更新之后的新位置.
		typename mesh_t::Point            pos(0.0,0.0,0.0);

		pos  = _m.point( _vh );

		_m.property( vp_pos_, _vh ) = pos;
	} // end of function smooth.

private: // data

	OpenMesh::VPropHandleT< typename mesh_t::Point > vp_pos_;//每一个点的更新之后的位置.
	OpenMesh::EPropHandleT< typename mesh_t::Point > ep_pos_;//每一个边上新增加的点的位置.
	OpenMesh::FPropHandleT< typename mesh_t::Point > fp_pos_;//增加了四边形之后每一个面都有一个中心点.
};

/*
// 上面那个继承的版本不知为什么不对, 所以这个写个叫完整的.
// 其实就是加入了attach, detach和two interface functions. 
template <typename MeshType, typename RealType = float>
class MidpointLinearSubdT 
{
public:
	typedef RealType                                real_t;
	typedef MeshType                                mesh_t;	

public:
	MidpointLinearSubdT(void) : attached_(NULL) 
	{ }

	~MidpointLinearSubdT()   { detach(); }

public: // Interface 1 
	bool operator () ( MeshType& _m, size_t _n )
	{    
		return prepare(_m) && subdivide( _m, _n ) && cleanup( _m );//对_m进行_n次细分操作
	} 

private:
	MeshType *attached_;
public: // Interface 2 
	bool attach( MeshType& _m )
	{
		if ( attached_ == &_m )
			return true;

		detach();
		if (prepare( _m ))
		{
			attached_ = &_m;
			return true;
		}
		return false;
	}

	/// Subdivide the attached \c _n times.
	/// \see SubdividerT(), attach(), detach()
	bool operator()( size_t _n )
	{
		return attached_ ? subdivide( *attached_, _n ) : false;
	}

	/// Detach an eventually attached mesh.
	/// \see SubdividerT(), attach(), operator()(size_t)
	void detach(void)
	{
		if ( attached_ )
		{
			cleanup( *attached_ );
			attached_ = NULL;
		}
	}
	// -------------------------
public:
	const char *name() const { return "Minpoint linear Subdivision"; }

protected:
	bool prepare( mesh_t& _m )
	{
		_m.add_property( vp_pos_ );
		_m.add_property( ep_pos_ );
		return true;
	}
	bool cleanup( mesh_t& _m )
	{
		_m.remove_property( vp_pos_ );
		_m.remove_property( ep_pos_ );
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
			// compute new positions for old vertices
			for ( vit  = _m.vertices_begin();
				vit != _m.vertices_end(); ++vit)
				smooth( _m, vit.handle() );//求出每一个旧的顶点被更新后的新坐标.


			// Compute position for new vertices and store them in the edge property
			for (eit=_m.edges_begin(); eit != _m.edges_end(); ++eit)
				compute_midpoint( _m, eit.handle() );//求出每一个边上先增加的点的坐标


			// Split each edge at midpoint and store precomputed positions (stored in
			// edge property ep_pos_) in the vertex property vp_pos_;

			// Attention! Creating new edges, hence make sure the loop ends correctly.
			e_end = _m.edges_end();
			for (eit=_m.edges_begin(); eit != e_end; ++eit)
				split_edge(_m, eit.handle() );


			// Commit changes in topology and reconsitute consistency

			// Attention! Creating new faces, hence make sure the loop ends correctly.
			f_end   = _m.faces_end();
			for (fit = _m.faces_begin(); fit != f_end; ++fit)
				split_face(_m, fit.handle() );

			// Commit changes in geometry
			for ( vit  = _m.vertices_begin();
				vit != _m.vertices_end(); ++vit)
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

		// Old and new Face
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

	void split_edge(mesh_t& _m, const typename mesh_t::EdgeHandle& _eh)
	{
		typename mesh_t::HalfedgeHandle
			heh     = _m.halfedge_handle(_eh, 0),
			opp_heh = _m.halfedge_handle(_eh, 1);

		typename mesh_t::HalfedgeHandle new_heh, opp_new_heh, t_heh;
		typename mesh_t::VertexHandle   vh;
		typename mesh_t::VertexHandle   vh1(_m.to_vertex_handle(heh));
		typename mesh_t::Point          zero(0,0,0);

		// new vertex
		vh                = _m.new_vertex( zero );

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
			{}
		}

		new_heh     = _m.new_edge(vh, vh1);
		opp_new_heh = _m.opposite_halfedge_handle(new_heh);
		_m.set_vertex_handle( heh, vh );

		_m.set_next_halfedge_handle(t_heh, opp_new_heh);
		_m.set_next_halfedge_handle(new_heh, _m.next_halfedge_handle(heh));
		_m.set_next_halfedge_handle(heh, new_heh);
		_m.set_next_halfedge_handle(opp_new_heh, opp_heh);

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

	void compute_midpoint(mesh_t& _m, const typename mesh_t::EdgeHandle& _eh)
	{

		typename mesh_t::HalfedgeHandle heh, opp_heh;

		heh      = _m.halfedge_handle( _eh, 0);
		opp_heh  = _m.halfedge_handle( _eh, 1);

		typename mesh_t::Point
			pos(_m.point(_m.to_vertex_handle(heh)));

		pos += ( _m.point(_m.to_vertex_handle(opp_heh)) );

		pos *= 0.5;
		_m.property( ep_pos_, _eh ) = pos;

	}

	void smooth(mesh_t& _m, const typename mesh_t::VertexHandle& _vh)
	{//求出这个顶点被更新之后的新位置.
		typename mesh_t::Point            pos(0.0,0.0,0.0);

		pos  = _m.point( _vh );

		_m.property( vp_pos_, _vh ) = pos;
	} // end of function smooth.

private: // data

	OpenMesh::VPropHandleT< typename mesh_t::Point > vp_pos_;//每一个点的更新之后的位置.
	OpenMesh::EPropHandleT< typename mesh_t::Point > ep_pos_;//每一个边上新增加的点的位置.
};*/

// usage function, 对_m进行_n次细分.
template <class MeshT>
inline void midpoint_linear_subdi(MeshT &_m, int _n) {

	typedef DGP::MidpointLinearSubdT<MeshT, double>  MidS;
	MidS subdivide_obj;

	//std::cout << "Subdivide " << _n  << " times with '" << subdivide_obj.name() << "'\n";

	// use interface 1
	subdivide_obj(_m, _n);
	// use interface 2
	//subdivide_obj.attach(_m); subdivide_obj(_n);  subdivide_obj.detach();

	// compute face & vertex normals
	_m.update_normals();
}
} // end of namespace DGP
#endif //dgp_midpointlinearsub_h