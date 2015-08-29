
#ifndef hoppe94loopsubt_h
#define hoppe94loopsubt_h

#include "OpenMeshAll.h"
#include "VertexFeatureType.h"

namespace DGP {
// ����һ������hoppe 97 Loop Subdivision��ģ����
// 1. ���ԶԺ��б߽��ģ�ͽ���ϸ�֣������ڱ߽��ķ��ߺ����ű߽�ߵ����������ܾ�û�����. 2009-01.

template <typename MeshType, typename RealType = double>
class Hoppe94LoopSubT : private OpenMesh::Utils::Noncopyable  
{
public:
	typedef RealType                                real_t;
	typedef MeshType                                mesh_t;	
	
public:
	Hoppe94LoopSubT(void) : attached_(NULL), vp_type_(NULL)//�������Կ�����һ��operator, ���԰󶨲�ͬ��mesh.
	{ }

	~Hoppe94LoopSubT()   { detach(); }

public: // Interface 1  ��Ϊ�������׺�interface 2��ɻ���, ����ȥ��.
	//bool operator () ( MeshType& _m, size_t _n )
	//{   
		//return prepare(_m) && subdivide( _m, _n ) && cleanup( _m );//��_m����_n��ϸ�ֲ���
	//} 

public: // Interface 2 , �岽��: attach, set feature type, operator(subdivise or/and calc the limit pos), get feature type, and last detach.
	//����ӿڱ�������attached_����Ҫ�Լ�add_property(*vp_type_), ������һ��һ��������������, ������ﻹ��Ҫremove_property().
	// �����е�mesh����Ҫidentify_all_sharp_features(...)��ȷ����ʼ��������.
	bool attach( MeshType& _m, OpenMesh::VPropHandleT<v_feature_type> & _vp_type) {
		if ( attached_ == &_m && vp_type_ == &_vp_type)
			return true;

		detach();
		if (prepare( _m )) {
			attached_ = &_m;
			vp_type_ = &_vp_type;
			return true;
		}
		return false;
	} 

	////���Ǿɵ�����, ����Ҳadd_property(*vp_type_),����һ��һ��������������, ������ﻹ��Ҫremove_property().
	////���Ǻ����ҷ���ԭ����������Լ���property,����ֻ��Ҫȡ�Ǹ�property���þͿ�����, ���������ظ�һ����Դ.
	//void set_edge_feature(typename mesh_t::EdgeHandle _eh, bool _t) { // 
	//	attached_->status(_eh).set_feature(_t);
	//}
	//void set_vertex_feature(typename mesh_t::VertexHandle _vh, v_feature_type _t) {
	//	attached_->property(*vp_type_, _vh) = _t;
	//} 

	/// Subdivide the attached \c _n times.
	/// \see SubdividerT(), attach(), detach()
	bool operator()( size_t _n )
	{
		return (attached_ && vp_type_) ? subdivide( *attached_, _n ) : false;
	}

	OpenMesh::Vec3d limit_pos(typename mesh_t::VertexHandle _vh) {
		OpenMesh::Vec3d limit_pos = attached_->point(_vh);

		v_feature_type vtype = attached_->property(*vp_type_, _vh);
		if (vtype == CORNER_VFT) { // croner vertex still unchanged

		} else if (vtype == CREASE_VFT && attached_->is_boundary(_vh) == false) {
			std::vector<typename mesh_t::VertexHandle> vec; //����crease edge����Ӧ�Ķ���.
			for (typename mesh_t::VOHIter voh_it(*attached_, _vh); voh_it; ++voh_it) {
				if (attached_->status(attached_->edge_handle(voh_it)).feature()) 
					vec.push_back(attached_->to_vertex_handle(voh_it.handle()));				
			}
			if (vec.size() != 2) std::cout << "Error: Hoppe94LoopSubT::limit_pos: ����Ӧ��ֻ������crease edge�Ŷ�.\n" << vec.size() << "; " << vec[0] << ", " << vec[1] << ".\n";
			limit_pos = attached_->point(_vh) * 0.666667 + (attached_->point(vec[0])+attached_->point(vec[1])) * 0.166667;			

		} else if (attached_->is_boundary(_vh)) {                   //
			typename mesh_t::HalfedgeHandle heh = attached_->halfedge_handle(_vh); // 
			if (attached_->is_boundary(heh) == false) { std::cout << "Error: limit_pos::set_limit_position: ���Ҫ���Ǳ߽�Ŷ���.\n"; }
			typename mesh_t::VertexHandle v0 = attached_->to_vertex_handle(heh);
			typename mesh_t::VertexHandle v1 = attached_->from_vertex_handle(attached_->prev_halfedge_handle(heh));
			limit_pos = attached_->point(_vh) * 0.666667 + (attached_->point(v0)+attached_->point(v1)) * 0.166667;

		} else if (vtype == SMOOTH_VFT || vtype == DART_VFT) {
			int k = attached_->valence(_vh);
			double inv_v = 1.0/double(k); double t = (3.0 + 2.0 * cos( 2.0 * M_PI * inv_v) ); // double beta = inv_v * (40.0 - t * t)/64.0;
			double k_alpha = 1.0 / ( 24.0/(40.0-t*t) + 1); double alpha = inv_v * k_alpha;

			typename mesh_t::Point limit_pos = OpenMesh::Vec3d(0, 0, 0);
			for (typename mesh_t::VertexVertexIter vv_it(*attached_, _vh); vv_it; ++vv_it) {
				limit_pos += attached_->point(vv_it.handle()) * alpha; 
			}
			limit_pos += attached_->point(_vh) * (1.0 - k_alpha);

		} else std::cout << "Error: Hoppe94LoopSubT::limit_pos: û���������͵Ķ�������Ŷ��1.\n"; 

		return limit_pos;
	}
	// set all the vertice to it's limit position, to get the limit mesh.
	// Ӧ���������operator��, ���޸�����ģ�͵�.
	bool set_limit_position() {
		for (TriMesh::VIter v_it(attached_->vertices_begin()), v_end(attached_->vertices_end()); v_it != v_end; ++v_it) {
			// calculate the limit position for all points
			attached_->set_point(v_it.handle(), limit_pos(v_it.handle()));
		}
		return true;
	}
	//cjren 09-05-12�˺������������, ֻ��smooth/dart������������, ��Ϊsharp feature���м������õķ����.
	// ����crease edge���ߵķ��಻ͬ, ��ôһ��crease vertex���ĸ�normal�ź�����?  ��ʵӦ��������,�ֱ��ڲ�ͬ��patch.
	// ��flat shading���ܿ�����, ������ÿһƬһƬ����Ⱦ���ܿ���features.
	OpenMesh::Vec3d limit_nor(typename mesh_t::VertexHandle _vh) {
		OpenMesh::Vec3d limit_nor(attached_->normal(_vh));
		OpenMesh::Vec3d p0(attached_->point(_vh));

		v_feature_type vtype = attached_->property(*vp_type_, _vh);
		if (vtype == CORNER_VFT) { // croner vertex still unchanged
			int k = attached_->valence(_vh); // k >= 3; then has k tengent vectors, are coplaner.
			std::vector<typename mesh_t::VHandle> vec; //to record the feature edges.
			typename mesh_t::HHandle h0 = attached_->voh_iter(_vh);
			if (attached_->status(attached_->edge_handle(h0)).feature() == true) {
				vec.push_back(attached_->to_vertex_handle(h0));
			}			
			typename mesh_t::HHandle h = attached_->opposite_halfedge_handle(attached_->prev_halfedge_handle(h0));
			for (; h != h0; ) { //��ôtroublesome��Ϊ�˱���˳��.
				if (attached_->status(attached_->edge_handle(h)).feature() == true) {
					vec.push_back(attached_->to_vertex_handle(h));
				}
				h = attached_->opposite_halfedge_handle(attached_->prev_halfedge_handle(h));
			} 
			if (vec.size() < 3) std::cerr << "Error: Hoppe94LoopSubT::limit_nor(): should >=3.\n";
			
			// every two neighoring tangent vectors 's cross product is one normal dir.
			std::vector<OpenMesh::Vec3d> nor_array(vec.size(), OpenMesh::Vec3d(0,0,0));
			for (int i = 0; i < vec.size(); ++i) { //cross product of tangent vectors.
				nor_array[i] = cross(attached_->point(vec[i]) - p0, attached_->point(vec[(i+1)%vec.size()]) - p0).normalize();
				limit_nor += nor_array[i];
			} 
			limit_nor /= nor_array.size(); limit_nor.normalize();
		} else if (vtype == CREASE_VFT && attached_->is_boundary(_vh) == false) {
			// ����ʵ��,��regular-crease edge��û�����⣬������non-regular crease�¾�Ч��������.
			// crease vertex has 3 tangent vec, one is along the edge, two are across the edge.
			 //crease vertex ����Ӧ������feature edge.
			std::vector<typename mesh_t::HHandle> vec; //crease edge����Ӧ����������.
			typename mesh_t::HHandle h0 = attached_->voh_iter(_vh);
			if (attached_->status(attached_->edge_handle(h0)).feature() == true) {
				vec.push_back(h0);
			}
			typename mesh_t::HHandle h = attached_->opposite_halfedge_handle(attached_->prev_halfedge_handle(h0));
			for (; h != h0; ) { //��ôtroublesome��Ϊ�˱���˳��.
				if (attached_->status(attached_->edge_handle(h)).feature() == true) {
					vec.push_back(h);
				}
				h = attached_->opposite_halfedge_handle(attached_->prev_halfedge_handle(h));
			}
			if (vec.size() != 2) { std::cout << "Error: Hoppe94LoopSubT::limit_nor(): ����Ӧ��ֻ������crease edge�Ŷ�.\n";
				std::cout << vec.size() << ": "; std::cout << vec[0] << ", "; std::cout << vec[1] << ".\n";
			}
			OpenMesh::Vec3d t0 = attached_->point(attached_->to_vertex_handle(vec[0])) - attached_->point(attached_->to_vertex_handle(vec[1]));
			
			int k = attached_->valence(_vh);

			std::vector<typename mesh_t::HHandle> harray;//������crease edge ����ߵ�tangent vector.
			for (typename mesh_t::HHandle h = vec[0]; h != vec[1]; ) {
				harray.push_back(h); 
				h = attached_->opposite_halfedge_handle(attached_->prev_halfedge_handle(h));
			} 
			harray.push_back(vec[1]); 
			int left_n = harray.size();
			OpenMesh::Vec3d t1(0, 0, 0);
			if (k == 6 && harray.size() == 4) {//regular crease edge.
				t1 = (-2.0) * p0 + (-1.0)*(attached_->point(attached_->to_vertex_handle(harray[0])) + attached_->point(attached_->to_vertex_handle(harray[3]))) 
					+ 2.0 * (attached_->point(attached_->to_vertex_handle(harray[1])) + attached_->point(attached_->to_vertex_handle(harray[2])));
			} else if (harray.size() == 2) {
				t1 = (-2.0) * p0 + attached_->point(attached_->to_vertex_handle(harray[0])) 
					+ attached_->point(attached_->to_vertex_handle(harray[1]));				
			} else if (harray.size() == 3) {
				t1 = (-1.0) * p0 + attached_->point(attached_->to_vertex_handle(harray[1]));	
			} else {// non-regular crease edge.
				double  theta = M_PI / (harray.size() - 1); 
				for (int i = 1; i <= harray.size() - 2; ++i) {
					t1 += (2.0*cos(theta) - 2)*(sin(i*theta)) * attached_->point(attached_->to_vertex_handle(harray[i])); 
				}
				t1 += sin(theta) * 
					(attached_->point(attached_->to_vertex_handle(harray[0])) + attached_->point(attached_->to_vertex_handle(harray[harray.size()-1]))); 
				//if (k != 6 && harray.size() == 4) 
				//std::cout << "n1: " << harray.size() << ", " << "t1: " << t1 << ". ";
			}
			

			harray.clear(); //������crease edge���ұߵ�tangent vector
			for (typename mesh_t::HHandle h = vec[0]; h != vec[1]; ) {
				harray.push_back(h); 
				h = attached_->next_halfedge_handle(attached_->opposite_halfedge_handle(h));
			} 
			harray.push_back(vec[1]); 
			int right_n = harray.size();
			if (left_n + right_n != k + 2) { std::cout << "Error: limit_nor():ö�ٰ��ʱ������.\n";}
			OpenMesh::Vec3d t2(0, 0, 0);
			if (k == 6 && harray.size() == 4) {//regular crease edge.
				t2 = (-2.0) * p0 + (-1.0)*(attached_->point(attached_->to_vertex_handle(harray[0])) + attached_->point(attached_->to_vertex_handle(harray[3]))) 
					+ 2.0 * (attached_->point(attached_->to_vertex_handle(harray[1])) + attached_->point(attached_->to_vertex_handle(harray[2])));
			} else if (harray.size() == 2) {
				t2 = (-2.0) * p0 + attached_->point(attached_->to_vertex_handle(harray[0])) 
					+ attached_->point(attached_->to_vertex_handle(harray[1]));		
			} else if (harray.size() == 3) {
				t2 = (-1.0) * p0 + attached_->point(attached_->to_vertex_handle(harray[1]));	
			} else {// non-regular crease edge.
				double  theta = M_PI * (harray.size() - 1);
				for (int i = 1; i <= harray.size() - 2; ++i) {
					t2 += (2.0*cos(theta) - 2)*(sin(i*theta)) * attached_->point(attached_->to_vertex_handle(harray[i]));
				}
				t2 += sin(theta) * 
					(attached_->point(attached_->to_vertex_handle(harray[0])) + attached_->point(attached_->to_vertex_handle(harray[harray.size()-1])));
				//std::cout << "n2: " << harray.size() << ", " << "t2: " << t1 << ". ";
			}
			OpenMesh::Vec3d n1 = cross(t0, t1).normalize(), n2 = cross(t2, t0).normalize();
			attached_->property(vp_t0_, _vh) = t0; 
			attached_->property(vp_t1_, _vh) = t1; 
			attached_->property(vp_t2_, _vh) = t2; 
			attached_->property(vp_n1_, _vh) = n1; 
			attached_->property(vp_n2_, _vh) = n2; 
			limit_nor = ((n1 + n2) * 0.5).normalize();
			
			
		} 
		//else if (attached_->is_boundary(_vh)) {                   //
		//	typename mesh_t::HalfedgeHandle heh = attached_->halfedge_handle(_vh); // 
		//	if (attached_->is_boundary(heh) == false) { std::cout << "Error: limit_nor()::set_limit_position: ���Ҫ���Ǳ߽�Ŷ���.\n"; }
		//	typename mesh_t::HHandle prevh = attached_->prev_halfedge_handle(heh);
		//	OpenMesh::Vec3d t0 = (attached_->point(attached_->to_vertex_handle(heh)) - attached_->point(attached_->from_vertex_handle(prevh)));
		//	
		//	int k = attached_->valence(_vh);
		//	std::vector<typename mesh_t::HHandle> vec(2); //crease edge����Ӧ����������.
		//	vec[0] = heh; vec[1] = attached_->opposite_halfedge_handle(prevh);
		//	std::vector<typename mesh_t::HHandle> harray;
		//	harray.clear(); //������crease edge���ұߵ�tangent vector
		//	for (typename mesh_t::HHandle h = vec[0]; h != vec[1]; ) {
		//		harray.push_back(h); 
		//		h = attached_->next_halfedge_handle(attached_->opposite_halfedge_handle(h));
		//	} 
		//	harray.push_back(vec[1]); 
		//	OpenMesh::Vec3d t2(0, 0, 0);
		//	if (k == 4 && harray.size() == 4) {//regular crease edge.
		//		t2 = (-2.0) * p0 + (-1.0)*(attached_->point(attached_->to_vertex_handle(harray[0])) + attached_->point(attached_->to_vertex_handle(harray[3]))) 
		//			+ 2.0 * (attached_->point(attached_->to_vertex_handle(harray[1])) + attached_->point(attached_->to_vertex_handle(harray[2])));
		//	} else if (harray.size() == 3) {
		//		 
		//	} else {// non-regular crease edge.
		//		double  theta = M_PI / (harray.size() - 1);
		//		for (int i = 1; i <= harray.size() - 2; ++i) {
		//			t2 += (2.0*cos(theta) - 2)*(sin(i*theta)) * attached_->point(attached_->to_vertex_handle(harray[i]));
		//		}
		//		t2 += sin(theta) * 
		//			(attached_->point(attached_->to_vertex_handle(harray[0])) + attached_->point(attached_->to_vertex_handle(harray[harray.size()-1])));

		//	}
		//	limit_nor = cross(t2, t0).normalize();

		//} 
		else if (vtype == SMOOTH_VFT || vtype == DART_VFT) {
			limit_nor = cloop_limit_nor(*attached_, _vh);
		} else std::cout << "Error: Hoppe94LoopSubT::limit_nor(): û���������͵Ķ�������Ŷ��1.\n"; 
		
		return limit_nor;
	}
	// 
	bool edge_feature(typename mesh_t::EdgeHandle _eh) {
		return attached_->status(_eh).feature();
	}
	v_feature_type vertex_feature(typename mesh_t::VertexHandle _vh) {
		return attached_->property(*vp_type_, _vh);
	}
	/// Detach an eventually attached mesh. �����������Ĳ������֮��ſ��Ե���detach()������.
	/// \see SubdividerT(), attach(), operator()(size_t)
	void detach(void)
	{
		if ( attached_ )
		{
			cleanup( *attached_ );
			attached_ = NULL;
			vp_type_ = NULL;
		}
	}
	// -------------------------
public:
	const char *name() const { return "Hoppe94 Loop Subdivision"; }

protected:
	bool prepare( mesh_t& _m )
	{
		_m.add_property( vp_pos_ );
		_m.add_property( ep_pos_ ); 
		_m.add_property( vp_t0_ ); 
		_m.add_property( vp_t1_ ); 
		_m.add_property( vp_t2_ ); 
		_m.add_property( vp_n1_ ); 
		_m.add_property( vp_n2_ ); 
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
			for (vit  = _m.vertices_begin(); vit != _m.vertices_end(); ++vit)
				smooth( _m, vit.handle() );//���ÿһ���ɵĶ��㱻���º��������.
			
			// Compute position for new vertices and store them in the edge property
			for (eit=_m.edges_begin(); eit != _m.edges_end(); ++eit)
				compute_midpoint( _m, eit.handle() );//���ÿһ�����������ӵĵ������

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
			for ( vit  = _m.vertices_begin(); vit != _m.vertices_end(); ++vit)
				_m.set_point(vit, _m.property( vp_pos_, vit ) );

		} // end of for

		_m.update_normals();
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

		// ����Ҳ���Ҽ����, ��Ϊ�¼���Ķ���vh����(new_heh, opp_new_heh)���������ɵ��µı߶���Ҫָ������
		if (_m.status(_eh).feature()) { // �ղŽ��з��ѵ��������crease edge
			_m.property(*vp_type_, vh) = CREASE_VFT;
			_m.status(_m.edge_handle(new_heh)).set_feature(true);
		} 
		else if (_m.is_boundary(_eh)) { 
			//������identify_all_sharp_features()�����н��߽�Ҳ��Ϊһ��feature edge�Ļ�,
			//��ô�߽����¼�����е�Ҳ��Ӧ��crease vertex.
			//������������crease����Χû������crease edge, ��Ӧ����������boundary edge. 
			//����û��ȷ�������Ƿ�����г��������ôӰ��. 2009-01.
			_m.property(*vp_type_, vh) = CREASE_VFT; 
		}
		else {
			_m.property(*vp_type_, vh) = SMOOTH_VFT;
			_m.status(_m.edge_handle(new_heh)).set_feature(false);
		}
	}

private: // geometry helper

	void compute_midpoint(mesh_t& _m, const typename mesh_t::EdgeHandle& _eh)
	{

		if (_m.status(_eh).feature() || _m.is_boundary(_eh)) {
			TriMesh::VertexHandle v0 = _m.to_vertex_handle(_m.halfedge_handle(_eh, 0));
			TriMesh::VertexHandle v1 = _m.to_vertex_handle(_m.halfedge_handle(_eh, 1));
			_m.property(ep_pos_, _eh) = (_m.point(v0) + _m.point(v1)) * 0.5;
		
		} else { // smooth edge, new vertex = 3/8 * v0 + 3/8 * v1 + 1/8 * v2 + 1/8 * v3
			TriMesh::VertexHandle v0 = _m.to_vertex_handle(_m.halfedge_handle(_eh, 0));
			TriMesh::VertexHandle v1 = _m.to_vertex_handle(_m.halfedge_handle(_eh, 1));
			TriMesh::VertexHandle v2 = _m.to_vertex_handle(_m.next_halfedge_handle(_m.halfedge_handle(_eh, 0)));
			TriMesh::VertexHandle v3 = _m.to_vertex_handle(_m.next_halfedge_handle(_m.halfedge_handle(_eh, 1)));

			_m.property(ep_pos_, _eh) = (_m.point(v0) + _m.point(v1)) * 0.375 + (_m.point(v2) + _m.point(v3)) * 0.125;
		}
	}


	void smooth(mesh_t& _m, const typename mesh_t::VertexHandle& _vh)
	{//���������㱻����֮�����λ��.
		typename mesh_t::Point            pos(0.0,0.0,0.0);

		v_feature_type vtype = _m.property(*vp_type_, _vh);
		if (vtype == CORNER_VFT) {

			pos = _m.point(_vh);
		} else if (vtype == CREASE_VFT && _m.is_boundary(_vh) == false) {

			std::vector<TriMesh::VertexHandle> vec; //����crease edge����Ӧ�Ķ���.
			for (TriMesh::VOHIter voh_it(_m, _vh); voh_it; ++voh_it) {
				if (_m.status(_m.edge_handle(voh_it)).feature()) 
					vec.push_back(_m.to_vertex_handle(voh_it.handle()));				
			}
			if (vec.size() != 2) std::cout << "Error: Hoppe94LoopSubT::smooth: ����Ӧ��ֻ������crease edge�Ŷ�, " 
				<< vec.size() << "; " << vec[0] << ", " << vec[1] << ".\n";
			pos = _m.point(_vh) * 0.75 + _m.point(vec[0]) * 0.125 + _m.point(vec[1]) * 0.125;
		} else if (_m.is_boundary(_vh)) {                   //

			TriMesh::HalfedgeHandle heh = _m.halfedge_handle(_vh); //de
			if (_m.is_boundary(heh) == false) { std::cout << "Error: Hoppe94LoopSubT::smooth: ���Ҫ���Ǳ߽�Ŷ���.\n"; }
			TriMesh::VertexHandle v0 = _m.to_vertex_handle(heh);
			TriMesh::VertexHandle v1 = _m.from_vertex_handle(_m.prev_halfedge_handle(heh));
			pos = _m.point(_vh) * 0.75
				+ _m.point(v0) * 0.125 + _m.point(v1) * 0.125;//��crease vertex�ĸ���ģʽ��һ����

		} else if (vtype == SMOOTH_VFT || vtype == DART_VFT) {
			int k = _m.valence(_vh);
			double beta = 0, alpha = 0;
			if (k == 3) { beta = 0.1875; alpha = 0.5625; } // 3.0/16.0 = 0.1875,   9.0/16.0 = 0.5625
			else if (k == 6) { beta = 0.0625; alpha = 0.375; } // 1.0/16.0 = 0.0625,    3.0/8.0 = 0.375
			else { 
				double inv_v = 1.0/double(k); double t = (3.0 + 2.0 * cos( 2.0 * M_PI * inv_v) );
				alpha  = (40.0 - t * t)/64.0; beta = inv_v * alpha;
			}
			pos = _m.point(_vh)*(1.0 - alpha); // 

			for (TriMesh::VVIter vv_it(_m, _vh); vv_it; ++vv_it) {
				pos += _m.point(vv_it.handle()) * beta;
			}
		} else std::cout << "Error: Hoppe94LoopSubT::smooth: û���������͵Ķ�������Ŷ�.\n";

		_m.property( vp_pos_, _vh ) = pos;
	} // end of function smooth.


private: // data

	MeshType *attached_; 
	OpenMesh::VPropHandleT<v_feature_type> *vp_type_;// 
	
	OpenMesh::VPropHandleT< typename mesh_t::Point > vp_pos_;//ÿһ����ĸ���֮���λ��.
	OpenMesh::EPropHandleT< typename mesh_t::Point > ep_pos_;//ÿһ�����������ӵĵ��λ��.

	OpenMesh::VPropHandleT< typename mesh_t::Point > vp_t0_;//ÿһ����ĸ���֮���λ��.
	OpenMesh::VPropHandleT< typename mesh_t::Point > vp_t1_;//ÿһ����ĸ���֮���λ��.
	OpenMesh::VPropHandleT< typename mesh_t::Point > vp_t2_;//ÿһ����ĸ���֮���λ��.
	OpenMesh::VPropHandleT< typename mesh_t::Point > vp_n1_;//ÿһ����ĸ���֮���λ��.
	OpenMesh::VPropHandleT< typename mesh_t::Point > vp_n2_;//ÿһ����ĸ���֮���λ��.

};

// ������crease edge _crease_h��ʱ�뷽��ö�ٵ�һ�������������������һ��along crease edge, һ��across crease edge.
// crease vertex���ķ������ cross(_t_long, _t_cross).normalize();ע�������.
// Reference H.Hoppe94. Piecewise smooth surface reconstruction. ��ֱ�ӿ������ʵ��MeshOp.C.
template <typename TriMeshT>
inline bool tangent_vectors_at_crease_vertex(const TriMeshT &_m, const OpenMesh::VPropHandleT<v_feature_type> &_vp_type, //input
											 const typename TriMeshT::HHandle &_crease_h, //input
											 OpenMesh::Vec3d &_t_along, OpenMesh::Vec3d &_t_across, //output
											 std::vector<typename TriMeshT::HHandle> &_smooth_edges_array)  // output
{
	if (_m.property(_vp_type, _m.from_vertex_handle(_crease_h)) != DGP::CREASE_VFT) {
		std::cerr << "Error: �����ǵĶ��㲻��crease����.\n"; 
		return false;
	} 
	typename TriMeshT::VHandle vp = _m.from_vertex_handle(_crease_h);
	typename TriMeshT::HHandle h1 = _crease_h, hn;
	_smooth_edges_array.clear(); // smooth edges between two crease edges.
	for (hn = _m.opposite_halfedge_handle(_m.prev_halfedge_handle(h1)); _m.status(_m.edge_handle(hn)).feature() == false; 
		hn = _m.opposite_halfedge_handle(_m.prev_halfedge_handle(hn)))
	{
		_smooth_edges_array.push_back(hn); 
	}
	_t_along = _m.point(_m.to_vertex_handle(hn)) - _m.point(_m.to_vertex_handle(h1));
	// ���������������Ϊ: h1, smooth_edges_array[], hn. 
	if (_smooth_edges_array.size() == 0) {
		_t_across = _m.point(vp) * 2.0 - (_m.point(_m.to_vertex_handle(hn)) + _m.point(_m.to_vertex_handle(h1)));
	}
	else if (_smooth_edges_array.size() == 1) {
		_t_across = _m.point(vp) - _m.point(_m.to_vertex_handle(_smooth_edges_array[0])); 
	} 
	else if (_smooth_edges_array.size() == 2 && _m.valence(vp) == 6) { // regular crease vertex
		_t_across = _m.point(vp) + 0.5*(_m.point(_m.to_vertex_handle(hn)) + _m.point(_m.to_vertex_handle(h1)))
			- (_m.point(_m.to_vertex_handle(_smooth_edges_array[0])) + _m.point(_m.to_vertex_handle(_smooth_edges_array[1])));
	} 
	else { // non-regular crease vertex
		double theta = M_PI/(_smooth_edges_array.size() + 2 - 1);
		double w1 = 1./(2.-2.*cos(theta));
		_t_across = (_m.point(_m.to_vertex_handle(hn)) + _m.point(_m.to_vertex_handle(h1))) * w1;
		for (int i = 0; i < _smooth_edges_array.size(); i++)
			_t_across -= _m.point(_m.to_vertex_handle(_smooth_edges_array[i])) *sin((i+1)*theta)/sin(theta); 
	}
	return true;
} // end of function tangent_vectors_at_crease_vertex
} //end of namespace DGP
#endif hoppe94loopsubt_h