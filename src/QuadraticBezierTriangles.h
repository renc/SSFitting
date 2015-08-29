#ifndef dgp_quadratic_bezier_triangles_h
#define dgp_quadratic_bezier_triangles_h 

#include "DGPcommon.h"

namespace DGP {

	/**
	�ο�����:
	Tamy Boubekeur, etc. QAS: Real-time Quadratic Approximation of Subdivision Surfaces. 2007.
	Tamy Boubekeur, etc. Approximation of Subdivision Surfaces for Interactive Applications. 2007.	
	����ָ����CPU�µ�ʵ��, ��Ȼû��GPU��ʵ��.
	*/
	class QuadraticBezierTriangles {
	public:
		QuadraticBezierTriangles() {			 
		}
		~ QuadraticBezierTriangles() {}

		void create_geometry(const TriMesh &_m, OpenMesh::FPropHandleT<TriMesh::HHandle> & fp_hh, 
			OpenMesh::VPropHandleT<DGP::v_feature_type> & vp_type_of_original_mesh) 
		{
			TriMesh m;
			DGP::mesh_copy(m, _m);
			DGP::cloop_subdi(m, 1);
			DGP::cloop_limit_surf(m); 

			f_pijk_ = std::vector<std::vector<OpenMesh::Vec3d> > (_m.n_faces(), std::vector<OpenMesh::Vec3d> (6, OpenMesh::Vec3d(0, 0, 0)));
			f_nijk_ = std::vector<std::vector<OpenMesh::Vec3d> > (_m.n_faces(), std::vector<OpenMesh::Vec3d> (6, OpenMesh::Vec3d(0, 0, 0)));
			f_pijk_sf_ = std::vector<std::vector<OpenMesh::Vec3d> > (_m.n_faces(), std::vector<OpenMesh::Vec3d> (6, OpenMesh::Vec3d(0, 0, 0)));
			f_nijk_sf_ = std::vector<std::vector<OpenMesh::Vec3d> > (_m.n_faces(), std::vector<OpenMesh::Vec3d> (6, OpenMesh::Vec3d(0, 0, 0)));

			TriMesh msf; // mesh with sharp features
			DGP::mesh_copy(msf, _m); // mesh_.status(eh).feature()������copy��ȥ��.
			//OpenMesh::VPropHandleT<DGP::v_feature_type> vp_type_of_original_mesh;
			//if (_m.get_property_handle(vp_type_of_original_mesh, "name:vp_type") == false) {
			//	std::cerr << "Error: QAS,ԭʼģ����Ӧ�������\"name:vp_type\"��property����¼ÿ�������������.\n";  
			//} // �Ӳ���������ģ����Ҫ�������е�property����������, һ��������������, ����ֱ���ڲ���������pass-by-reference.
			OpenMesh::VPropHandleT<DGP::v_feature_type> vp_type;
			msf.add_property(vp_type); 
			for (TriMesh::VIter v_it(msf.vertices_begin()), v_end(msf.vertices_end()); v_it != v_end; ++v_it) {
				msf.property(vp_type, v_it) = _m.property(vp_type_of_original_mesh, v_it);
			}


			DGP::Hoppe94LoopSubT<TriMesh> subdi_obj;
			subdi_obj.attach(msf, vp_type); 
			subdi_obj(1);  

			//OpenMesh::FPropHandleT<TriMesh::HHandle> fp_hh;
			//if (_m.get_property_handle(fp_hh, "name:fp_hh") ) {
			if (fp_hh.is_valid() ) {
				for (TriMesh::ConstFaceIter f_it(_m.faces_begin()), f_end(_m.faces_end()); f_it != f_end; ++f_it) {
					TriMesh::HHandle h = _m.property(fp_hh, f_it);//_m����f_it������ȷ��ö��3������˳��İ��. 
					//TriMesh::HHandle hsf = msf.property(hp_hh_sf, h);
					// ����û�п��Ǵ��߽������
					if (_m.from_vertex_handle(h) != m.from_vertex_handle(h)) {
						h = DGP::backward3_halfedge_handle(m, h);
					}
					if (_m.from_vertex_handle(_m.property(fp_hh, f_it)) == m.from_vertex_handle(h)) { 
						TriMesh::HHandle h_tmp = h;
						OpenMesh::Vec3d p0 = m.point(m.from_vertex_handle(h));
						OpenMesh::Vec3d n0 = DGP::cloop_limit_nor(m, m.from_vertex_handle(h));//m.normal(m.from_vertex_handle(h));//
						
						OpenMesh::Vec3d p0e = m.point(m.to_vertex_handle(h));
						OpenMesh::Vec3d n0e = DGP::cloop_limit_nor(m, m.to_vertex_handle(h));//m.normal(m.to_vertex_handle(h));//
						
						h = m.next_halfedge_handle(h); // p0e -> p2e
						OpenMesh::Vec3d p2e = m.point(m.to_vertex_handle(h));
						OpenMesh::Vec3d n2e = DGP::cloop_limit_nor(m, m.to_vertex_handle(h));//m.normal(m.to_vertex_handle(h));//
						
						h = m.next_halfedge_handle(m.opposite_halfedge_handle(h));
						OpenMesh::Vec3d p1e = m.point(m.to_vertex_handle(h));
						OpenMesh::Vec3d n1e = DGP::cloop_limit_nor(m, m.to_vertex_handle(h));//m.normal(m.to_vertex_handle(h));//
						
						OpenMesh::Vec3d p1 = m.point(_m.to_vertex_handle(_m.property(fp_hh, f_it)));
						OpenMesh::Vec3d n1 = DGP::cloop_limit_nor(m, _m.to_vertex_handle(_m.property(fp_hh, f_it)));
						//m.normal(_m.to_vertex_handle(_m.property(fp_hh, f_it)));//
						OpenMesh::Vec3d p2 = m.point(_m.to_vertex_handle(_m.next_halfedge_handle(_m.property(fp_hh, f_it))));
						OpenMesh::Vec3d n2 = DGP::cloop_limit_nor(m, _m.to_vertex_handle(_m.next_halfedge_handle(_m.property(fp_hh, f_it))));
						//m.normal(_m.to_vertex_handle(_m.next_halfedge_handle(_m.property(fp_hh, f_it))));//

						OpenMesh::Vec3d pijk[6] = { p0, p1, p2, 0.5*(4.0*p0e-p0-p1), 0.5*(4.0*p1e-p1-p2), 0.5*(4.0*p2e-p0-p2) };//control points.
						OpenMesh::Vec3d nijk[6] = { n0, n1, n2, 0.5*(4.0*n0e-n0-n1), 0.5*(4.0*n1e-n1-n2), 0.5*(4.0*n2e-n0-n2) };//control normals.
						f_pijk_[f_it.handle().idx()] = std::vector<OpenMesh::Vec3d>(pijk, pijk+6);
						f_nijk_[f_it.handle().idx()] = std::vector<OpenMesh::Vec3d>(nijk, nijk+6);
	
						// ----------------------------------------------------------
						h = h_tmp; // for sharp feature.	 
						p0 = subdi_obj.limit_pos(msf.from_vertex_handle(h));
						n0 = subdi_obj.limit_nor(msf.from_vertex_handle(h));// ���ַ������ǲ��Ե�, ��Ϊһ��crease vertex�ڲ�ͬ���������в�ͬ�ķ���.

						p0e = subdi_obj.limit_pos(msf.to_vertex_handle(h));
						n0e = subdi_obj.limit_nor(msf.to_vertex_handle(h));// //

						h = msf.next_halfedge_handle(h); // p0e -> p2e
						p2e = subdi_obj.limit_pos(msf.to_vertex_handle(h));
						n2e = subdi_obj.limit_nor(msf.to_vertex_handle(h));// //

						h = msf.next_halfedge_handle(msf.opposite_halfedge_handle(h));
						p1e = subdi_obj.limit_pos(msf.to_vertex_handle(h));
						n1e = subdi_obj.limit_nor(msf.to_vertex_handle(h));// ;//

						p1 = subdi_obj.limit_pos(_m.to_vertex_handle(_m.property(fp_hh, f_it)));
						n1 = subdi_obj.limit_nor(_m.to_vertex_handle(_m.property(fp_hh, f_it)));
						
						p2 = subdi_obj.limit_pos(_m.to_vertex_handle(_m.next_halfedge_handle(_m.property(fp_hh, f_it))));
						n2 = subdi_obj.limit_nor(_m.to_vertex_handle(_m.next_halfedge_handle(_m.property(fp_hh, f_it))));
						
						OpenMesh::Vec3d pijk_sf[6] = { p0, p1, p2, 0.5*(4.0*p0e-p0-p1), 0.5*(4.0*p1e-p1-p2), 0.5*(4.0*p2e-p0-p2) };//control points.
						OpenMesh::Vec3d nijk_sf[6] = { n0, n1, n2, 0.5*(4.0*n0e-n0-n1), 0.5*(4.0*n1e-n1-n2), 0.5*(4.0*n2e-n0-n2) };//control normals.
						f_pijk_sf_[f_it.handle().idx()] = std::vector<OpenMesh::Vec3d>(pijk_sf, pijk_sf+6);
						f_nijk_sf_[f_it.handle().idx()] = std::vector<OpenMesh::Vec3d>(nijk_sf, nijk_sf+6);
					} else { std::cout << "Error: ";}
				}
			}
			else std::cerr << "Error: QAS,ԭʼģ����Ӧ�������\"name:fp_hh\"��property����¼ÿ������������.\n";

			subdi_obj.detach(); 
			msf.remove_property(vp_type); 
		}

		// ���LOD level�µ�������Ⱦ�Ķ���,����, and index����.
		void evaluate_bezier_points(int _level) {

			if (_level < 0) return;
			int gNumRenderVertices = 0;
			int gNumRenderIndices = 0;
			ver_lod_render_.clear();
			nor_lod_render_.clear();
			idx_lod_render_.clear();
			ver_lod_render_sf_.clear();
			nor_lod_render_sf_.clear();

			for (int fh = 0, size_end = f_pijk_.size(); fh < size_end; ++fh) {
				double r = 1, s = 0, t = 1 - r - s; // barycentric coordinate, ��ʼ��(1,0,0)��������.
				double inc = 1.0 / (_level+1);//barycentric coor�ĵ���step side.
				int numVerticesBeforeAdd = gNumRenderVertices; //Used below for indices, gNumRenderVertices��ʼʱΪ0.
				for (int j = 0; j < (_level+2); ++j) //For each row of vertices, ��LODΪ_levelʱ��, ϸ��֮��������ζ�����_level+2��. 
				{
					s = inc*j;
					t = 1.0 - s - r;

					for (int k = 0; k < j + 1; ++k) //For each vertex in that row
					{
						// Compute new vertex 
						ver_lod_render_.push_back(Pos(fh, r, s, t));
						ver_lod_render_sf_.push_back(Pos(fh, r, s, t, true)); // sharp features

						// ��ÿһ������ĵ�ķ���.
						nor_lod_render_.push_back(Nor(fh, r, s, t).normalize());
						nor_lod_render_sf_.push_back(Nor(fh, r, s, t, true).normalize());

						gNumRenderVertices++;

						//printf("%f %f %f\n", s, r, t);
						s -= inc;
						t += inc;
					}

					r -= inc;
				} // end of each level

				// Generate new indices
				int x = numVerticesBeforeAdd + 0;
				int y = numVerticesBeforeAdd + 1;
				for (int j=0; j<(_level+1); j++) //For each row of triangles
				{
					int indices[3];

					indices[0] = x;
					indices[1] = y;

					// Add all polygons in this row 
					for (int k=0; k<2*j+1; k++)
					{
						//Get new index
						if (k & 0x1)
							indices[2] = ++x;
						else
							indices[2] = ++y;

						//printf("%d %d %d\n", indices[0], indices[1], indices[2]);

						//Add triangle
						idx_lod_render_.push_back(indices[0]);
						idx_lod_render_.push_back(indices[1]);
						idx_lod_render_.push_back(indices[2]); 
						gNumRenderIndices += 3;

						//Swap indices for proper ordering
						if (k & 0x1)
							indices[0] = indices[2];
						else
							indices[1] = indices[2];
					}

					x++;
					y++;
				}
			} // end of each face
			if (ver_lod_render_.size() != nor_lod_render_.size()) 
				std::cerr << "Error: ver_lod_render_.size " << ver_lod_render_.size() << " !=  nor_lod_render_ size " << nor_lod_render_.size() << ".\n";
			if (ver_lod_render_.size() != gNumRenderVertices) 
				std::cerr << "Error: ver_lod_render_.size " << ver_lod_render_.size() << " !=  gNumRenderVertices " << gNumRenderVertices << ".\n";
			if (idx_lod_render_.size() != gNumRenderIndices) 
				std::cerr << "Error: idx_lod_render_.size " << idx_lod_render_.size() << " !=  gNumRenderIndices " << gNumRenderIndices << ".\n";
			std::cout << "QAS level: " << _level << " ver_lod_render_.size():" 
				<< ver_lod_render_.size() << ", idx_lod_render_.size():" << idx_lod_render_.size() << ".\n";
		}

		unsigned int get_ver_lod_render_size() { return ver_lod_render_.size(); }
		unsigned int get_idx_lod_render_size() { return idx_lod_render_.size(); }
		OpenMesh::Vec3d * get_ver_lod_render_data(bool _sf = false) { 
			if (_sf == false ) return &ver_lod_render_[0]; 
			else return &ver_lod_render_sf_[0]; 
		}
		OpenMesh::Vec3d * get_nor_lod_render_data(bool _sf = false) { 
			if (_sf == false) return &nor_lod_render_[0]; 
			else return &nor_lod_render_sf_[0];
		}
		unsigned int * get_idx_lod_render_data() { return &idx_lod_render_[0]; }
		const std::vector<OpenMesh::Vec3d> & f_pijk(const TriMesh::FaceHandle _fh, bool _sf) const {
			if (_sf == false) return f_pijk_[_fh.idx()];
			else return f_pijk_sf_[_fh.idx()];
		}
		OpenMesh::Vec3d pos(TriMesh::FHandle _fh, double _u, double _v, double _w, bool _sf) {
			return Pos(_fh.idx(), _u, _v, _w, _sf);
		}
		OpenMesh::Vec3d dPos_du(TriMesh::FHandle _fh, double _u, double _v, double _w, bool _sf) { // ��u������
			int fh = _fh.idx();
			if (_sf == false) return f_pijk_[fh][0]*2*_u + 
				f_pijk_[fh][3]*2*_v + f_pijk_[fh][5]*2*(1-2*_u-_v) + 
				OpenMesh::Vec3d(0,0,0) - f_pijk_[fh][4]*2*_v - f_pijk_[fh][2]*2*(1-_u-_v);
			else return f_pijk_sf_[fh][0]*2*_u + 
				f_pijk_sf_[fh][3]*2*_v + f_pijk_sf_[fh][5]*2*(1-2*_u-_v) + 
				OpenMesh::Vec3d(0,0,0) - f_pijk_sf_[fh][4]*2*_v - f_pijk_sf_[fh][2]*2*(1-_u-_v);
		}
		OpenMesh::Vec3d dPos_dv(TriMesh::FHandle _fh, double _u, double _v, double _w, bool _sf) {//��v������
			int fh = _fh.idx();
			if (_sf == false) return OpenMesh::Vec3d(0,0,0) + 
						f_pijk_[fh][3]*2*_u - f_pijk_[fh][5]*2*_u + 
				f_pijk_[fh][1]*2*_v + f_pijk_[fh][4]*2*(1-_u-2*_v) - f_pijk_[fh][2]*2*(1-_u-_v);
			else OpenMesh::Vec3d(0,0,0) + 
						f_pijk_sf_[fh][3]*2*_u - f_pijk_sf_[fh][5]*2*_u + 
				f_pijk_sf_[fh][1]*2*_v + f_pijk_sf_[fh][4]*2*(1-_u-2*_v) - f_pijk_sf_[fh][2]*2*(1-_u-_v);
		}
		OpenMesh::Vec3d dPos_dw(TriMesh::FHandle _fh, double _u, double _v, double _w, bool _sf) {//��w������
			int fh = _fh.idx();
			if (_sf == false) return f_pijk_[fh][0]*(-2)*(1-_v-_w) + 
						f_pijk_[fh][3]*(-2)*_v + f_pijk_[fh][5]*2*(1-_v-2*_w) + 
				OpenMesh::Vec3d(0,0,0) + f_pijk_[fh][4]*2*_v + f_pijk_[fh][2]*2*_w;
			else return f_pijk_sf_[fh][0]*(-2)*(1-_v-_w) + 
						f_pijk_sf_[fh][3]*(-2)*_v + f_pijk_sf_[fh][5]*2*(1-_v-2*_w) + 
				OpenMesh::Vec3d(0,0,0) + f_pijk_sf_[fh][4]*2*_v + f_pijk_sf_[fh][2]*2*_w;
		}
	private:
		std::vector<std::vector<OpenMesh::Vec3d> > f_pijk_; // ÿһ�����6��������, p200, p020, p002, p110, p011, p101 
		std::vector<std::vector<OpenMesh::Vec3d> > f_nijk_;
		
		//��Ҫ����ǰֱ��ȷ��u = 1 - v - w
		// ֱ����6��pos/nor������2��bezier triangle�Ĺ�ʽ��ֵ.
		OpenMesh::Vec3d Pos(int fh, double u, double v, double w, bool _sf = false) {
			if (_sf == false) return f_pijk_[fh][0]*u*u + 
				f_pijk_[fh][3]*2*u*v + f_pijk_[fh][5]*2*u*w + 
				f_pijk_[fh][1]*v*v + f_pijk_[fh][4]*2*v*w + f_pijk_[fh][2]*w*w;
			else return f_pijk_sf_[fh][0]*u*u + 
				f_pijk_sf_[fh][3]*2*u*v + f_pijk_sf_[fh][5]*2*u*w + 
				f_pijk_sf_[fh][1]*v*v + f_pijk_sf_[fh][4]*2*v*w + f_pijk_sf_[fh][2]*w*w;
		}
		OpenMesh::Vec3d Nor(int fh, double u, double v, double w, bool _sf = false) {
			if (_sf == false ) return f_nijk_[fh][0]*u*u + 
				f_nijk_[fh][3]*2*u*v + f_nijk_[fh][5]*2*u*w + 
				f_nijk_[fh][1]*v*v + f_nijk_[fh][4]*2*v*w + f_nijk_[fh][2]*w*w;
			else return f_nijk_sf_[fh][0]*u*u + 
				f_nijk_sf_[fh][3]*2*u*v + f_nijk_sf_[fh][5]*2*u*w + 
				f_nijk_sf_[fh][1]*v*v + f_nijk_sf_[fh][4]*2*v*w + f_nijk_sf_[fh][2]*w*w;
		}

		// ���������Ⱦ��
		std::vector<OpenMesh::Vec3d> ver_lod_render_;
		std::vector<OpenMesh::Vec3d> nor_lod_render_;
		std::vector<unsigned int> idx_lod_render_; //no matter consider sharp features or not, will be same.

		// add the supprot to sharp features(Hoppe94 Loop subdivision surfaces); 
		std::vector< std::vector<OpenMesh::Vec3d > > f_pijk_sf_;// 6 control points of one sharp feature patch.
		std::vector< std::vector<OpenMesh::Vec3d > > f_nijk_sf_;// 6 control normals of one sharp feature patch.
		std::vector<OpenMesh::Vec3d> ver_lod_render_sf_;
		std::vector<OpenMesh::Vec3d> nor_lod_render_sf_;
	}; // end of class 

} // end of namespace DGP

#endif // dgp_quadratic_bezier_triangles_h