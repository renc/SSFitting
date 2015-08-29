// 1. ��Ӧcurvature�ļ����ǲο�CourseExamples������.

#ifndef dgp_common_laplaciant_h
#define dgp_common_laplaciant_h



#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))

namespace DGP {

	template <class MeshT, class RealT = double>
	class LaplacianT {
	public:
		LaplacianT(MeshT &_mesh):mesh_(_mesh) {
			mesh_.add_property(eweight_);
			mesh_.add_property(varea_);
			mesh_.add_property(vweight_);
			mesh_.add_property(vlaplace_vector_);
			mesh_.add_property(vmean_curvature_);
			mesh_.add_property(vmean_curvature_normal_); 

			mesh_.request_vertex_colors();
		}
		~LaplacianT() {
			mesh_.remove_property(eweight_);
			mesh_.remove_property(varea_);
			mesh_.remove_property(vweight_);
			mesh_.remove_property(vlaplace_vector_);
			mesh_.remove_property(vmean_curvature_);
			mesh_.remove_property(vmean_curvature_normal_); 
		}
	
		void calc_wights() {
			std::cout << "-";
			//��һ����, ��ߵ�Ȩ��, Ҳ�����������ص�cotangent weight.
			MeshT::EdgeIter          e_it, e_end(mesh_.edges_end());
			MeshT::HalfedgeHandle    h0, h1, h2;
			MeshT::VertexHandle      v0, v1;
			MeshT::Point             p0, p1, p2, d0, d1;
			RealT					w, b(0.99);

			for (e_it=mesh_.edges_begin(); e_it!=e_end; ++e_it)
			{
				w  = 0.0; // �����������wӦ��������ߵ�cotangent weight
				// Fig 8."Surface Parameterization A Tutorial and Survey_FH05_M.S.Floater and K.Hormann" 
				//      v0
				//   /  | \
				// v2   |  v2
				//   \  |  /
				//    \v1 /     
				h0 = mesh_.halfedge_handle(e_it.handle(), 0);
				v0 = mesh_.to_vertex_handle(h0);
				p0 = mesh_.point(v0);

				h1 = mesh_.halfedge_handle(e_it.handle(), 1);
				v1 = mesh_.to_vertex_handle(h1);
				p1 = mesh_.point(v1);

				if (!mesh_.is_boundary(h0))
				{
					h2 = mesh_.next_halfedge_handle(h0);
					p2 = mesh_.point(mesh_.to_vertex_handle(h2));
					d0 = (p0 - p2).normalize();
					d1 = (p1 - p2).normalize();
					w += 1.0 / tan(acos(max(-b, 
						min(b, (d0|d1)))));//�����������,acos�����������ļн�,1/tan = cot 
				}

				if (!mesh_.is_boundary(h1))
				{
					h2 = mesh_.next_halfedge_handle(h1);
					p2 = mesh_.point(mesh_.to_vertex_handle(h2));
					d0 = (p0 - p2).normalize();
					d1 = (p1 - p2).normalize();
					w += 1.0 / tan(acos(max(-b, min(b, (d0|d1)))));
				}

				// force weights to be non-negative for higher robustness
				w = max(w, (double)0.0);

				mesh_.property(eweight_, e_it) = w;
			} // end of for every edge

			// �ڶ�����,��ÿһ����������ܵ����
			MeshT::VertexIter        v_it, v_end(mesh_.vertices_end());
			MeshT::VertexFaceIter    vf_it;
			MeshT::FaceVertexIter    fv_it;
			RealT area;
			// compute the area associated to a vertex. ����һ�������������, ��һ������ֵ
			// The simplest such area is the barycentric area, 
			// given by 1/3 the area of all incident triangles. ����ʹ�õľ����ⷽ��.
			// A better measure of the area is the Voronoi area associated to the vertex.
			for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
			{
				area = 0.0;

				for (vf_it=mesh_.vf_iter(v_it); vf_it; ++vf_it)
				{
					fv_it = mesh_.fv_iter(vf_it);

					const MeshT::Point& P = mesh_.point(fv_it);  ++fv_it;
					const MeshT::Point& Q = mesh_.point(fv_it);  ++fv_it;
					const MeshT::Point& R = mesh_.point(fv_it);

					area += ((Q-P)%(R-P)).norm() * 0.5f * 0.3333f;//���������������֮һ
				}

				mesh_.property(varea_, v_it) = (fabs(area)>FLT_MIN ? 1.0 / (2.0 * area) : 0.0);
			}
			// ��������
			MeshT::HalfedgeHandle    h;
			MeshT::EdgeHandle        e;
			MeshT::VertexVertexIter  vv_it;
			MeshT::Point             laplace(0.0, 0.0, 0.0);
			RealT					 ww = 0.0;

			for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
			{	// ��ʾ�˴�һ��������ô��������ܵİ�ߺͱ�
				laplace = MeshT::Point(0,0,0);
				ww = 0.0;

				if (!mesh_.is_boundary(v_it.handle()))
				{
					for (vv_it=mesh_.vv_iter(v_it); vv_it; ++vv_it)
					{
						h = vv_it.current_halfedge_handle();
						e = mesh_.edge_handle(h);

						laplace += (double)mesh_.property(eweight_, e) * (mesh_.point(vv_it) - mesh_.point(v_it));
						ww += mesh_.property(eweight_, e);
					}
					mesh_.property(vweight_, v_it) = ww;
					mesh_.property(vlaplace_vector_, v_it) = laplace / ww; // laplace vector; 
					mesh_.property(vmean_curvature_normal_, v_it) = laplace * mesh_.property(varea_, v_it.handle());
					mesh_.property(vmean_curvature_, v_it) = mesh_.property(vmean_curvature_normal_, v_it).norm();//����������Ϊ���������ֵ, �ⶫ��������mean curvature.
					
				}
			}
		}

		void color_coding() //����mesh_.property(vmean_curvature_, v_it)����ɫ.
		{

			MeshT::VertexIter  v_it, v_end(mesh_.vertices_end());
			MeshT::Scalar      curv, min_curv(FLT_MAX), max_curv(-FLT_MAX);
			OpenMesh::Vec3f	col(255, 255, 255);

			// put all curvature values into one array
			std::vector<MeshT::Scalar> curv_values;
			curv_values.reserve(mesh_.n_vertices());
			for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it) {
				curv_values.push_back(mesh_.property(vmean_curvature_, v_it));
			}

			// discard upper and lower 5%
			unsigned int n = curv_values.size()-1;
			unsigned int i = n / 20;
			std::sort(curv_values.begin(), curv_values.end());
			min_curv = curv_values[i];
			max_curv = curv_values[n-1-i];

			// define uniform color intervalls [v0,v1,v2,v3,v4]
			MeshT::Scalar v0, v1, v2, v3, v4;
			v0 = min_curv + 0.0/4.0 * (max_curv - min_curv);
			v1 = min_curv + 1.0/4.0 * (max_curv - min_curv);
			v2 = min_curv + 2.0/4.0 * (max_curv - min_curv);
			v3 = min_curv + 3.0/4.0 * (max_curv - min_curv);
			v4 = min_curv + 4.0/4.0 * (max_curv - min_curv); 

			// map curvatures to colors
			for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
			{
				curv = mesh_.property(vmean_curvature_, v_it); 
				col = OpenMesh::Vec3f(252,255,255);  
				unsigned char u;

				if (curv < v0)
				{
					col = OpenMesh::Vec3f(0, 0, 255);
				}
				else if (curv > v4) 
				{
					col = OpenMesh::Vec3f(255, 0, 0);
				}

				else if (curv <= v2) 
				{
					if (curv <= v1) // [v0, v1]
					{
						u = (unsigned char) (255.0 * (curv - v0) / (v1 - v0));
						col = OpenMesh::Vec3f(0, u, 255);
					}      
					else // ]v1, v2]
					{
						u = (unsigned char) (255.0 * (curv - v1) / (v2 - v1));
						col = OpenMesh::Vec3f(0, 255, 255-u);
					}
				}
				else 
				{
					if (curv <= v3) // ]v2, v3]
					{
						u = (unsigned char) (255.0 * (curv - v2) / (v3 - v2));
						col = OpenMesh::Vec3f(u, 255, 0);
					}
					else // ]v3, v4]
					{
						u = (unsigned char) (255.0 * (curv - v3) / (v4 - v3));
						col = OpenMesh::Vec3f(255, 255-u, 0);
					}
				} 
				mesh_.set_color(v_it, (OpenMesh::Vec3uc)col);//ǰ����request_vertex_colors()
			} // end of for each vertex.
		} // end of color_coding().

		// const menber function,
		// the const object , a pointer or reference to a const object may be used to call only const member functions.
		RealT edge_weight(typename MeshT::HalfedgeHandle _hh, bool _is_normalized = false) const { 
			if (_is_normalized == false) {
				return mesh_.property(eweight_, mesh_.edge_handle(_hh)); // wij = cotAij + cotBij
			} else {
				return mesh_.property(eweight_, mesh_.edge_handle(_hh)) / mesh_.property(vweight_, mesh_.from_vertex_handle(_hh));
			}			
		}
		RealT vertex_area(typename MeshT::VertexHandle _vh) const {
			return mesh_.property(varea_, _vh);
		}
		RealT vertex_weight(typename MeshT::VertexHandle _vh, bool _is_normalized = false) const {
			return _is_normalized == true ? 1 : mesh_.property(weight, _vh);
		}
		typename MeshT::Point vertex_laplace_vector(typename MeshT::VertexHandle _vh) {
			return mesh_.property(vlaplace_vector_, _vh);
		}
		RealT vertex_mean_curvature(typename MeshT::VertexHandle _vh) const {
			return mesh_.property(vmean_curvature_, _vh);
		}
		typename MeshT::Point vertex_mean_curvature_normal(typename MeshT::VertexHandle _vh) const {
			return mesh_.property(vmean_curvature_normal_, _vh);
		}
	private:
		MeshT &mesh_;
		OpenMesh::EPropHandleT<RealT> eweight_;//���������ı�������contangent weight, wij = weight(eij) = weight(vi, vj)
		OpenMesh::VPropHandleT<RealT> varea_;  //�������ܵ�����ĵ��� 1/Ai
		OpenMesh::VPropHandleT<RealT> vweight_;//���� wi = sum(wij)
		OpenMesh::VPropHandleT<OpenMesh::VectorT<RealT, 3> > vlaplace_vector_;//��������laplace vector����ֵ
		OpenMesh::VPropHandleT<OpenMesh::VectorT<RealT, 3> > vmean_curvature_normal_;//��������mean curvature normal����ֵ
		OpenMesh::VPropHandleT<RealT> vmean_curvature_;//��������mean curvature����ֵ 
	};
} // end of namespace DGP

#endif //dgp_common_laplaciant_h