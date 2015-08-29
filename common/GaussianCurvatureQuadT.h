
#ifndef _dgp_common_guassian_curvature_quad_h
#define _dgp_common_guassian_curvature_quad_h

#include "OpenMeshMeshType.h"

namespace DGP {
// 实现的内容是: 
// Angle deficit approximation of Gaussian curvature and its convergence over
// quadrilateral meshes, Dan Liu, Guoliang Xu, 2007.
// K(pi) = 4/A(pi) * [2PI - sum_n((theta_j + alpha_j + beta_j) * 0.5)], 
// for vertex i, which valence is n.
// Notes: 1. quadrilateral mesh. 2. no boundary!!
template <class MeshT>
class GaussianCurvatureQuadT {
public:
	GaussianCurvatureQuadT(MeshT &_m): mesh_(_m) {
		mesh_.add_property(vp_cur_); 
	}
	~GaussianCurvatureQuadT() {
		mesh_.remove_property(vp_cur_);
	}

	bool calc_wights() {
		for (MeshT::ConstFaceIter f_it(mesh_.faces_begin()), f_end(mesh_.faces_end());
			f_it != f_end; ++f_it) {
				 
		} 
		for (MeshT::ConstVertexIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end());
			v_it != v_end; ++v_it) {
				//v3<----------- v2
				//|     h2     /|
				//| h3          | h1
				//|             |
				//|/___________\|
				//v0     h0     v1
				double area = 0;
				double angle = 0; // in radian.
				for (MeshT::ConstVertexOHalfedgeIter voh_it(mesh_, v_it); voh_it; ++voh_it) {
					MeshT::HHandle h0 = voh_it.handle();
					MeshT::HHandle h1 = mesh_.next_halfedge_handle(h0);
					MeshT::HHandle h2 = mesh_.next_halfedge_handle(h1);
					MeshT::HHandle h3 = mesh_.next_halfedge_handle(h2);
					if (h0 != mesh_.next_halfedge_handle(h3)) {
						std::cerr << "Error: maybe there's boundary??\n";
						return false; 
					} 
					MeshT::VHandle v0 = v_it.handle();
					MeshT::VHandle v1 = mesh_.to_vertex_handle(h0); 
					MeshT::VHandle v2 = mesh_.to_vertex_handle(h1); 
					MeshT::VHandle v3 = mesh_.to_vertex_handle(h2); 
					area += cross(mesh_.point(v2)-mesh_.point(v1), mesh_.point(v1)-mesh_.point(v0)).norm() * 0.5;
					area += cross(mesh_.point(v0)-mesh_.point(v3), mesh_.point(v3)-mesh_.point(v2)).norm() * 0.5;

					angle += acos(dot((mesh_.point(v3)-mesh_.point(v0)).normalize(),
						(mesh_.point(v1)-mesh_.point(v0)).normalize())); // angle v3v0v1,
					angle += acos(dot((mesh_.point(v3)-mesh_.point(v0)).normalize(),
						(mesh_.point(v2)-mesh_.point(v0)).normalize())); // angle v3v0v2,
					angle += acos(dot((mesh_.point(v2)-mesh_.point(v0)).normalize(),
						(mesh_.point(v1)-mesh_.point(v0)).normalize())); // angle v2v0v1,

				}
				mesh_.property(vp_cur_, v_it ) = 4.0 / area * (2*M_PI - angle * 0.5);
			
		}
		return true;
	} // end of fuc calc_weights();

	void color_coding() {
		MeshT::VertexIter  v_it, v_end(mesh_.vertices_end());
		MeshT::Scalar      curv, min_curv(FLT_MAX), max_curv(-FLT_MAX);
		OpenMesh::Vec3f	col(255, 255, 255);

		// put all curvature values into one array
		std::vector<MeshT::Scalar> curv_values;
		curv_values.reserve(mesh_.n_vertices());
		for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it) {
			curv_values.push_back(mesh_.property(vp_cur_, v_it));
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
			curv = mesh_.property(vp_cur_, v_it); 
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
			mesh_.set_color(v_it, (OpenMesh::Vec3uc)col);
		} // end of for each vertex. 
	} // end of fuc color_coding().

private:
	MeshT &mesh_;//
	OpenMesh::VPropHandleT<double> vp_cur_;//the area of vertex
	 

}; // end of class GaussianCurvatureQuadT

} // end of namespace DGP
#endif //_dgp_common_guassian_curvature_quad_h