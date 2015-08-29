
// include 
//#include "../../../src/CourseExamples/01-MeshViewer/gl.hh"
#include "DecimationModel.h" 
#include <OpenMesh/Tools/Utils/Timer.hh> //OpenMesh::Utils::Timer  timer; timer.start(); ...timer.stop(); std::cerr << "done (" << timer.as_string() << ")\n";

//#include "../../../src/CourseExamples/04-Fairing/TaucsSolver.hh" // renc 

#include <algorithm> // for max_element algorithm

#include <assert.h>
#include <fstream>

//#include "../../../../taucs/taucsaddon.h" //��һ��ʹ�÷���, ���ɵ�exe�ļ���. // renc 

#include "lib_IDSSsrc/bnd_loop.h"  //renc


// ----------------
void DecimationModel::try_bounding_process() {
	std::cout << "\n";
	std::cout << "========��ʼ����bounding    ===========\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.
	
	if (upperbounding_mesh_.n_vertices())  upperbounding_mesh_.clear(); //if (upperbounding_mesh_.empty() == false) upperbounding_mesh_.clear();//renc 20150829
	if (lowerbounding_mesh_.n_vertices()) lowerbounding_mesh_.clear(); //if (lowerbounding_mesh_.empty() == false) lowerbounding_mesh_.clear();
	upperbounding_mesh_ = controlmesh_forbounding_;
	lowerbounding_mesh_ = controlmesh_forbounding_;

	for (TriMesh::FaceIter f_it(controlmesh_forbounding_.faces_begin()), f_end(controlmesh_forbounding_.faces_end()); f_it != f_end; ++f_it) {
		// ����ÿһ����, ���ҳ��ײ���6���Ǹ�����
		TriMesh::VertexHandle v0, v1, v2;
		int valence_v0 = 6; 
		TriMesh::HalfedgeHandle h0 = controlmesh_forbounding_.halfedge_handle(f_it);
		for (int i = 0; i < 3; ++i) {
			if (controlmesh_forbounding_.valence(controlmesh_forbounding_.to_vertex_handle(h0)) != 6) {
				v0 = controlmesh_forbounding_.to_vertex_handle(h0); 
				valence_v0 = controlmesh_forbounding_.valence(v0);
				break;
			}
			h0 = controlmesh_forbounding_.next_halfedge_handle(h0);
		}
		if (v0.is_valid() == false) {
			v0 = controlmesh_forbounding_.to_vertex_handle(h0);
			valence_v0 = 6;
		}
		TriMesh::HHandle h1 = controlmesh_forbounding_.next_halfedge_handle(h0), h2 = controlmesh_forbounding_.prev_halfedge_handle(h0);
		v1 = controlmesh_forbounding_.to_vertex_handle(h1); v2 = controlmesh_forbounding_.from_vertex_handle(h0);
		if (controlmesh_forbounding_.valence(v1) != 6 && controlmesh_forbounding_.valence(v2) != 6) {
			std::cout << "Error: v1, v2 valence should be 6.\n";
		}
		
		// ����v0,v1,v2�������һ����ĵ���.
		int k = valence_v0 + 6;//k����
		std::vector<TriMesh::VertexHandle> vec(k);//vec[0] ...vec[n+5];
		vec[0] = v0; vec[1] = v1; vec[2] = v2;
		TriMesh::HalfedgeHandle h3 = controlmesh_forbounding_.opposite_halfedge_handle(controlmesh_forbounding_.prev_halfedge_handle(
			controlmesh_forbounding_.opposite_halfedge_handle(controlmesh_forbounding_.prev_halfedge_handle(controlmesh_forbounding_.opposite_halfedge_handle(controlmesh_forbounding_.next_halfedge_handle(h0))))));
		vec[3] = controlmesh_forbounding_.to_vertex_handle(h3);
		TriMesh::HalfedgeHandle h4 = controlmesh_forbounding_.opposite_halfedge_handle(controlmesh_forbounding_.prev_halfedge_handle(h3));
		vec[4] = controlmesh_forbounding_.to_vertex_handle(h4);
		vec[5] = controlmesh_forbounding_.to_vertex_handle(controlmesh_forbounding_.opposite_halfedge_handle(controlmesh_forbounding_.prev_halfedge_handle(h4)));
		assert(vec[5] == controlmesh_forbounding_.to_vertex_handle(controlmesh_forbounding_.next_halfedge_handle(controlmesh_forbounding_.opposite_halfedge_handle(h2))));
		TriMesh::HHandle h6 = controlmesh_forbounding_.opposite_halfedge_handle(controlmesh_forbounding_.prev_halfedge_handle(
			controlmesh_forbounding_.opposite_halfedge_handle(controlmesh_forbounding_.prev_halfedge_handle(controlmesh_forbounding_.opposite_halfedge_handle(h2)))));
		vec[6] = controlmesh_forbounding_.to_vertex_handle(h6);
		vec[7] = controlmesh_forbounding_.to_vertex_handle(controlmesh_forbounding_.opposite_halfedge_handle(controlmesh_forbounding_.prev_halfedge_handle(h6)));
		TriMesh::HalfedgeHandle h8 = controlmesh_forbounding_.opposite_halfedge_handle(controlmesh_forbounding_.prev_halfedge_handle(controlmesh_forbounding_.opposite_halfedge_handle(h0)));
		TriMesh::HalfedgeHandle hn5 = controlmesh_forbounding_.next_halfedge_handle(controlmesh_forbounding_.opposite_halfedge_handle(h1)); // hn+5
		int i = 8;
		for (TriMesh::HHandle h = h8; h != hn5; h = controlmesh_forbounding_.opposite_halfedge_handle(controlmesh_forbounding_.prev_halfedge_handle(h))) {
			vec[i++] = controlmesh_forbounding_.to_vertex_handle(h);
		}
		vec[i] = controlmesh_forbounding_.to_vertex_handle(hn5);
		assert(i == (valence_v0 + 5));

		double *coeff = new double[k*3]; 
		double *upper = new double[9]; //v0, v1, v2�ĸ���x,y,z����.
		double *lower = new double[9];
		for (int i = 0, end = k; i < k; ++i) { //k����, 0, 1, 2�ֱ�x, y, z��������.
			coeff[i*3+0] = controlmesh_forbounding_.point(vec[i])[0]; 
			coeff[i*3+1] = controlmesh_forbounding_.point(vec[i])[1]; 
			coeff[i*3+2] = controlmesh_forbounding_.point(vec[i])[2]; 
		} 
		for (int i = 0; i < 9; ++i) {
			upper[i] = 0; lower[i] = 0;
		}
		for(int m = 0; m < 3; ++m) { //m:0,1,2�ֱ��Ӧxyz��������
			BoundLoop(valence_v0, &coeff[m], 3, &upper[m], &lower[m], 3);
		} 
		upperbounding_mesh_.set_point(v0, OpenMesh::Vec3d(upper[0], upper[1], upper[2]));
		upperbounding_mesh_.set_point(v1, OpenMesh::Vec3d(upper[3], upper[4], upper[5]));
		upperbounding_mesh_.set_point(v2, OpenMesh::Vec3d(upper[6], upper[7], upper[8]));
		lowerbounding_mesh_.set_point(v0, OpenMesh::Vec3d(lower[0], lower[1], lower[2]));
		lowerbounding_mesh_.set_point(v1, OpenMesh::Vec3d(lower[3], lower[4], lower[5]));
		lowerbounding_mesh_.set_point(v2, OpenMesh::Vec3d(lower[6], lower[7], lower[8]));
				
	}
	std::cout << "========��������bounding    ===========\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.

}
