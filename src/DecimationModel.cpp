// ***************************************************************
//  
//  DecimationModel 
//  Copyright (C) 2007 - by rencanjiang All Rights Reserved
//  rencanjiang@163.com	
//	-------------------------------------------------------------
//  
//  version: 1.0
//  data: 12/10/2007
// ***************************************************************
//  
//  File description:
// ***************************************************************

// ============================================================================
// include 
#include "DecimationModel.h" 
//#include "../../../src/CourseExamples/04-Fairing/TaucsSolver.hh"

#include "GUIFrame.h" 

#include <algorithm> // for max_element algorithm
#include <assert.h>
#include <fstream>
 

// ============================================================================
// declare the static variable
DecimationModel *DecimationModel::currentDV = NULL;

DecimationModel::DecimationModel()
: MeshModel()
{ 
	vertex_optimal_placement = false; 
	
	meshtype_toberendered_ = ORIGINAL;

	is_view_fitting_error_of_mesh2_ = false;
	is_view_fitting_error_of_sm_ = false;


	mesh_.add_property(vp_eval_pos_); //为后面的evaluation做prepare.
	mesh_.add_property(vp_eval_);
	mesh_.add_property(vp_closest_pos_);
	mesh_.add_property(vp_closest_);
	mesh_.add_property(vp_project_pos_);
	mesh_.add_property(vp_project_);
	for (TriMesh::VIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++v_it) {
		mesh_.property(vp_eval_, v_it) = false;
		mesh_.property(vp_closest_, v_it) = false;
	}
	mesh_.add_property(vp_fitting_error_);

	fitting_terminate_ = false;
	DecimationModel::currentDV = this; 
}

// ============================================================================
// 简化部分
// 平滑部分
// 中点采样部分
// 
// ============================================================================
 
void DecimationModel::draw() {
	
	if (meshtype_toberendered_ == ORIGINAL) {
	
		
		MeshModel::draw();
	} else if (meshtype_toberendered_ == SIMPLIFIED || meshtype_toberendered_ == RESAMPLE_MIDPOINT) {
		// 第一步, 要是做了参数化的话先显示所有顶点, 包括那些被deleted但是参数化到base complex面上的点.
		for (TriMesh::VIter v_it(mesh2_.vertices_begin()), v_end(mesh2_.vertices_end()); v_it != v_end; ++v_it) {
			if ( v_it.handle() == test_vertex_||v_it.handle().idx() == 69000 ) {//|| v_it.handle().idx() == 1493
				glColor3f(0, 1, 0); glPointSize(10); 
				glBegin(GL_POINTS);
				GL::glVertex(mesh2_.point(v_it));
				glEnd();glPointSize(1);
			}
			if (v_it.handle() == test_vertex1_ || v_it.handle().idx() == 56656 || v_it.handle().idx() == 66696 ||v_it.handle().idx() == 203340) {//
				glColor3f(0, 1, 1); glPointSize(10); 
				glBegin(GL_POINTS);
				GL::glVertex(mesh2_.point(v_it));
				glEnd();glPointSize(1);
			}
		}
		for (TriMesh::VIter v_it(initial_control_mesh_.vertices_begin()), v_end(initial_control_mesh_.vertices_end()); v_it != v_end; ++v_it) {
			if (v_it.handle() == test_vertex_ ) { //显示新分裂的顶点
				glColor3f(0, 1, 0); glPointSize(10); 
				glBegin(GL_POINTS);
				GL::glVertex(initial_control_mesh_.point(v_it));
				glEnd();glPointSize(1);
			}
		}

	 
			//int fandisk_arr[] = { //这些是一部分特征边
			//	9, 269, 265,
			//	261, 257, 254,
			//	250, 246, 241,
			//	237, 233, 229,
			//	225, 221, 217,
			//	214, 210, 206,
			//	//202, 995, 9762

			//	9922, 9925, 9928,
			//	9931, 9934, 9939,
			//	9943, 9947, 9951,
			//	9954, 9958, 9962,
			//	9966, 9970, 9974,
			//	9978, 9982, 9986,
			//	9990, 9994, 9998,
			//	10002, 10006, 10011,
			//	10014, 784, 788

			//}; 
			//std::vector<int> fandisk(fandisk_arr, fandisk_arr+18+27);
			//for (int i = 0; i < 45; ++i) {
			//	glColor3f(1, 0, 1); // 显示特征边
			//	glLineWidth(3.5);
			//	glBegin(GL_LINES);
			//	GL::glVertex(mesh2_.point(mesh2_.to_vertex_handle(mesh2_.halfedge_handle(mesh2_.edge_handle(fandisk_arr[i]), 0))));
			//	GL::glVertex(mesh2_.point(mesh2_.to_vertex_handle(mesh2_.halfedge_handle(mesh2_.edge_handle(fandisk_arr[i]), 1))));
			//	glEnd();
			//	glLineWidth(1.0);

			//}
		 
		if (do_parameterization_) {
			//(1) 先显示所有顶点, 包括那些被deleted但是参数化到base complex面上的点.
			// 但是后来我觉得这可以省略, 因为下面显示mesh_collapsed_时候就会显示这些被参数化的点了.
			glDisable(GL_LIGHTING); 
			
			for (TriMesh::EIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end()); e_it != e_end; ++e_it) {
				if (mesh_.status(e_it).deleted() == false) {
					if (mesh_.status(e_it).feature()) {
						glColor3f(1, 0, 1); // 显示特征边
						glLineWidth(3.5);
						glBegin(GL_LINES);
						GL::glVertex(mesh_.point(mesh_.to_vertex_handle(mesh_.halfedge_handle(e_it.handle(), 0))));
						GL::glVertex(mesh_.point(mesh_.to_vertex_handle(mesh_.halfedge_handle(e_it.handle(), 1))));
						glEnd();
						glLineWidth(1.0);/**/

					}
					if (meshtype_toberendered_ == RESAMPLE_MIDPOINT) {
						if (mesh_.property(empl_, e_it.handle()) == true) { //已经求出中点的了.
							glColor3f(0, 1, 0);//黄色 显示那些采样点
							glPointSize(4);
							//glBegin(GL_POINTS);
							//GL::glVertex(mesh_.property(emp_, e_it.handle()));
							//glEnd();
							glPointSize(1);
							glPushMatrix();
							glTranslatef( mesh_.property(emp_, e_it.handle())[0], mesh_.property(emp_, e_it.handle())[1], mesh_.property(emp_, e_it.handle())[2] );
							glutSolidSphere(bb.diagonal_lenght*0.003, 20, 20);
							glPopMatrix();
						} else { // 还没有求出中点的
							glColor3f(0.f, 1.f, 0.f);//绿色 显示那些没有对中点采样成功的边
							glLineWidth(6.5);
							glBegin(GL_LINES);
							GL::glVertex(mesh_.point(mesh_.to_vertex_handle(mesh_.halfedge_handle(e_it, 0))));
							GL::glVertex(mesh_.point(mesh_.to_vertex_handle(mesh_.halfedge_handle(e_it, 1))));
							glEnd();
							glLineWidth(1.0);
						} 
					} // 如果是中点采样的话显示采样点或不成功的边.				
				}
			}
			
			// (4) 显示被参数化的网格, 作为参考.用hiden line线框模式显示
			if (mesh_collapsed_render_or_not_) {
				DGP::hidderline_render(mesh_collapsed_, OpenMesh::Vec3f(0.0, 0.0, 1.0));				
			}
			 
			/////////////////////////////
			//draw_mode_ = HIDDENLINE;//LINEANDSOLIDFLAT; //NONE; // //这里设置了的话用户就不可以选择了. 
		} //end of if (do_parameterization_).
		// 第二步, 正常显示.
		MeshModel::draw();
	} else if (meshtype_toberendered_ == REFINED_SIMPLIFIED) {	
		//glBegin(GL_POINTS);
		//GL::glColor(OpenMesh::Vec3f(1, 0, 0));
		//glPointSize(3.0);
		//for (TriMesh::FIter f_it(mesh_.faces_begin()), f_end(mesh_.faces_end()); f_it != f_end; ++f_it) {
		//	for (std::vector<TriMesh::VHandle>::const_iterator it(mesh_.property(fvset_, f_it.handle()).begin()), it_end(mesh_.property(fvset_, f_it.handle()).end()); 
		//		it != it_end; ++it) {
		//			if (32 == mesh_.property(vf_of_rsm_, *it).idx())
		//			GL::glVertex(refined_simplified_mesh_.point(mesh_.property(vf0_of_rsm_, *it)) * mesh_.property(vbc_of_rsm_, *it)[0]
		//			+ refined_simplified_mesh_.point(mesh_.property(vf1_of_rsm_, *it)) * mesh_.property(vbc_of_rsm_, *it)[1]
		//			+ refined_simplified_mesh_.point(mesh_.property(vf2_of_rsm_, *it)) * mesh_.property(vbc_of_rsm_, *it)[2]); 

		//	}
		//}glPointSize(1);
		//glEnd();
		MeshModel::draw();
	} else if (meshtype_toberendered_ == INITIAL_CONTROL ||meshtype_toberendered_ == LOOP_ONE
		|| meshtype_toberendered_ == LIMIT) {
			
		if (upperbounding_mesh_.n_vertices()) { 
		//if (upperbounding_mesh_.empty() == false) {// renc: 20150829
			DGP::hidderline_render(upperbounding_mesh_, OpenMesh::Vec3f(0, 0, 1));
			DGP::hidderline_render(lowerbounding_mesh_, OpenMesh::Vec3f(0, 1, 0));			
		}
		glBegin(GL_LINES);
		glLineWidth(3);
		glColor3f(1, 0, 0);
		for (int i = 0; i < laplacian_.size(); ++i) {
			TriMesh::VHandle vh = refined_simplified_mesh_.vertex_handle(i);
			GL::glVertex(refined_simplified_mesh_.point(vh));
			GL::glVertex(refined_simplified_mesh_.point(vh) + laplacian_[i]); 
		}
		glColor3f(0, 1, 0);
		for (int i = 0; i < laplacian_rcm_.size(); ++i) {
			TriMesh::VHandle vh = refined_simplified_mesh_.vertex_handle(i);
			GL::glVertex(refined_simplified_mesh_.point(vh));
			GL::glVertex(refined_simplified_mesh_.point(vh) + laplacian_rcm_[i]); 
		}
		glEnd();
		// 下面画出采样点到evaluation坐标的连线, 也画出到closest坐标的连线, 可视化比较是否离采样点近了.
		glLineWidth(4);
		glBegin(GL_LINES);
		GL::glColor(OpenMesh::Vec3f(0, 1, 0));
		for (TriMesh::VIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++v_it) {
			if (mesh_.status(v_it).deleted()) { //原始网格的点到其在极限曲面上的位置的连线
				if (std::find(vsplit_array_.begin(), vsplit_array_.end(), v_it.handle()) != vsplit_array_.end())
				if (mesh_.property(vp_eval_, v_it)) {
					GL::glVertex(mesh_.property(vp_eval_pos_, v_it.handle())); 
					GL::glVertex(mesh2_.point(v_it));
				}
			}						 
		} 
		GL::glColor(OpenMesh::Vec3f(1, 1, 0));//最近点
		for (TriMesh::VIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++v_it) {
			if (mesh_.status(v_it).deleted()) {
				if (std::find(vsplit_array_.begin(), vsplit_array_.end(), v_it.handle()) != vsplit_array_.end())
				if (mesh_.property(vp_closest_, v_it)) {
					GL::glVertex(mesh_.property(vp_closest_pos_, v_it)); 
					GL::glVertex(mesh2_.point(v_it)); 
				}
			}						 
		} 
		GL::glColor(OpenMesh::Vec3f(1, 0, 0));//project point
		for (TriMesh::VIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++v_it) {
			if (mesh_.status(v_it).deleted()) {
				if (std::find(vsplit_array_.begin(), vsplit_array_.end(), v_it.handle()) != vsplit_array_.end())
				if (mesh_.property(vp_project_, v_it)) {
					GL::glVertex(mesh_.property(vp_project_pos_, v_it)); 
					GL::glVertex(mesh2_.point(v_it)); 
				}
			}						 
		} 
		glEnd();
		glLineWidth(1);
		 
		if (is_view_fitting_error_of_mesh2_) {
			glDisable(GL_LIGHTING);
			glShadeModel(GL_SMOOTH);

			// 被取代的原因详见MeshViewer::draw(_draw_mode)函数
			glEnableClientState(GL_VERTEX_ARRAY);
			glEnableClientState(GL_NORMAL_ARRAY);
			glEnableClientState(GL_COLOR_ARRAY);
			GL::glVertexPointer(mesh2_.points());
			GL::glNormalPointer(mesh2_.vertex_normals());
			GL::glColorPointer(mesh2_.vertex_colors());

			glDrawElements(GL_TRIANGLES, mesh2_indices_for_render_.size(), GL_UNSIGNED_INT, &mesh2_indices_for_render_[0]);

			glDisableClientState(GL_VERTEX_ARRAY);
			glDisableClientState(GL_NORMAL_ARRAY);
			glDisableClientState(GL_COLOR_ARRAY);
		}  
		if (is_view_fitting_error_of_sm_) {
			glDisable(GL_LIGHTING);
			glShadeModel(GL_SMOOTH);

			glBegin(GL_TRIANGLES);
			for (TriMesh::FIter f_it(simplified_mesh_.faces_begin()), f_end(simplified_mesh_.faces_end()); f_it!=f_end; ++f_it) { // 这种方式应该学学啊
				TriMesh::FVIter fv_it = simplified_mesh_.fv_iter(f_it.handle()); 
				GL::glColor(simplified_mesh_.property(vp_color_fitting_error_sm_, fv_it));
				GL::glVertex(simplified_mesh_.point(fv_it));
				++fv_it;
				GL::glColor(simplified_mesh_.property(vp_color_fitting_error_sm_, fv_it));
				GL::glVertex(simplified_mesh_.point(fv_it));
				++fv_it;
				GL::glColor(simplified_mesh_.property(vp_color_fitting_error_sm_, fv_it));
				GL::glVertex(simplified_mesh_.point(fv_it));
			}
			glEnd();
		} 
		MeshModel::draw(); 
	} else { std::cout << "Error: there isn't this render type.\n"; }

	//draw_mesh2_ =  true;//false;
	if (draw_mesh2_ == true) { // 是否显示原始网格, 为了方便比较.
		DGP::hidderline_render(mesh2_, OpenMesh::Vec3f(1, 0, 0));
	}
}
void DecimationModel::pick(float _p0, float _p1, float _p2) {
	OpenMesh::Vec3d p(_p0, _p1, _p2);
	//// find closest edge
	Mesh::EHandle     eh;// the same as Mesh::VertexHandle vh;
	Mesh::Scalar      d, dmin(FLT_MAX);
	for (TriMesh::EIter e_it=mesh_.edges_begin(), e_end=mesh_.edges_end(); e_it!=e_end; ++e_it)
	{
		d = ((mesh_.point(mesh_.to_vertex_handle(mesh_.halfedge_handle(e_it, 0)))
			+ mesh_.point(mesh_.to_vertex_handle(mesh_.halfedge_handle(e_it, 1)))) * 0.5 - p).sqrnorm();//还没有开方的距离
		if (d < dmin) { dmin = d; eh = e_it.handle(); }
	}
	std::cout << "Pick edge " << eh.idx() << "\n";
}
void DecimationModel::update() { //每一次openmesh后都update,此时不一定要do parameterization
	do_parameterization_ = false;  

	refined_simplified_mesh_.clear();
	MeshModel::update();
}

// ============================================================================