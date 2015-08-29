#include "DecimationModel.h" 

//#include "../../../src/CourseExamples/04-Fairing/TaucsSolver.hh"

#include <algorithm> // for max_element algorithm
#include <map>
#include <assert.h>
#include <fstream>
//#include "../../../../taucs/taucsaddon.h" //��һ��ʹ�÷���, ���ɵ�exe�ļ���.
#include "taucs/taucsaddon.h" //��һ��ʹ�÷���, ���ɵ�exe�ļ���.
#include "EquationSolver_lgq/equationsolver.h"

// -------------- ���ﲻֱ��ʹ��MidpointLinearSubdT���ԭ����split_edge_rsm��������Ҫһ���ر�Ĵ���.
void DecimationModel::split_face_rsm(TriMesh::FaceHandle& _fh)
{
	TriMesh::HalfedgeHandle
		heh1(refined_simplified_mesh_.halfedge_handle(_fh)),
		heh2(refined_simplified_mesh_.next_halfedge_handle(refined_simplified_mesh_.next_halfedge_handle(heh1))),
		heh3(refined_simplified_mesh_.next_halfedge_handle(refined_simplified_mesh_.next_halfedge_handle(heh2)));

	TriMesh::FHandle f0, f1, f2, f3;
	// Cutting off every corner of the 6_gon
	f1 = corner_cutting_rsm(heh1 );
	f2 = corner_cutting_rsm(heh2 );
	f3 = corner_cutting_rsm(heh3 );
	f0 = _fh;

	// --------------
	TriMesh::VHandle v1 = refined_simplified_mesh_.from_vertex_handle(heh1);
	TriMesh::VHandle v2 = refined_simplified_mesh_.from_vertex_handle(heh2);
	TriMesh::VHandle v3 = refined_simplified_mesh_.from_vertex_handle(heh3);
	TriMesh::VHandle v4 = refined_simplified_mesh_.to_vertex_handle(heh1);
	TriMesh::VHandle v5 = refined_simplified_mesh_.to_vertex_handle(heh2);
	TriMesh::VHandle v6 = refined_simplified_mesh_.to_vertex_handle(heh3);
	TriMesh::VHandle v1old = voldhandles_[vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v1)).point_handle()];
	TriMesh::VHandle v2old = voldhandles_[vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v2)).point_handle()];
	TriMesh::VHandle v3old = voldhandles_[vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v3)).point_handle()];

	TriMesh::FHandle f0old(simplified_mesh_.property(fp_fh_sm_, f0)); 

	OpenMesh::Vec3d vbc_tmp(0, 0, 0);
	for (std::vector<TriMesh::VHandle>::const_iterator it(mesh_.property(fvset_, f0old).begin()), it_end(mesh_.property(fvset_, f0old).end()); it != it_end; ++it) {
		if (mesh_.status(*it).deleted() == false) std::cout << "Error: split_face_rsm, *it deleted.\n";
		if (mesh_.property(vf_, *it) != f0old) std::cout << "Error: split_face_rsm, vf_ of *it.\n";
		TriMesh::VHandle tmp = mesh_.property(vf0_, *it);
		if (tmp == v1old) vbc_tmp[0] = mesh_.property(vbc_, *it)[0];
		else if (tmp == v2old) vbc_tmp[1] = mesh_.property(vbc_, *it)[0];
		else if (tmp == v3old) vbc_tmp[2] = mesh_.property(vbc_, *it)[0];
		else std::cout << "Error: 0.\n";
		tmp = mesh_.property(vf1_, *it);
		if (tmp == v1old) vbc_tmp[0] = mesh_.property(vbc_, *it)[1];
		else if (tmp == v2old) vbc_tmp[1] = mesh_.property(vbc_, *it)[1];
		else if (tmp == v3old) vbc_tmp[2] = mesh_.property(vbc_, *it)[1];
		else std::cout << "Error: 1.\n";
		tmp = mesh_.property(vf2_, *it);
		if (tmp == v1old) vbc_tmp[0] = mesh_.property(vbc_, *it)[2];
		else if (tmp == v2old) vbc_tmp[1] = mesh_.property(vbc_, *it)[2];
		else if (tmp == v3old) vbc_tmp[2] = mesh_.property(vbc_, *it)[2];
		else std::cout << "Error: 2.\n";
		// ���˻�õ�vbc_tmp�������(v1, v2, v3)��˵��.
		if (vbc_tmp[1] > 0.5) {
			mesh_.property(vf_of_rsm_, *it) = f2;
			mesh_.property(vf0_of_rsm_, *it) = v4; mesh_.property(vf1_of_rsm_, *it) = v2; mesh_.property(vf2_of_rsm_, *it) = v5;
			mesh_.property(vbc_of_rsm_, *it) = OpenMesh::Vec3d(1 - (vbc_tmp[1]-0.5) * 2 - vbc_tmp[2] * 2, (vbc_tmp[1]-0.5) * 2, vbc_tmp[2] * 2);

			refined_simplified_mesh_.property(fvset_rsm_, f2).push_back(*it);
		} else if (vbc_tmp[2] > 0.5) {
			mesh_.property(vf_of_rsm_, *it) = f3;
			mesh_.property(vf0_of_rsm_, *it) = v6; mesh_.property(vf1_of_rsm_, *it) = v5; mesh_.property(vf2_of_rsm_, *it) = v3;
			mesh_.property(vbc_of_rsm_, *it) = OpenMesh::Vec3d(1 - vbc_tmp[1] * 2 - (vbc_tmp[2]-0.5) * 2, vbc_tmp[1] * 2, (vbc_tmp[2] - 0.5) * 2);

			refined_simplified_mesh_.property(fvset_rsm_, f3).push_back(*it); 
		} else if (vbc_tmp[1] + vbc_tmp[2] <= 0.5) {
			mesh_.property(vf_of_rsm_, *it) = f1;
			mesh_.property(vf0_of_rsm_, *it) = v1; mesh_.property(vf1_of_rsm_, *it) = v4; mesh_.property(vf2_of_rsm_, *it) = v6;
			mesh_.property(vbc_of_rsm_, *it) = OpenMesh::Vec3d(1 - vbc_tmp[1] * 2 - vbc_tmp[2] * 2, vbc_tmp[1] * 2, vbc_tmp[2] * 2);

			refined_simplified_mesh_.property(fvset_rsm_, f1).push_back(*it);
		} else {
			mesh_.property(vf_of_rsm_, *it) = f0;
			mesh_.property(vf0_of_rsm_, *it) = v5; mesh_.property(vf1_of_rsm_, *it) = v6; mesh_.property(vf2_of_rsm_, *it) = v4;
			mesh_.property(vbc_of_rsm_, *it) = OpenMesh::Vec3d(1 - (0.5 -vbc_tmp[1]) * 2 - (0.5 - vbc_tmp[2]) * 2, (0.5 -vbc_tmp[1]) * 2, (0.5 - vbc_tmp[2]) * 2);

			refined_simplified_mesh_.property(fvset_rsm_, f0).push_back(*it);    //���������û�д��.
		}
	} //����refined_simplified_mesh_��ÿһ���涼������һЩ�������������ϵ�ԭʼ�㼯.mesh_�ϵı�ɾ���ĵ�Ҳ��refined_simplified_mesh_�����Լ��Ĳ���������.

	if (refined_simplified_mesh_.property(fvset_rsm_, f0).size() + refined_simplified_mesh_.property(fvset_rsm_, f1).size()
		+ refined_simplified_mesh_.property(fvset_rsm_, f2).size() + refined_simplified_mesh_.property(fvset_rsm_, f3).size() 
		!= mesh_.property(fvset_, f0old).size() )
		std::cout << "Error: size.\n";
}
TriMesh::FHandle DecimationModel::corner_cutting_rsm(TriMesh::HalfedgeHandle& _he)
{
	// Define Halfedge Handles
	TriMesh::HalfedgeHandle
		heh1(_he),
		heh5(heh1),
		heh6(refined_simplified_mesh_.next_halfedge_handle(heh1));

	// Cycle around the polygon to find correct Halfedge
	for (; refined_simplified_mesh_.next_halfedge_handle(refined_simplified_mesh_.next_halfedge_handle(heh5)) != heh1;
		heh5 = refined_simplified_mesh_.next_halfedge_handle(heh5))
	{}

	TriMesh::VertexHandle
		vh1 = refined_simplified_mesh_.to_vertex_handle(heh1),
		vh2 = refined_simplified_mesh_.to_vertex_handle(heh5);

	TriMesh::HalfedgeHandle
		heh2(refined_simplified_mesh_.next_halfedge_handle(heh5)),
		heh3(refined_simplified_mesh_.new_edge( vh1, vh2)),
		heh4(refined_simplified_mesh_.opposite_halfedge_handle(heh3));

	// Old and new Face
	TriMesh::FaceHandle     fh_old(refined_simplified_mesh_.face_handle(heh6));
	TriMesh::FaceHandle     fh_new(refined_simplified_mesh_.new_face());


	// Re-Set Handles around old Face
	refined_simplified_mesh_.set_next_halfedge_handle(heh4, heh6);
	refined_simplified_mesh_.set_next_halfedge_handle(heh5, heh4);

	refined_simplified_mesh_.set_face_handle(heh4, fh_old);
	refined_simplified_mesh_.set_face_handle(heh5, fh_old);
	refined_simplified_mesh_.set_face_handle(heh6, fh_old);
	refined_simplified_mesh_.set_halfedge_handle(fh_old, heh4);

	// Re-Set Handles around new Face
	refined_simplified_mesh_.set_next_halfedge_handle(heh1, heh3);
	refined_simplified_mesh_.set_next_halfedge_handle(heh3, heh2);

	refined_simplified_mesh_.set_face_handle(heh1, fh_new);
	refined_simplified_mesh_.set_face_handle(heh2, fh_new);
	refined_simplified_mesh_.set_face_handle(heh3, fh_new);

	refined_simplified_mesh_.set_halfedge_handle(fh_new, heh1);

	return fh_new;
}


void DecimationModel::split_edge_rsm(TriMesh::EdgeHandle& _eh)
{
	TriMesh::HalfedgeHandle
		heh     = refined_simplified_mesh_.halfedge_handle(_eh, 0),
		opp_heh = refined_simplified_mesh_.halfedge_handle(_eh, 1);

	TriMesh::HalfedgeHandle new_heh, opp_new_heh, t_heh;
	TriMesh::VertexHandle   vh;
	TriMesh::VertexHandle   vh1(refined_simplified_mesh_.to_vertex_handle(heh));

	// new vertex, ���ǸĶ��ĵط�, ���ӵ��е���������������ʱ�����������. 
	TriMesh::HalfedgeHandle h0 = simplified_mesh_.halfedge_handle(_eh, 0);  
	TriMesh::VHandle v0 = simplified_mesh_.to_vertex_handle(h0);           //����ߵ�������.
	TriMesh::VHandle v1 = simplified_mesh_.from_vertex_handle(h0); 
	int p0 = vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v0)).point_handle();
	int p1 = vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v1)).point_handle();
	TriMesh::VHandle old_v0 = voldhandles_[p0], old_v1 = voldhandles_[p1]; //��������Ķ�Ӧ�ľɵ�.
	if (mesh_.status(old_v0).deleted()) std::cout << "Error: old v0 deleted.\n";
	if (mesh_.status(old_v1).deleted()) std::cout << "Error: old v1 deleted.\n";
	TriMesh::EHandle old_e = mesh_.edge_handle(mesh_.find_halfedge(old_v0, old_v1)); 
	if (mesh_.property(empl_, old_e) == false ) std::cout << "Error: simplified_mesh_��Ӧ��mesh_�ϻ��б�û���е����, ���ǲ��÷�����.\n";
	vh = refined_simplified_mesh_.new_vertex(mesh_.property(emp_, old_e)); 
	/*std::cout << " " << mesh_.property(emp_, old_e) << " ";*/
	// Re-link mesh entities //Ŀ���Ƕ�λt_heh
	if (refined_simplified_mesh_.is_boundary(_eh))
	{
		for (t_heh = heh;
			refined_simplified_mesh_.next_halfedge_handle(t_heh) != opp_heh;
			t_heh = refined_simplified_mesh_.opposite_halfedge_handle(refined_simplified_mesh_.next_halfedge_handle(t_heh)))
		{}
	}
	else
	{
		for (t_heh = refined_simplified_mesh_.next_halfedge_handle(opp_heh);
			refined_simplified_mesh_.next_halfedge_handle(t_heh) != opp_heh;
			t_heh = refined_simplified_mesh_.next_halfedge_handle(t_heh) )
		{} 
	}

	new_heh     = refined_simplified_mesh_.new_edge(vh, vh1);
	opp_new_heh = refined_simplified_mesh_.opposite_halfedge_handle(new_heh);
	refined_simplified_mesh_.set_vertex_handle( heh, vh );

	refined_simplified_mesh_.set_next_halfedge_handle(t_heh, opp_new_heh);
	refined_simplified_mesh_.set_next_halfedge_handle(new_heh, refined_simplified_mesh_.next_halfedge_handle(heh));
	refined_simplified_mesh_.set_next_halfedge_handle(heh, new_heh);
	refined_simplified_mesh_.set_next_halfedge_handle(opp_new_heh, opp_heh);

	if (refined_simplified_mesh_.face_handle(opp_heh).is_valid())
	{
		refined_simplified_mesh_.set_face_handle(opp_new_heh, refined_simplified_mesh_.face_handle(opp_heh));
		refined_simplified_mesh_.set_halfedge_handle(refined_simplified_mesh_.face_handle(opp_new_heh), opp_new_heh);
	}

	refined_simplified_mesh_.set_face_handle( new_heh, refined_simplified_mesh_.face_handle(heh) );
	refined_simplified_mesh_.set_halfedge_handle( vh, new_heh);
	refined_simplified_mesh_.set_halfedge_handle( refined_simplified_mesh_.face_handle(heh), heh );
	refined_simplified_mesh_.set_halfedge_handle( vh1, opp_new_heh );

	// Never forget this, when playing with the topology
	refined_simplified_mesh_.adjust_outgoing_halfedge( vh );
	refined_simplified_mesh_.adjust_outgoing_halfedge( vh1 );

	// ����Ҳ���Ҽ����, ��Ϊ�¼���Ķ���vh����(new_heh, opp_new_heh)���������ɵ��µı߶���Ҫָ������
	// ����û�д���߽������.
	if (refined_simplified_mesh_.status(_eh).feature()) { // �ղŽ��з��ѵ��������crease edge
		refined_simplified_mesh_.property(vp_type_rsm_, vh) = DGP::CREASE_VFT;
		refined_simplified_mesh_.status(refined_simplified_mesh_.edge_handle(new_heh)).set_feature(true);
	} else {
		refined_simplified_mesh_.property(vp_type_rsm_, vh) = DGP::SMOOTH_VFT;
		refined_simplified_mesh_.status(refined_simplified_mesh_.edge_handle(new_heh)).set_feature(false);
	} 
	Xm_vertex_.push_back(vh);//�����ĵ���뵽����
}

void DecimationModel::midpointsubdivision_1_process() { //�򵥵��е�ϸ��, �ɵ㲻�ö�, �߼����µ�Ϳ�����.
	std::cout << "\n";
	std::cout << "========��ʼ1���е�ϸ�ֹ���=============\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.
	crease_edge_2_mid_pos_.clear();

	if (midpoint_sample_succeed_ == false) {
		std::cout << "Error: �е������û����ȫ�ɹ�, ��ʱ���ܼ�������.\n";

	} else {
		refined_simplified_mesh_ = simplified_mesh_;//���Ƽ򻯵�����, ����ÿ������������Լ���Щ����crease edge����Щ����. �����������ͬ��.
		refined_simplified_mesh_.add_property(vp_type_rsm_);
		mesh_.add_property(vf_of_rsm_);//ԭʼ�����ϵĶ���Ҳ���refined_simplified_mesh_�в�������Ϣ.
		mesh_.add_property(vf0_of_rsm_); mesh_.add_property(vf1_of_rsm_); mesh_.add_property(vf2_of_rsm_); mesh_.add_property(vbc_of_rsm_);
		refined_simplified_mesh_.add_property(fvset_rsm_);

		Xm_vertex_.clear();
		es_vec_.clear();
		int n_corner = 0; // 
		for (TriMesh::VIter v_it(simplified_mesh_.vertices_begin()), v_end(simplified_mesh_.vertices_end()); v_it != v_end; ++v_it) {
			refined_simplified_mesh_.property(vp_type_rsm_, v_it) = simplified_mesh_.property(vp_type_sm_, v_it);

			Xm_vertex_.push_back(v_it.handle());//��ʵ�Ҿ�����������simplified_mesh_.vertex_handle(0��mesh_.n_vertices()-1)���ϸ񰴴��������б���
			
		}
		//assert(Xm_vertex_.size() == simplified_mesh_.n_vertices() - n_corner);

		// ����ģ��OpenMesh/tools/subdivider/uniform/loopt.hh�е���1-to-4����.
		for (TriMesh::EIter e_it = refined_simplified_mesh_.edges_begin(), e_end(refined_simplified_mesh_.edges_end()); 
			e_it != e_end; ++e_it) { //�����Ǹ�e_endһ��Ҫ������, ��������İ��Ҳ�μӷ��ѵ�.
			refined_simplified_mesh_.status(e_it).set_feature(simplified_mesh_.status(e_it).feature()); // �������crease edge����smooth edge

			es_vec_.push_back(e_it.handle());

			split_edge_rsm(e_it.handle());//�����Ķ���Ҳ���뵽Xm_vertex_����ĺ�벿��.

			if (refined_simplified_mesh_.status(e_it).feature()) 
				crease_edge_2_mid_pos_[e_it.handle()] 
				= refined_simplified_mesh_.point(Xm_vertex_[Xm_vertex_.size()-1]);//�ռ�����Ǹ�����.
		}
		int i = 0; //To confirm the vertices list.
		for (TriMesh::VIter v_it(refined_simplified_mesh_.vertices_begin()), v_end(refined_simplified_mesh_.vertices_end()); v_it != v_end; ++v_it) {
			//if (refined_simplified_mesh_.property(vp_type_rsm_, v_it) != DGP::CORNER_VFT) 
			if (v_it.handle() != Xm_vertex_[i++]) std::cout << "Error: vertex not match.\n"; 
			//std::cout << v_it.handle().idx() << ": " << refined_simplified_mesh_.point(v_it).norm() << "\n";
		}
		//assert(Xm_vertex_.size() == refined_simplified_mesh_.n_vertices() - n_corner);

		for (TriMesh::FIter f_it = refined_simplified_mesh_.faces_begin(), f_end(refined_simplified_mesh_.faces_end()); f_it != f_end; ++f_it)
			split_face_rsm(f_it.handle());

		// compute face & vertex normals
		refined_simplified_mesh_.update_normals();
		std::cout << "������õ����� " << refined_simplified_mesh_.n_vertices() << " vertices, " 
			<< refined_simplified_mesh_.n_edges() << " edges, " << refined_simplified_mesh_.n_faces()    << " faces\n";

		// ------------
		int fe = 0;
		for (TriMesh::EIter e_it(refined_simplified_mesh_.edges_begin()), e_end(refined_simplified_mesh_.edges_end()); e_it != e_end; ++e_it) {
			if (refined_simplified_mesh_.status(e_it).feature()) ++fe;
		}
		int sv = 0, dv = 0, crv = 0, cov = 0, bv = 0;
		for (TriMesh::VIter v_it = refined_simplified_mesh_.vertices_begin(), v_end = refined_simplified_mesh_.vertices_end(); v_it != v_end; ++ v_it) {
			if (refined_simplified_mesh_.property(vp_type_rsm_, v_it) == DGP::SMOOTH_VFT) sv++;
			if (refined_simplified_mesh_.property(vp_type_rsm_, v_it) == DGP::DART_VFT) dv++;
			if (refined_simplified_mesh_.property(vp_type_rsm_, v_it) == DGP::CREASE_VFT) crv++;
			if (refined_simplified_mesh_.property(vp_type_rsm_, v_it) == DGP::CORNER_VFT) cov++;
			if (refined_simplified_mesh_.is_boundary(v_it) == true) bv++; 
		}
		std::cout << "������: " << fe << ".\n";
		std::cout << "������: smooth: " << sv << ", dart " << dv << ", crease: " << crv << ", corner: " << cov << ", boundary: " << bv << ".\n";
		//assert(n_corner == cov);
		n_corv_of_sm_ = cov;//����ǵ���Ŀ,�����Ŀ�ǲ����.both simplified_mesh_ and resimplified_simplifed_mesh_ have the same corners

		set_mesh_toberendered(&refined_simplified_mesh_, REFINED_SIMPLIFIED);		
	}
	std::cout << "========����1���е�ϸ�ֹ���=============\n";
}
// ============================================================================
bool linear_equation_solover(const std::vector<std::map<int, taucsType> > & Amn_data, std::vector<OpenMesh::Vec3d> &Xn, 
							 const std::vector<OpenMesh::Vec3d> &Bm_data, const int m, const int n,
							 const std::string &_str)
{
	// Amn Xn = Bm.
	if (Amn_data.size() != n) { std::cout << "Error: linear_equation_solover() ��������.\n"; return false;}
	if (Bm_data.size() != m)  { std::cout << "Error: linear_equation_solover() ��������.\n"; return false;}

	taucs_ccs_matrix * Amn = CreateTaucsMatrixFromColumns(Amn_data, m, TAUCS_DOUBLE);
	taucs_ccs_matrix * Amn_t = MatrixTranspose(Amn); 
	taucs_ccs_matrix *Amn_tAmn = Mul2NonSymmMatSymmResult(Amn_t, Amn); 
	double *Bm = new double[m*3];
	for (int i = 0; i < m; ++i) {
		Bm[i] = Bm_data[i][0];
		Bm[i+m] = Bm_data[i][1];
		Bm[i+m+m] = Bm_data[i][2];
	}
	double *Amn_tBm = new double[n*3]; //notice here is n not m, for Amn_t is n*m.
	MulNonSymmMatrixVector(Amn_t, Bm, Amn_tBm);                // ÿ����һά, 3ά������3��. (n m)*(m 1)=(n 1)
	MulNonSymmMatrixVector(Amn_t, Bm + m, Amn_tBm + n);    
	MulNonSymmMatrixVector(Amn_t, Bm + 2*m, Amn_tBm + 2*n);

	double *Xn_temp = new double[n*3];

	char* solver_options[] = { "taucs.factor.LLT=true", NULL };
	int error_code;
	bool return_code;
	error_code = taucs_linsolve(Amn_tAmn, NULL, 3, Xn_temp, Amn_tBm, solver_options, NULL);//���2�Ǳ�ʾx��A_tb���Ƕ�ά��, �ؼ��������һ������A_tA
	if (error_code != TAUCS_SUCCESS) {
		std::cerr << "Solver failed (" << error_code << ")."  << _str << std::endl;
		return_code = false;
	} else {
		std::cout << "Equations solved successfully." << _str << "\n"; 
		Xn = std::vector<OpenMesh::Vec3d>(n, OpenMesh::Vec3d(0,0,0));
		for (int i = 0; i < n; ++i) {
			Xn[i][0] = Xn_temp[i];	
			Xn[i][1] = Xn_temp[i+n];	
			Xn[i][2] = Xn_temp[i+2*n];	
		}
		return_code = true;
	}
	taucs_ccs_free(Amn); taucs_ccs_free(Amn_t);	taucs_ccs_free(Amn_tAmn);  
	delete [] Bm;
	delete [] Amn_tBm;
	delete [] Xn_temp;
	return return_code;
}
void DecimationModel::calc_crease_first(int m, int n) {
	//����simplified_mesh_�ж��ٸ�crease vertex, ÿһ��crease vertex��Ӧ��һ������.
	std::vector<TriMesh::VertexHandle> crease_vertex_sm;
	crease_vertex_sm.clear();
	std::map<TriMesh::VHandle, int> crease_vertex_idx_sm; //��������crease vertex�������е���.
	for (TriMesh::VertexIter v_it(simplified_mesh_.vertices_begin()), v_end(simplified_mesh_.vertices_end());
		v_it != v_end; ++ v_it) {
			if (simplified_mesh_.property(vp_type_sm_, v_it) == DGP::CREASE_VFT) {
				crease_vertex_sm.push_back(v_it.handle());	 
				crease_vertex_idx_sm[v_it.handle()] = crease_vertex_sm.size() - 1;
			}
	}
	//����simplified_mesh_�ж��ٸ�crease edge, ÿһ����refined_simplified_mesh_����һ��������, �ֶ�Ӧ��һ������
	std::vector<TriMesh::EdgeHandle> crease_edge_sm;
	crease_edge_sm.clear();
	for (TriMesh::EdgeIter e_it(simplified_mesh_.edges_begin()), e_end(simplified_mesh_.edges_end());
		e_it != e_end; ++e_it) {
			if (simplified_mesh_.status(e_it).feature()) {
				TriMesh::VHandle v0 = simplified_mesh_.to_vertex_handle(simplified_mesh_.halfedge_handle(e_it, 0));
				TriMesh::VHandle v1 = simplified_mesh_.to_vertex_handle(simplified_mesh_.halfedge_handle(e_it, 1));
				//��Ҫ���˵㶼��corner��. ��Ҫ���˵㶼��dart��.
				if ((simplified_mesh_.property(vp_type_sm_, v0) == DGP::CORNER_VFT ||  simplified_mesh_.property(vp_type_sm_, v0) == DGP::DART_VFT)
					&& (simplified_mesh_.property(vp_type_sm_, v1) == DGP::CORNER_VFT ||  simplified_mesh_.property(vp_type_sm_, v1) == DGP::DART_VFT)) {

				} else crease_edge_sm.push_back(e_it.handle()); 
			}
	}
	//if (crease_edge_sm.size() != crease_edge_2_mid_pos_.size()) { //�����ǶԵģ�����������ȥ�������˶���corner��crease��
	//	std::cout << "Error: simplified_mesh_�ϵ���������Ӧ��һֱ��.\n";
	//}
	// ���������Ϣ����������.
	if (crease_vertex_sm.size() > 0) {
		int c = crease_vertex_sm.size();   // num of column.
		int r = c + crease_edge_sm.size(); // num of row
		std::vector<std::map<int, taucsType> > Crc_data(c);//Crc[col][row]
		double *Xc;
		double *Qr;
		Xc = new double[c * 3];
		Qr = new double[r * 3];
		// ����Crc�����ǰc��. �����һ��, ����simplified_mesh_��initial_control_mesh_�е�һ��crease vertex.
		for (int i = 0; i < c; ++i) {
			// �� i ��
			TriMesh::VHandle vh = crease_vertex_sm[i];
			// this vh as a limit position at right hand side.
			OpenMesh::Vec3d limit_pos = refined_simplified_mesh_.point(vh) * 48; // 48? check the notice.
			Qr[i] = limit_pos[0]; Qr[i + r] = limit_pos[1]; Qr[i + r+r] = limit_pos[2]; 

			// for the left hand side
			// find two crease edges along this vh.
			std::vector<TriMesh::VHandle> temp; //�ҳ����crease vertex���������ڵ�crease edge
			for (TriMesh::VOHIter it(simplified_mesh_, vh); it; ++it) {
				if (simplified_mesh_.status(simplified_mesh_.edge_handle(it.handle())).feature() == true) {
					temp.push_back(simplified_mesh_.to_vertex_handle(it.handle()));
				} 
			}
			if (temp.size() != 2) std::cout << "Error: crease vertex should has only two crease edges.\n";
			// temp[0,1] can only be dart, corner or crease.
			if (simplified_mesh_.property(vp_type_sm_, temp[0]) == DGP::CREASE_VFT) {
				int idx = crease_vertex_idx_sm[temp[0]];
				Crc_data[idx][i] = 8;
			} else { // dart or corner
				OpenMesh::Vec3d cor = refined_simplified_mesh_.point(temp[0]) * 8;
				Qr[i] -= cor[0]; Qr[i + r] -= cor[1]; Qr[i + r+r] -= cor[2];   
			}
			if (simplified_mesh_.property(vp_type_sm_, temp[1]) == DGP::CREASE_VFT) {
				int idx = crease_vertex_idx_sm[temp[1]];
				Crc_data[idx][i] = 8;
			} else { // dart or corner
				OpenMesh::Vec3d cor = refined_simplified_mesh_.point(temp[1]) * 8;
				Qr[i] -= cor[0]; Qr[i + r] -= cor[1]; Qr[i + r+r] -= cor[2];   
			}
			// �Խ����ϵ�Ԫ��
			Crc_data[i][i] = 32;
		} // end of i < c.
		// ����Crc����ĺ�r-c��, crease_edge_sm.size()��ô���, һ��crease edgeһ������.
		for (int i = 0; i < r-c; ++i) {
			// �� c+i ��
			int j = c+i; // row j
			TriMesh::EHandle eh = crease_edge_sm[i]; //�����Ѿ���֤�˲������˵㶼��corner��.
			// for right hand side.
			OpenMesh::Vec3d limit_pos = crease_edge_2_mid_pos_[eh] * 48;
			Qr[j] = limit_pos[0]; Qr[j + r] = limit_pos[1]; Qr[j + r+r] = limit_pos[2]; 
			// for left hand side.
			TriMesh::HHandle heh = simplified_mesh_.halfedge_handle(eh, 0);
			TriMesh::VHandle vl = simplified_mesh_.to_vertex_handle(heh);//left vertex 
			TriMesh::VHandle vr = simplified_mesh_.from_vertex_handle(heh); //right vertex 
			// vl and vr������Ҫô��crease or (dart/corner), ����������
			// 1. vl and vr are crease.
			if (simplified_mesh_.property(vp_type_sm_, vl) == DGP::CREASE_VFT
				&& simplified_mesh_.property(vp_type_sm_, vr) == DGP::CREASE_VFT) {
					// ϵ������Ӧ����vll, vl, vr, vrr: 1, 23, 23, 1.(sum=48)
					// �ȴ���vl���.
					std::vector<TriMesh::VHandle> temp; //�ҳ����crease vertex���������ڵ�crease edge
					for (TriMesh::VOHIter it(simplified_mesh_, vl); it; ++it) {
						if (simplified_mesh_.status(simplified_mesh_.edge_handle(it.handle())).feature() == true) {
							temp.push_back(simplified_mesh_.to_vertex_handle(it.handle()));
						} 
					}
					if (temp.size() != 2) std::cout << "Error: crease vertex should has only two crease edges.\n";
					TriMesh::VHandle vll;
					if (temp[0] == vr) vll = temp[1];
					else vll = temp[0];
					// �ٴ���vr���.
					temp.clear();
					for (TriMesh::VOHIter it(simplified_mesh_, vr); it; ++it) {
						if (simplified_mesh_.status(simplified_mesh_.edge_handle(it.handle())).feature() == true) {
							temp.push_back(simplified_mesh_.to_vertex_handle(it.handle()));
						} 
					}
					if (temp.size() != 2) std::cout << "Error: crease vertex should has only two crease edges.\n";
					TriMesh::VHandle vrr;
					if (temp[0] == vl) vrr = temp[1];
					else vrr = temp[0];
					// ϵ������ 1, 23, 23, 1
					int vl_idx = crease_vertex_idx_sm[vl];
					int vr_idx = crease_vertex_idx_sm[vr];
					Crc_data[vl_idx][j] = 23;
					Crc_data[vr_idx][j] = 23;

					if (simplified_mesh_.property(vp_type_sm_, vll) == DGP::CREASE_VFT) {
						int vll_idx = crease_vertex_idx_sm[vll];
						Crc_data[vll_idx][j] = 1;
					} else { // dart or corner
						OpenMesh::Vec3d cor = refined_simplified_mesh_.point(vll) * 1.0;
						Qr[j] -= cor[0]; Qr[j + r] -= cor[1]; Qr[j + r+r] -= cor[2]; 
					}
					if (simplified_mesh_.property(vp_type_sm_, vrr) == DGP::CREASE_VFT) {
						int vrr_idx = crease_vertex_idx_sm[vrr];
						Crc_data[vrr_idx][j] = 1;
					} else { // dart or corner
						OpenMesh::Vec3d cor = refined_simplified_mesh_.point(vrr) * 1.0;
						Qr[j] -= cor[0]; Qr[j + r] -= cor[1]; Qr[j + r+r] -= cor[2]; 
					}
			}
			// 2. vl is crease and vr is not crease(maybe dart or corner).
			else if (simplified_mesh_.property(vp_type_sm_, vl) == DGP::CREASE_VFT
				&& simplified_mesh_.property(vp_type_sm_, vr) != DGP::CREASE_VFT) {
					// ϵ������Ӧ���� vll, vl, vr : 1, 22, 25.
					// �ȴ���vl���.
					std::vector<TriMesh::VHandle> temp; //�ҳ����crease vertex���������ڵ�crease edge
					for (TriMesh::VOHIter it(simplified_mesh_, vl); it; ++it) {
						if (simplified_mesh_.status(simplified_mesh_.edge_handle(it.handle())).feature() == true) {
							temp.push_back(simplified_mesh_.to_vertex_handle(it.handle()));
						} 
					}
					if (temp.size() != 2) std::cout << "Error: crease vertex should has only two crease edges.\n";
					TriMesh::VHandle vll;
					if (temp[0] == vr) vll = temp[1];
					else vll = temp[0];
					// ϵ������ vll, vl, vr : 1, 22, 25.
					int vl_idx = crease_vertex_idx_sm[vl]; 
					Crc_data[vl_idx][j] = 22; 
					OpenMesh::Vec3d cor = refined_simplified_mesh_.point(vr) * 25;
					Qr[j] -= cor[0]; Qr[j + r] -= cor[1]; Qr[j + r+r] -= cor[2]; 

					if (simplified_mesh_.property(vp_type_sm_, vll) == DGP::CREASE_VFT) {
						int vll_idx = crease_vertex_idx_sm[vll];
						Crc_data[vll_idx][j] = 1;
					} else { // dart or corner
						OpenMesh::Vec3d cor = refined_simplified_mesh_.point(vll) * 1.0;
						Qr[j] -= cor[0]; Qr[j + r] -= cor[1]; Qr[j + r+r] -= cor[2]; 
					} 
			}
			// 3. vl is not, and vr is crease.
			else if (simplified_mesh_.property(vp_type_sm_, vl) != DGP::CREASE_VFT
				&& simplified_mesh_.property(vp_type_sm_, vr) == DGP::CREASE_VFT) {
					// ϵ������Ӧ����vl, vr, vrr: 25, 22, 1.(sum=48)
					// �ȴ���vr���.
					std::vector<TriMesh::VHandle > temp;
					for (TriMesh::VOHIter it(simplified_mesh_, vr); it; ++it) {
						if (simplified_mesh_.status(simplified_mesh_.edge_handle(it.handle())).feature() == true) {
							temp.push_back(simplified_mesh_.to_vertex_handle(it.handle()));
						} 
					}
					if (temp.size() != 2) std::cout << "Error: crease vertex should has only two crease edges.\n";
					TriMesh::VHandle vrr;
					if (temp[0] == vl) vrr = temp[1];
					else vrr = temp[0];
					// ϵ������ vl, vr, vrr: 25, 22, 1
					int vr_idx = crease_vertex_idx_sm[vr];
					Crc_data[vr_idx][j] = 22;
					OpenMesh::Vec3d cor = refined_simplified_mesh_.point(vl) * 25;
					Qr[j] -= cor[0]; Qr[j + r] -= cor[1]; Qr[j + r+r] -= cor[2]; 

					if (simplified_mesh_.property(vp_type_sm_, vrr) == DGP::CREASE_VFT) {
						int vrr_idx = crease_vertex_idx_sm[vrr];
						Crc_data[vrr_idx][j] = 1;
					} else { // dart or corner
						OpenMesh::Vec3d cor = refined_simplified_mesh_.point(vrr) * 1.0;
						Qr[j] -= cor[0]; Qr[j + r] -= cor[1]; Qr[j + r+r] -= cor[2]; 
					}
			}
			else {
				std::cout << "Error: ��Ӧ���е���������İ�." << simplified_mesh_.property(vp_type_sm_, vl) << ", " 
					<< simplified_mesh_.property(vp_type_sm_, vr) << ".\n";
			}
		} // fill in the Crc_data
		// �ⷽ�� Crc Xc = Qr.
		taucs_ccs_matrix *Crc = CreateTaucsMatrixFromColumns(Crc_data, r, TAUCS_DOUBLE);//PrintTaucsCCSMatrix(Bmn);
		taucs_ccs_matrix *Crc_t = MatrixTranspose(Crc); //PrintTaucsCCSMatrix(Cmn_t);//n*m //
		taucs_ccs_matrix *Crc_tCrc = Mul2NonSymmMatSymmResult(Crc_t, Crc);  //PrintTaucsCCSMatrix(Cmn_tCmn);//n*n

		double *Crc_tQr = new double[c * 3];//��������n�е�, ��֮ǰŪ����m, �˷������timeȥcheck error.
		MulNonSymmMatrixVector(Crc_t, Qr, Crc_tQr);                // ÿ����һά, 3ά������3��. (n m)*(m 1)=(n 1)
		MulNonSymmMatrixVector(Crc_t, Qr + r, Crc_tQr + c);    
		MulNonSymmMatrixVector(Crc_t, Qr + r + r, Crc_tQr + c + c);

		// call TAUCS to solve the system
		char* solver_options[] = { "taucs.factor.LLT=true", NULL };
		int error_code;
		error_code = taucs_linsolve(Crc_tCrc, NULL, 3, Xc, Crc_tQr, solver_options, NULL);//���2�Ǳ�ʾx��A_tb���Ƕ�ά��, �ؼ��������һ������A_tA
		if (error_code != TAUCS_SUCCESS) {
			std::cerr << "Solver failed (" << error_code << ")." << std::endl;
		} else {
			std::cout << "Equations solved successfully. Crc Xc = Qr for crease.\n";
			for (int i = 0; i < c; ++i) { //std::cout << Xn[i] << ", " << Xn[i+n] << ", " <<  Xn[i+n+n] << ".\t";
				initial_control_mesh_.set_point(crease_vertex_sm[i], 
					OpenMesh::Vec3d(Xc[i], Xc[i+c], Xc[i+c+c]));
			}
		} 
		taucs_ccs_free(Crc); taucs_ccs_free(Crc_t); taucs_ccs_free(Crc_tCrc);
		delete[] Xc; delete [] Qr; delete [] Crc_tQr;
	} //2009-04-30������2��, �������crease vertex��Ч��Ҫ��!!!
}
void DecimationModel::calc_AmmBmn(std::vector<std::map<int, taucsType> > &Amm_data, std::vector<std::map<int, taucsType> > &Bmn_data) {
	const int m = Xm_vertex_.size(), n = simplified_mesh_.n_vertices();
	if (m != refined_simplified_mesh_.n_vertices()) std::cout << "Error: ��Ӧ��һ�µ�, ������������.\n";
	for (int i = 0; i < n; ++i) { //һ��һ�е���д, ����дBmn_data��ǰn��, Ҳ���Ǹ�����n���ɵ������
		//TriMesh::VertexHandle vh(simplified_mesh_.vertex_handle(i)); //���濪ʼ��дBmn_data�еĵ�i�е�, ��Ӧ��Xm�еĵ�i������(vh)
		//if (vh != Xm_vertex_[i]) std::cout << "Error: Xm_vertex_[i] should be equal to simplified_mesh_.vertex_handle(i), when i < n.\n";
		//std::cout << "i " << i << ", " << simplified_mesh_.property(vp_type_sm_, vh) << ". ";
		TriMesh::VertexHandle vh(Xm_vertex_[i]);
		if (simplified_mesh_.property(vp_type_sm_, vh) == DGP::CORNER_VFT) {
			Bmn_data[vh.idx()][i] = 1;
		} else 
			if (simplified_mesh_.property(vp_type_sm_, vh) == DGP::CREASE_VFT) {
				std::vector<TriMesh::VertexHandle> vec; //����crease edge����Ӧ�Ķ���.
				for (TriMesh::VOHIter voh_it(simplified_mesh_, vh); voh_it; ++voh_it) {
					if (simplified_mesh_.status(simplified_mesh_.edge_handle(voh_it.handle())).feature()) 
						vec.push_back(simplified_mesh_.to_vertex_handle(voh_it.handle()));				
				}
				if (vec.size() != 2) std::cout << "Error: ����Ӧ��ֻ������crease edge�Ŷ�.\n" << vec.size() << "; " << vec[0] << ", " << vec[1] << ".\n";
				//std::cout << vec.size() << ", " << vec[0] << ", " << vec[1] << ".\n";
				Bmn_data[vh.idx()][i] = 0.75; Bmn_data[vec[0].idx()][i] = 0.125; Bmn_data[vec[1].idx()][i] = 0.125; //
			} else if (simplified_mesh_.is_boundary(vh)) {                   //
				TriMesh::HalfedgeHandle heh = simplified_mesh_.halfedge_handle(vh); //de
				if (simplified_mesh_.is_boundary(heh) == false) { std::cout << "Error: ���Ҫ���Ǳ߽�Ŷ���.\n"; }
				TriMesh::VertexHandle v0 = simplified_mesh_.to_vertex_handle(heh);
				TriMesh::VertexHandle v1 = simplified_mesh_.from_vertex_handle(simplified_mesh_.prev_halfedge_handle(heh));
				Bmn_data[vh.idx()][i] = 0.75; Bmn_data[v0.idx()][i] = 0.125; Bmn_data[v1.idx()][i] = 0.125; //��crease vertex�ĸ���ģʽ��һ����
			} else if (simplified_mesh_.property(vp_type_sm_, vh) == DGP::SMOOTH_VFT || simplified_mesh_.property(vp_type_sm_, vh) == DGP::DART_VFT) {
				int k = simplified_mesh_.valence(vh);
				double beta = 0, alpha = 0;
				if (k == 3) { beta = 0.1875; alpha = 0.5625; } // 3.0/16.0 = 0.1875,   9.0/16.0 = 0.5625
				else if (k == 6) { beta = 0.0625; alpha = 0.375; } // 1.0/16.0 = 0.0625,    3.0/8.0 = 0.375
				else { 
					double inv_v = 1.0/double(k); double t = (3.0 + 2.0 * cos( 2.0 * M_PI * inv_v) );
					alpha  = (40.0 - t * t)/64.0; beta = inv_v * alpha;
				}
				Bmn_data[vh.idx()][i] = 1.0 - alpha; // ��Ӧ��Xm_vertex_[i]���i������, �ǵ�i��, �����ǵڼ��о�����Χ���ھӶ������(�����Ƕ�������)
				for (TriMesh::VertexOHalfedgeIter voh_it(simplified_mesh_, vh); voh_it; ++voh_it) {
					Bmn_data[simplified_mesh_.to_vertex_handle(voh_it).idx()][i] = beta;
				}
			} else std::cout << "Error: û���������͵Ķ�������Ŷ��1.\n";		
	} 
	for (int i = n; i < m; ++i) { // ��дBmn_data�ĵ�n��m-1��, ��es_vec_.size()��ô����
		TriMesh::EdgeHandle eh = es_vec_[i - n];
		if (simplified_mesh_.status(eh).feature() || simplified_mesh_.is_boundary(eh)) {
			TriMesh::VertexHandle v0 = simplified_mesh_.to_vertex_handle(simplified_mesh_.halfedge_handle(eh, 0));
			TriMesh::VertexHandle v1 = simplified_mesh_.to_vertex_handle(simplified_mesh_.halfedge_handle(eh, 1));
			Bmn_data[v0.idx()][i] = 0.5; Bmn_data[v1.idx()][i] = 0.5; 

		} else { // smooth edge, new vertex = 3/8 * v0 + 3/8 * v1 + 1/8 * v2 + 1/8 * v3
			TriMesh::VertexHandle v0 = simplified_mesh_.to_vertex_handle(simplified_mesh_.halfedge_handle(eh, 0));
			TriMesh::VertexHandle v1 = simplified_mesh_.to_vertex_handle(simplified_mesh_.halfedge_handle(eh, 1));
			TriMesh::VertexHandle v2 = simplified_mesh_.to_vertex_handle(simplified_mesh_.next_halfedge_handle(simplified_mesh_.halfedge_handle(eh, 0)));
			TriMesh::VertexHandle v3 = simplified_mesh_.to_vertex_handle(simplified_mesh_.next_halfedge_handle(simplified_mesh_.halfedge_handle(eh, 1)));//std::cout << v0 << ", " << v1 << v2 << ", " << v3 << ";";
			Bmn_data[v0.idx()][i] = 0.375; Bmn_data[v1.idx()][i] = 0.375; Bmn_data[v2.idx()][i] = 0.125; Bmn_data[v3.idx()][i] = 0.125;  
		}
	}  
	for (int i = 0; i < m; ++i) { //һ��һ�е���д, ���m��ָm��, ��i�ж�Ӧ��Xm�еĵ�i������
		TriMesh::VertexHandle vh = Xm_vertex_[i]; //����refined simplified mesh�����ӹ�ϵ�Լ��������, ���refined simplified mesh���е��limit position.
		//if (vh != refined_simplified_mesh_.vertex_handle(i)) std::cout << "Error: Xm�����m������Ӧ�ú�refined simplified mesh�еĶ����Ӧ��.\n";

		if (refined_simplified_mesh_.property(vp_type_rsm_, vh) == DGP::CORNER_VFT) {
			Amm_data[vh.idx()][i] = 1;
		} else 
			if (refined_simplified_mesh_.property(vp_type_rsm_, vh) == DGP::CREASE_VFT) {
				std::vector<TriMesh::VertexHandle> vec; //����crease edge����Ӧ�Ķ���.
				for (TriMesh::VOHIter voh_it(refined_simplified_mesh_, vh); voh_it; ++voh_it) {
					if (refined_simplified_mesh_.status(refined_simplified_mesh_.edge_handle(voh_it.handle())).feature() == true) 
						vec.push_back(refined_simplified_mesh_.to_vertex_handle(voh_it.handle()));				
				}
				if (vec.size() != 2) std::cout << "Error: ����Ӧ��ֻ������crease edge�Ŷ�(2).\n";
				Amm_data[vh.idx()][i] = 0.666667; Amm_data[vec[0].idx()][i] = 0.166667; Amm_data[vec[1].idx()][i] = 0.166667; //
			} else if (refined_simplified_mesh_.is_boundary(vh)) {
				TriMesh::HalfedgeHandle heh = refined_simplified_mesh_.halfedge_handle(vh); //de
				if (refined_simplified_mesh_.is_boundary(heh) == false) { std::cout << "Error: ���Ҫ���Ǳ߽�Ŷ���2.\n"; }
				TriMesh::VertexHandle v0 = refined_simplified_mesh_.to_vertex_handle(heh);
				TriMesh::VertexHandle v1 = refined_simplified_mesh_.from_vertex_handle(refined_simplified_mesh_.prev_halfedge_handle(heh));
				Amm_data[vh.idx()][i] = 0.666667; Amm_data[v0.idx()][i] = 0.166667; Amm_data[v1.idx()][i] = 0.166667; //

			} else if (refined_simplified_mesh_.property(vp_type_rsm_, vh) == DGP::SMOOTH_VFT || refined_simplified_mesh_.property(vp_type_rsm_, vh) == DGP::DART_VFT) {
				int k = refined_simplified_mesh_.valence(vh);
				double inv_v = 1.0/double(k); double t = (3.0 + 2.0 * cos( 2.0 * M_PI * inv_v) ); // double beta = inv_v * (40.0 - t * t)/64.0;
				double k_alpha = 1.0 / ( 24.0/(40.0-t*t) + 1); double alpha = inv_v * k_alpha;
				for (TriMesh::VertexVertexIter vv_it(refined_simplified_mesh_, vh); vv_it; ++vv_it) {
					Amm_data[vv_it.handle().idx()][i] = alpha; 
				}
				Amm_data[vh.idx()][i] = 1.0 - k_alpha; 
			} else std::cout << "Error: û���������͵Ķ�������Ŷ��2.\n";	
	}
}
// ����ⷽ������innitial control mesh����������ѡ��һ.
void DecimationModel::create_initialcontrolmesh_process() {
	std::cout << "\n"; //������ԭʼ����ֱ�����Cmn*Xn=Qm����, ��û������ͬ����ķ������.
	std::cout << "========��ʼ��������������============\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.

	if (Xm_vertex_.size() != (simplified_mesh_.n_vertices() + es_vec_.size())) std::cout << "Error: �����ϲ�ƥ��1.\n";
	if (Xm_vertex_.size() != refined_simplified_mesh_.n_vertices()) std::cout << "Error: �����ϲ�ƥ��2.\n";
	//for (int i = 0, end = refined_simplified_mesh_.n_vertices(); i < end; ++ i) 
	// if (refined_simplified_mesh_.property(vp_type_rsm_, refined_simplified_mesh_.vertex_handle(i)) != DGP::CORNER_VFT) 
	//   assert(Xm_vertex_[i] == refined_simplified_mesh_.vertex_handle());

	initial_control_mesh_ = simplified_mesh_;
	initial_control_mesh_.add_property(vp_type_icm_); 
	initial_control_mesh_.add_property(vp_closestp_icm_);      
	initial_control_mesh_.add_property(vp_distances_);

	for (TriMesh::VertexIter v_it = simplified_mesh_.vertices_begin(), v_end = simplified_mesh_.vertices_end(); v_it != v_end; ++v_it) { //���Ƶ�����
		initial_control_mesh_.property(vp_type_icm_, v_it) = simplified_mesh_.property(vp_type_sm_, v_it);
	}
	for (TriMesh::EdgeIter e_it = simplified_mesh_.edges_begin(), e_end = simplified_mesh_.edges_end(); e_it != e_end; ++e_it) {//���Ʊ�����
		initial_control_mesh_.status(e_it).set_feature(simplified_mesh_.status(e_it.handle()).feature());
	} //std::cout << "1. ";
	// ���Bmn_data, 
	int m = Xm_vertex_.size(), n = simplified_mesh_.n_vertices();//m���� n����, mҲ�Ƿ������з��̵ĸ���.
	std::vector<std::map<int, taucsType> > Bmn_data(n);//n����. �����Bmn_data[�к�][�к�]
	//std::cout << "2, " << m << ", " << n << ", " << n_corv_of_sm_ << ", " << simplified_mesh_.n_edges() << "=" << es_vec_.size() << ".\n";
	 
	// ����Amm_data
	std::vector<std::map<int, taucsType> > Amm_data(m);//���m��ָm��
	calc_AmmBmn(Amm_data, Bmn_data);

	taucs_ccs_matrix *Bmn = CreateTaucsMatrixFromColumns(Bmn_data, m, TAUCS_DOUBLE);//PrintTaucsCCSMatrix(Bmn);
	taucs_ccs_matrix *Amm = CreateTaucsMatrixFromColumns(Amm_data, m, TAUCS_DOUBLE);//PrintTaucsCCSMatrix(Amm);	
	//taucs_ccs_matrix *Bmn_t = MatrixTranspose(Bmn); //PrintTaucsCCSMatrix(Bmn_t);// n*m
	//taucs_ccs_matrix *Amm_t = MatrixTranspose(Amm); //PrintTaucsCCSMatrix(Amm_t);// m*m
	taucs_ccs_matrix *Cmn = Mul2NonSymmetricMatrices(Amm, Bmn); //PrintTaucsCCSMatrix(Cmn); 
	taucs_ccs_matrix *Cmn_t = MatrixTranspose(Cmn); //PrintTaucsCCSMatrix(Cmn_t);//n*m //
	//taucs_ccs_matrix *Cmn_t2 = Mul2NonSymmetricMatrices(Bmn_t, Amm_t);//
	//EqualTest(Cmn_t, Cmn_t2);//ͨ����.
	taucs_ccs_matrix *Cmn_tCmn = Mul2NonSymmMatSymmResult(Cmn_t, Cmn);  //PrintTaucsCCSMatrix(Cmn_tCmn);//n*n
	//taucs_ccs_matrix *Cmn_tCmn2 = Mul2NonSymmMatSymmResult(Mul2NonSymmetricMatrices(Bmn_t, Amm_t), Mul2NonSymmetricMatrices(Amm, Bmn));
	//EqualTest(Cmn_tCmn, Cmn_tCmn2);//ͨ����.
	//taucs_ccs_matrix *Amm_tAmm = Mul2NonSymmMatSymmResult(Amm_t, Amm); //m*m
	//taucs_ccs_matrix *Bmn_tBmn = Mul2NonSymmMatSymmResult(Bmn_t, Bmn); //n*n
 
	//����right hand side, Qm. ����ȥ�ƽ�refined simplified mesh.
	double *Qm = new double[m * 3];
	for (int i = 0; i < m; ++i) {
		TriMesh::Point p = refined_simplified_mesh_.point(Xm_vertex_[i]);//std::cout << p << "\t";
		Qm[i] = p[0]; Qm[i+m] = p[1]; Qm[i+m+m] = p[2];
	}
	double *Cmn_tQm = new double[n * 3];//��������n�е�, ��֮ǰŪ����m, �˷������timeȥcheck error.
	MulNonSymmMatrixVector(Cmn_t, Qm, Cmn_tQm);                // ÿ����һά, 3ά������3��. (n m)*(m 1)=(n 1)
	MulNonSymmMatrixVector(Cmn_t, Qm + m, Cmn_tQm + n);    
	MulNonSymmMatrixVector(Cmn_t, Qm + m + m, Cmn_tQm + n + n);

	//��Ҫ����δ֪Xn
	double *Xn = new double[n * 3];

	// call TAUCS to solve the system
	char* solver_options[] = { "taucs.factor.LLT=true", NULL };
	int error_code;
	error_code = taucs_linsolve(Cmn_tCmn, NULL, 3, Xn, Cmn_tQm, solver_options, NULL);//���2�Ǳ�ʾx��A_tb���Ƕ�ά��, �ؼ��������һ������A_tA
	if (error_code != TAUCS_SUCCESS) {
		std::cerr << "Solver failed (" << error_code << ")." << std::endl;
	} else {
		std::cout << "Equations solved successfully.\n";
		for (int i = 0; i < n; ++i) { //std::cout << Xn[i] << ", " << Xn[i+n] << ", " <<  Xn[i+n+n] << ".\t";
			initial_control_mesh_.set_point(Xm_vertex_[i], OpenMesh::Vec3d(Xn[i], Xn[i+n], Xn[i+n+n]));
		}
	}/**/

	/*// �����������, ����Xm, ����Xn
	double *Amm_tQm = new double[m * 3];//right hand side
	MulNonSymmMatrixVector(Amm_t, Qm, Amm_tQm);                // ÿ����һά, 3ά������3��.
	MulNonSymmMatrixVector(Amm_t, Qm + m, Amm_tQm + m);
	MulNonSymmMatrixVector(Amm_t, Qm + m + m, Amm_tQm + m + m);
	double *Xm = new double[m * 3];//result
	error_code = taucs_linsolve(Amm_tAmm, NULL, 3, Xm, Amm_tQm, solver_options, NULL);//���2�Ǳ�ʾx��A_tb���Ƕ�ά��, �ؼ��������һ������A_tA
	if (error_code != TAUCS_SUCCESS) {
	std::cerr << "Solver failed (" << error_code << ")." << std::endl;
	} else {
	std::cout << "Equations solved.\n";
	for (int i = 0; i < m; ++i) { //std::cout << Xm[i] << ", " << Xm[i+m] << ", " <<  Xm[i+m+m] << ".\t";
	refined_simplified_mesh_.set_point(Xm_vertex_[i], OpenMesh::Vec3d(Xm[i], Xm[i+m], Xm[i+m+m]));
	}
	}
	double *Bmn_tXm = new double[n * 3];//right hand side
	MulNonSymmMatrixVector(Bmn_t, Xm, Bmn_tXm);                // ÿ����һά, 3ά������3��.
	MulNonSymmMatrixVector(Bmn_t, Xm + m, Bmn_tXm + n);
	MulNonSymmMatrixVector(Bmn_t, Xm + m + m, Bmn_tXm + n + n);
	Xn = new double[n * 3];//result
	error_code = taucs_linsolve(Bmn_tBmn, NULL, 3, Xn, Bmn_tXm, solver_options, NULL);//���2�Ǳ�ʾx��A_tb���Ƕ�ά��, �ؼ��������һ������A_tA
	if (error_code != TAUCS_SUCCESS) {
	std::cerr << "Solver failed (" << error_code << ")." << std::endl;
	} else {
	std::cout << "Equations solved successfully.\n";
	for (int i = 0; i < n; ++i) { //std::cout << Xn[i] << ", " << Xn[i+n] << ", " <<  Xn[i+n+n] << ".\t";
	initial_control_mesh_.set_point(Xm_vertex_[i], OpenMesh::Vec3d(Xn[i], Xn[i+n], Xn[i+n+n]));
	}
	}*/


	// free up the memory
	taucs_ccs_free(Bmn); taucs_ccs_free(Amm);
	taucs_ccs_free(Cmn); taucs_ccs_free(Cmn_t); taucs_ccs_free(Cmn_tCmn);
	delete [] Qm;
	delete [] Cmn_tQm;
	delete [] Xn;

	//std::cout << initial_control_mesh_.n_vertices() << ", \n";//" << 
	set_mesh_toberendered(&initial_control_mesh_, INITIAL_CONTROL);
	std::cout << "========������������������============\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.
}
void DecimationModel::create_initialcontrolmesh_process2() {
	// ����汾�����̷ֽ⣬���Ƿֽ�ò�����.
	std::cout << "\n";
	std::cout << "========��ʼ��������������v2==========\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=. 

	if (Xm_vertex_.size() != (simplified_mesh_.n_vertices() + es_vec_.size())) {
		std::cout << "Error: �����ϲ�ƥ��1." << Xm_vertex_.size() << " != " << simplified_mesh_.n_vertices()
			<< " + " << es_vec_.size() << ". " << refined_simplified_mesh_.n_vertices() << ".\n";
	}
	if (Xm_vertex_.size() != refined_simplified_mesh_.n_vertices()) std::cout << "Error: �����ϲ�ƥ��2.\n";
	std::vector<TriMesh::VHandle> smooth_dart_vertex_rsm; smooth_dart_vertex_rsm.clear();
	int idx_Xm =0;
	for (TriMesh::VIter v_it(refined_simplified_mesh_.vertices_begin()), v_end(refined_simplified_mesh_.vertices_end());
		v_it != v_end; ++v_it) {
			if (refined_simplified_mesh_.property(vp_type_rsm_, v_it) == DGP::DART_VFT
				|| refined_simplified_mesh_.property(vp_type_rsm_, v_it) == DGP::SMOOTH_VFT) {
					smooth_dart_vertex_rsm.push_back(v_it.handle()); 
			}
			if (v_it.handle() != Xm_vertex_[idx_Xm++]) std::cout << "Error: this should be the same.\n";
	}
	const int m = Xm_vertex_.size(), n = simplified_mesh_.n_vertices();//m���� n����, mҲ�Ƿ������з��̵ĸ���.
  
	// To figure out the geometry of the initial_control_mesh_,
	// be ware that, the topology of the icm is the same as simplified_mesh_. 
	initial_control_mesh_ = simplified_mesh_;
	initial_control_mesh_.add_property(vp_type_icm_); 
	initial_control_mesh_.add_property(vp_closestp_icm_);
	initial_control_mesh_.add_property(vp_distances_);

	for (TriMesh::VertexIter v_it = simplified_mesh_.vertices_begin(), v_end = simplified_mesh_.vertices_end(); v_it != v_end; ++v_it) { //���Ƶ�����
		initial_control_mesh_.property(vp_type_icm_, v_it) = simplified_mesh_.property(vp_type_sm_, v_it);
	}
	for (TriMesh::EdgeIter e_it = simplified_mesh_.edges_begin(), e_end = simplified_mesh_.edges_end(); e_it != e_end; ++e_it) {//���Ʊ�����
		initial_control_mesh_.status(e_it).set_feature(simplified_mesh_.status(e_it.handle()).feature());
	} //std::cout << "1. ";

	//����simplified_mesh_�ж��ٸ�crease vertex, ÿһ��crease vertex��Ӧ��һ������.
	std::vector<TriMesh::VertexHandle> crease_vertex_sm;
	crease_vertex_sm.clear();
	std::map<TriMesh::VHandle, int> crease_vertex_idx_sm; //��������crease vertex�������е���.
	for (TriMesh::VertexIter v_it(simplified_mesh_.vertices_begin()), v_end(simplified_mesh_.vertices_end());
		v_it != v_end; ++ v_it) {
			if (simplified_mesh_.property(vp_type_sm_, v_it) == DGP::CREASE_VFT) {
				crease_vertex_sm.push_back(v_it.handle());	 
				crease_vertex_idx_sm[v_it.handle()] = crease_vertex_sm.size() - 1;
			}
	}
	//����simplified_mesh_�ж��ٸ�crease edge, ÿһ����refined_simplified_mesh_����һ��������, �ֶ�Ӧ��һ������
	std::vector<TriMesh::EdgeHandle> crease_edge_sm;
	crease_edge_sm.clear();
	for (TriMesh::EdgeIter e_it(simplified_mesh_.edges_begin()), e_end(simplified_mesh_.edges_end());
		e_it != e_end; ++e_it) {
			if (simplified_mesh_.status(e_it).feature()) {
				TriMesh::VHandle v0 = simplified_mesh_.to_vertex_handle(simplified_mesh_.halfedge_handle(e_it, 0));
				TriMesh::VHandle v1 = simplified_mesh_.to_vertex_handle(simplified_mesh_.halfedge_handle(e_it, 1));
				//��Ҫ���˵㶼��corner��. ��Ҫ���˵㶼��dart��.
				if ((simplified_mesh_.property(vp_type_sm_, v0) == DGP::CORNER_VFT ||  simplified_mesh_.property(vp_type_sm_, v0) == DGP::DART_VFT)
					&& (simplified_mesh_.property(vp_type_sm_, v1) == DGP::CORNER_VFT ||  simplified_mesh_.property(vp_type_sm_, v1) == DGP::DART_VFT)) {

				} else crease_edge_sm.push_back(e_it.handle()); 
			}
	}
	//if (crease_edge_sm.size() != crease_edge_2_mid_pos_.size()) { //�����ǶԵģ�����������ȥ�������˶���corner��crease��
	//	std::cout << "Error: simplified_mesh_�ϵ���������Ӧ��һֱ��.\n";
	//}
	// ���������Ϣ����������.
	if (crease_vertex_sm.size() > 0) {
		int c = crease_vertex_sm.size();   // num of column.
		int r = c + crease_edge_sm.size(); // num of row
		std::vector<std::map<int, taucsType> > Crc_data(c);//Crc[col][row]
		double *Xc;
		double *Qr;
		Xc = new double[c * 3];
		Qr = new double[r * 3];
		// ����Crc�����ǰc��. �����һ��, ����simplified_mesh_��initial_control_mesh_�е�һ��crease vertex.
		for (int i = 0; i < c; ++i) {
			// �� i ��
			TriMesh::VHandle vh = crease_vertex_sm[i];
			// this vh as a limit position at right hand side.
			OpenMesh::Vec3d limit_pos = refined_simplified_mesh_.point(vh) * 48; // 48? check the notice.
			Qr[i] = limit_pos[0]; Qr[i + r] = limit_pos[1]; Qr[i + r+r] = limit_pos[2]; 

			// for the left hand side
			// find two crease edges along this vh.
			std::vector<TriMesh::VHandle> temp; //�ҳ����crease vertex���������ڵ�crease edge
			for (TriMesh::VOHIter it(simplified_mesh_, vh); it; ++it) {
				if (simplified_mesh_.status(simplified_mesh_.edge_handle(it.handle())).feature() == true) {
					temp.push_back(simplified_mesh_.to_vertex_handle(it.handle()));
				} 
			}
			if (temp.size() != 2) std::cout << "Error: crease vertex should has only two crease edges.\n";
			// temp[0,1] can only be dart, corner or crease.
			if (simplified_mesh_.property(vp_type_sm_, temp[0]) == DGP::CREASE_VFT) {
				int idx = crease_vertex_idx_sm[temp[0]];
				Crc_data[idx][i] = 8;
			} else { // dart or corner
				OpenMesh::Vec3d cor = refined_simplified_mesh_.point(temp[0]) * 8;
				Qr[i] -= cor[0]; Qr[i + r] -= cor[1]; Qr[i + r+r] -= cor[2];   
			}
			if (simplified_mesh_.property(vp_type_sm_, temp[1]) == DGP::CREASE_VFT) {
				int idx = crease_vertex_idx_sm[temp[1]];
				Crc_data[idx][i] = 8;
			} else { // dart or corner
				OpenMesh::Vec3d cor = refined_simplified_mesh_.point(temp[1]) * 8;
				Qr[i] -= cor[0]; Qr[i + r] -= cor[1]; Qr[i + r+r] -= cor[2];   
			}
			// �Խ����ϵ�Ԫ��
			Crc_data[i][i] = 32;
		} // end of i < c.
		// ����Crc����ĺ�r-c��, crease_edge_sm.size()��ô���, һ��crease edgeһ������.
		for (int i = 0; i < r-c; ++i) {
			// �� c+i ��
			int j = c+i; // row j
			TriMesh::EHandle eh = crease_edge_sm[i]; //�����Ѿ���֤�˲������˵㶼��corner��.
			// for right hand side.
			OpenMesh::Vec3d limit_pos = crease_edge_2_mid_pos_[eh] * 48;
			Qr[j] = limit_pos[0]; Qr[j + r] = limit_pos[1]; Qr[j + r+r] = limit_pos[2]; 
			// for left hand side.
			TriMesh::HHandle heh = simplified_mesh_.halfedge_handle(eh, 0);
			TriMesh::VHandle vl = simplified_mesh_.to_vertex_handle(heh);//left vertex 
			TriMesh::VHandle vr = simplified_mesh_.from_vertex_handle(heh); //right vertex 
			// vl and vr������Ҫô��crease or (dart/corner), ����������
			// 1. vl and vr are crease.
			if (simplified_mesh_.property(vp_type_sm_, vl) == DGP::CREASE_VFT
				&& simplified_mesh_.property(vp_type_sm_, vr) == DGP::CREASE_VFT) {
					// ϵ������Ӧ����vll, vl, vr, vrr: 1, 23, 23, 1.(sum=48)
					// �ȴ���vl���.
					std::vector<TriMesh::VHandle> temp; //�ҳ����crease vertex���������ڵ�crease edge
					for (TriMesh::VOHIter it(simplified_mesh_, vl); it; ++it) {
						if (simplified_mesh_.status(simplified_mesh_.edge_handle(it.handle())).feature() == true) {
							temp.push_back(simplified_mesh_.to_vertex_handle(it.handle()));
						} 
					}
					if (temp.size() != 2) std::cout << "Error: crease vertex should has only two crease edges.\n";
					TriMesh::VHandle vll;
					if (temp[0] == vr) vll = temp[1];
					else vll = temp[0];
					// �ٴ���vr���.
					temp.clear();
					for (TriMesh::VOHIter it(simplified_mesh_, vr); it; ++it) {
						if (simplified_mesh_.status(simplified_mesh_.edge_handle(it.handle())).feature() == true) {
							temp.push_back(simplified_mesh_.to_vertex_handle(it.handle()));
						} 
					}
					if (temp.size() != 2) std::cout << "Error: crease vertex should has only two crease edges.\n";
					TriMesh::VHandle vrr;
					if (temp[0] == vl) vrr = temp[1];
					else vrr = temp[0];
					// ϵ������ 1, 23, 23, 1
					int vl_idx = crease_vertex_idx_sm[vl];
					int vr_idx = crease_vertex_idx_sm[vr];
					Crc_data[vl_idx][j] = 23;
					Crc_data[vr_idx][j] = 23;

					if (simplified_mesh_.property(vp_type_sm_, vll) == DGP::CREASE_VFT) {
						int vll_idx = crease_vertex_idx_sm[vll];
						Crc_data[vll_idx][j] = 1;
					} else { // dart or corner
						OpenMesh::Vec3d cor = refined_simplified_mesh_.point(vll) * 1.0;
						Qr[j] -= cor[0]; Qr[j + r] -= cor[1]; Qr[j + r+r] -= cor[2]; 
					}
					if (simplified_mesh_.property(vp_type_sm_, vrr) == DGP::CREASE_VFT) {
						int vrr_idx = crease_vertex_idx_sm[vrr];
						Crc_data[vrr_idx][j] = 1;
					} else { // dart or corner
						OpenMesh::Vec3d cor = refined_simplified_mesh_.point(vrr) * 1.0;
						Qr[j] -= cor[0]; Qr[j + r] -= cor[1]; Qr[j + r+r] -= cor[2]; 
					}
			}
			// 2. vl is crease and vr is not crease(maybe dart or corner).
			else if (simplified_mesh_.property(vp_type_sm_, vl) == DGP::CREASE_VFT
				&& simplified_mesh_.property(vp_type_sm_, vr) != DGP::CREASE_VFT) {
					// ϵ������Ӧ���� vll, vl, vr : 1, 22, 25.
					// �ȴ���vl���.
					std::vector<TriMesh::VHandle> temp; //�ҳ����crease vertex���������ڵ�crease edge
					for (TriMesh::VOHIter it(simplified_mesh_, vl); it; ++it) {
						if (simplified_mesh_.status(simplified_mesh_.edge_handle(it.handle())).feature() == true) {
							temp.push_back(simplified_mesh_.to_vertex_handle(it.handle()));
						} 
					}
					if (temp.size() != 2) std::cout << "Error: crease vertex should has only two crease edges.\n";
					TriMesh::VHandle vll;
					if (temp[0] == vr) vll = temp[1];
					else vll = temp[0];
					// ϵ������ vll, vl, vr : 1, 22, 25.
					int vl_idx = crease_vertex_idx_sm[vl]; 
					Crc_data[vl_idx][j] = 22; 
					OpenMesh::Vec3d cor = refined_simplified_mesh_.point(vr) * 25;
					Qr[j] -= cor[0]; Qr[j + r] -= cor[1]; Qr[j + r+r] -= cor[2]; 

					if (simplified_mesh_.property(vp_type_sm_, vll) == DGP::CREASE_VFT) {
						int vll_idx = crease_vertex_idx_sm[vll];
						Crc_data[vll_idx][j] = 1;
					} else { // dart or corner
						OpenMesh::Vec3d cor = refined_simplified_mesh_.point(vll) * 1.0;
						Qr[j] -= cor[0]; Qr[j + r] -= cor[1]; Qr[j + r+r] -= cor[2]; 
					} 
			}
			// 3. vl is not, and vr is crease.
			else if (simplified_mesh_.property(vp_type_sm_, vl) != DGP::CREASE_VFT
				&& simplified_mesh_.property(vp_type_sm_, vr) == DGP::CREASE_VFT) {
					// ϵ������Ӧ����vl, vr, vrr: 25, 22, 1.(sum=48)
					// �ȴ���vr���.
					std::vector<TriMesh::VHandle > temp;
					for (TriMesh::VOHIter it(simplified_mesh_, vr); it; ++it) {
						if (simplified_mesh_.status(simplified_mesh_.edge_handle(it.handle())).feature() == true) {
							temp.push_back(simplified_mesh_.to_vertex_handle(it.handle()));
						} 
					}
					if (temp.size() != 2) std::cout << "Error: crease vertex should has only two crease edges.\n";
					TriMesh::VHandle vrr;
					if (temp[0] == vl) vrr = temp[1];
					else vrr = temp[0];
					// ϵ������ vl, vr, vrr: 25, 22, 1
					int vr_idx = crease_vertex_idx_sm[vr];
					Crc_data[vr_idx][j] = 22;
					OpenMesh::Vec3d cor = refined_simplified_mesh_.point(vl) * 25;
					Qr[j] -= cor[0]; Qr[j + r] -= cor[1]; Qr[j + r+r] -= cor[2]; 

					if (simplified_mesh_.property(vp_type_sm_, vrr) == DGP::CREASE_VFT) {
						int vrr_idx = crease_vertex_idx_sm[vrr];
						Crc_data[vrr_idx][j] = 1;
					} else { // dart or corner
						OpenMesh::Vec3d cor = refined_simplified_mesh_.point(vrr) * 1.0;
						Qr[j] -= cor[0]; Qr[j + r] -= cor[1]; Qr[j + r+r] -= cor[2]; 
					}
			}
			else {
				std::cout << "Error: ��Ӧ���е���������İ�." << simplified_mesh_.property(vp_type_sm_, vl) << ", " 
					<< simplified_mesh_.property(vp_type_sm_, vr) << ".\n";
			}
		} // fill in the Crc_data
		// �ⷽ�� Crc Xc = Qr.
		taucs_ccs_matrix *Crc = CreateTaucsMatrixFromColumns(Crc_data, r, TAUCS_DOUBLE);//PrintTaucsCCSMatrix(Bmn);
		taucs_ccs_matrix *Crc_t = MatrixTranspose(Crc); //PrintTaucsCCSMatrix(Cmn_t);//n*m //
		taucs_ccs_matrix *Crc_tCrc = Mul2NonSymmMatSymmResult(Crc_t, Crc);  //PrintTaucsCCSMatrix(Cmn_tCmn);//n*n

		double *Crc_tQr = new double[c * 3];//��������n�е�, ��֮ǰŪ����m, �˷������timeȥcheck error.
		MulNonSymmMatrixVector(Crc_t, Qr, Crc_tQr);                // ÿ����һά, 3ά������3��. (n m)*(m 1)=(n 1)
		MulNonSymmMatrixVector(Crc_t, Qr + r, Crc_tQr + c);    
		MulNonSymmMatrixVector(Crc_t, Qr + r + r, Crc_tQr + c + c);

		// call TAUCS to solve the system
		char* solver_options[] = { "taucs.factor.LLT=true", NULL };
		int error_code;
		error_code = taucs_linsolve(Crc_tCrc, NULL, 3, Xc, Crc_tQr, solver_options, NULL);//���2�Ǳ�ʾx��A_tb���Ƕ�ά��, �ؼ��������һ������A_tA
		if (error_code != TAUCS_SUCCESS) {
			std::cerr << "Solver failed (" << error_code << ")." << std::endl;
		} else {
			std::cout << "Equations solved successfully.\n";
			for (int i = 0; i < c; ++i) { //std::cout << Xn[i] << ", " << Xn[i+n] << ", " <<  Xn[i+n+n] << ".\t";
				initial_control_mesh_.set_point(crease_vertex_sm[i], 
					OpenMesh::Vec3d(Xc[i], Xc[i+c], Xc[i+c+c]));
			}
		} 
	} //2009-04-30������2��, �������crease vertex��Ч��Ҫ��!!!

	// �����smooth and dart vertex ----------
	// 2009-04-30�Թ�������2����smooth�Ĵ���copy������, Ч�������񷽷�2��������.
	// �����������������AmmBmn, ֮����ȡ������Ҫ��ϵ���ŵ������ұ�.
	std::vector<std::map<int, taucsType> > Bmn_data(n);//n����. �����Bmn_data[�к�][�к�]
	//std::cout << "2, " << m << ", " << n << ", " << n_corv_of_sm_ << ", " << simplified_mesh_.n_edges() << "=" << es_vec_.size() << ".\n";
	for (int i = 0; i < n; ++i) { //һ��һ�е���д, ����дBmn_data��ǰn��, Ҳ���Ǹ�����n���ɵ������
		//TriMesh::VertexHandle vh(simplified_mesh_.vertex_handle(i)); //���濪ʼ��дBmn_data�еĵ�i�е�, ��Ӧ��Xm�еĵ�i������(vh)
		//if (vh != Xm_vertex_[i]) std::cout << "Error: Xm_vertex_[i] should be equal to simplified_mesh_.vertex_handle(i), when i < n.\n";
		//std::cout << "i " << i << ", " << simplified_mesh_.property(vp_type_sm_, vh) << ". ";
		TriMesh::VertexHandle vh(Xm_vertex_[i]);
		if (simplified_mesh_.property(vp_type_sm_, vh) == DGP::CORNER_VFT) {
			Bmn_data[vh.idx()][i] = 1;
		} else 
			if (simplified_mesh_.property(vp_type_sm_, vh) == DGP::CREASE_VFT) {
				std::vector<TriMesh::VertexHandle> vec; //����crease edge����Ӧ�Ķ���.
				for (TriMesh::VOHIter voh_it(simplified_mesh_, vh); voh_it; ++voh_it) {
					if (simplified_mesh_.status(simplified_mesh_.edge_handle(voh_it.handle())).feature()) 
						vec.push_back(simplified_mesh_.to_vertex_handle(voh_it.handle()));				
				}
				if (vec.size() != 2) std::cout << "Error: ����Ӧ��ֻ������crease edge�Ŷ�.\n" << vec.size() << "; " << vec[0] << ", " << vec[1] << ".\n";
				//std::cout << vec.size() << ", " << vec[0] << ", " << vec[1] << ".\n";
				Bmn_data[vh.idx()][i] = 0.75; Bmn_data[vec[0].idx()][i] = 0.125; Bmn_data[vec[1].idx()][i] = 0.125; //
			} else if (simplified_mesh_.is_boundary(vh)) {                   //
				TriMesh::HalfedgeHandle heh = simplified_mesh_.halfedge_handle(vh); //de
				if (simplified_mesh_.is_boundary(heh) == false) { std::cout << "Error: ���Ҫ���Ǳ߽�Ŷ���.\n"; }
				TriMesh::VertexHandle v0 = simplified_mesh_.to_vertex_handle(heh);
				TriMesh::VertexHandle v1 = simplified_mesh_.from_vertex_handle(simplified_mesh_.prev_halfedge_handle(heh));
				Bmn_data[vh.idx()][i] = 0.75; Bmn_data[v0.idx()][i] = 0.125; Bmn_data[v1.idx()][i] = 0.125; //��crease vertex�ĸ���ģʽ��һ����
			} else if (simplified_mesh_.property(vp_type_sm_, vh) == DGP::SMOOTH_VFT || simplified_mesh_.property(vp_type_sm_, vh) == DGP::DART_VFT) {
				int k = simplified_mesh_.valence(vh);
				double beta = 0, alpha = 0;
				if (k == 3) { beta = 0.1875; alpha = 0.5625; } // 3.0/16.0 = 0.1875,   9.0/16.0 = 0.5625
				else if (k == 6) { beta = 0.0625; alpha = 0.375; } // 1.0/16.0 = 0.0625,    3.0/8.0 = 0.375
				else { 
					double inv_v = 1.0/double(k); double t = (3.0 + 2.0 * cos( 2.0 * M_PI * inv_v) );
					alpha  = (40.0 - t * t)/64.0; beta = inv_v * alpha;
				}
				Bmn_data[vh.idx()][i] = 1.0 - alpha; // ��Ӧ��Xm_vertex_[i]���i������, �ǵ�i��, �����ǵڼ��о�����Χ���ھӶ������(�����Ƕ�������)
				for (TriMesh::VertexOHalfedgeIter voh_it(simplified_mesh_, vh); voh_it; ++voh_it) {
					Bmn_data[simplified_mesh_.to_vertex_handle(voh_it).idx()][i] = beta;
				}
			} else std::cout << "Error: û���������͵Ķ�������Ŷ��1.\n";		
	} 
	for (int i = n; i < m; ++i) { // ��дBmn_data�ĵ�n��m-1��, ��es_vec_.size()��ô����
		TriMesh::EdgeHandle eh = es_vec_[i - n];
		if (simplified_mesh_.status(eh).feature() || simplified_mesh_.is_boundary(eh)) {
			TriMesh::VertexHandle v0 = simplified_mesh_.to_vertex_handle(simplified_mesh_.halfedge_handle(eh, 0));
			TriMesh::VertexHandle v1 = simplified_mesh_.to_vertex_handle(simplified_mesh_.halfedge_handle(eh, 1));
			Bmn_data[v0.idx()][i] = 0.5; Bmn_data[v1.idx()][i] = 0.5; 

		} else { // smooth edge, new vertex = 3/8 * v0 + 3/8 * v1 + 1/8 * v2 + 1/8 * v3
			TriMesh::VertexHandle v0 = simplified_mesh_.to_vertex_handle(simplified_mesh_.halfedge_handle(eh, 0));
			TriMesh::VertexHandle v1 = simplified_mesh_.to_vertex_handle(simplified_mesh_.halfedge_handle(eh, 1));
			TriMesh::VertexHandle v2 = simplified_mesh_.to_vertex_handle(simplified_mesh_.next_halfedge_handle(simplified_mesh_.halfedge_handle(eh, 0)));
			TriMesh::VertexHandle v3 = simplified_mesh_.to_vertex_handle(simplified_mesh_.next_halfedge_handle(simplified_mesh_.halfedge_handle(eh, 1)));//std::cout << v0 << ", " << v1 << v2 << ", " << v3 << ";";
			Bmn_data[v0.idx()][i] = 0.375; Bmn_data[v1.idx()][i] = 0.375; Bmn_data[v2.idx()][i] = 0.125; Bmn_data[v3.idx()][i] = 0.125;  
		}
	}  
	// ����Amm_data
	std::vector<std::map<int, taucsType> > Amm_data(m);//���m��ָm��
	for (int i = 0; i < m; ++i) { //һ��һ�е���д, ���m��ָm��, ��i�ж�Ӧ��Xm�еĵ�i������
		TriMesh::VertexHandle vh = Xm_vertex_[i]; //����refined simplified mesh�����ӹ�ϵ�Լ��������, ���refined simplified mesh���е��limit position.
		//if (vh != refined_simplified_mesh_.vertex_handle(i)) std::cout << "Error: Xm�����m������Ӧ�ú�refined simplified mesh�еĶ����Ӧ��.\n";

		if (refined_simplified_mesh_.property(vp_type_rsm_, vh) == DGP::CORNER_VFT) {
			Amm_data[vh.idx()][i] = 1;
		} else 
			if (refined_simplified_mesh_.property(vp_type_rsm_, vh) == DGP::CREASE_VFT) {
				std::vector<TriMesh::VertexHandle> vec; //����crease edge����Ӧ�Ķ���.
				for (TriMesh::VOHIter voh_it(refined_simplified_mesh_, vh); voh_it; ++voh_it) {
					if (refined_simplified_mesh_.status(refined_simplified_mesh_.edge_handle(voh_it.handle())).feature() == true) 
						vec.push_back(refined_simplified_mesh_.to_vertex_handle(voh_it.handle()));				
				}
				if (vec.size() != 2) std::cout << "Error: ����Ӧ��ֻ������crease edge�Ŷ�(2).\n";
				Amm_data[vh.idx()][i] = 0.666667; Amm_data[vec[0].idx()][i] = 0.166667; Amm_data[vec[1].idx()][i] = 0.166667; //
			} else if (refined_simplified_mesh_.is_boundary(vh)) {
				TriMesh::HalfedgeHandle heh = refined_simplified_mesh_.halfedge_handle(vh); //de
				if (refined_simplified_mesh_.is_boundary(heh) == false) { std::cout << "Error: ���Ҫ���Ǳ߽�Ŷ���2.\n"; }
				TriMesh::VertexHandle v0 = refined_simplified_mesh_.to_vertex_handle(heh);
				TriMesh::VertexHandle v1 = refined_simplified_mesh_.from_vertex_handle(refined_simplified_mesh_.prev_halfedge_handle(heh));
				Amm_data[vh.idx()][i] = 0.666667; Amm_data[v0.idx()][i] = 0.166667; Amm_data[v1.idx()][i] = 0.166667; //

			} else if (refined_simplified_mesh_.property(vp_type_rsm_, vh) == DGP::SMOOTH_VFT || refined_simplified_mesh_.property(vp_type_rsm_, vh) == DGP::DART_VFT) {
				int k = refined_simplified_mesh_.valence(vh);
				double inv_v = 1.0/double(k); double t = (3.0 + 2.0 * cos( 2.0 * M_PI * inv_v) ); // double beta = inv_v * (40.0 - t * t)/64.0;
				double k_alpha = 1.0 / ( 24.0/(40.0-t*t) + 1); double alpha = inv_v * k_alpha;
				for (TriMesh::VertexVertexIter vv_it(refined_simplified_mesh_, vh); vv_it; ++vv_it) {
					Amm_data[vv_it.handle().idx()][i] = alpha; 
				}
				Amm_data[vh.idx()][i] = 1.0 - k_alpha; 
			} else std::cout << "Error: û���������͵Ķ�������Ŷ��2.\n";	
	} 

	taucs_ccs_matrix *Bmn = CreateTaucsMatrixFromColumns(Bmn_data, m, TAUCS_DOUBLE);//PrintTaucsCCSMatrix(Bmn);
	taucs_ccs_matrix *Amm = CreateTaucsMatrixFromColumns(Amm_data, m, TAUCS_DOUBLE);//PrintTaucsCCSMatrix(Amm);	
	//taucs_ccs_matrix *Bmn_t = MatrixTranspose(Bmn); //PrintTaucsCCSMatrix(Bmn_t);// n*m
	//taucs_ccs_matrix *Amm_t = MatrixTranspose(Amm); //PrintTaucsCCSMatrix(Amm_t);// m*m
	taucs_ccs_matrix *Cmn = Mul2NonSymmetricMatrices(Amm, Bmn); //PrintTaucsCCSMatrix(Cmn); 
	if (Cmn->n != n) { std::cout << "Error: ����Ӧ����һ����.\n"; }
	// ����ķ�����ʽΪCmn => ����ʽ����Cmn_data -> ȥ��Cmn__data => ����ʽ���� Cmn__arr
	// -> delete row Cm_n_arr => ����ʽ����Cm_n__data.=> taucs_ccs Cmn_ 
	std::vector<std::map<int, taucsType> > Cmn_data(Cmn->n);
	CreateColumnsFromTaucsMatrix(Cmn, Cmn_data, m);

	double *Qm = new double[m * 3];
	for (int i = 0; i < m; ++i) {
		TriMesh::Point p = refined_simplified_mesh_.point(Xm_vertex_[i]);//std::cout << p << "\t";
		Qm[i] = p[0]; Qm[i+m] = p[1]; Qm[i+m+m] = p[2];
	}

	std::vector<TriMesh::VertexHandle> crease_corner_vertex_sm, smooth_dart_vertex_sm;
	crease_corner_vertex_sm.clear(); smooth_dart_vertex_sm.clear();
	for (TriMesh::VertexIter v_it(simplified_mesh_.vertices_begin()), v_end(simplified_mesh_.vertices_end());
		v_it != v_end; ++ v_it) {
			if (simplified_mesh_.property(vp_type_sm_, v_it) == DGP::CREASE_VFT || 
				simplified_mesh_.property(vp_type_sm_, v_it) == DGP::CORNER_VFT) {
					crease_corner_vertex_sm.push_back(v_it.handle());	  
			} else 
				smooth_dart_vertex_sm.push_back(v_it.handle());
	}
	if (n-crease_corner_vertex_sm.size() != smooth_dart_vertex_sm.size()) {
		std::cout << "Error: n-crease_corner_vertex_sm.size() != smooth_dart_vertex_sm.size().\n"; 
	}
	if (smooth_dart_vertex_sm.size() > 0) {
		// ȥ�У���crease��corner���Ӧ����ȥ����
		std::map<int,taucsType>::const_iterator rit;
		for (int i = 0; i < crease_corner_vertex_sm.size(); ++i) {
			int c = crease_corner_vertex_sm[i].idx();
			for (rit = Cmn_data[c].begin();rit!= Cmn_data[c].end();++rit) { // how many rows in this column.
				int r = rit->first;
				double a = rit->second;
				OpenMesh::Vec3d p = initial_control_mesh_.point(crease_corner_vertex_sm[i]) * a; 
				Qm[r] -= p[0]; Qm[r + m] -= p[1]; Qm[r + m+m] -= p[2]; // this is the key. 

				Cmn_data[c][r] = 0; // no matter, will be deleted later.
			}
		} 

		std::vector<std::map<int, taucsType> > Cmn__data(smooth_dart_vertex_sm.size()); // num of column.
		for (int i = 0; i < smooth_dart_vertex_sm.size(); ++i) {
			Cmn__data[i] = Cmn_data[smooth_dart_vertex_sm[i].idx()];
		}

		// ת�ɾ�����ʽ����ȥ��, 
		std::vector<std::vector<double> > Cmn__arr(m, std::vector<double>(smooth_dart_vertex_sm.size(), 0));
		for (int i = 0; i < smooth_dart_vertex_sm.size(); ++i) {
			for (std::map<int, taucsType>::iterator it = Cmn__data[i].begin(), end = Cmn__data[i].end(); 
				it != end; ++it) {
					Cmn__arr[it->first][i] = it->second;  
			}
		}
		// ������Ч��
		std::vector<int > valid_row, invalid_row; // valid row of Qm, if this row has one non-zero.
		for (int i = 0; i < m; ++i) {
			bool all_zero = true;
			for (int j = 0; j < smooth_dart_vertex_sm.size(); ++j) {
				if (Cmn__arr[i][j] != 0) {
					all_zero = false;
					valid_row.push_back(i); 
					break;
				}
			}
			if (all_zero) invalid_row.push_back(i);
		} // valid_row.size() < m;
		if (valid_row.size() + invalid_row.size() != m) std::cout << "Error: ��������.\n";
		// ȥ��, in row array form
		std::vector<std::vector<double> > Cm_n__arr(valid_row.size(), std::vector<double>(smooth_dart_vertex_sm.size(), 0));
		for (int i = 0; i < valid_row.size(); ++i)
			Cm_n__arr[i] = Cmn__arr[valid_row[i]];
		//�������Է���valid_row.size() != smooth_dart_vertex_rsm.size() 
		//����fandiskģ�� valid_row.size() < smooth_dart_vertex_rsm.size() 
		//����bunnyģ��   valid_row.size() > smooth_dart_vertex_rsm.size() 

		// only change from row array to the column array form
		std::vector<std::map<int, taucsType> > Cm_n__data(smooth_dart_vertex_sm.size());
		for (int i = 0; i < smooth_dart_vertex_sm.size(); ++i) {
			for (int j = 0; j < valid_row.size(); ++j) { //row no.
				if (Cm_n__arr[j][i] != 0) {
					Cm_n__data[i][j] = Cm_n__arr[j][i]; 
				}
			} 
		}
		double *Qm_valid = new double[valid_row.size() * 3];
		for (int i = 0; i < valid_row.size(); ++i) {
			Qm_valid[i] = Qm[valid_row[i]];
			Qm_valid[i + valid_row.size()] = Qm[valid_row[i] + m];
			Qm_valid[i + valid_row.size()+ valid_row.size()] = Qm[valid_row[i] + m + m];
		}
		taucs_ccs_matrix *Cmn_ = CreateTaucsMatrixFromColumns(Cm_n__data, valid_row.size(), TAUCS_DOUBLE); // valid_row.size()*smooth_dart_vertex_sm.size()
		taucs_ccs_matrix *Cmn_t = MatrixTranspose(Cmn_); 
		taucs_ccs_matrix *Cmn_tCmn = Mul2NonSymmMatSymmResult(Cmn_t, Cmn_);   
		double *Cmn_tQm = new double[smooth_dart_vertex_sm.size() * 3];//��������n�е�, ��֮ǰŪ����m, �˷������timeȥcheck error.
		MulNonSymmMatrixVector(Cmn_t, Qm_valid, Cmn_tQm);                // ÿ����һά, 3ά������3��. (n m)*(m 1)=(n 1)
		MulNonSymmMatrixVector(Cmn_t, Qm_valid + valid_row.size(), Cmn_tQm + smooth_dart_vertex_sm.size());    
		MulNonSymmMatrixVector(Cmn_t, Qm_valid + valid_row.size() + valid_row.size(), Cmn_tQm + smooth_dart_vertex_sm.size() + smooth_dart_vertex_sm.size());

		//��Ҫ����δ֪Xn
		double *Xn = new double[smooth_dart_vertex_sm.size() * 3];

		// call TAUCS to solve the systemd
		char* solver_options[] = { "taucs.factor.LLT=true", NULL };
		int error_code;
		error_code = taucs_linsolve(Cmn_tCmn, NULL, 3, Xn, Cmn_tQm, solver_options, NULL);//���2�Ǳ�ʾx��A_tb���Ƕ�ά��, �ؼ��������һ������A_tA
		if (error_code != TAUCS_SUCCESS) {
			std::cerr << "Solver failed (" << error_code << ")." << std::endl;
		} else {
			std::cout << "Equations solved successfully.\n"; 
			for (int i = 0; i < smooth_dart_vertex_sm.size(); ++i)
				initial_control_mesh_.set_point(smooth_dart_vertex_sm[i], OpenMesh::Vec3d(Xn[i], Xn[i+smooth_dart_vertex_sm.size()], Xn[i+smooth_dart_vertex_sm.size()+smooth_dart_vertex_sm.size()]));		
		}
	} // there're no smooth / dart vertices.


	set_mesh_toberendered(&initial_control_mesh_, INITIAL_CONTROL);
	//set_mesh_toberendered(&refined_simplified_mesh_, REFINED_SIMPLIFIED);
	std::cout << "========������������������============\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.
}
void DecimationModel::create_initialcontrolmesh_process3() {
	// ����汾�����̷ֽ⣬
	std::cout << "\n";
	std::cout << "========��ʼ��������������v3==========\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.

	if (Xm_vertex_.size() != (simplified_mesh_.n_vertices() + es_vec_.size())) {
		std::cout << "Error: �����ϲ�ƥ��1." << Xm_vertex_.size() << " != " << simplified_mesh_.n_vertices()
			<< " + " << es_vec_.size() << ". " << refined_simplified_mesh_.n_vertices() << ".\n";
	}
	if (Xm_vertex_.size() != refined_simplified_mesh_.n_vertices()) std::cout << "Error: �����ϲ�ƥ��2.\n";
	std::vector<TriMesh::VHandle> smooth_dart_vertex_rsm; smooth_dart_vertex_rsm.clear();
	int idx_Xm =0;
	for (TriMesh::VIter v_it(refined_simplified_mesh_.vertices_begin()), v_end(refined_simplified_mesh_.vertices_end());
		v_it != v_end; ++v_it) {
			if (refined_simplified_mesh_.property(vp_type_rsm_, v_it) == DGP::DART_VFT
				|| refined_simplified_mesh_.property(vp_type_rsm_, v_it) == DGP::SMOOTH_VFT) {
					smooth_dart_vertex_rsm.push_back(v_it.handle()); 
			}
			if (v_it.handle() != Xm_vertex_[idx_Xm++]) std::cout << "Error: this should be the same.\n";
	}
	const int m = Xm_vertex_.size(), n = simplified_mesh_.n_vertices();//m���� n����, mҲ�Ƿ������з��̵ĸ���.
	
	std::vector<OpenMesh::Vec3d> laplacian_Qm(m, OpenMesh::Vec3d(0,0,0)); // �������rsm��ÿһ�������laplacian����.
	for (int i = 0; i < refined_simplified_mesh_.n_vertices(); ++i) {
		// every vertex, no matter what type it is, calc it's laplacian.
		TriMesh::VHandle vh = refined_simplified_mesh_.vertex_handle(i);
		OpenMesh::Vec3d p(0,0,0);
		double weight = -1.0 / refined_simplified_mesh_.valence(vh);
		for (TriMesh::VVIter vv_it(refined_simplified_mesh_, vh); vv_it; ++vv_it) {
			p += refined_simplified_mesh_.point(vv_it); 
		}
		laplacian_Qm[i] = refined_simplified_mesh_.point(vh) + p * weight; 
	}
	laplacian_ = laplacian_Qm;
	std::vector< std::map< int, taucsType > > lmm_data(m);     // in column
	std::vector< std::map<int , taucsType > > lmm_data_row(m); // in row 
	idx_Xm = 0;
	for (TriMesh::VIter v_it(refined_simplified_mesh_.vertices_begin()), v_end(refined_simplified_mesh_.vertices_end());
		v_it != v_end; ++v_it) {
			if (v_it.handle().idx() != idx_Xm++) std::cout << "Error: idx_Xm should ..\n";
			double weight = -1.0 / refined_simplified_mesh_.valence(v_it);
			for (TriMesh::VVIter vv_it(refined_simplified_mesh_, v_it); vv_it; ++vv_it) {
				//if ((refined_simplified_mesh_.property(vp_type_rsm_, vv_it) == DGP::CREASE_VFT) 
				//	|| (refined_simplified_mesh_.property(vp_type_rsm_, vv_it) == DGP::CORNER_VFT)) 
				//{
				//	laplacian_Qm[v_it.handle().idx()] -= refined_simplified_mesh_.point(vv_it) * weight;
				//} else { //vv_it is smooth or dart.
				lmm_data[vv_it.handle().idx()][v_it.handle().idx()] = weight; 
				lmm_data_row[v_it.handle().idx()][vv_it.handle().idx()] = weight; 

				//}				
			} 
			lmm_data[v_it.handle().idx()][v_it.handle().idx()] = 1;
			lmm_data_row[v_it.handle().idx()][v_it.handle().idx()] = 1; 

	}
	if (laplacian_Qm.size() != m) std::cout << "Error: laplacian_Qm.size() shoule equal to m.\n";

	 
	// To figure out the geometry of the initial_control_mesh_,
	// be ware that, the topology of the icm is the same as simplified_mesh_. 
	initial_control_mesh_ = simplified_mesh_;
	initial_control_mesh_.add_property(vp_type_icm_); 
	initial_control_mesh_.add_property(vp_closestp_icm_);
	initial_control_mesh_.add_property(vp_distances_);

	for (TriMesh::VertexIter v_it = simplified_mesh_.vertices_begin(), v_end = simplified_mesh_.vertices_end(); v_it != v_end; ++v_it) { //���Ƶ�����
		initial_control_mesh_.property(vp_type_icm_, v_it) = simplified_mesh_.property(vp_type_sm_, v_it);
	}
	for (TriMesh::EdgeIter e_it = simplified_mesh_.edges_begin(), e_end = simplified_mesh_.edges_end(); e_it != e_end; ++e_it) {//���Ʊ�����
		initial_control_mesh_.status(e_it).set_feature(simplified_mesh_.status(e_it.handle()).feature());
	} //std::cout << "1. ";

	calc_crease_first(m, n); // ����simplified_mesh_��ÿһ��crease vertex ÿһ��crease edge����Ӧ��һ������, ����crease vertex.	

	// �����smooth and dart vertex ----------
	// 2009-04-30�Թ�������2����smooth�Ĵ���copy������, Ч�������񷽷�2��������.
	// �����������������AmmBmn, ֮����ȡ������Ҫ��ϵ���ŵ������ұ�.
	std::vector<std::map<int, taucsType> > Bmn_data(n);//n����. �����Bmn_data[�к�][�к�]
	std::vector<std::map<int, taucsType> > Amm_data(m);// ����Amm_data//���m��ָm��
	calc_AmmBmn(Amm_data, Bmn_data); 

	taucs_ccs_matrix *Bmn = CreateTaucsMatrixFromColumns(Bmn_data, m, TAUCS_DOUBLE);//PrintTaucsCCSMatrix(Bmn);
	taucs_ccs_matrix *Amm = CreateTaucsMatrixFromColumns(Amm_data, m, TAUCS_DOUBLE);//PrintTaucsCCSMatrix(Amm);	
	//taucs_ccs_matrix *Bmn_t = MatrixTranspose(Bmn); //PrintTaucsCCSMatrix(Bmn_t);// n*m
	//taucs_ccs_matrix *Amm_t = MatrixTranspose(Amm); //PrintTaucsCCSMatrix(Amm_t);// m*m
	taucs_ccs_matrix *Cmn = Mul2NonSymmetricMatrices(Amm, Bmn); //PrintTaucsCCSMatrix(Cmn); 
	if (Cmn->n != n) { std::cout << "Error: ����Ӧ����һ����.\n"; }
	// ����ķ�����ʽΪCmn => ����ʽ����Cmn_data -> ȥ��Cmn__data => ����ʽ���� Cmn__arr
	// -> delete row Cm_n_arr => ����ʽ����Cm_n__data.=> taucs_ccs Cmn_ 
	std::vector<std::map<int, taucsType> > Cmn_data(Cmn->n); //������
	CreateColumnsFromTaucsMatrix(Cmn, Cmn_data, m);  
	std::vector<std::map<int, taucsType> > Cmn_data_row(m); //������.
	for (int i = 0; i < Cmn->n; ++i) { //��Cmn_data[i]��
		for (std::map<int, taucsType>::iterator it = Cmn_data[i].begin(), it_end=Cmn_data[i].end();
			it != it_end; ++it) {
				Cmn_data_row[it->first][i] = it->second;
		}
	}

	double *Qm = new double[m * 3];
	std::vector<OpenMesh::Vec3d> Qm_data(m, OpenMesh::Vec3d(0,0,0));
	for (int i = 0; i < m; ++i) {
		TriMesh::Point p = refined_simplified_mesh_.point(Xm_vertex_[i]);//std::cout << p << "\t";
		Qm[i] = p[0]; Qm[i+m] = p[1]; Qm[i+m+m] = p[2];

		Qm_data[i] = p;
	}

	//// �ⷨ1.
	//taucs_ccs_matrix *Cmn_t = MatrixTranspose(Cmn); //PrintTaucsCCSMatrix(Cmn_t);//n*m //
	//taucs_ccs_matrix *Cmn_tCmn = Mul2NonSymmMatSymmResult(Cmn_t, Cmn);  //PrintTaucsCCSMatrix(Cmn_tCmn);//n*n
	//double *Cmn_tQm = new double[n * 3];//��������n�е�, ��֮ǰŪ����m, �˷������timeȥcheck error.
	//MulNonSymmMatrixVector(Cmn_t, Qm, Cmn_tQm);                // ÿ����һά, 3ά������3��. (n m)*(m 1)=(n 1)
	//MulNonSymmMatrixVector(Cmn_t, Qm + m, Cmn_tQm + n);    
	//MulNonSymmMatrixVector(Cmn_t, Qm + m + m, Cmn_tQm + n + n);
	//
	//double *Xn = new double[n * 3];//��Ҫ����δ֪Xn
	//// call TAUCS to solve the system
	//char* solver_options[] = { "taucs.factor.LLT=true", NULL };
	//int error_code;
	//error_code = taucs_linsolve(Cmn_tCmn, NULL, 3, Xn, Cmn_tQm, solver_options, NULL);//���2�Ǳ�ʾx��A_tb���Ƕ�ά��, �ؼ��������һ������A_tA
	//if (error_code != TAUCS_SUCCESS) {
	//	std::cerr << "Solver failed (" << error_code << ")." << std::endl;
	//} else {
	//	std::cout << "Equations solved successfully. CmnXn=Qm\n";
	//	for (int i = 0; i < n; ++i) { //std::cout << Xn[i] << ", " << Xn[i+n] << ", " <<  Xn[i+n+n] << ".\t";
	//		initial_control_mesh_.set_point(Xm_vertex_[i], OpenMesh::Vec3d(Xn[i], Xn[i+n], Xn[i+n+n]));
	//	}
	//}

	// �ⷨ2.
	//std::vector<OpenMesh::Vec3d> Xn_data(n, OpenMesh::Vec3d(0,0,0));
	//if (linear_equation_solover(Cmn_data, Xn_data, Qm_data, m, n, "Cmn_data Xn_data = Qm_data")) {
	//	for (int i = 0; i < n; ++i) { //std::cout << Xn[i] << ", " << Xn[i+n] << ", " <<  Xn[i+n+n] << ".\t";
	//		initial_control_mesh_.set_point(Xm_vertex_[i], Xn_data[i]);
	//	}
	//}

	//laplacian_rcm_ = std::vector<OpenMesh::Vec3d> (m, OpenMesh::Vec3d(0,0,0));
	//for (int i = 0; i < m; ++i) {
	//	for (std::map<int, taucsType>::iterator it = lmm_data_row[i].begin(), it_end=lmm_data_row[i].end();
	//		it != it_end; ++it) {
	//			laplacian_rcm_[i] += OpenMesh::Vec3d(Qm[it->first],Qm[it->first+m],Qm[it->first+m+m]) * it->second;
	//	}
	//} //����lmm�����Ƿ�����, ����Ƿ����rsm�ĵ���laplacian����һ��, ����ǶԵ�.

	std::vector<TriMesh::VertexHandle> crease_corner_vertex_sm, smooth_dart_vertex_sm;
	crease_corner_vertex_sm.clear(); smooth_dart_vertex_sm.clear();
	for (TriMesh::VertexIter v_it(simplified_mesh_.vertices_begin()), v_end(simplified_mesh_.vertices_end());
		v_it != v_end; ++ v_it) {
			if (simplified_mesh_.property(vp_type_sm_, v_it) == DGP::CREASE_VFT || 
				simplified_mesh_.property(vp_type_sm_, v_it) == DGP::CORNER_VFT) {
				crease_corner_vertex_sm.push_back(v_it.handle());	  
			} else 
				smooth_dart_vertex_sm.push_back(v_it.handle());
	}
	if (n-crease_corner_vertex_sm.size() != smooth_dart_vertex_sm.size()) {
		std::cout << "Error: n-crease_corner_vertex_sm.size() != smooth_dart_vertex_sm.size().\n"; 
	}
	if (smooth_dart_vertex_sm.size() > 0) {
		// ȥ�У���crease��corner���Ӧ����ȥ���� 
		for (int i = 0; i < crease_corner_vertex_sm.size(); ++i) {
			int c = crease_corner_vertex_sm[i].idx();
			for (std::map<int,taucsType>::const_iterator rit = Cmn_data[c].begin();rit!= Cmn_data[c].end();++rit) { // how many rows in this column.
				int r = rit->first;
				double a = rit->second;
				OpenMesh::Vec3d p = initial_control_mesh_.point(crease_corner_vertex_sm[i]) * a; 
				Qm[r] -= p[0]; Qm[r + m] -= p[1]; Qm[r + m+m] -= p[2]; // this is the key. 

				Cmn_data[c][r] = 0; // no matter, will be deleted later.
			}
		} 

		std::vector<std::map<int, taucsType> > Cmn__data(smooth_dart_vertex_sm.size()); // num of column.
		std::vector<std::map<int, taucsType> > Cmn__data_row(m); // num of row.
		for (int i = 0; i < smooth_dart_vertex_sm.size(); ++i) {
			Cmn__data[i] = Cmn_data[smooth_dart_vertex_sm[i].idx()];

			for (std::map<int, taucsType>::iterator it=Cmn__data[i].begin(), it_end=Cmn__data[i].end();
				it != it_end; ++it) 
				Cmn__data_row[it->first][i] = it->second;
		}

		// ת�ɾ�����ʽ����ȥ��, 
		std::vector<std::vector<double> > Cmn__arr(m, std::vector<double>(smooth_dart_vertex_sm.size(), 0));
		for (int i = 0; i < smooth_dart_vertex_sm.size(); ++i) {
			for (std::map<int, taucsType>::iterator it = Cmn__data[i].begin(), end = Cmn__data[i].end(); 
				it != end; ++it) {
					Cmn__arr[it->first][i] = it->second;  
			}
		}
		// ������Ч��
		std::vector<int > valid_row, invalid_row; // valid row of Qm, if this row has one non-zero.
		for (int i = 0; i < m; ++i) {
			bool all_zero = true;
			for (int j = 0; j < smooth_dart_vertex_sm.size(); ++j) {
				if (Cmn__arr[i][j] != 0) {
					all_zero = false;
					valid_row.push_back(i); 
					break;
				}
			}
			if (all_zero) invalid_row.push_back(i);
		} // valid_row.size() < m;
		if (valid_row.size() + invalid_row.size() != m) std::cout << "Error: ��������.\n";
		std::vector<std::vector<double> > Cm_n__arr(valid_row.size(), std::vector<double>(smooth_dart_vertex_sm.size(), 0));
		// ȥ��
		for (int i = 0; i < valid_row.size(); ++i)
			Cm_n__arr[i] = Cmn__arr[valid_row[i]];
		//�������Է���valid_row.size() != smooth_dart_vertex_rsm.size() 
		//����fandiskģ�� valid_row.size() < smooth_dart_vertex_rsm.size() 
		//����bunnyģ��   valid_row.size() > smooth_dart_vertex_rsm.size() 

		// only change from row array to the column array form
		std::vector<std::map<int, taucsType> > Cm_n__data(smooth_dart_vertex_sm.size());
		for (int i = 0; i < smooth_dart_vertex_sm.size(); ++i) {
			for (int j = 0; j < valid_row.size(); ++j) { //row no.
				if (Cm_n__arr[j][i] != 0) {
					Cm_n__data[i][j] = Cm_n__arr[j][i]; 
				}
			} 
		}
		double *Qm_valid = new double[valid_row.size() * 3];
		std::vector<OpenMesh::Vec3d> Qm_valid_data(valid_row.size(), OpenMesh::Vec3d(0,0,0));
		for (int i = 0; i < valid_row.size(); ++i) {
			Qm_valid[i] = Qm[valid_row[i]];
			Qm_valid[i + valid_row.size()] = Qm[valid_row[i] + m];
			Qm_valid[i + valid_row.size()+ valid_row.size()] = Qm[valid_row[i] + m + m];

			Qm_valid_data[i] = OpenMesh::Vec3d(Qm_valid[i], Qm_valid[i + valid_row.size()], Qm_valid[i + valid_row.size()+ valid_row.size()]);
		}

		// try to add the laplacian. 
		taucs_ccs_matrix *LCmn = Mul2NonSymmetricMatrices(CreateTaucsMatrixFromColumns(lmm_data, m, TAUCS_DOUBLE), Cmn);
		std::vector<std::map<int, taucsType> > LCmn_data(n);
		CreateColumnsFromTaucsMatrix(LCmn, LCmn_data, m); 
		std::vector<std::map<int, taucsType> > LCmn_data_row(m); //������.
		for (int i = 0; i < n; ++i) { //��Cmn_data[i]��
			for (std::map<int, taucsType>::iterator it = LCmn_data[i].begin(), it_end=LCmn_data[i].end();
				it != it_end; ++it) {
					LCmn_data_row[it->first][i] = it->second;
			}
		}
		double *LQm = new double[m*3];
		std::vector<OpenMesh::Vec3d> LQm_data(m, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < m; ++i) {
			LQm[i] = laplacian_Qm[i][0]; 
			LQm[i+m] = laplacian_Qm[i][1]; 
			LQm[i+m+m] = laplacian_Qm[i][2]; 

			LQm_data[i] = laplacian_Qm[i]; 
		} 

		std::vector<OpenMesh::Vec3d> Xn_data(n, OpenMesh::Vec3d(0,0,0));
		if (linear_equation_solover(LCmn_data, Xn_data, LQm_data, m, n, "LCmn*Xn=LQm")) {
			for (int i = 0; i < m; ++i) {
				OpenMesh::Vec3d p(0,0,0);
				for (std::map<int, taucsType>::iterator it = LCmn_data_row[i].begin(), it_end=LCmn_data_row[i].end();
					it != it_end; ++it) {
						p += Xn_data[it->first] * it->second;
				} 
				/*std::cout << (p - LQm_data[i]).norm() << ", ";*/
			}
		}
		
		// ȥ�У���crease��corner���Ӧ����ȥ���� 
		for (int i = 0; i < crease_corner_vertex_sm.size(); ++i) {
			int c = crease_corner_vertex_sm[i].idx();
			for (std::map<int,taucsType>::const_iterator rit = LCmn_data[c].begin();rit!= LCmn_data[c].end();++rit) { // how many rows in this column.
				int r = rit->first;
				double a = rit->second;
				OpenMesh::Vec3d p = initial_control_mesh_.point(crease_corner_vertex_sm[i]) * a; 
				LQm[r] -= p[0]; LQm[r + m] -= p[1]; LQm[r + m+m] -= p[2]; // this is the key. 
				LQm_data[r] -= p;

				LCmn_data[c][r] = 0; // no matter, will be deleted later.
			}
		} 
		std::vector<std::map<int, taucsType> > LCmn__data(smooth_dart_vertex_sm.size()); // num of column.
		std::vector<std::map<int, taucsType> > LCmn__data_row(m); // num of row.
		for (int i = 0; i < smooth_dart_vertex_sm.size(); ++i) {
			LCmn__data[i] = LCmn_data[smooth_dart_vertex_sm[i].idx()];

			for (std::map<int, taucsType>::iterator it=LCmn__data[i].begin(), it_end=LCmn__data[i].end();
				it != it_end; ++it) 
				LCmn__data_row[it->first][i] = it->second;
		}

		std::vector<std::map<int, taucsType> > D2mn__data(smooth_dart_vertex_sm.size());
		for (int c = 0; c < smooth_dart_vertex_sm.size(); ++c) {
			for (std::map<int, taucsType>::iterator it=Cmn__data[c].begin(), it_end=Cmn__data[c].end();
				it != it_end; ++it) 
				D2mn__data[c][it->first] = it->second;
		}
		for (int c = 0; c < smooth_dart_vertex_sm.size(); ++c) {
			for (std::map<int, taucsType>::iterator it=LCmn__data[c].begin(), it_end=LCmn__data[c].end();
				it != it_end; ++it) 
				D2mn__data[c][it->first+m] = it->second;
		}
		double *Q2m = new double[m*2*3];
		for (int i = 0; i < m; ++i) {
			Q2m[i] = Qm[i];
			Q2m[i+2*m] = Qm[i+m];
			Q2m[i+4*m] = Qm[i+m+m];
		}
		for (int i = m; i < 2*m; ++i) {
			Q2m[i] = LQm_data[i-m][0];
			Q2m[i+2*m] = LQm_data[i-m][1];
			Q2m[i+4*m] = LQm_data[i-m][2];
		}
		
		std::vector<OpenMesh::Vec3d> Q2m_data(2*m, OpenMesh::Vec3d(0,0,0));
		for (int i = 0; i < 2*m; ++i) {
			Q2m_data[i][0] =  Q2m[i];
			Q2m_data[i][1] =  Q2m[i+2*m];
			Q2m_data[i][2] =  Q2m[i+4*m];
		}
		std::vector<OpenMesh::Vec3d> Xn__data(smooth_dart_vertex_sm.size(), OpenMesh::Vec3d(0,0,0));
		if (linear_equation_solover(D2mn__data, Xn__data, Q2m_data, 2*m, smooth_dart_vertex_sm.size(), "D2mn_data*Xn__data=Q2m_data")) {
		//if (linear_equation_solover(LCmn__data, Xn__data, laplacian_Qm, m, smooth_dart_vertex_sm.size())) {
				for (int i =0; i < smooth_dart_vertex_sm.size(); ++i) {
				initial_control_mesh_.set_point(smooth_dart_vertex_sm[i], Xn__data[i]);
			}
		}

		
		
		//taucs_ccs_matrix *D2mn_ = CreateTaucsMatrixFromColumns(D2mn__data, 2*m, TAUCS_DOUBLE);
		//taucs_ccs_matrix *D2mn_t = MatrixTranspose(D2mn_); 
		//taucs_ccs_matrix *D2mn_tD2mn = Mul2NonSymmMatSymmResult(D2mn_t, D2mn_); 
		//double *D2mn_tQ2m = new double[smooth_dart_vertex_sm.size()*3];
		//MulNonSymmMatrixVector(D2mn_t, Q2m, D2mn_tQ2m);                // ÿ����һά, 3ά������3��. (n m)*(m 1)=(n 1)
		//MulNonSymmMatrixVector(D2mn_t, Q2m + 2*m, D2mn_tQ2m + smooth_dart_vertex_sm.size());    
		//MulNonSymmMatrixVector(D2mn_t, Q2m + 4*m, D2mn_tQ2m + 2*smooth_dart_vertex_sm.size());

		//��Ҫ����δ֪Xn
		//double *Xn = new double[smooth_dart_vertex_sm.size() * 3];
		// call TAUCS to solve the systemd
		//char* solver_options[] = { "taucs.factor.LLT=true", NULL };
		//int error_code;
		//error_code = taucs_linsolve(D2mn_tD2mn, NULL, 3, Xn, D2mn_tQ2m, solver_options, NULL);//���2�Ǳ�ʾx��A_tb���Ƕ�ά��, �ؼ��������һ������A_tA
		//if (error_code != TAUCS_SUCCESS) {
		//	std::cerr << "Solver failed (" << error_code << ")." << std::endl;
		//} else {
		//	std::cout << "Equations solved successfully.\n"; 
		//	for (int i = 0; i < smooth_dart_vertex_sm.size(); ++i)
		//		initial_control_mesh_.set_point(smooth_dart_vertex_sm[i], OpenMesh::Vec3d(Xn[i], Xn[i+smooth_dart_vertex_sm.size()], Xn[i+smooth_dart_vertex_sm.size()+smooth_dart_vertex_sm.size()]));		
		//}

	} // there're no smooth / dart vertices.
	

	set_mesh_toberendered(&initial_control_mesh_, INITIAL_CONTROL); 
	std::cout << "========������������������============\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.
}
// ----------------------
void DecimationModel::loop_subdivision_process() {
	std::cout << "\n";
	std::cout << "========��ʼ����һ��loopϸ��============\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.

	// do the thing here. 
	DGP::Hoppe94LoopSubT<TriMesh> subdi_obj;
	subdi_obj.attach(initial_control_mesh_, vp_type_icm_);
	subdi_obj(1);
	subdi_obj.detach();

	initial_control_mesh_.update_normals();
	set_mesh_toberendered(&initial_control_mesh_, LOOP_ONE);
	 
	// ��controlmesh_forbounding_ ���ǵ���"���һ��loop subd�õ���mesh", 
	// ��Ϊ����һ����"���һ��loop subd�õ���mesh"��limit surfʱ, ��ϣ��������bounding.
	//if (controlmesh_forbounding_.empty() == false) controlmesh_forbounding_.clear(); // renc 20150829
	if (controlmesh_forbounding_.n_vertices()) controlmesh_forbounding_.clear(); 
	controlmesh_forbounding_ = initial_control_mesh_; 

	// -------
	std::cout << "ϸ�ֻ�õ����� " << initial_control_mesh_.n_vertices() << " vertices, " << initial_control_mesh_.n_edges() << " edges, " << initial_control_mesh_.n_faces()    << " faces\n";

	int fe = 0;
	for (TriMesh::EIter e_it(initial_control_mesh_.edges_begin()), e_end(initial_control_mesh_.edges_end()); e_it != e_end; ++e_it) {
		if (initial_control_mesh_.status(e_it).feature()) ++fe;
	}
	int sv = 0, dv = 0, crv = 0, cov = 0, bv = 0;
	for (TriMesh::VIter v_it = initial_control_mesh_.vertices_begin(), v_end = initial_control_mesh_.vertices_end(); v_it != v_end; ++ v_it) {
		if (initial_control_mesh_.property(vp_type_icm_, v_it) == DGP::SMOOTH_VFT) sv++;
		if (initial_control_mesh_.property(vp_type_icm_, v_it) == DGP::DART_VFT) dv++;
		if (initial_control_mesh_.property(vp_type_icm_, v_it) == DGP::CREASE_VFT) crv++;
		if (initial_control_mesh_.property(vp_type_icm_, v_it) == DGP::CORNER_VFT) cov++;
		if (initial_control_mesh_.is_boundary(v_it) == true) bv++; 
	}
	std::cout << "������: " << fe << ".\n";
	std::cout << "������: smooth: " << sv << ", dart " << dv << ", crease: " << crv << ", corner: " << cov << ", boundary: " << bv << ".\n";

	std::cout << "========��������һ��loopϸ��============\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.
}
// ---------------
void DecimationModel::create_limit_surface_process()
{
	std::cout << "\n";
	std::cout << "========��ʼ����limit surface===========\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.

	DGP::Hoppe94LoopSubT<TriMesh> subdi_obj;
	subdi_obj.attach(initial_control_mesh_, vp_type_icm_);
	subdi_obj.set_limit_position();
	subdi_obj.detach();

	initial_control_mesh_.update_normals();
	set_mesh_toberendered(&initial_control_mesh_, LIMIT);
	std::cout << "========��������limit surface===========\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.
}
// ---------------
bool DecimationModel::fitting_process() { //ֻ��Ϊ�˷������.
	while (fitting_terminate_ == false) {
		resample_edge_midpoint_process();     //������Щ��������Ϊ��refinementʱ����Ҫ��������,���Բ���.add_property(...);
		midpointsubdivision_1_process();
		if (feature_considered_) {
			create_initialcontrolmesh_process3();
		} else create_initialcontrolmesh_process2();
		loop_subdivision_process();
		//create_limit_surface_process();

		evaluation_process();
		refine_process();

		// ���綼�Ѿ����������ж�����, �Ǿ���ֹ, �õ�һ����С��������. 2009-07-07
		if (simplified_mesh_.n_vertices() == mesh2_.n_vertices()) fitting_terminate_ = true;
	}
	std::cout << "Info: Fitting Terminated.\n";
	
	return true;
}