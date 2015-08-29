#include "DecimationModel.h" 

// ---------------
void DecimationModel::find_closestp_icm() {
	std::cout << "����find_closestp_icm().\n";
	//����limit surface��m��sample points, ��original surface����m����Ӧ��closest points��Ϊfoot point
	// ������Ϊ�˱�֤ refined_simplified_mesh_��initial_control_mesh_������ͬ��index����, ����Է����Ժ�ļ���.
	int n = simplified_mesh_.n_vertices();
	int m = initial_control_mesh_.n_vertices(); // 
	if (m - n != simplified_mesh_.n_edges()) std::cout << "Error: m-n == simplified_mesh_.n_edges().\n";
	if (refined_simplified_mesh_.n_vertices() != m) std::cout << "Error: refined simplified mesh��limit surface�Ķ�����Ӧ����ͬ.\n";
	//std::cout << "1,"; // for test
	for (int i = 0; i < n; ++i) {
		TriMesh::VertexHandle vh = initial_control_mesh_.vertex_handle(i);//limit surface��ǰn������.
		if (vh != simplified_mesh_.vertex_handle(i) || vh != refined_simplified_mesh_.vertex_handle(i)) std::cout << "Error: ǰn������Ӧ����ƥ���.\n";
		TriMesh::VertexHandle footpoint = voldhandles_[vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, vh)).point_handle()]; 
		double distance = (initial_control_mesh_.point(vh) - mesh_.point(footpoint)).norm();
		//std::cout << "distance " << distance << ", ";
		bool stop = false;
		while (stop == false) {
			stop = true;
			for (TriMesh::VVIter vv_it(mesh2_, footpoint); vv_it; ++vv_it) {
				if ((mesh2_.point(vv_it) - initial_control_mesh_.point(vh)).norm() < distance) {
					footpoint = vv_it.handle(); 
					distance = (mesh2_.point(vv_it) - initial_control_mesh_.point(vh)).norm(); //std::cout << "distance " << distance << ",\t";
					stop = false;
				}
			}
		}
		initial_control_mesh_.property(vp_closestp_icm_, vh) = footpoint;
		initial_control_mesh_.property(vp_distances_, vh) = distance;
	} // i = n
	// ���濼����
	for (int i = 0; i < m -n; ++i) { //simplified_mesh_.edge_handle(i) ��Ӧrefined_simplified_mesh_.vertex_handle(n + i)
		TriMesh::VertexHandle vh = initial_control_mesh_.vertex_handle(n + i); //vertex: n ... m-1. һ��ϸ��ʱ��m-n���߼�����µ�. 
		if (vh != refined_simplified_mesh_.vertex_handle(n + i)) std::cout << "Error: n...m-1����Ӧ��ƥ��ŶԵ�.\n";
		TriMesh::HalfedgeHandle h0 = simplified_mesh_.halfedge_handle(simplified_mesh_.edge_handle(i), 0);
		TriMesh::VHandle v0 = simplified_mesh_.to_vertex_handle(h0);
		TriMesh::VHandle v1 = simplified_mesh_.from_vertex_handle(h0); 
		int p0 = vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v0)).point_handle();
		int p1 = vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v1)).point_handle();
		TriMesh::VHandle old_v0 = voldhandles_[p0], old_v1 = voldhandles_[p1]; 
		TriMesh::EHandle old_e = mesh_.edge_handle(mesh_.find_halfedge(old_v0, old_v1)); 
		if (mesh_.property(empl_, old_e) == false ) std::cout << "Error: simplified_mesh_��Ӧ��mesh_�ϻ��б�û���е����.\n";
		TriMesh::VertexHandle footpoint = mesh_.property(emp_closest_vh_, old_e); // -1
		//std::cout << "8," << footpoint << ", ";

		double distance = (initial_control_mesh_.point(vh) - mesh_.point(footpoint)).norm();//std::cout << "initial distance " << distance << ",\t";
		bool stop = false;
		while (stop==false) {
			stop = true;
			for (TriMesh::VVIter vv_it(mesh2_, footpoint); vv_it; ++vv_it) {
				if ((mesh2_.point(vv_it) - initial_control_mesh_.point(vh)).norm() < distance) {
					footpoint = vv_it.handle(); 
					distance = (mesh2_.point(vv_it) - initial_control_mesh_.point(vh)).norm(); //std::cout << "distance " << distance << ",\t";
					stop = false;
				}
			}
		}
		//std::cout << "7,";
		initial_control_mesh_.property(vp_closestp_icm_, vh) = footpoint;
		initial_control_mesh_.property(vp_distances_, vh) = distance;		
	}
	std::cout << "�뿪find_closestp_icm().\n";
}
bool DecimationModel::test_is_leaf_of_vertex_hierarchy(TriMesh::VHandle _vh_in_simplified_mesh) {
	// �����ɷ���, �Ѿ���ԭ�Ӷ���(��Ӧ��vertex hierarchy�ϵ�Ҷ�ӽڵ�)ʱ����true.
	TriMesh::VHandle v_to = _vh_in_simplified_mesh;  
	int v_to_node_handle = simplified_mesh_.property(vp_node_handle_sm_, v_to);
	//TriMesh::VertexHandle v_to_old = voldhandles_[vhierarchy_.node(v_to_node_handle).point_handle()];  
	// for test, ע��, ����ʱ��ע����v_to_node_handle ������-1��, -1��ʾ�������vertex handle��ÿһ��������Ӧ��Node����.
	// �����������Ӧ��Node����Ҫactive��
	if (vhierarchy_.node(v_to_node_handle).is_active() == false ) { std::cout << "Error: Testing: active node�ſ���split.\n"; }
	int lc_node_handle = vhierarchy_.node(v_to_node_handle).lchild_node_handle();
	//int rc_node_handle = vhierarchy_.node(v_to_node_handle).rchild_node_handle();
	if (lc_node_handle != -1) { // rc_node_handle != -1 //����ڵ㲻��Ҷ�ӽڵ�,��������split
		return false;
	} else {
		std::cout << "Error: Testing: ����Ҷ�ӽڵ�, not alow to split.\n";
		return true;//��Ҷ�ӽڵ�.
	}
}
bool DecimationModel::test_is_triangle_flipping_after_vsplit(TriMesh::VHandle _vh_in_simplified_mesh) {
	// �����ɷ���, ����,���Ѻ���γ�triangle flipping, ����true.
	TriMesh::VHandle v_to = _vh_in_simplified_mesh;  
	int v_to_node_handle = simplified_mesh_.property(vp_node_handle_sm_, v_to);
	// for test, ע��, ����ʱ��ע����v_to_node_handle ������-1��, -1��ʾ�������vertex handle��ÿһ��������Ӧ��Node����.
	// �����������Ӧ��Node����Ҫactive��
	if (vhierarchy_.node(v_to_node_handle).is_active() == false ) { std::cout << "Error: Testing: active node�ſ���split.\n"; }
	int lc_node_handle = vhierarchy_.node(v_to_node_handle).lchild_node_handle();
	int rc_node_handle = vhierarchy_.node(v_to_node_handle).rchild_node_handle();
	if (lc_node_handle == -1) {
		std::cout << "Error: test_is_triangle_flipping_after_vsplit()���ﲻ����Ҷ�ӷ��ѵ�.\n";//ǰ��Ӧ���Ѿ�������Թ����˰�.
		return true;
	}

	// lc_node_handle != -1 // rc_node_handle != -1 //����ڵ㲻��Ҷ�ӽڵ�,��������split
	//vsplit(to)->(from, to), Ҳ����vspit(v_to)->(v0, v_to). vsplit(v1)->(v0, v1)

	TriMesh::Point v_from_point = vpoints_[vhierarchy_.node(lc_node_handle).point_handle()]; 
	TriMesh::VHandle v_from = simplified_mesh_.add_vertex(v_from_point);//���������from��
	TriMesh::VHandle v_l = TriMesh::VertexHandle(-1);//��ʾinvalid vertex handle
	TriMesh::VHandle v_r = TriMesh::VertexHandle(-1);//�����ǰsplitʱ���active left/right vertex
	int fc0 = vhierarchy_.node(v_to_node_handle).fund_cut_node_handle0();//std::cout << "fc0 n h: " << fc0 << "\n";//for test, fundamental cut vertex vl^.
	if (fc0 != -1) { // fc0 == -1��ʾv_l = -1, �Ǳ߽����
		DGP::VHierarchyNode fund_cut_node = vhierarchy_.node(fc0);
		while (fund_cut_node.is_active() == false) {  //std::cout << ": " << voldhandles_[fund_cut_node.point_handle()];
			fund_cut_node = vhierarchy_.node(fund_cut_node.parent_node_handle());
		}
		v_l = fund_cut_node.vertex_handle(); //std::cout << ": " << voldhandles_[fund_cut_node.point_handle()];
	}
	int fc1 = vhierarchy_.node(v_to_node_handle).fund_cut_node_handle1();//std::cout << "fc1 n h: " << fc1 << "\n";//for test, fundamental cut vertex vr^.
	if (fc1 != -1) {
		DGP::VHierarchyNode fund_cut_node = vhierarchy_.node(fc1); //std::cout << fund_cut_node.self_node_handle() << "\n";
		while (fund_cut_node.is_active() == false) {  //std::cout << ":: " << voldhandles_[fund_cut_node.point_handle()];
			fund_cut_node = vhierarchy_.node(fund_cut_node.parent_node_handle());
		}  
		v_r = fund_cut_node.vertex_handle();  //std::cout << ":: " << voldhandles_[fund_cut_node.point_handle()];
	}
	// �������ҳ�simplified_mesh_���Ѷ���ʱ���������Ϣ.

	simplified_mesh_.vertex_split(v_from, v_to, v_l, v_r);//std::cout << "0k.\n";
	//std::cout << "vplit simplified_mesh_: " << v_from << ", " << v_to << ", " << v_l << ", " << v_r << "; funt cut ver: " << fc0 << ", " << fc1 << ".\n";//for test//

	// �����Ƿ�����vsplit֮���triangle filpping.
	bool triangle_flipped = false;
	for (TriMesh::VertexFaceIter vf_it(simplified_mesh_, v_from); vf_it; ++vf_it) {
		simplified_mesh_.set_normal(vf_it.handle(), simplified_mesh_.calc_face_normal(vf_it.handle()));
	}
	for (TriMesh::VertexFaceIter vf_it(simplified_mesh_, v_from); vf_it; ++vf_it) {
		for (TriMesh::VertexFaceIter vff_it(simplified_mesh_, v_from); vff_it; ++vff_it) {
			TriMesh::Normal f0 = (simplified_mesh_.normal(vf_it));
			TriMesh::Normal fi = (simplified_mesh_.normal(vff_it));
			if (dot(f0, fi)/(f0.norm() * fi.norm()) < -0.866) {
				std::cout << "Error: Testing����, triangle flipped after this vsplit. v_from's valence=" << simplified_mesh_.valence(v_from) << ".\n";
				triangle_flipped = true; break;
			}
		}
		if (triangle_flipped) break;
	} std::cout << "dadk ";
	simplified_mesh_.collapse(simplified_mesh_.find_halfedge(v_from, v_to)); std::cout << "xx ";
	simplified_mesh_.garbage_collection(); std::cout << "aa ";//һ��Ҫ��һ������û��������ȥ�ǵ�.
	simplified_mesh_.update_normals();
	if (triangle_flipped) {
		return true;
	} else {
		// ��������γ������η�ת, �ͷ���false.
		return false;
	}	
}
bool DecimationModel::vsplit_refine(TriMesh::VertexHandle v_to) {
	int v_to_node_handle = simplified_mesh_.property(vp_node_handle_sm_, v_to);
	TriMesh::VertexHandle v_to_old = voldhandles_[vhierarchy_.node(v_to_node_handle).point_handle()];
	std::cout << "Info:��Ҫ���ѵ���simplified_mesh_.v_to: " << v_to /*<< ", v_to_node_handle: " << v_to_node_handle */
		<< ". mesh_.v_to_old: " << v_to_old << ".\n";
	test_vertex_ = v_to_old;//
	// for test, ע��, ����ʱ��ע����v_to_node_handle ������-1��, -1��ʾ�������vertex handle��ÿһ��������Ӧ��Node����.
	// �����������Ӧ��Node����Ҫactive��
	if (vhierarchy_.node(v_to_node_handle).is_active() == false ) { std::cout << "Error: active node�ſ���split.\n"; }
	int lc_node_handle = vhierarchy_.node(v_to_node_handle).lchild_node_handle();
	int rc_node_handle = vhierarchy_.node(v_to_node_handle).rchild_node_handle();
	if (lc_node_handle != -1) { // rc_node_handle != -1 //����ڵ㲻��Ҷ�ӽڵ�,��������split
		//vsplit(to)->(from, to), Ҳ����vspit(v_to)->(v0, v_to). vsplit(v1)->(v0, v1)
		if (v_to_old != voldhandles_[vhierarchy_.node(rc_node_handle).point_handle()]) std::cout << "Error: ������Ӧ����ָ��ͬһ���ɶ����.\n";
		TriMesh::VHandle v_from_old = voldhandles_[vhierarchy_.node(lc_node_handle).point_handle()];
		test_vertex1_ = v_from_old;
		std::cout << "mesh_.vsplit: to->(from, to) = (" << v_from_old << ", " << v_to_old << "), type: (" 
			<< mesh_.property(vp_type_, v_from_old) << ", " << mesh_.property(vp_type_, v_to_old) << ").\n";//for test

		TriMesh::Point v_from_point = vpoints_[vhierarchy_.node(lc_node_handle).point_handle()]; 
		TriMesh::VHandle v_from = simplified_mesh_.add_vertex(v_from_point);//���������from��
		TriMesh::VHandle v_l = TriMesh::VertexHandle(-1);//��ʾinvalid vertex handle
		TriMesh::VHandle v_r = TriMesh::VertexHandle(-1);//�����ǰsplitʱ���active left/right vertex
		int fc0 = vhierarchy_.node(v_to_node_handle).fund_cut_node_handle0();//std::cout << "fc0 n h: " << fc0 << "\n";//for test, fundamental cut vertex vl^.
		if (fc0 != -1) { // fc0 == -1��ʾv_l = -1, �Ǳ߽����
			DGP::VHierarchyNode fund_cut_node = vhierarchy_.node(fc0);
			while (fund_cut_node.is_active() == false) {  //std::cout << ": " << voldhandles_[fund_cut_node.point_handle()];
				fund_cut_node = vhierarchy_.node(fund_cut_node.parent_node_handle());
			}
			v_l = fund_cut_node.vertex_handle(); //std::cout << ": " << voldhandles_[fund_cut_node.point_handle()];
		}
		int fc1 = vhierarchy_.node(v_to_node_handle).fund_cut_node_handle1();//std::cout << "fc1 n h: " << fc1 << "\n";//for test, fundamental cut vertex vr^.
		if (fc1 != -1) {
			DGP::VHierarchyNode fund_cut_node = vhierarchy_.node(fc1); //std::cout << fund_cut_node.self_node_handle() << "\n";
			while (fund_cut_node.is_active() == false) {  //std::cout << ":: " << voldhandles_[fund_cut_node.point_handle()];
				fund_cut_node = vhierarchy_.node(fund_cut_node.parent_node_handle());
			}  
			v_r = fund_cut_node.vertex_handle();  //std::cout << ":: " << voldhandles_[fund_cut_node.point_handle()];
		}
		TriMesh::VHandle v_l_old;
		if (v_l.is_valid()) {
			v_l_old = voldhandles_[vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v_l)).point_handle()];
			if (mesh_.status(v_l_old).deleted()) std::cout << "Error: v_l_old should be not deleted.\n";
		}
		TriMesh::VHandle v_r_old;
		if (v_r.is_valid()) {
			v_r_old = voldhandles_[vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v_r)).point_handle()];
			if (mesh_.status(v_r_old).deleted()) std::cout << "Error: v_r_old should be not deleted.\n";
		}
		std::cout << "mesh_.(vl, vr) = (" << v_l_old << ", " << v_r_old << "); type=(" << mesh_.property(vp_type_, v_l_old) << ", " << mesh_.property(vp_type_, v_r_old) << ").\n";
		// �������ҳ�simplified_mesh_��mesh_���Ѷ���ʱ���������Ϣ.

		simplified_mesh_.vertex_split(v_from, v_to, v_l, v_r);//std::cout << "0k.\n";
		std::cout << "vplit simplified_mesh_: " << v_from << ", " << v_to << ", " << v_l << ", " << v_r << "; funt cut ver: " << fc0 << ", " << fc1 << ".\n";//for test//

		// �����Ƿ�����vsplit֮���triangle filpping.
		bool triangle_flipped = false;
		for (TriMesh::VertexFaceIter vf_it(simplified_mesh_, v_from); vf_it; ++vf_it) {
			simplified_mesh_.set_normal(vf_it.handle(), simplified_mesh_.calc_face_normal(vf_it.handle()));
		}
		for (TriMesh::VertexFaceIter vf_it(simplified_mesh_, v_from); vf_it; ++vf_it) {
			for (TriMesh::VertexFaceIter vff_it(simplified_mesh_, v_from); vff_it; ++vff_it) {
				TriMesh::Normal f0 = (simplified_mesh_.normal(vf_it));
				TriMesh::Normal fi = (simplified_mesh_.normal(vff_it));
				if (dot(f0, fi)/(f0.norm() * fi.norm()) < -0.866) {
					std::cout << "Error: ����, triangle flipped after this vsplit. v_from's valence=" << simplified_mesh_.valence(v_from) << ".\n";
					triangle_flipped = true; break;
				}
			}
			if (triangle_flipped) break;
		}
		if (triangle_flipped) {
			//simplified_mesh_.collapse(simplified_mesh_.find_halfedge(v_from, v_to)); 
			//return false;
		} /**/
		// ���û�з�ת, ����
		vhierarchy_.node(v_to_node_handle).set_active(false);
		vhierarchy_.node(lc_node_handle).set_vertex_handle(v_from);	simplified_mesh_.property(vp_node_handle_sm_, v_from) = lc_node_handle;
		vhierarchy_.node(lc_node_handle).set_active(true);
		vhierarchy_.node(rc_node_handle).set_vertex_handle(v_to);	simplified_mesh_.property(vp_node_handle_sm_, v_to) = rc_node_handle;
		vhierarchy_.node(rc_node_handle).set_active(true);	//

		// simplified_mesh_����֮��, mesh_Ҳ����Ӧ�ķ���.
		if (mesh_.status(v_from_old).deleted() == false) std::cout << "Error: v_from_old should be deleted.\n";// ��ʱ��v_from_old����֮ǰɾȥ�ĵ�.
		//std::cout << "v_from_old, vbc_: " << mesh_.property(vbc_, v_from_old) << ". point: " << mesh_.point(v_from_old) << ".\n";
		if (mesh_.status(v_to_old).deleted()) std::cout << "Error: v_to_old should be not deleted.\n";
		mesh_.vertex_split(v_from_old, v_to_old, v_l_old, v_r_old);
		std::cout << "vplit mesh_: " << v_from_old << ", " << v_to_old << ", " << v_l_old << ", " << v_r_old << ".\n";//for test//

		mesh_.status(v_from_old).set_deleted(false); //���¼���v_from_old����.
		//if (mesh_.status(v_from_old).deleted() == true) std::cout << "Error: v_from_old should be not deleted, now.\n";
		// v_from ����one-ring�ı߶���Ҫ�������е�, ��Ҫ����mesh_�²���
		// ������һ��Ҫ������: �������������߶���û�б仯��, Ҳ����Ҫ���²���.
		for (TriMesh::VertexEdgeIter ve_it(mesh_, v_from_old); ve_it; ++ve_it) { //��Щ���������ɵı�.
			mesh_.property(empl_, ve_it) = false; 
		}
		midpoint_sample_succeed_ = false; // v_from_old������, �������ڱ߶���û���е������.

		if (mesh_.is_boundary(v_from_old) == false) { // v_to_old -> (v_from_old, v_to_old) �ڲ����
			// �����ɵĶ���v_from��v_from_old�Ķ���������ʵ��Ӧ�����õ�.
			if (mesh_.property(vp_type_, v_to_old) != simplified_mesh_.property(vp_type_sm_, v_to)) std::cout << "Error: v_to's should be = v_to_old's type.\n";

			TriMesh::HalfedgeHandle cehe0, cehe1;
			// ��v_to/v_to_old�ĵ�����Ϊdart, crease, cornerʱ��, ��v_from/v_from_old����crease��, ��Ϊdart/corner�����Ϊfrom�������������.
			if (mesh_.property(vp_type_, v_from_old) == DGP::CREASE_VFT) 
			{	//����v_from_old������crease��, ��ôһͷһβ�Ķ���crease edges�Ŷ�
				std::cout << "����������: v_to_old's vp_type_: " << mesh_.property(vp_type_, v_to_old) << ", v_from_old's vp_type_: " << mesh_.property(vp_type_, v_from_old) << ".\n";
				//����simplified_mesh_�ϵĵ�ͱߵ���������.
				int n_crease_edges = 0; //�����ȴ���һ�������.
				for (TriMesh::VertexEdgeIter ve_it(simplified_mesh_, v_from); ve_it; ++ve_it) {
					if (simplified_mesh_.status(ve_it.handle()).feature())  ++ n_crease_edges;
				}
				if (n_crease_edges != 1) { //������Ӧ����һ���ߵ�,����(v_feature and v_from), ������ǵĻ��ͳ����������账����������
					if (mesh_.status(mesh_.edge_handle(mesh_.find_halfedge(v_l_old, v_to_old))).feature() == false 
						&& mesh_.status(mesh_.edge_handle(mesh_.find_halfedge(v_r_old, v_to_old))).feature() == false) { //���ǲ�����ֵ������
							std::cout << "Error: ����Ӧ����ֻ����1�������߲ŶԵ�.\n";
					} else if (mesh_.status(mesh_.edge_handle(mesh_.find_halfedge(v_l_old, v_to_old))).feature()) {
						simplified_mesh_.status(simplified_mesh_.edge_handle(simplified_mesh_.find_halfedge(v_from, v_l))).set_feature(true);
						mesh_.status(mesh_.edge_handle(mesh_.find_halfedge(v_from_old, v_l_old))).set_feature(true);
						simplified_mesh_.status(simplified_mesh_.edge_handle(simplified_mesh_.find_halfedge(v_l, v_to))).set_feature(false);
						mesh_.status(mesh_.edge_handle(mesh_.find_halfedge(v_l_old, v_to_old))).set_feature(false);

						for (int i = 0; i < mesh_.property(hep_heset_, mesh_.find_halfedge(v_to_old, v_l_old)).size(); ++i) {
							mesh_.property(hep_heset_, mesh_.find_halfedge(v_from_old, v_l_old)).push_back(mesh_.property(hep_heset_, mesh_.find_halfedge(v_to_old, v_l_old))[i]);
							TriMesh::VHandle vh = mesh2_.to_vertex_handle(mesh_.property(hep_heset_, mesh_.find_halfedge(v_to_old, v_l_old))[i]);
							if (vh != v_l_old) mesh_.property(vp_feature_edge_, vh) = mesh_.edge_handle(mesh_.find_halfedge(v_from_old, v_l_old));
						} 
						mesh_.property(hep_heset_, mesh_.find_halfedge(v_to_old, v_l_old)).clear();
						for (int i = 0; i < mesh_.property(hep_heset_, mesh_.find_halfedge(v_l_old, v_to_old)).size(); ++i) {
							mesh_.property(hep_heset_, mesh_.find_halfedge(v_l_old, v_from_old)).push_back(mesh_.property(hep_heset_, mesh_.find_halfedge(v_l_old, v_to_old))[i]);
						}
						mesh_.property(hep_heset_, mesh_.find_halfedge(v_l_old, v_to_old)).clear();
					} else if (mesh_.status(mesh_.edge_handle(mesh_.find_halfedge(v_r_old, v_to_old))).feature()) {
						simplified_mesh_.status(simplified_mesh_.edge_handle(simplified_mesh_.find_halfedge(v_from, v_r))).set_feature(true);
						mesh_.status(mesh_.edge_handle(mesh_.find_halfedge(v_from_old, v_r_old))).set_feature(true);
						simplified_mesh_.status(simplified_mesh_.edge_handle(simplified_mesh_.find_halfedge(v_r, v_to))).set_feature(false);
						mesh_.status(mesh_.edge_handle(mesh_.find_halfedge(v_r_old, v_to_old))).set_feature(false);

						for (int i = 0; i < mesh_.property(hep_heset_, mesh_.find_halfedge(v_to_old, v_r_old)).size(); ++i) {
							mesh_.property(hep_heset_, mesh_.find_halfedge(v_from_old, v_r_old)).push_back(mesh_.property(hep_heset_, mesh_.find_halfedge(v_to_old, v_r_old))[i]);
							TriMesh::VHandle vh = mesh2_.to_vertex_handle(mesh_.property(hep_heset_, mesh_.find_halfedge(v_to_old, v_r_old))[i]);
							if (vh != v_r_old) mesh_.property(vp_feature_edge_, vh) = mesh_.edge_handle(mesh_.find_halfedge(v_from_old, v_r_old));
						} 
						mesh_.property(hep_heset_, mesh_.find_halfedge(v_to_old, v_r_old)).clear();
						for (int i = 0; i < mesh_.property(hep_heset_, mesh_.find_halfedge(v_r_old, v_to_old)).size(); ++i) {
							mesh_.property(hep_heset_, mesh_.find_halfedge(v_r_old, v_from_old)).push_back(mesh_.property(hep_heset_, mesh_.find_halfedge(v_r_old, v_to_old))[i]);
						}
						mesh_.property(hep_heset_, mesh_.find_halfedge(v_r_old, v_to_old)).clear();
					}						
				} /*
				  int n_crease_vs = 0; // for test
				  for (TriMesh::VertexVertexIter vv_it(simplified_mesh_, v_from); vv_it; ++vv_it) {
				  if (simplified_mesh_.property(vp_type_sm_, vv_it.handle()) == DGP::CREASE_VFT) {
				  ++ n_crease_vs;
				  std::cout << " CV " << voldhandles_[vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, vv_it.handle())).point_handle()];
				  } else std::cout << " SV " << voldhandles_[vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, vv_it.handle())).point_handle()];
				  }
				  std::cout << "n_crease_vs: " << n_crease_vs << "; " << simplified_mesh_.valence(v_from) << ".\n";
				  std::cout << "n_crease_edges: " << n_crease_edges << ".\n";*/

				// ����v_from��crease����.
				simplified_mesh_.status(simplified_mesh_.edge_handle(simplified_mesh_.find_halfedge(v_from, v_to))).set_feature(true);
				n_crease_edges = 2; // ������simplified_mesh_��Ҳ������crease edges��.��Ϊ��ʱ��������crease edge�ϳ�һ����.
				simplified_mesh_.property(vp_type_sm_, v_from) = DGP::CREASE_VFT;

				// ��Ȼ������������, �������������߾Ϳ����ȷ����������������������㲢���е���.
				for (TriMesh::VertexOHalfedgeIter voh_it(mesh_, v_from_old); voh_it; ++voh_it) {
					if (mesh_.status(mesh_.edge_handle(voh_it)).feature()) {
						cehe1 = voh_it; break;
					}
				} // ���ҳ�cehe1�����ǰ�������������ø������ɵ�������. 
				cehe0 = mesh_.find_halfedge(v_from_old, v_to_old);
				mesh_.status(mesh_.edge_handle(cehe0)).set_feature(true);				
				TriMesh::HHandle ocehe0 = mesh_.opposite_halfedge_handle(cehe0), ocehe1 = mesh_.opposite_halfedge_handle(cehe1); //std::cout << "here1.\n";
				int hep_heset_size = mesh_.property(hep_heset_, ocehe1).size(); //std::cout << "here2.\n"; // for test.
				if (mesh_.property(hep_heset_, cehe1).size() != hep_heset_size) std::cout << "Error: ocehe1��cehe1��hep_heset_��СӦ����һ���ŶԵ�(�����е�����).\n";
				// ��ʱcehe1 and ocehe1���ϵ���������ԭ����û�з���ʱ������е�, ����ſ�ʼ��Щ������ķ���.
				//std::cout << "hep_heset_size " << hep_heset_size << ".\n";
				int i = 0; // ����ocehe1�ϵ�hep_heset_  
				std::vector<TriMesh::HalfedgeHandle> tmp;
				for (; i < mesh_.property(hep_heset_, ocehe1).size(); ++i) {
					if (mesh2_.to_vertex_handle(mesh_.property(hep_heset_, ocehe1)[i]) != v_from_old) {
						tmp.push_back(mesh_.property(hep_heset_, ocehe1)[i]);
					} else break;
				} 
				tmp.push_back(mesh_.property(hep_heset_, ocehe1)[i]);// 0 - i ԭ����
				for (++i; i < mesh_.property(hep_heset_, ocehe1).size(); ++i) {
					mesh_.property(hep_heset_, cehe0).push_back(mesh_.property(hep_heset_, ocehe1)[i]);
				} 
				mesh_.property(hep_heset_, ocehe1).clear();
				for (std::vector<TriMesh::HHandle>::iterator it(tmp.begin()), end(tmp.end()); it != end; ++it) {
					mesh_.property(hep_heset_, ocehe1).push_back(*it);
				}
				if ((mesh_.property(hep_heset_, ocehe1).size() + mesh_.property(hep_heset_, cehe0).size()) != hep_heset_size) std::cout << "Error: cehe0, ocehe1��hep_heset_�������.\n";
				//std::cout << mesh_.property(hep_heset_, ocehe1).size() << ", " <<  mesh_.property(hep_heset_, cehe0).size() << ".\n";
				//����_hh���ϵĵ�Ĳ�����
				for (int i = 0; i < mesh_.property(hep_heset_, cehe0).size() -1; ++i) {
					mesh_.property(vp_feature_edge_, mesh2_.to_vertex_handle(mesh_.property(hep_heset_, cehe0)[i])) = mesh_.edge_handle(cehe0);
				}
				for (int i = 0; i < mesh_.property(hep_heset_, ocehe1).size() -1; ++i) {
					mesh_.property(vp_feature_edge_, mesh2_.to_vertex_handle(mesh_.property(hep_heset_, ocehe1)[i])) = mesh_.edge_handle(ocehe1);
				}

				i = 0; // ����cehe1�ϵ�hep_heset_ //std::cout << "cehe1's ��ʼsize: " << mesh_.property(hep_heset_, cehe1).size() << ".\n";
				for (; i < hep_heset_size; ++i) {
					if (mesh2_.to_vertex_handle(mesh_.property(hep_heset_, cehe1)[i]) != v_from_old) {
						mesh_.property(hep_heset_, ocehe0).push_back(mesh_.property(hep_heset_, cehe1)[i]);
					} else break;
				}
				mesh_.property(hep_heset_, ocehe0).push_back(mesh_.property(hep_heset_, cehe1)[i]); //std::cout << "i: " << i << ", ocehe0.size " << mesh_.property(hep_heset_, ocehe0).size() << ".\n";
				tmp.clear();
				for (++i; i < hep_heset_size; ++i) {
					tmp.push_back(mesh_.property(hep_heset_, cehe1)[i]);
				}
				mesh_.property(hep_heset_, cehe1).clear(); //std::cout << "tmp's size " << tmp.size() << ".\n";
				for (std::vector<TriMesh::HHandle>::iterator it(tmp.begin()), end(tmp.end()); it != end; ++it) {
					mesh_.property(hep_heset_, cehe1).push_back(*it);
				} //std::cout << "cehe1's ����size " << mesh_.property(hep_heset_, cehe1).size() << ".\n";
				if ((mesh_.property(hep_heset_, ocehe0).size() + mesh_.property(hep_heset_, cehe1).size()) != hep_heset_size) std::cout << "Error: ocehe0, cehe1��hep_heset_�������.\n";
				for (int i = 0; i < mesh_.property(hep_heset_, ocehe0).size() -1; ++i) {
					mesh_.property(vp_feature_edge_, mesh2_.to_vertex_handle(mesh_.property(hep_heset_, ocehe0)[i])) = mesh_.edge_handle(ocehe0);
				}
				for (int i = 0; i < mesh_.property(hep_heset_, cehe1).size() -1; ++i) {
					mesh_.property(vp_feature_edge_, mesh2_.to_vertex_handle(mesh_.property(hep_heset_, cehe1)[i])) = mesh_.edge_handle(cehe1);
				}
				if ( mesh_.property(hep_heset_, cehe0).size() != mesh_.property(hep_heset_, ocehe0).size()) std::cout << "Error: cehe0, ocehe0��hep_heset_Ӧ�õȴ�Ŷ�.\n";
				if ( mesh_.property(hep_heset_, cehe1).size() != mesh_.property(hep_heset_, ocehe1).size()) std::cout << "Error: cehe1, ocehe1��hep_heset_Ӧ�õȴ�Ŷ�.\n";
				//std::cout << " " << mesh_.property(hep_heset_, cehe0).size() << ", " << mesh_.property(hep_heset_, cehe1).size() << ".\n";
				if (mesh2_.to_vertex_handle(mesh_.property(hep_heset_, cehe0)[mesh_.property(hep_heset_, cehe0).size()-1]) != v_to_old) { std::cout << "Error: cehe0.\n";}
				if (mesh2_.to_vertex_handle(mesh_.property(hep_heset_, ocehe0)[mesh_.property(hep_heset_, ocehe0).size()-1]) != v_from_old) { std::cout << "Error: ocehe0.\n";}
				if (mesh2_.to_vertex_handle(mesh_.property(hep_heset_, ocehe1)[mesh_.property(hep_heset_, ocehe1).size()-1]) != v_from_old) { std::cout << "Error: ocehe1.\n";}
				// ��������������߷���ʱ���������Ϣ�ط���.
				// -----------------------------------------------------------------------------------
				//�Ըշ��ѵ�������crease edge���е�
				std::vector<TriMesh::HHandle> cehes; cehes.push_back(cehe0); cehes.push_back(cehe1);
				for (int ii = 0; ii < 2; ++ii) { 
					TriMesh::HalfedgeHandle cehe = cehes[ii];
					double len = 0;
					for (int i = 0; i < mesh_.property(hep_heset_, cehe).size(); ++i) {
						TriMesh::HHandle h = mesh_.property(hep_heset_, cehe)[i];
						len += (mesh2_.point(mesh2_.from_vertex_handle(h)) - mesh2_.point(mesh2_.to_vertex_handle(h))).norm();
					}
					double half_len = len / 2.0, len_tmp = 0;
					for (int i = 0; i < mesh_.property(hep_heset_, cehe).size(); ++i) {// �������crease edge���е�
						TriMesh::HHandle h = mesh_.property(hep_heset_, cehe)[i];
						len_tmp += (mesh2_.point(mesh2_.from_vertex_handle(h)) - mesh2_.point(mesh2_.to_vertex_handle(h))).norm();
						if (len_tmp >= half_len) {
							mesh_.property(emp_, mesh_.edge_handle(cehe)) 
								= mesh2_.point(mesh2_.to_vertex_handle(h))
								+ (mesh2_.point(mesh2_.from_vertex_handle(h)) - mesh2_.point(mesh2_.to_vertex_handle(h))).normalize() * (len_tmp - half_len);
							mesh_.property(empl_, mesh_.edge_handle(cehe)) = true;
							mesh_.property(emp_closest_vh_, mesh_.edge_handle(cehe)) = mesh2_.to_vertex_handle(h);
							break;
						}
					}
				}				
			} else { 
				std::cout << "����smooth��: v_to_old's vp_type_: " << mesh_.property(vp_type_, v_to_old) << ", v_from_old's vp_type_: " << mesh_.property(vp_type_, v_from_old) << ".\n";
				simplified_mesh_.property(vp_type_sm_, v_from) = DGP::SMOOTH_VFT;
				if (mesh_.property(vp_type_, v_from_old) != DGP::SMOOTH_VFT) std::cout << "Error: mesh_'s v_from_old's type Ӧ���Ǳ������вŶ��.\n";
				//if (v_from_old.idx() == 2207 && v_to_old.idx() == 2195) {
				//	std::cout << mesh_.valence(TriMesh::VHandle(2004)) << " ?= " << simplified_mesh_.valence(v_from) <<  ".\n";
				//	for (TriMesh::VVIter vv_it(mesh_, TriMesh::VHandle(2004)); vv_it; ++vv_it) {
				//		std::cout << vv_it.handle() << ": " << mesh_.property(vp_type_, vv_it) << ".\t";
				//	} std::cout << ".\n";
				//}
			}// end of if-else v_from_old��crease��

			// ------------------------------------------------------------------------------------------------------
			// ��������˵�ķ����Լ��������͵��ж�, �������v_from_oldһ����������Ĳ�������Ϣ.

			// ɾȥv_from_old������ԭ���Ǹ��������ϵļ�¼
			TriMesh::FHandle vf_v_from_old = mesh_.property(vf_, v_from_old);
			std::vector<TriMesh::VHandle>::iterator to_be_earsed 
				= remove((mesh_.property(fvset_, vf_v_from_old)).begin(), (mesh_.property(fvset_, vf_v_from_old)).end(), v_from_old);
			(mesh_.property(fvset_, vf_v_from_old)).erase(to_be_earsed, (mesh_.property(fvset_, vf_v_from_old)).end());

			// ѹƽv_from_old��һ����, �������û�д������
			TriMesh::HalfedgeHandle h0 = mesh_.find_halfedge(v_from_old, v_to_old);
			if (mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h0)) != v_l_old) std::cout << "Error: ѹƽv_from_oldʱ, h0.\n";
			std::vector<TriMesh::HalfedgeHandle> loop_h;
			std::vector<TriMesh::VertexHandle> loop; 
			int feature_vertex_index = -1;
			TriMesh::HHandle hh = h0;
			do {
				loop_h.push_back(hh); 
				loop.push_back(mesh_.to_vertex_handle(hh));
				if (mesh_.status(mesh_.edge_handle(hh)).feature()) feature_vertex_index = loop_h.size() - 1;
				if (mesh_.status(mesh_.to_vertex_handle(hh))  .deleted()) std::cout << "Error: v_from_old��һ����㶼����deleted��.\n";
				hh = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(hh));
			} while (hh != h0);
			if (loop_h.size() != mesh_.valence(v_from_old)) std::cout << "Error: loop_h's size.\n";

			int n = loop.size();
			/* double cir_length = 0; // ������λԲ�Ĳ�������.
			for (int i = 0; i < n; ++i) {
			cir_length += (mesh2_.point(loop[i]) - mesh2_.point(loop[(i+1)%n])).norm();
			}
			double l = 0, angle = 0, wij = 0;
			double sum_wij = 0.0;
			OpenMesh::Vec2d sum_wh(0, 0, 0);
			for (int i = 0; i < n; ++i) { //����vi��one-ring�ϵ���n���ڽӵ�
			// fix the boundary/one-ring vertices, 
			angle = l / cir_length * (2.0*M_PI); //std::cout << "ang: " << angle << std::endl;
			mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(0.5*cos(angle)+0.5, 0.5*sin(angle)+0.5);
			//�����λԲ����(0.5, 0.5)ΪԲ��, �뾶��0.5. Բ��Ҫ�Ƕ���(0, 0)�Ļ�, ���籾����(0, 0.5)�������������Ϊ(e, 0.5),e��һ����С��С��ֵ.
			l += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();

			// ���vi��vj�����ϵ�ϵ��, ע������ʹ�õ���mean value coordiante, ������harmonic map�е�cotangent weight
			TriMesh::VertexHandle v0 = v_from_old; // = mesh_.from_vertex_handle(loop_h[i]);
			TriMesh::VertexHandle v1 = loop[i];// =mesh_.to_vertex_handle(loop_h[i]);
			TriMesh::VertexHandle v2 = loop[(i+1)%n];//= mesh_.to_vertex_handle(mesh_.next_halfedge_handle(loop_h[i]));//������������ȷ��ǰ���ǷǱ߽����
			TriMesh::VertexHandle v3 = loop[(i+n-1)%n];//mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(loop_h[i])));
			OpenMesh::Vec3d v1v0 = mesh2_.point(v1) - mesh2_.point(v0);	// ֮ǰ������ֵ�-1
			OpenMesh::Vec3d v2v0 = mesh2_.point(v2) - mesh2_.point(v0);
			OpenMesh::Vec3d v3v0 = mesh2_.point(v3) - mesh2_.point(v0);
			double v1v0v2_angle = acos(dot(v1v0, v2v0) / (v1v0.norm() * v2v0.norm()));//ʸ��v1v0��v2v0�ļн�.
			double v3v0v1_angle = acos(dot(v3v0, v1v0) / (v3v0.norm() * v1v0.norm()));//ʸ��v3v0��v1v0�ļн�.
			wij = (tan(v1v0v2_angle / 2) + tan(v3v0v1_angle / 2)) / v1v0.norm();
			if (wij < 0) { std::cout << "Error: wij is negative.\n"; }

			sum_wij += wij;
			sum_wh += wij * mesh_.property(vuv_, loop[i]);
			}
			mesh_.property(vuv_, v_from_old) = sum_wh * (1.0/sum_wij);//�����vi�����ڵ�λԲ�ϵĲ���������.
			*/
			//����ѹ�ȷ���ʹ����Щ���λ�����������, ��Ϊ������������from�㲻������������ֱ����.
			double rou_angle = 0.0;  //round angle
			std::vector<double> vec_angle;
			for (int i=0 ; i<loop.size(); ++i) {//����ܵ�angle
				OpenMesh::Vec3d d1 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).normalize();
				OpenMesh::Vec3d d2 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).normalize();
				rou_angle += acos(dot(d1, d2));
				vec_angle.push_back(acos(dot(d1, d2)));
			}  //std::cout << "angle " << rou_angle << ".\n";//std::cout << "len: " << length << std::endl; //for test	
			if (vec_angle.size() != n) { std::cout << "Error: vec_angle��loop����������һ��.\n"; }
			double angle_scale_ratio = 2 * M_PI / rou_angle; //���ű���
			double temp_sum_angle = 0.0, l = 0;
			for (int i = 0; i < loop.size(); ++i) {
				temp_sum_angle += (vec_angle[i] * angle_scale_ratio);
				l = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).norm(); 
				l *= angle_scale_ratio;
				mesh_.property(vuv_, loop[(i+1)%n]) = OpenMesh::Vec2d(l*cos(temp_sum_angle), l*sin(temp_sum_angle));
			}
			//std::cout << feature_vertex_index ;//			
			/*// ��������ֻ����Ҫ�ڳ�ʼ������ʱ����Ҫ�������, ����һ���ѹƽӦ�ò���Ҫ��.
			if (feature_vertex_index != -1 && feature_vertex_index < n) {
			double len = 0;
			for (int i = 0; i < mesh_.property(hep_heset_, loop_h[0]).size(); ++i) {
			TriMesh::HHandle h = mesh_.property(hep_heset_, loop_h[0])[i];
			len += (mesh2_.point(mesh2_.from_vertex_handle(h)) - mesh2_.point(mesh2_.to_vertex_handle(h))).norm();
			}
			double tmp = len;
			for (int i = 0; i < mesh_.property(hep_heset_, loop_h[feature_vertex_index]).size(); ++i) {
			TriMesh::HHandle h = mesh_.property(hep_heset_, loop_h[feature_vertex_index])[i];
			len += (mesh2_.point(mesh2_.from_vertex_handle(h)) - mesh2_.point(mesh2_.to_vertex_handle(h))).norm();
			}
			tmp = tmp / len;
			mesh_.property(vuv_, v_from_old) = mesh_.property(vuv_, v_to_old) + (mesh_.property(vuv_, loop[feature_vertex_index]) - mesh_.property(vuv_, v_to_old)) * tmp;
			} else {*/
			mesh_.property(vuv_, v_from_old) = OpenMesh::Vec2d(0, 0); //���һ��������� 
			//}

			// ȡ��ԭ��1��n-2������Ĳ�������, ���������ѹƽ�õ����²�����, �����������������
			std::vector<TriMesh::VertexHandle> free_vertices;
			int free_vertices_size = 0;
			for (int i = 1; i <= n-2; ++i) { //loop_h[0 and n-1]ָ����������γɵ�.
				TriMesh::FaceHandle fh = mesh_.face_handle(mesh_.next_halfedge_handle(loop_h[i]));
				free_vertices_size += mesh_.property(fvset_, fh).size();
				//std::cout << mesh_.status(fh).deleted() << ", " << mesh_.property(fvset_, fh).size() << ".\n";

				for (std::vector<TriMesh::VHandle>::iterator it(mesh_.property(fvset_, fh).begin()), it_end(mesh_.property(fvset_, fh).end()); it != it_end; ++it) {
					free_vertices.push_back(*it);
					//if (mesh_.property(vp_type_, *it) != DGP::CREASE_VFT) {
					mesh_.property(vuv_, *it) = mesh_.property(vuv_, mesh_.property(vf0_, *it)) * (mesh_.property(vbc_, *it))[0]
					+ mesh_.property(vuv_, mesh_.property(vf1_, *it)) * (mesh_.property(vbc_, *it))[1]
					+ mesh_.property(vuv_,mesh_.property(vf2_, *it) ) * (mesh_.property(vbc_, *it))[2];
					//} else { // crease vertex��mesh_.property(vuv_, *it)ֵ�������ر�����.
					// ������ôif-else����,��Ϊ, һ�������߶�������������, ��Щ���ϵ������㶼���ܼ��뵽free_vertices������, ��������ȴû�и�����Щ�����Ĳ���������
					//} // end of if *it �Ƿ���crease vertex, �����*itҪô��smoothҪô����crease.
					// ����*it�Ƿ������㶼���뵽free_vertices�������Ա�����ͳһ�����䵽������.

				}
				mesh_.property(fvset_, fh).clear();
			}

			if (free_vertices.size() != free_vertices_size) std::cout << "Error: free_vertices's size.\n";
			// free_vertices������Ҳ��������Щ�����߷���ʱ���������, ��Щ������Ĳ�����������ر�����, ��ʵǰ��.
			if (mesh_.property(vp_type_, v_from_old) == DGP::CREASE_VFT) {
				if (!cehe0.is_valid() || !cehe1.is_valid()) std::cout << "Error: ��Ȼv_from_old��������, ��cehe0, cehe1�͸����������.\n";
				double len_cehe0 = 0; 
				for (int i = 0; i < mesh_.property(hep_heset_, cehe0).size(); ++i) {
					len_cehe0 += 
						(mesh2_.point(mesh2_.to_vertex_handle(mesh_.property(hep_heset_, cehe0)[i])) - mesh2_.point(mesh2_.from_vertex_handle(mesh_.property(hep_heset_, cehe0)[i]))).norm();
				} //std::cout << len_cehe0 << ".\n";
				double tmp = 0;
				for (int i = 0; i < mesh_.property(hep_heset_, cehe0).size(); ++i) {
					tmp += 
						(mesh2_.point(mesh2_.to_vertex_handle(mesh_.property(hep_heset_, cehe0)[i])) - mesh2_.point(mesh2_.from_vertex_handle(mesh_.property(hep_heset_, cehe0)[i]))).norm();
					mesh_.property(vuv_, mesh2_.to_vertex_handle(mesh_.property(hep_heset_, cehe0)[i])) 
						= mesh_.property(vuv_, v_from_old) + (mesh_.property(vuv_, v_to_old) - mesh_.property(vuv_, v_from_old)) * (tmp / len_cehe0);
				}
				double len_cehe1 = 0;
				for (int i = 0; i < mesh_.property(hep_heset_, cehe1).size(); ++i) {
					len_cehe1 += 
						(mesh2_.point(mesh2_.to_vertex_handle(mesh_.property(hep_heset_, cehe1)[i])) - mesh2_.point(mesh2_.from_vertex_handle(mesh_.property(hep_heset_, cehe1)[i]))).norm();
				} //std::cout << len_cehe1 << ".\n";
				tmp = 0;
				for (int i = 0; i < mesh_.property(hep_heset_, cehe1).size(); ++i) {
					tmp += 
						(mesh2_.point(mesh2_.to_vertex_handle(mesh_.property(hep_heset_, cehe1)[i])) - mesh2_.point(mesh2_.from_vertex_handle(mesh_.property(hep_heset_, cehe1)[i]))).norm();
					mesh_.property(vuv_, mesh2_.to_vertex_handle(mesh_.property(hep_heset_, cehe1)[i])) 
						= mesh_.property(vuv_, v_from_old) + (mesh_.property(vuv_, mesh_.to_vertex_handle(cehe1)) - mesh_.property(vuv_, v_from_old)) * (tmp / len_cehe1);
				}
			} 

			// �ж�free_vertices����ĵ������ĸ�������.
			OpenMesh::VPropHandleT<bool> is_face_located;
			mesh_.add_property(is_face_located);
			for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
				mesh_.property(is_face_located, *freev_it) = false;
			}
			for (int i = 0; i < n; ++ i) { // ѭ��n����, ��n�����Ǽ���to->(from, to)֮���γɵ�
				// �ж���������Ҫ���¶�λ��free vertices�Ƿ������triangle(v_from_old, loop[i], loop[i+1])֮��
				for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
					if (mesh_.property(is_face_located, *freev_it) == false) {

						OpenMesh::Vec2d pp = mesh_.property(vuv_, *freev_it);//���free vertex�Ĳ�������
						OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, v_from_old), mesh_.property(vuv_, loop[i]), mesh_.property(vuv_, loop[(i+1)%n]), pp);
						//std::cout << bc << std::endl;//for test

						if (DGP::is_valid_barycentric_coordinate(bc)) { //Ҫ�����bc��������������>=0

							mesh_.property(is_face_located, *freev_it) = true;

							// For interior vertex, its' barycentric coordinates is bc according to the triangle(0, i, i+1).
							//std::cout << *freev_it << ": " << v_from_old << "," << loop[i] << "," << loop[i+1] << ": " << bc << "\n"; // for test

							TriMesh::FaceHandle fh = mesh_.face_handle(loop_h[i]);
							if (fh.is_valid() == false) std::cerr << "Error: face invalid." << std::endl;//Only for assert

							// �������о���ÿһ�β�������Ӧ�û�õĻ��Ǹ��ĵ���Ϣ.
							mesh_.property(vf_,  *freev_it) = fh;
							mesh_.property(vf0_, *freev_it) = v_from_old; mesh_.property(vf1_, *freev_it) = loop[i]; mesh_.property(vf2_, *freev_it) = loop[(i+1)%n];
							mesh_.property(vbc_, *freev_it) = bc;

							(mesh_.property(fvset_, fh)).push_back(*freev_it);

							//��deletedҲ���Ǳ�parameterized�Ķ����3d����.
							mesh_collapsed_.set_point(*freev_it, mesh_.point(v_from_old) * bc[0] + mesh_.point(loop[i]) * bc[1] + mesh_.point(loop[(i+1)%n]) * bc[2]);////
						}
					} // ���free vertex�����������
				} // һ����free vertex �����������			
			} // ���е�free vertices ������n������
			// ��mesh_.property(vp_type_, v_from_old) == DGP::CREASE_VFTʱ��, free_vertices�������������������ߵ�(����)��Ĳ���������(bc0, 0, bc2)
			for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
				if (mesh_.property(is_face_located, *freev_it) == false) { std::cout << "Error: free_vertices ��û�ж�λ��.\n"; }
			}
			mesh_.remove_property(is_face_located);
			mesh_collapsed_.set_point(v_from_old, mesh_.point(v_from_old));
			free_vertices_size = 0; //����ȷ��һ���Ƿ������֤
			for (int i = 0; i < n; ++i) {
				free_vertices_size += mesh_.property(fvset_, mesh_.face_handle(loop_h[i])).size();
			}
			if (free_vertices_size != free_vertices.size()) std::cout << "Error: free_vertices's size �����¶�Ϊ�������ʹ���.\n";


			// �����γɵ�n���ߵ��е�, ��Щ�е��Ǳ����������, ����û������ȥ.
			flatten_vertex_one_ring_resample_edges_midpoint(v_from_old);
			int n_tmp_count = 0;//�������水��ѹƽһ����ķ��������ĳɹ���.
			for (int i = 0; i < n; ++i) {
				if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])))  ++ n_tmp_count;
			}
			std::cout << "Info: flatten vertex one ring resample " << n_tmp_count << "/" << mesh_.valence(v_from_old) << ".\n";

			for (int i = 0; i < n; ++i) {
				// ���ж������mesh_.edge_handle(loop_h[i])�Ƿ���ԭ�ӱ���, ���Ǻܿ��ܵ�.
				if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) { 
					TriMesh::VHandle v0 = mesh_.to_vertex_handle(loop_h[i]);
					for (TriMesh::VVIter vv_it(mesh2_, v_from_old); vv_it; ++ vv_it) {
						if (vv_it.handle() == v0) { 
							mesh_.property(emp_, mesh_.edge_handle(loop_h[i])) = (mesh_.point(v0) + mesh_.point(v_from_old)) / 2.0;
							mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) = true;
							mesh_.property(emp_closest_vh_, mesh_.edge_handle(loop_h[i])) = v0;
						}
					}
				}

				if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) { // ����߻�û���е����, ����Ĳ�����ֻ��ɹ���.
					flatten_face_resample_3edge_midpoint(mesh_.face_handle(loop_h[i]));
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) { 
						//flatten_face_resample_3edge_midpoint(mesh_.face_handle(loop_h[(i-1)%n]));
					}
					std::cout << "resample loop_h[" << i << "]: " << mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) << ".\n";//���0��ʾǰ��Ĳ�����û�гɹ�.

					TriMesh::FHandle fhi = mesh_.face_handle(loop_h[i]), fhi_1 = mesh_.face_handle(loop_h[(i+n-1)%n]);
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) {
						for (std::vector<TriMesh::VHandle>::iterator it(mesh_.property(fvset_, fhi).begin()), end(mesh_.property(fvset_, fhi).end()); it != end; ++it) {
							TriMesh::VertexHandle vh = *it; //��fhi�ϵĲ�������
							for (TriMesh::VertexOHalfedgeIter voh_it(mesh2_, vh); voh_it; ++voh_it) {
								TriMesh::VHandle v0 = vh, v1 = mesh2_.to_vertex_handle(voh_it);
								if (mesh_.property(vf_, v1) == fhi_1) {
									TriMesh::VHandle v2 = mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(voh_it));
									// ������ʵ������v2����ܲ�������������������.��������û�б�deleted�Ļ����ǲ��������������ϵ����.
									if (mesh_.property(vf_, v2) == fhi|| mesh_.property(vf_, v2) == fhi_1) {
										OpenMesh::Vec2d pp = (mesh_.property(vuv_, loop[i]) + mesh_.property(vuv_, v_from_old)) * 0.5;//����ߵ��е�Ĳ�������
										OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, v0), mesh_.property(vuv_, v1), mesh_.property(vuv_, v2), pp);
										//std::cout << bc << std::endl;//for test

										//if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //Ҫ�����bc��������������>=0
										if (DGP::is_valid_barycentric_coordinate(bc)) { //Ҫ�����bc��������������>=0
											mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) = true; std::cout << "���������㶨�˱��е����.\n";
											mesh_.property(emp_, mesh_.edge_handle(loop_h[i])) = mesh2_.point(v0) * bc[0] + mesh2_.point(v1) * bc[1] + mesh2_.point(v2) * bc[2];
											mesh_.property(emp_closest_vh_, mesh_.edge_handle(loop_h[i])) = v0;
										} 
									} else { // end of if . v2 is inside the valid region
										//std::cout << "Error: free vertex's neighbor v2 isn't in fhi-1.\n"; 
									}
								} else { // end of if.����
									//std::cout << "Error: free vertex's neighbor isn't in fhi-1.\n"; 
								}
							}
						} // end of for. ѭ��fhi�ϵ����в�������.
					}
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) { std::cout << "Waring: ����е����ֻ�����������.\n"; 
					use_closepoints(mesh_.edge_handle(loop_h[i]));
					}

					// ��ʵ���з�����ʵ������е����������������¶��Ѿ��ɹ���, ����ֻ��Ԥ������.
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) {
						// ����������ʲô�취��ҪŪ���е���, ��������
						std::cout << "�������� fhi.fvset_.size: " << mesh_.property(fvset_, fhi).size() << ".\n";
						std::vector<TriMesh::VertexHandle> feature_vers;
						if (mesh_.status(mesh_.edge_handle(loop_h[(i+1)%n])).feature()) {
							std::vector<TriMesh::HHandle> heset = mesh_.property(hep_heset_, loop_h[(i+1)%n]);
							for (int i = 0; i < heset.size(); ++i) {
								feature_vers.push_back(mesh2_.to_vertex_handle(heset[i]));
							}
						}  
						if (mesh_.status(mesh_.edge_handle(mesh_.next_halfedge_handle(loop_h[i]))).feature()) {
							std::vector<TriMesh::HHandle> heset = mesh_.property(hep_heset_, mesh_.next_halfedge_handle(loop_h[i]));
							for (int i = 0; i < heset.size(); ++i) {
								feature_vers.push_back(mesh2_.to_vertex_handle(heset[i]));
							}
						} std::cout << "feature_vers's size: " << feature_vers.size() << ". ";
						for (int j = 0; j < feature_vers.size(); ++j) { std::cout << "j " << j << ", "; //
						for (TriMesh::VOHIter voh_it(mesh2_, feature_vers[j]); voh_it; ++voh_it) {// && mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false
							TriMesh::VHandle v0 = feature_vers[j], v1 = mesh2_.to_vertex_handle(voh_it);
							if (mesh_.status(v1).deleted()) {
								if (mesh_.property(vf_, v1) == fhi_1) {
									std::cout << "a.\t";
									TriMesh::VHandle v2 = mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(voh_it));
									if (v2.is_valid() == false) std::cout << "Error: v2 isn't valid.\t";
									if (mesh_.property(vp_type_, v2) == DGP::CREASE_VFT) { std::cout << "aa\t";
									if (mesh_.property(vp_feature_edge_, v2) == mesh_.edge_handle(loop_h[(i+1)%n])
										|| mesh_.property(vp_feature_edge_, v2) == mesh_.edge_handle(mesh_.next_halfedge_handle(loop_h[i]))) {
											mesh_.property(vuv_, v2) = mesh_.property(vuv_, mesh_.property(vf0_, v2)) * (mesh_.property(vbc_, v2))[0]
											+ mesh_.property(vuv_, mesh_.property(vf1_, v2)) * (mesh_.property(vbc_, v2))[1]
											+ mesh_.property(vuv_, mesh_.property(vf0_, v2)) * (mesh_.property(vbc_, v2))[2];std::cout << "aaa\t";
									}
									} else { // ����Ӧ���������v2�Ĳ����������.
										std::cout << "v2 isn't crease, " << mesh_.property(vp_type_, v2) << ".\t";//v2's type == smooth
										mesh_.property(vuv_, v2) = mesh_.property(vuv_, mesh_.property(vf0_, v2)) * (mesh_.property(vbc_, v2))[0]
										+ mesh_.property(vuv_, mesh_.property(vf1_, v2)) * (mesh_.property(vbc_, v2))[1]
										+ mesh_.property(vuv_, mesh_.property(vf0_, v2)) * (mesh_.property(vbc_, v2))[2];std::cout << "bbb\t";
									}
									std::cout <<  mesh_.property(vuv_, v2) << ".\n";
									OpenMesh::Vec2d pp = (mesh_.property(vuv_, loop[i]) + mesh_.property(vuv_, v_from_old)) * 0.5;//����ߵ��е�Ĳ�������
									OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, v0), mesh_.property(vuv_, v1), mesh_.property(vuv_, v2), pp);
									std::cout << bc << ".\n";//for test

									if (fabs(bc[0]) < 1.0e-10) bc[0] = 0; if (fabs(bc[1]) < 1.0e-10) bc[1] = 0; if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;

									//if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //Ҫ�����bc��������������>=0
									if ((bc[0] >= 0) && (bc[1] >= 0) && (bc[2] >= 0)) { //Ҫ�����bc��������������>=0
										mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) = true; std::cout << "True." << mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) << ".\t";
										mesh_.property(emp_, mesh_.edge_handle(loop_h[i])) = mesh2_.point(v0) * bc[0] + mesh2_.point(v1) * bc[1] + mesh2_.point(v2) * bc[2];
										mesh_.property(emp_closest_vh_, mesh_.edge_handle(loop_h[i])) = v0;
									} 

								} 
							} else { // mesh2_.to_vertex_handle(voh_it) ���Ǳ�deleted�Ĳ�������,
								if (mesh2_.to_vertex_handle(voh_it) == v_from_old) {}
							} 
						} 
						} // end of for. feature_vers.size()
					} // end of if .������ʹ�õ��Ǳ�������.

					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) { std::cout << "Error: ����, �����µı�û������е����.\n"; }
				} // end of if�˱߻�û�е����
			} //end of for. n����, ��n���±��е����

			// simplified_mesh_ and mesh_�������γɵ���֮���Ӧ��.�����е�ϸ��ʱ���õ�.
			mesh_.property(fp_fh_, mesh_.face_handle(loop_h[0])) = simplified_mesh_.face_handle(simplified_mesh_.find_halfedge(v_from, v_to));
			mesh_.property(fp_fh_, mesh_.face_handle(loop_h[n-1])) = simplified_mesh_.face_handle(simplified_mesh_.find_halfedge(v_to, v_from));
			simplified_mesh_.property(fp_fh_sm_, simplified_mesh_.face_handle(simplified_mesh_.find_halfedge(v_from, v_to))) = mesh_.face_handle(loop_h[0]);
			simplified_mesh_.property(fp_fh_sm_, simplified_mesh_.face_handle(simplified_mesh_.find_halfedge(v_to, v_from))) = mesh_.face_handle(loop_h[n-1]);
		} else { /*//  mesh_.is_boundary(v_from_old) == true; // v_to_old -> (v_from_old, v_to_old) �߽����
				 // v_l_old������v_r_old��is_valid() == false
				 if (mesh_.is_boundary(v_to_old) == false ) std::cout << "Error: v_to_old��v_from_old���Ǳ߽��ŶԵ��.\n";
				 if (!simplified_mesh_.is_boundary(v_to) || !simplified_mesh_.is_boundary(v_from)) std::cout << "Error: v_to��v_from���Ǳ߽��ŶԵ��.\n";
				 std::cout << "���ѱ߽綥��: " << ".\n";
				 if (v_l_old.is_valid() == false) { // ����� 
				 TriMesh::HHandle h0 = mesh_.find_halfedge(v_from_old, v_to_old);	
				 if (mesh_.is_boundary(h0) == false) std::cout << "Error: ���(v_from_old, v_to_old)Ӧ���Ǳ߽�Ŷ�.\n"; 
				 std::vector<TriMesh::VertexHandle> loop;
				 std::vector<TriMesh::HalfedgeHandle> loop_h;
				 TriMesh::HHandle h = h0;
				 do {
				 loop_h.push_back(h);
				 loop.push_back(mesh_.to_vertex_handle(h));
				 h = mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(h));
				 } while (h != h0);
				 int n = loop.size();
				 if (n != mesh_.valence(v_from_old)) std::cout << "Error: v_from_old's valence should equal to the loop's size.\n";

				 double rou_angle = 0.0;  //round angle
				 std::vector<double> vec_angle;
				 for (int i=1 ; i<= n-1; ++i) {//����ܵ�angle
				 OpenMesh::Vec3d d1 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).normalize();
				 OpenMesh::Vec3d d2 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i-1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).normalize();
				 rou_angle += acos(dot(d1, d2));
				 vec_angle.push_back(acos(dot(d1, d2)));
				 } 
				 if (vec_angle.size() != n -1) { std::cout << "Error: vec_angle��loop����������һ��.\n"; }
				 double angle_scale_ratio = M_PI / rou_angle; //���ű���
				 double temp_sum_angle = 0.0, l = 0;
				 for (int i = 1; i <= n-1; ++i) { //����loop[1...n-1]�Ĳ���������
				 temp_sum_angle += (vec_angle[i-1] * angle_scale_ratio);
				 l = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).norm(); 
				 l *= angle_scale_ratio;
				 mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(l*cos(temp_sum_angle), l*sin(temp_sum_angle));
				 }
				 mesh_.property(vuv_, loop[0]) = OpenMesh::Vec2d((mesh_.point(loop[0]) - mesh_.point(v_from_old)).norm() * angle_scale_ratio, 0);
				 mesh_.property(vuv_, v_from_old) = OpenMesh::Vec2d(0, 0); //���һ���������

				 std::vector<TriMesh::VHandle> free_vertices;
				 int free_vertices_size = 0;
				 for (int i = 2; i <= n-1; ++i) { //
				 TriMesh::FaceHandle fh = mesh_.face_handle(mesh_.next_halfedge_handle(loop_h[i]));
				 free_vertices_size += mesh_.property(fvset_, fh).size();
				 for (std::vector<TriMesh::VHandle>::iterator it(mesh_.property(fvset_, fh).begin()), it_end(mesh_.property(fvset_, fh).end()); it != it_end; ++it) {
				 if (*it != v_from_old) { 
				 free_vertices.push_back(*it);
				 TriMesh::VertexHandle vf0 = mesh_.property(vf0_, *it); TriMesh::VertexHandle vf1 = mesh_.property(vf1_, *it); TriMesh::VertexHandle vf2 = mesh_.property(vf2_, *it);  
				 mesh_.property(vuv_, *it) = mesh_.property(vuv_, vf0) * (mesh_.property(vbc_, *it))[0] + mesh_.property(vuv_, vf1) * (mesh_.property(vbc_, *it))[1]
				 + mesh_.property(vuv_, vf2) * (mesh_.property(vbc_, *it))[2];
				 }
					}
					mesh_.property(fvset_, fh).clear();
					}
					free_vertices_size -= 1;//��ȥ�Ǹ�v_from_old
					if (free_vertices.size() != free_vertices_size) std::cout << "Error: free_vertices's size.\n";
					// �ж�free_vertices����ĵ������ĸ�������.
					OpenMesh::VPropHandleT<bool> is_face_located;
					mesh_.add_property(is_face_located);
					for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
					mesh_.property(is_face_located, *freev_it) = false;
					}
					for (int i = 0; i <= n-2; ++ i) { 
					for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
					if (mesh_.property(is_face_located, *freev_it) == false) {
					OpenMesh::Vec2d pp = mesh_.property(vuv_, *freev_it);//���free vertex�Ĳ�������
					OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, v_from_old), mesh_.property(vuv_, loop[i]), mesh_.property(vuv_, loop[i+1]), pp);
					if (fabs(bc[0]) < 1.0e-10) bc[0] = 0;  if (fabs(bc[1]) < 1.0e-10) bc[1] = 0;  if (fabs(bc[2]) < 1.0e-10) bc[2] = 0; 

					if ((bc[0] >= 0) && (bc[1] >= 0) && (bc[2] >= 0)) { //Ҫ�����bc��������������>=0
					mesh_.property(is_face_located, *freev_it) = true;
					TriMesh::FaceHandle fh = mesh_.face_handle(loop_h[i+1]);
					if (fh.is_valid() == false) std::cerr << "Error: face invalid." << std::endl;//Only for assert

					// �������о���ÿһ�β�������Ӧ�û�õĻ��Ǹ��ĵ���Ϣ.
					mesh_.property(vf_,  *freev_it) = fh;
					mesh_.property(vf0_, *freev_it) = v_from_old; mesh_.property(vf1_, *freev_it) = loop[i]; mesh_.property(vf2_, *freev_it) = loop[i+1];
					mesh_.property(vbc_, *freev_it) = bc;

					(mesh_.property(fvset_, fh)).push_back(*freev_it);

					TriMesh::Point vfp0 = mesh_.point(v_from_old); TriMesh::Point vfp1 = mesh_.point(loop[i]); TriMesh::Point vfp2 = mesh_.point(loop[i+1]);
					TriMesh::Point vp = vfp0 * bc[0] + vfp1 * bc[1] + vfp2 * bc[2];//��deletedҲ���Ǳ�parameterized�Ķ����3d����.
					mesh_collapsed_.set_point(*freev_it, vp);////
					}
					} // ���free vertex�����������
					} // һ����free vertex �����������			
					} // ���е�free vertices ������n������
					for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
					if (mesh_.property(is_face_located, *freev_it) == false) { std::cout << "Error: free_vertices ��û�ж�λ��.\n"; }
					}
					mesh_.remove_property(is_face_located);
					mesh_collapsed_.set_point(v_from_old, mesh_.point(v_from_old));
					free_vertices_size = 0; //����ȷ��һ���Ƿ������֤
					for (int i = 1; i <= n-1; ++i) {
					free_vertices_size += mesh_.property(fvset_, mesh_.face_handle(loop_h[i])).size();
					}
					if (free_vertices_size != free_vertices.size()) std::cout << "Error: free_vertices's size �����¶�Ϊ�������ʹ���.\n";
					// ��ߵ��е����
					resample_boundary_edge_midpoint(mesh_, loop_h[0], v_from_old, v_to_old);
					resample_boundary_edge_midpoint(mesh_, mesh_.opposite_halfedge_handle(loop_h[n-1]), loop[n-1], v_from_old);
					for (int i = 1; i <= n-2; ++i) {
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) { // ����߻�û���е����, ����Ĳ�����ֻ��ɹ���.
					TriMesh::FHandle fhi = mesh_.face_handle(loop_h[i]), fhi_1 = mesh_.face_handle(loop_h[(i+1)%n]);
					for (std::vector<TriMesh::VHandle>::iterator it(mesh_.property(fvset_, fhi).begin()), end(mesh_.property(fvset_, fhi).end()); it != end; ++it) {
					TriMesh::VertexHandle vh = *it; //��fhi�ϵĲ�������
					for (TriMesh::VertexOHalfedgeIter voh_it(mesh2_, vh); voh_it; ++voh_it) {
					TriMesh::VHandle v0 = vh, v1 = mesh2_.to_vertex_handle(voh_it), v2 = mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(voh_it));
					// ������ʵ������v2����ܲ�������������������.��������û�б�deleted�Ļ����ǲ��������������ϵ����.
					if (mesh_.property(vf_, v2) == fhi|| mesh_.property(vf_, v2) == fhi_1) {
					OpenMesh::Vec2d pp = (mesh_.property(vuv_, loop[i]) + mesh_.property(vuv_, v_from_old)) * 0.5;//����ߵ��е�Ĳ�������
					OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, v0), mesh_.property(vuv_, v1), mesh_.property(vuv_, v2), pp);
					if (fabs(bc[0]) < 1.0e-10) bc[0] = 0; if (fabs(bc[1]) < 1.0e-10) bc[1] = 0; if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;

					if ((bc[0] >= 0) && (bc[1] >= 0) && (bc[2] >= 0)) { //Ҫ�����bc��������������>=0
					mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) = true;
					mesh_.property(emp_, mesh_.edge_handle(loop_h[i])) = mesh2_.point(v0) * bc[0] + mesh2_.point(v1) * bc[1] + mesh2_.point(v2) * bc[2];
					mesh_.property(emp_closest_vh_, mesh_.edge_handle(loop_h[i])) = v0;
					} 
					} else { // end of if . v2 is inside the valid region
					std::cout << "Error: free vertex's neighbor v2 isn't in fhi-1.\n"; 
					}
					}
					} // end of for. ѭ��fhi�ϵ����в�������.
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) {
					std::cout << "Error: ��û���е�����ɹ�.\n";
					}
					}
					}

					} else if (v_r_old.is_valid() == false) { //�����
					TriMesh::HHandle h0 = mesh_.find_halfedge(v_from_old, v_to_old);	
					if (mesh_.is_boundary(mesh_.opposite_halfedge_handle(h0)) == false) std::cout << "Error: ���(v_from_old, v_to_old)�Ķ԰��Ӧ���Ǳ߽�Ŷ�.\n"; 
					std::vector<TriMesh::VertexHandle> loop;
					std::vector<TriMesh::HalfedgeHandle> loop_h;
					TriMesh::HHandle h = h0;
					do {
					loop_h.push_back(h);
					loop.push_back(mesh_.to_vertex_handle(h));
					h = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(h));
					} while (h != h0);
					int n = loop.size();
					if (n != mesh_.valence(v_from_old)) std::cout << "Error: v_from_old's valence should equal to the loop's size.\n";

					double rou_angle = 0.0;  //round angle
					std::vector<double> vec_angle;
					for (int i=1 ; i<= n-1; ++i) {//����ܵ�angle
					OpenMesh::Vec3d d1 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).normalize();
					OpenMesh::Vec3d d2 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i-1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).normalize();
					rou_angle += acos(dot(d1, d2));
					vec_angle.push_back(acos(dot(d1, d2)));
					} 
					if (vec_angle.size() != n -1) { std::cout << "Error: vec_angle��loop����������һ��.\n"; }
					double angle_scale_ratio = M_PI / rou_angle; //���ű���
					double temp_sum_angle = 0.0, l = 0;
					for (int i = 1; i <= n-1; ++i) { //����loop[1...n-1]�Ĳ���������
					temp_sum_angle += (vec_angle[i-1] * angle_scale_ratio);
					l = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).norm(); 
					l *= angle_scale_ratio;
					mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(l*cos(temp_sum_angle), l*sin(temp_sum_angle));
					}
					mesh_.property(vuv_, loop[0]) = OpenMesh::Vec2d((mesh_.point(loop[0]) - mesh_.point(v_from_old)).norm() * angle_scale_ratio, 0);
					mesh_.property(vuv_, v_from_old) = OpenMesh::Vec2d(0, 0); //���һ���������

					std::vector<TriMesh::VHandle> free_vertices;
					int free_vertices_size = 0;
					for (int i = 1; i <= n-2; ++i) { //
					TriMesh::FaceHandle fh = mesh_.face_handle(mesh_.next_halfedge_handle(loop_h[i]));
					free_vertices_size += mesh_.property(fvset_, fh).size();
					for (std::vector<TriMesh::VHandle>::iterator it(mesh_.property(fvset_, fh).begin()), it_end(mesh_.property(fvset_, fh).end()); it != it_end; ++it) {
					if (*it != v_from_old) { 
					free_vertices.push_back(*it);
					TriMesh::VertexHandle vf0 = mesh_.property(vf0_, *it); TriMesh::VertexHandle vf1 = mesh_.property(vf1_, *it); TriMesh::VertexHandle vf2 = mesh_.property(vf2_, *it);  
					mesh_.property(vuv_, *it) = mesh_.property(vuv_, vf0) * (mesh_.property(vbc_, *it))[0] + mesh_.property(vuv_, vf1) * (mesh_.property(vbc_, *it))[1]
					+ mesh_.property(vuv_, vf2) * (mesh_.property(vbc_, *it))[2];
					}
					}
					mesh_.property(fvset_, fh).clear();
					}
					free_vertices_size -= 1;//��ȥ�Ǹ�v_from_old
					if (free_vertices.size() != free_vertices_size) std::cout << "Error: free_vertices's size.\n";
					// �ж�free_vertices����ĵ������ĸ�������.
					OpenMesh::VPropHandleT<bool> is_face_located;
					mesh_.add_property(is_face_located);
					for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
					mesh_.property(is_face_located, *freev_it) = false;
					}
					for (int i = 0; i <= n-2; ++ i) { 
					for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
					if (mesh_.property(is_face_located, *freev_it) == false) {
					OpenMesh::Vec2d pp = mesh_.property(vuv_, *freev_it);//���free vertex�Ĳ�������
					OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, v_from_old), mesh_.property(vuv_, loop[i]), mesh_.property(vuv_, loop[i+1]), pp);
					if (fabs(bc[0]) < 1.0e-10) bc[0] = 0;  if (fabs(bc[1]) < 1.0e-10) bc[1] = 0;  if (fabs(bc[2]) < 1.0e-10) bc[2] = 0; 

					if ((bc[0] >= 0) && (bc[1] >= 0) && (bc[2] >= 0)) { //Ҫ�����bc��������������>=0
					mesh_.property(is_face_located, *freev_it) = true;
					TriMesh::FaceHandle fh = mesh_.face_handle(loop_h[i]);
					if (fh.is_valid() == false) std::cerr << "Error: face invalid." << std::endl;//Only for assert

					// �������о���ÿһ�β�������Ӧ�û�õĻ��Ǹ��ĵ���Ϣ.
					mesh_.property(vf_,  *freev_it) = fh;
					mesh_.property(vf0_, *freev_it) = v_from_old; mesh_.property(vf1_, *freev_it) = loop[i]; mesh_.property(vf2_, *freev_it) = loop[i+1];
					mesh_.property(vbc_, *freev_it) = bc;

					(mesh_.property(fvset_, fh)).push_back(*freev_it);

					TriMesh::Point vfp0 = mesh_.point(v_from_old); TriMesh::Point vfp1 = mesh_.point(loop[i]); TriMesh::Point vfp2 = mesh_.point(loop[i+1]);
					TriMesh::Point vp = vfp0 * bc[0] + vfp1 * bc[1] + vfp2 * bc[2];//��deletedҲ���Ǳ�parameterized�Ķ����3d����.
					mesh_collapsed_.set_point(*freev_it, vp);////
					}
					} // ���free vertex�����������
					} // һ����free vertex �����������			
					} // ���е�free vertices ������n������
					for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
					if (mesh_.property(is_face_located, *freev_it) == false) { std::cout << "Error: free_vertices ��û�ж�λ��.\n"; }
					}
					mesh_.remove_property(is_face_located);
					mesh_collapsed_.set_point(v_from_old, mesh_.point(v_from_old));
					free_vertices_size = 0; //����ȷ��һ���Ƿ������֤
					for (int i = 0; i <= n-2; ++i) {
					free_vertices_size += mesh_.property(fvset_, mesh_.face_handle(loop_h[i])).size();
					}
					if (free_vertices_size != free_vertices.size()) std::cout << "Error: free_vertices's size �����¶�Ϊ�������ʹ���.\n";
					// ��ߵ��е����
					resample_boundary_edge_midpoint(mesh_, mesh_.opposite_halfedge_handle(loop_h[0]), v_to_old, v_from_old);
					resample_boundary_edge_midpoint(mesh_, loop_h[n-1], v_from_old, loop[n-1]);
					for (int i = 1; i <= n-2; ++i) {
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) { // ����߻�û���е����, ����Ĳ�����ֻ��ɹ���.
					TriMesh::FHandle fhi = mesh_.face_handle(loop_h[i]), fhi_1 = mesh_.face_handle(loop_h[i-1]);
					for (std::vector<TriMesh::VHandle>::iterator it(mesh_.property(fvset_, fhi).begin()), end(mesh_.property(fvset_, fhi).end()); it != end; ++it) {
					TriMesh::VertexHandle vh = *it; //��fhi�ϵĲ�������
					for (TriMesh::VertexOHalfedgeIter voh_it(mesh2_, vh); voh_it; ++voh_it) {
					TriMesh::VHandle v0 = vh, v1 = mesh2_.to_vertex_handle(voh_it), v2 = mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(voh_it));
					// ������ʵ������v2����ܲ�������������������.��������û�б�deleted�Ļ����ǲ��������������ϵ����.
					if (mesh_.property(vf_, v2) == fhi|| mesh_.property(vf_, v2) == fhi_1) {
					OpenMesh::Vec2d pp = (mesh_.property(vuv_, loop[i]) + mesh_.property(vuv_, v_from_old)) * 0.5;//����ߵ��е�Ĳ�������
					OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, v0), mesh_.property(vuv_, v1), mesh_.property(vuv_, v2), pp);
					if (fabs(bc[0]) < 1.0e-10) bc[0] = 0; if (fabs(bc[1]) < 1.0e-10) bc[1] = 0; if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;

					if ((bc[0] >= 0) && (bc[1] >= 0) && (bc[2] >= 0)) { //Ҫ�����bc��������������>=0
					mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) = true;
					mesh_.property(emp_, mesh_.edge_handle(loop_h[i])) = mesh2_.point(v0) * bc[0] + mesh2_.point(v1) * bc[1] + mesh2_.point(v2) * bc[2];
					mesh_.property(emp_closest_vh_, mesh_.edge_handle(loop_h[i])) = v0;
					} 
					} else { // end of if . v2 is inside the valid region
					std::cout << "Error: free vertex's neighbor v2 isn't in fhi-1.\n"; 
					}
					}
					} // end of for. ѭ��fhi�ϵ����в�������.
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) {
					std::cout << "Error: ��û���е�����ɹ�.\n";
					}
					}
					}			 
					}*/
		} // end of if (mesh_.is_boundary(v_from_old).v_to_old -> (v_from_old, v_to_old) �ķ��ڲ��ͱ߽�
		// �߽������������ܴ���, ��������:simplified_mesh_ and mesh_�������γɵ���֮�仹û�ж�Ӧ��.�����е�ϸ��ʱ���õ�.

		midpoint_sample_succeed_ = true;//�����±߶��е�����ɹ���
		for (TriMesh::VertexEdgeIter ve_it(mesh_, v_from_old); ve_it; ++ve_it) { // ����v_from_old��һ�����Ƿ��ȴ���е�����ɹ���.
			if (mesh_.property(empl_, ve_it.handle()) == false) midpoint_sample_succeed_ = false;
		}
	} else{ //end of if lc_node_handle != -1 �ɷ���
		std::cout << "Info: vsplit_refine(v_to): lc_node_handle == -1. �Ѿ���Ҷ�ӽڵ���.\n";
	}
	return true;
} // end of function vsplit_refine.
bool DecimationModel::error_driven_selective_refine() {
	std::cout << "����error_driven_selective_refine().\n";
	// ����һ ��icm�������ҷ��ѵĶ���. icm����һ��loopϸ��֮��, ���ж��㶼����ԭʼģ���ϵ��������.
	// ���ⷽ��Ӧ�ò���. Ӧ�ô�ԭʼģ����������������Ȼ���simplified_mesh_���������.
	find_closestp_icm();

	double largest_err = 0;
	TriMesh::VertexHandle vh_largest_err; 
	int i = 0;
	for (TriMesh::VIter v_it(initial_control_mesh_.vertices_begin()), v_end(initial_control_mesh_.vertices_end()); v_it != v_end; ++v_it) {
		// initial_control_mesh_����һ��loopϸ��֮���refined_simplified_mesh_�Ķ����Ӧ.
		if (v_it.handle() != refined_simplified_mesh_.vertex_handle(i)) std::cout << "Error: ������mesh��indexӦ����match��" <<  ".\n"; 
		++ i;
		if (initial_control_mesh_.property(vp_distances_,v_it) > largest_err) {
			largest_err = initial_control_mesh_.property(vp_distances_,v_it);
			vh_largest_err =v_it.handle();
		}
	}
	
	
	std::cout << "largest fitting error: " << largest_err << ".\n";
	std::cout << "simplified mesh's vertices num: " << simplified_mesh_.n_vertices() << ", vh_largest_err: " << vh_largest_err.idx() << ".\n";

	// --------------------------------------------------------
	// ��simplified mesh���ҿ��ܱ����ѵĶ���.
	TriMesh::VertexHandle v0, v1, v2, v3;
	if (vh_largest_err.idx() >= simplified_mesh_.n_vertices()) { //��ʱvh_largest_err�������е�ϸ��ʱ�����ɵ�, ��Ҫ����simplified_mesh_�ϵĵ�
		std::vector<TriMesh::VertexHandle> vec;	//��vh_largest_err��һ���򶥵�, �������ҽ�����������simplified_mesh_�ϵĵ�.
		for (TriMesh::VVIter vv_it(refined_simplified_mesh_, vh_largest_err); vv_it; ++vv_it) {
			if (vv_it.handle().idx() < simplified_mesh_.n_vertices()) {
				vec.push_back(vv_it.handle()); //std::cout << vv_it.handle() << ".\n";				
			}
		}		//std::cout << "vec.size(): " << vec.size() << ".\n";//==2 //for test
		v0 = vec[0]; v1 = vec[1];
		TriMesh::HHandle h0 = simplified_mesh_.find_halfedge(vec[0], vec[1]);
		TriMesh::HHandle h1 = simplified_mesh_.opposite_halfedge_handle(h0);
		if (simplified_mesh_.is_boundary(h0) == false) {
			vec.push_back(simplified_mesh_.to_vertex_handle(simplified_mesh_.next_halfedge_handle(h0)));
			v2 = vec[2];
		}
		if (simplified_mesh_.is_boundary(h1) == false) {
			vec.push_back(simplified_mesh_.to_vertex_handle(simplified_mesh_.next_halfedge_handle(h1)));
			v3 = vec[3];
		}
	} else { // ��֤v0������simplified_mesh_�ϵĵ�
		v0 = vh_largest_err;
	} 
	// -----------------------------------------------------------------------------
	TriMesh::VertexHandle v_to = v0;// v_to�����Ӧ��v_to_node_handle����ڵ�Node������.
	if (v1.is_valid()) { // ��Ҫ��v0, v1��ѡ�������һ��������
		if (vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v0)).is_active() == false) std::cout << "Error: v0's node should be active.\n";
		if (vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v1)).is_active() == false) std::cout << "Error: v1's node should be active.\n";
		TriMesh::VertexHandle v0_old = voldhandles_[vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v0)).point_handle()];
		TriMesh::VertexHandle v1_old = voldhandles_[vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v1)).point_handle()];
		TriMesh::FaceHandle f0_old = mesh_.face_handle(mesh_.find_halfedge(v0_old, v1_old));
		TriMesh::FaceHandle f1_old = mesh_.face_handle(mesh_.find_halfedge(v1_old, v0_old));
		int v0_lc = vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v0)).lchild_node_handle();
		int v1_lc = vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v1)).lchild_node_handle();
		if (v0_lc != -1) {
			int lc_point_handle = vhierarchy_.node(v0_lc).point_handle();
			TriMesh::VHandle lc_old_vh = voldhandles_[lc_point_handle];    // v0���ѵĻ��ͷ��ѳ������.�������ʱ����������������.
			if (mesh_.status(lc_old_vh).deleted() == false) std::cout << "Error: ��ʱ��v0����ڵ���mesh_���Ǳ�deleted�ŶԵ�.\n"; 
			if (mesh_.property(vf_, lc_old_vh) == f0_old|| mesh_.property(vf_, lc_old_vh) == f1_old) v_to = v0;
		}
		if (v1_lc != -1) {
			TriMesh::VHandle lc_old_vh = voldhandles_[vhierarchy_.node(v1_lc).point_handle()];
			if (mesh_.status(lc_old_vh).deleted() == false) std::cout << "Error: ��ʱ��vl����ڵ���mesh_���Ǳ�deleted�ŶԵ�.\n";
			if (mesh_.property(vf_, lc_old_vh) == f0_old|| mesh_.property(vf_, lc_old_vh) == f1_old) v_to = v1;
		}
		if (v0_lc == -1) v_to = v1; // v0�Ѿ���Ҷ�Ӳ��ɷ�����, �Ǿ�ѡ��ĵ��.
	}

	// ----------------------------
	// �����Ѿ�ѡ���v_to�Ƿ���Ҷ�ӻ��߷��ѻ���ɷ�ת.
	int n_change = 0;
	while (true) {
		if (test_is_leaf_of_vertex_hierarchy(v_to) || test_is_triangle_flipping_after_vsplit(v_to)) {
			// ��ʾ���v_toҪô�Ѿ���Ҷ����, // ��ʾ���ѻ���ɷ�ת, ����Ҫ����v_to������.
			std::cout << "bb";
			if (v_to == v0) { 
				if (n_change % 3 == 0) {
					if (v2.is_valid()) v_to = v2;
					else if (v3.is_valid()) v_to = v3; 
				}
				if (v1.is_valid()) v_to = v1;//
			} else {
				if (n_change % 3 == 0) {
					if (v2.is_valid()) v_to = v2;
					else if (v3.is_valid()) v_to = v3; 
				}
				v_to = v0;//
			}
			++ n_change;std::cout << "cc"; 
		} else break;
	}
	//
	std::cout << "---------------------------------------------------.\n";
	vsplit_refine(v_to); //�������, ���е����, ���²�������Ϣ.
	// 
	std::cout << "�뿪error_driven_selective_refine().\n";
	return true;
}
void DecimationModel::error_driven_selective_refine_process() {
	std::cout << "\n";
	std::cout << "========��ʼ���������refine============\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.
	if (error_driven_selective_refine() == false) std::cout << "Error: error driven refine.\n";
	set_mesh_toberendered(&simplified_mesh_, SIMPLIFIED);
	std::cout << "���ڵ�simplified_mesh_: v " << simplified_mesh_.n_vertices() << ", e " << simplified_mesh_.n_edges() << ", f " << simplified_mesh_.n_faces() << ".\n";
	std::cout << "========�������������refine============\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.
}
// ---------------
void DecimationModel::refine_process() {
	std::cout << "\n";
	std::cout << "========��ʼ���������refine============\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.
	// 2009-05-23, ����simplified mesh��ÿһ������������ȷ������.
	int n_split = 0; // for test. ÿ�η���10����
	//for (TriMesh::VIter v_it(simplified_mesh_.vertices_begin()), v_end(simplified_mesh_.vertices_end()); v_it != v_end; ++v_it) {
	//	if (simplified_mesh_.property(vp_fitting_error_sm_, v_it) > standard_error_sm_) {
	//		TriMesh::VHandle v_to = v_it.handle();
	//		if (test_is_leaf_of_vertex_hierarchy(v_to) || test_is_triangle_flipping_after_vsplit(v_to)) {
	//			continue;
	//		} else {
	//			vsplit_refine(v_to); //�������, ���е����, ���²�������Ϣ.
	//			++n_split;
	//		}
	//	}
	//} 
	
	std::vector<TriMesh::VHandle > array_tobe_split  = vsplit_array_sm_;//= vsplit_array_; //
	for (std::vector<TriMesh::VertexHandle>::const_iterator it(array_tobe_split.begin()), it_end(array_tobe_split.end()); it != it_end; ++it) {
		TriMesh::VHandle vh = *it;// use this to split.
		if (test_is_leaf_of_vertex_hierarchy(vh) || test_is_triangle_flipping_after_vsplit(vh)) {
			continue;
		} else {
			vsplit_refine(vh); //�������, ���е����, ���²�������Ϣ.
			++n_split;
		}
	}
	std::cout << "n_split: " << n_split << ".\n";

	set_mesh_toberendered(&simplified_mesh_, RESAMPLE_MIDPOINT);
	std::cout << "���ڵ�simplified_mesh_: v " << simplified_mesh_.n_vertices() << ", e " << simplified_mesh_.n_edges() << ", f " << simplified_mesh_.n_faces() << ".\n";
	std::cout << "========�������������refine============\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.

}
