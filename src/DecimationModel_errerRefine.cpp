#include "DecimationModel.h" 

// ---------------
void DecimationModel::find_closestp_icm() {
	std::cout << "进入find_closestp_icm().\n";
	//假设limit surface上m个sample points, 在original surface上找m个对应的closest points作为foot point
	// 下面是为了保证 refined_simplified_mesh_和initial_control_mesh_上有相同的index而已, 这可以方便以后的计算.
	int n = simplified_mesh_.n_vertices();
	int m = initial_control_mesh_.n_vertices(); // 
	if (m - n != simplified_mesh_.n_edges()) std::cout << "Error: m-n == simplified_mesh_.n_edges().\n";
	if (refined_simplified_mesh_.n_vertices() != m) std::cout << "Error: refined simplified mesh和limit surface的顶点数应该相同.\n";
	//std::cout << "1,"; // for test
	for (int i = 0; i < n; ++i) {
		TriMesh::VertexHandle vh = initial_control_mesh_.vertex_handle(i);//limit surface的前n个顶点.
		if (vh != simplified_mesh_.vertex_handle(i) || vh != refined_simplified_mesh_.vertex_handle(i)) std::cout << "Error: 前n个顶点应该是匹配的.\n";
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
	// 下面考虑中
	for (int i = 0; i < m -n; ++i) { //simplified_mesh_.edge_handle(i) 对应refined_simplified_mesh_.vertex_handle(n + i)
		TriMesh::VertexHandle vh = initial_control_mesh_.vertex_handle(n + i); //vertex: n ... m-1. 一次细分时候m-n条边加入的新点. 
		if (vh != refined_simplified_mesh_.vertex_handle(n + i)) std::cout << "Error: n...m-1顶点应该匹配才对的.\n";
		TriMesh::HalfedgeHandle h0 = simplified_mesh_.halfedge_handle(simplified_mesh_.edge_handle(i), 0);
		TriMesh::VHandle v0 = simplified_mesh_.to_vertex_handle(h0);
		TriMesh::VHandle v1 = simplified_mesh_.from_vertex_handle(h0); 
		int p0 = vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v0)).point_handle();
		int p1 = vhierarchy_.node(simplified_mesh_.property(vp_node_handle_sm_, v1)).point_handle();
		TriMesh::VHandle old_v0 = voldhandles_[p0], old_v1 = voldhandles_[p1]; 
		TriMesh::EHandle old_e = mesh_.edge_handle(mesh_.find_halfedge(old_v0, old_v1)); 
		if (mesh_.property(empl_, old_e) == false ) std::cout << "Error: simplified_mesh_对应的mesh_上还有边没有中点采样.\n";
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
	std::cout << "离开find_closestp_icm().\n";
}
bool DecimationModel::test_is_leaf_of_vertex_hierarchy(TriMesh::VHandle _vh_in_simplified_mesh) {
	// 当不可分裂, 已经是原子顶点(对应于vertex hierarchy上的叶子节点)时返回true.
	TriMesh::VHandle v_to = _vh_in_simplified_mesh;  
	int v_to_node_handle = simplified_mesh_.property(vp_node_handle_sm_, v_to);
	//TriMesh::VertexHandle v_to_old = voldhandles_[vhierarchy_.node(v_to_node_handle).point_handle()];  
	// for test, 注意, 检查的时候注意那v_to_node_handle 不能是-1的, -1表示这个顶点vertex handle还每一个和它对应的Node相连.
	// 这个顶点所对应的Node是需要active的
	if (vhierarchy_.node(v_to_node_handle).is_active() == false ) { std::cout << "Error: Testing: active node才可以split.\n"; }
	int lc_node_handle = vhierarchy_.node(v_to_node_handle).lchild_node_handle();
	//int rc_node_handle = vhierarchy_.node(v_to_node_handle).rchild_node_handle();
	if (lc_node_handle != -1) { // rc_node_handle != -1 //这个节点不是叶子节点,可以往下split
		return false;
	} else {
		std::cout << "Error: Testing: 遇到叶子节点, not alow to split.\n";
		return true;//是叶子节点.
	}
}
bool DecimationModel::test_is_triangle_flipping_after_vsplit(TriMesh::VHandle _vh_in_simplified_mesh) {
	// 当不可分裂, 例如,分裂后会形成triangle flipping, 返回true.
	TriMesh::VHandle v_to = _vh_in_simplified_mesh;  
	int v_to_node_handle = simplified_mesh_.property(vp_node_handle_sm_, v_to);
	// for test, 注意, 检查的时候注意那v_to_node_handle 不能是-1的, -1表示这个顶点vertex handle还每一个和它对应的Node相连.
	// 这个顶点所对应的Node是需要active的
	if (vhierarchy_.node(v_to_node_handle).is_active() == false ) { std::cout << "Error: Testing: active node才可以split.\n"; }
	int lc_node_handle = vhierarchy_.node(v_to_node_handle).lchild_node_handle();
	int rc_node_handle = vhierarchy_.node(v_to_node_handle).rchild_node_handle();
	if (lc_node_handle == -1) {
		std::cout << "Error: test_is_triangle_flipping_after_vsplit()这里不该有叶子分裂的.\n";//前面应该已经负责测试过的了啊.
		return true;
	}

	// lc_node_handle != -1 // rc_node_handle != -1 //这个节点不是叶子节点,可以往下split
	//vsplit(to)->(from, to), 也就是vspit(v_to)->(v0, v_to). vsplit(v1)->(v0, v1)

	TriMesh::Point v_from_point = vpoints_[vhierarchy_.node(lc_node_handle).point_handle()]; 
	TriMesh::VHandle v_from = simplified_mesh_.add_vertex(v_from_point);//这个新增的from点
	TriMesh::VHandle v_l = TriMesh::VertexHandle(-1);//表示invalid vertex handle
	TriMesh::VHandle v_r = TriMesh::VertexHandle(-1);//求出当前split时候的active left/right vertex
	int fc0 = vhierarchy_.node(v_to_node_handle).fund_cut_node_handle0();//std::cout << "fc0 n h: " << fc0 << "\n";//for test, fundamental cut vertex vl^.
	if (fc0 != -1) { // fc0 == -1表示v_l = -1, 是边界情况
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
	// 上面是找出simplified_mesh_分裂顶点时候所需的信息.

	simplified_mesh_.vertex_split(v_from, v_to, v_l, v_r);//std::cout << "0k.\n";
	//std::cout << "vplit simplified_mesh_: " << v_from << ", " << v_to << ", " << v_l << ", " << v_r << "; funt cut ver: " << fc0 << ", " << fc1 << ".\n";//for test//

	// 测试是否发生了vsplit之后的triangle filpping.
	bool triangle_flipped = false;
	for (TriMesh::VertexFaceIter vf_it(simplified_mesh_, v_from); vf_it; ++vf_it) {
		simplified_mesh_.set_normal(vf_it.handle(), simplified_mesh_.calc_face_normal(vf_it.handle()));
	}
	for (TriMesh::VertexFaceIter vf_it(simplified_mesh_, v_from); vf_it; ++vf_it) {
		for (TriMesh::VertexFaceIter vff_it(simplified_mesh_, v_from); vff_it; ++vff_it) {
			TriMesh::Normal f0 = (simplified_mesh_.normal(vf_it));
			TriMesh::Normal fi = (simplified_mesh_.normal(vff_it));
			if (dot(f0, fi)/(f0.norm() * fi.norm()) < -0.866) {
				std::cout << "Error: Testing严重, triangle flipped after this vsplit. v_from's valence=" << simplified_mesh_.valence(v_from) << ".\n";
				triangle_flipped = true; break;
			}
		}
		if (triangle_flipped) break;
	} std::cout << "dadk ";
	simplified_mesh_.collapse(simplified_mesh_.find_halfedge(v_from, v_to)); std::cout << "xx ";
	simplified_mesh_.garbage_collection(); std::cout << "aa ";//一定要这一步否则并没有真正除去那点.
	simplified_mesh_.update_normals();
	if (triangle_flipped) {
		return true;
	} else {
		// 如果不会形成三角形翻转, 就返回false.
		return false;
	}	
}
bool DecimationModel::vsplit_refine(TriMesh::VertexHandle v_to) {
	int v_to_node_handle = simplified_mesh_.property(vp_node_handle_sm_, v_to);
	TriMesh::VertexHandle v_to_old = voldhandles_[vhierarchy_.node(v_to_node_handle).point_handle()];
	std::cout << "Info:将要分裂的是simplified_mesh_.v_to: " << v_to /*<< ", v_to_node_handle: " << v_to_node_handle */
		<< ". mesh_.v_to_old: " << v_to_old << ".\n";
	test_vertex_ = v_to_old;//
	// for test, 注意, 检查的时候注意那v_to_node_handle 不能是-1的, -1表示这个顶点vertex handle还每一个和它对应的Node相连.
	// 这个顶点所对应的Node是需要active的
	if (vhierarchy_.node(v_to_node_handle).is_active() == false ) { std::cout << "Error: active node才可以split.\n"; }
	int lc_node_handle = vhierarchy_.node(v_to_node_handle).lchild_node_handle();
	int rc_node_handle = vhierarchy_.node(v_to_node_handle).rchild_node_handle();
	if (lc_node_handle != -1) { // rc_node_handle != -1 //这个节点不是叶子节点,可以往下split
		//vsplit(to)->(from, to), 也就是vspit(v_to)->(v0, v_to). vsplit(v1)->(v0, v1)
		if (v_to_old != voldhandles_[vhierarchy_.node(rc_node_handle).point_handle()]) std::cout << "Error: 这两个应该是指向同一个旧顶点的.\n";
		TriMesh::VHandle v_from_old = voldhandles_[vhierarchy_.node(lc_node_handle).point_handle()];
		test_vertex1_ = v_from_old;
		std::cout << "mesh_.vsplit: to->(from, to) = (" << v_from_old << ", " << v_to_old << "), type: (" 
			<< mesh_.property(vp_type_, v_from_old) << ", " << mesh_.property(vp_type_, v_to_old) << ").\n";//for test

		TriMesh::Point v_from_point = vpoints_[vhierarchy_.node(lc_node_handle).point_handle()]; 
		TriMesh::VHandle v_from = simplified_mesh_.add_vertex(v_from_point);//这个新增的from点
		TriMesh::VHandle v_l = TriMesh::VertexHandle(-1);//表示invalid vertex handle
		TriMesh::VHandle v_r = TriMesh::VertexHandle(-1);//求出当前split时候的active left/right vertex
		int fc0 = vhierarchy_.node(v_to_node_handle).fund_cut_node_handle0();//std::cout << "fc0 n h: " << fc0 << "\n";//for test, fundamental cut vertex vl^.
		if (fc0 != -1) { // fc0 == -1表示v_l = -1, 是边界情况
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
		// 上面是找出simplified_mesh_和mesh_分裂顶点时候所需的信息.

		simplified_mesh_.vertex_split(v_from, v_to, v_l, v_r);//std::cout << "0k.\n";
		std::cout << "vplit simplified_mesh_: " << v_from << ", " << v_to << ", " << v_l << ", " << v_r << "; funt cut ver: " << fc0 << ", " << fc1 << ".\n";//for test//

		// 测试是否发生了vsplit之后的triangle filpping.
		bool triangle_flipped = false;
		for (TriMesh::VertexFaceIter vf_it(simplified_mesh_, v_from); vf_it; ++vf_it) {
			simplified_mesh_.set_normal(vf_it.handle(), simplified_mesh_.calc_face_normal(vf_it.handle()));
		}
		for (TriMesh::VertexFaceIter vf_it(simplified_mesh_, v_from); vf_it; ++vf_it) {
			for (TriMesh::VertexFaceIter vff_it(simplified_mesh_, v_from); vff_it; ++vff_it) {
				TriMesh::Normal f0 = (simplified_mesh_.normal(vf_it));
				TriMesh::Normal fi = (simplified_mesh_.normal(vff_it));
				if (dot(f0, fi)/(f0.norm() * fi.norm()) < -0.866) {
					std::cout << "Error: 严重, triangle flipped after this vsplit. v_from's valence=" << simplified_mesh_.valence(v_from) << ".\n";
					triangle_flipped = true; break;
				}
			}
			if (triangle_flipped) break;
		}
		if (triangle_flipped) {
			//simplified_mesh_.collapse(simplified_mesh_.find_halfedge(v_from, v_to)); 
			//return false;
		} /**/
		// 如果没有翻转, 继续
		vhierarchy_.node(v_to_node_handle).set_active(false);
		vhierarchy_.node(lc_node_handle).set_vertex_handle(v_from);	simplified_mesh_.property(vp_node_handle_sm_, v_from) = lc_node_handle;
		vhierarchy_.node(lc_node_handle).set_active(true);
		vhierarchy_.node(rc_node_handle).set_vertex_handle(v_to);	simplified_mesh_.property(vp_node_handle_sm_, v_to) = rc_node_handle;
		vhierarchy_.node(rc_node_handle).set_active(true);	//

		// simplified_mesh_分裂之后, mesh_也做相应的分裂.
		if (mesh_.status(v_from_old).deleted() == false) std::cout << "Error: v_from_old should be deleted.\n";// 这时候v_from_old还是之前删去的点.
		//std::cout << "v_from_old, vbc_: " << mesh_.property(vbc_, v_from_old) << ". point: " << mesh_.point(v_from_old) << ".\n";
		if (mesh_.status(v_to_old).deleted()) std::cout << "Error: v_to_old should be not deleted.\n";
		mesh_.vertex_split(v_from_old, v_to_old, v_l_old, v_r_old);
		std::cout << "vplit mesh_: " << v_from_old << ", " << v_to_old << ", " << v_l_old << ", " << v_r_old << ".\n";//for test//

		mesh_.status(v_from_old).set_deleted(false); //从新激活v_from_old顶点.
		//if (mesh_.status(v_from_old).deleted() == true) std::cout << "Error: v_from_old should be not deleted, now.\n";
		// v_from 顶点one-ring的边都需要重新求中点, 主要是在mesh_下操作
		// 这里有一点要清晰的: 基网格上其它边都是没有变化的, 也不需要重新采样.
		for (TriMesh::VertexEdgeIter ve_it(mesh_, v_from_old); ve_it; ++ve_it) { //这些都是新生成的边.
			mesh_.property(empl_, ve_it) = false; 
		}
		midpoint_sample_succeed_ = false; // v_from_old被激活, 新生的邻边都是没有中点采样的.

		if (mesh_.is_boundary(v_from_old) == false) { // v_to_old -> (v_from_old, v_to_old) 内部情况
			// 新生成的顶点v_from和v_from_old的顶点类型其实都应该设置的.
			if (mesh_.property(vp_type_, v_to_old) != simplified_mesh_.property(vp_type_sm_, v_to)) std::cout << "Error: v_to's should be = v_to_old's type.\n";

			TriMesh::HalfedgeHandle cehe0, cehe1;
			// 当v_to/v_to_old的点类型为dart, crease, corner时候, 那v_from/v_from_old都是crease点, 因为dart/corner不会成为from点而被收缩掉的.
			if (mesh_.property(vp_type_, v_from_old) == DGP::CREASE_VFT) 
			{	//假如v_from_old本来是crease点, 那么一头一尾的都是crease edges才对
				std::cout << "分裂特征点: v_to_old's vp_type_: " << mesh_.property(vp_type_, v_to_old) << ", v_from_old's vp_type_: " << mesh_.property(vp_type_, v_from_old) << ".\n";
				//设置simplified_mesh_上的点和边的特征类型.
				int n_crease_edges = 0; //下面先处理一特殊情况.
				for (TriMesh::VertexEdgeIter ve_it(simplified_mesh_, v_from); ve_it; ++ve_it) {
					if (simplified_mesh_.status(ve_it.handle()).feature())  ++ n_crease_edges;
				}
				if (n_crease_edges != 1) { //本来是应该有一个边的,就是(v_feature and v_from), 如果不是的话就出现了下面需处理的特殊情况
					if (mesh_.status(mesh_.edge_handle(mesh_.find_halfedge(v_l_old, v_to_old))).feature() == false 
						&& mesh_.status(mesh_.edge_handle(mesh_.find_halfedge(v_r_old, v_to_old))).feature() == false) { //这是不会出现的情况吧
							std::cout << "Error: 现在应该先只有有1个特征边才对的.\n";
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

				// 设置v_from的crease类型.
				simplified_mesh_.status(simplified_mesh_.edge_handle(simplified_mesh_.find_halfedge(v_from, v_to))).set_feature(true);
				n_crease_edges = 2; // 现在在simplified_mesh_上也是两个crease edges了.因为简化时候是两个crease edge合成一个了.
				simplified_mesh_.property(vp_type_sm_, v_from) = DGP::CREASE_VFT;

				// 既然有两个特征边, 那这两个特征边就可以先分配参数化到其上面的特征点并求中点了.
				for (TriMesh::VertexOHalfedgeIter voh_it(mesh_, v_from_old); voh_it; ++voh_it) {
					if (mesh_.status(mesh_.edge_handle(voh_it)).feature()) {
						cehe1 = voh_it; break;
					}
				} // 先找出cehe1这个以前的特征边在设置刚新生成的特征边. 
				cehe0 = mesh_.find_halfedge(v_from_old, v_to_old);
				mesh_.status(mesh_.edge_handle(cehe0)).set_feature(true);				
				TriMesh::HHandle ocehe0 = mesh_.opposite_halfedge_handle(cehe0), ocehe1 = mesh_.opposite_halfedge_handle(cehe1); //std::cout << "here1.\n";
				int hep_heset_size = mesh_.property(hep_heset_, ocehe1).size(); //std::cout << "here2.\n"; // for test.
				if (mesh_.property(hep_heset_, cehe1).size() != hep_heset_size) std::cout << "Error: ocehe1和cehe1的hep_heset_大小应该是一样才对的(内容有点区别).\n";
				// 这时cehe1 and ocehe1边上的特征点是原来还没有分裂时候的所有点, 下面才开始这些特征点的分配.
				//std::cout << "hep_heset_size " << hep_heset_size << ".\n";
				int i = 0; // 分配ocehe1上的hep_heset_  
				std::vector<TriMesh::HalfedgeHandle> tmp;
				for (; i < mesh_.property(hep_heset_, ocehe1).size(); ++i) {
					if (mesh2_.to_vertex_handle(mesh_.property(hep_heset_, ocehe1)[i]) != v_from_old) {
						tmp.push_back(mesh_.property(hep_heset_, ocehe1)[i]);
					} else break;
				} 
				tmp.push_back(mesh_.property(hep_heset_, ocehe1)[i]);// 0 - i 原来的
				for (++i; i < mesh_.property(hep_heset_, ocehe1).size(); ++i) {
					mesh_.property(hep_heset_, cehe0).push_back(mesh_.property(hep_heset_, ocehe1)[i]);
				} 
				mesh_.property(hep_heset_, ocehe1).clear();
				for (std::vector<TriMesh::HHandle>::iterator it(tmp.begin()), end(tmp.end()); it != end; ++it) {
					mesh_.property(hep_heset_, ocehe1).push_back(*it);
				}
				if ((mesh_.property(hep_heset_, ocehe1).size() + mesh_.property(hep_heset_, cehe0).size()) != hep_heset_size) std::cout << "Error: cehe0, ocehe1的hep_heset_分配错了.\n";
				//std::cout << mesh_.property(hep_heset_, ocehe1).size() << ", " <<  mesh_.property(hep_heset_, cehe0).size() << ".\n";
				//更新_hh边上的点的参数边
				for (int i = 0; i < mesh_.property(hep_heset_, cehe0).size() -1; ++i) {
					mesh_.property(vp_feature_edge_, mesh2_.to_vertex_handle(mesh_.property(hep_heset_, cehe0)[i])) = mesh_.edge_handle(cehe0);
				}
				for (int i = 0; i < mesh_.property(hep_heset_, ocehe1).size() -1; ++i) {
					mesh_.property(vp_feature_edge_, mesh2_.to_vertex_handle(mesh_.property(hep_heset_, ocehe1)[i])) = mesh_.edge_handle(ocehe1);
				}

				i = 0; // 分配cehe1上的hep_heset_ //std::cout << "cehe1's 初始size: " << mesh_.property(hep_heset_, cehe1).size() << ".\n";
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
				} //std::cout << "cehe1's 后来size " << mesh_.property(hep_heset_, cehe1).size() << ".\n";
				if ((mesh_.property(hep_heset_, ocehe0).size() + mesh_.property(hep_heset_, cehe1).size()) != hep_heset_size) std::cout << "Error: ocehe0, cehe1的hep_heset_分配错了.\n";
				for (int i = 0; i < mesh_.property(hep_heset_, ocehe0).size() -1; ++i) {
					mesh_.property(vp_feature_edge_, mesh2_.to_vertex_handle(mesh_.property(hep_heset_, ocehe0)[i])) = mesh_.edge_handle(ocehe0);
				}
				for (int i = 0; i < mesh_.property(hep_heset_, cehe1).size() -1; ++i) {
					mesh_.property(vp_feature_edge_, mesh2_.to_vertex_handle(mesh_.property(hep_heset_, cehe1)[i])) = mesh_.edge_handle(cehe1);
				}
				if ( mesh_.property(hep_heset_, cehe0).size() != mesh_.property(hep_heset_, ocehe0).size()) std::cout << "Error: cehe0, ocehe0的hep_heset_应该等大才对.\n";
				if ( mesh_.property(hep_heset_, cehe1).size() != mesh_.property(hep_heset_, ocehe1).size()) std::cout << "Error: cehe1, ocehe1的hep_heset_应该等大才对.\n";
				//std::cout << " " << mesh_.property(hep_heset_, cehe0).size() << ", " << mesh_.property(hep_heset_, cehe1).size() << ".\n";
				if (mesh2_.to_vertex_handle(mesh_.property(hep_heset_, cehe0)[mesh_.property(hep_heset_, cehe0).size()-1]) != v_to_old) { std::cout << "Error: cehe0.\n";}
				if (mesh2_.to_vertex_handle(mesh_.property(hep_heset_, ocehe0)[mesh_.property(hep_heset_, ocehe0).size()-1]) != v_from_old) { std::cout << "Error: ocehe0.\n";}
				if (mesh2_.to_vertex_handle(mesh_.property(hep_heset_, ocehe1)[mesh_.property(hep_heset_, ocehe1).size()-1]) != v_from_old) { std::cout << "Error: ocehe1.\n";}
				// 至此完成了特征边分裂时候的特征信息重分配.
				// -----------------------------------------------------------------------------------
				//对刚分裂的这两个crease edge求中点
				std::vector<TriMesh::HHandle> cehes; cehes.push_back(cehe0); cehes.push_back(cehe1);
				for (int ii = 0; ii < 2; ++ii) { 
					TriMesh::HalfedgeHandle cehe = cehes[ii];
					double len = 0;
					for (int i = 0; i < mesh_.property(hep_heset_, cehe).size(); ++i) {
						TriMesh::HHandle h = mesh_.property(hep_heset_, cehe)[i];
						len += (mesh2_.point(mesh2_.from_vertex_handle(h)) - mesh2_.point(mesh2_.to_vertex_handle(h))).norm();
					}
					double half_len = len / 2.0, len_tmp = 0;
					for (int i = 0; i < mesh_.property(hep_heset_, cehe).size(); ++i) {// 计算这个crease edge的中点
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
				std::cout << "分裂smooth点: v_to_old's vp_type_: " << mesh_.property(vp_type_, v_to_old) << ", v_from_old's vp_type_: " << mesh_.property(vp_type_, v_from_old) << ".\n";
				simplified_mesh_.property(vp_type_sm_, v_from) = DGP::SMOOTH_VFT;
				if (mesh_.property(vp_type_, v_from_old) != DGP::SMOOTH_VFT) std::cout << "Error: mesh_'s v_from_old's type 应该是本来就有才对喔.\n";
				//if (v_from_old.idx() == 2207 && v_to_old.idx() == 2195) {
				//	std::cout << mesh_.valence(TriMesh::VHandle(2004)) << " ?= " << simplified_mesh_.valence(v_from) <<  ".\n";
				//	for (TriMesh::VVIter vv_it(mesh_, TriMesh::VHandle(2004)); vv_it; ++vv_it) {
				//		std::cout << vv_it.handle() << ": " << mesh_.property(vp_type_, vv_it) << ".\t";
				//	} std::cout << ".\n";
				//}
			}// end of if-else v_from_old是crease点

			// ------------------------------------------------------------------------------------------------------
			// 上面完成了点的分裂以及顶点类型的判定, 下面更新v_from_old一邻域面里面的参数化信息.

			// 删去v_from_old顶点在原来那个参数面上的记录
			TriMesh::FHandle vf_v_from_old = mesh_.property(vf_, v_from_old);
			std::vector<TriMesh::VHandle>::iterator to_be_earsed 
				= remove((mesh_.property(fvset_, vf_v_from_old)).begin(), (mesh_.property(fvset_, vf_v_from_old)).end(), v_from_old);
			(mesh_.property(fvset_, vf_v_from_old)).erase(to_be_earsed, (mesh_.property(fvset_, vf_v_from_old)).end());

			// 压平v_from_old的一邻域, 这里好像没有错误的了
			TriMesh::HalfedgeHandle h0 = mesh_.find_halfedge(v_from_old, v_to_old);
			if (mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h0)) != v_l_old) std::cout << "Error: 压平v_from_old时, h0.\n";
			std::vector<TriMesh::HalfedgeHandle> loop_h;
			std::vector<TriMesh::VertexHandle> loop; 
			int feature_vertex_index = -1;
			TriMesh::HHandle hh = h0;
			do {
				loop_h.push_back(hh); 
				loop.push_back(mesh_.to_vertex_handle(hh));
				if (mesh_.status(mesh_.edge_handle(hh)).feature()) feature_vertex_index = loop_h.size() - 1;
				if (mesh_.status(mesh_.to_vertex_handle(hh))  .deleted()) std::cout << "Error: v_from_old的一邻域点都不该deleted的.\n";
				hh = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(hh));
			} while (hh != h0);
			if (loop_h.size() != mesh_.valence(v_from_old)) std::cout << "Error: loop_h's size.\n";

			int n = loop.size();
			/* double cir_length = 0; // 建立单位圆的参数化域.
			for (int i = 0; i < n; ++i) {
			cir_length += (mesh2_.point(loop[i]) - mesh2_.point(loop[(i+1)%n])).norm();
			}
			double l = 0, angle = 0, wij = 0;
			double sum_wij = 0.0;
			OpenMesh::Vec2d sum_wh(0, 0, 0);
			for (int i = 0; i < n; ++i) { //顶点vi在one-ring上的有n个邻接点
			// fix the boundary/one-ring vertices, 
			angle = l / cir_length * (2.0*M_PI); //std::cout << "ang: " << angle << std::endl;
			mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(0.5*cos(angle)+0.5, 0.5*sin(angle)+0.5);
			//这个单位圆是以(0.5, 0.5)为圆心, 半径是0.5. 圆心要是定在(0, 0)的话, 例如本来是(0, 0.5)的坐标会机器误差为(e, 0.5),e是一个很小很小的值.
			l += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();

			// 求出vi和vj这半边上的系数, 注意这里使用的是mean value coordiante, 而不是harmonic map中的cotangent weight
			TriMesh::VertexHandle v0 = v_from_old; // = mesh_.from_vertex_handle(loop_h[i]);
			TriMesh::VertexHandle v1 = loop[i];// =mesh_.to_vertex_handle(loop_h[i]);
			TriMesh::VertexHandle v2 = loop[(i+1)%n];//= mesh_.to_vertex_handle(mesh_.next_halfedge_handle(loop_h[i]));//这以下两句正确的前提是非边界情况
			TriMesh::VertexHandle v3 = loop[(i+n-1)%n];//mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(loop_h[i])));
			OpenMesh::Vec3d v1v0 = mesh2_.point(v1) - mesh2_.point(v0);	// 之前这里出现的-1
			OpenMesh::Vec3d v2v0 = mesh2_.point(v2) - mesh2_.point(v0);
			OpenMesh::Vec3d v3v0 = mesh2_.point(v3) - mesh2_.point(v0);
			double v1v0v2_angle = acos(dot(v1v0, v2v0) / (v1v0.norm() * v2v0.norm()));//矢量v1v0和v2v0的夹角.
			double v3v0v1_angle = acos(dot(v3v0, v1v0) / (v3v0.norm() * v1v0.norm()));//矢量v3v0和v1v0的夹角.
			wij = (tan(v1v0v2_angle / 2) + tan(v3v0v1_angle / 2)) / v1v0.norm();
			if (wij < 0) { std::cout << "Error: wij is negative.\n"; }

			sum_wij += wij;
			sum_wh += wij * mesh_.property(vuv_, loop[i]);
			}
			mesh_.property(vuv_, v_from_old) = sum_wh * (1.0/sum_wij);//这就是vi顶点在单位圆上的参数化坐标.
			*/
			//这种压迫方法使得有些点错位到别的面上了, 因为分裂特征边上from点不再两个特征边直线上.
			double rou_angle = 0.0;  //round angle
			std::vector<double> vec_angle;
			for (int i=0 ; i<loop.size(); ++i) {//求出总的angle
				OpenMesh::Vec3d d1 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).normalize();
				OpenMesh::Vec3d d2 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).normalize();
				rou_angle += acos(dot(d1, d2));
				vec_angle.push_back(acos(dot(d1, d2)));
			}  //std::cout << "angle " << rou_angle << ".\n";//std::cout << "len: " << length << std::endl; //for test	
			if (vec_angle.size() != n) { std::cout << "Error: vec_angle和loop数组容量不一致.\n"; }
			double angle_scale_ratio = 2 * M_PI / rou_angle; //缩放比率
			double temp_sum_angle = 0.0, l = 0;
			for (int i = 0; i < loop.size(); ++i) {
				temp_sum_angle += (vec_angle[i] * angle_scale_ratio);
				l = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).norm(); 
				l *= angle_scale_ratio;
				mesh_.property(vuv_, loop[(i+1)%n]) = OpenMesh::Vec2d(l*cos(temp_sum_angle), l*sin(temp_sum_angle));
			}
			//std::cout << feature_vertex_index ;//			
			/*// 后来想想只是需要在初始参数化时候需要考虑这个, 现在一般的压平应该不需要吧.
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
			mesh_.property(vuv_, v_from_old) = OpenMesh::Vec2d(0, 0); //完成一邻域参数化 
			//}

			// 取出原来1到n-2面里面的参数化点, 根据上面的压平得到的新参数域, 重新求出参数化坐标
			std::vector<TriMesh::VertexHandle> free_vertices;
			int free_vertices_size = 0;
			for (int i = 1; i <= n-2; ++i) { //loop_h[0 and n-1]指向的面是新形成的.
				TriMesh::FaceHandle fh = mesh_.face_handle(mesh_.next_halfedge_handle(loop_h[i]));
				free_vertices_size += mesh_.property(fvset_, fh).size();
				//std::cout << mesh_.status(fh).deleted() << ", " << mesh_.property(fvset_, fh).size() << ".\n";

				for (std::vector<TriMesh::VHandle>::iterator it(mesh_.property(fvset_, fh).begin()), it_end(mesh_.property(fvset_, fh).end()); it != it_end; ++it) {
					free_vertices.push_back(*it);
					//if (mesh_.property(vp_type_, *it) != DGP::CREASE_VFT) {
					mesh_.property(vuv_, *it) = mesh_.property(vuv_, mesh_.property(vf0_, *it)) * (mesh_.property(vbc_, *it))[0]
					+ mesh_.property(vuv_, mesh_.property(vf1_, *it)) * (mesh_.property(vbc_, *it))[1]
					+ mesh_.property(vuv_,mesh_.property(vf2_, *it) ) * (mesh_.property(vbc_, *it))[2];
					//} else { // crease vertex其mesh_.property(vuv_, *it)值在下面特别的求出.
					// 不能这么if-else做的,因为, 一邻域的外边都可能是特征边, 这些边上的特征点都可能加入到free_vertices数组上, 但是下面却没有给出这些特殊点的参数化坐标
					//} // end of if *it 是否是crease vertex, 这里的*it要么是smooth要么就是crease.
					// 无论*it是否特征点都加入到free_vertices数组中以便下面统一求其落到哪面上.

				}
				mesh_.property(fvset_, fh).clear();
			}

			if (free_vertices.size() != free_vertices_size) std::cout << "Error: free_vertices's size.\n";
			// free_vertices数组中也包含了那些特征边分裂时候的特征点, 这些特征点的参数坐标可以特别地求出, 其实前面.
			if (mesh_.property(vp_type_, v_from_old) == DGP::CREASE_VFT) {
				if (!cehe0.is_valid() || !cehe1.is_valid()) std::cout << "Error: 既然v_from_old是特征点, 那cehe0, cehe1就该有意义的了.\n";
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

			// 判断free_vertices数组的点落在哪个新面中.
			OpenMesh::VPropHandleT<bool> is_face_located;
			mesh_.add_property(is_face_located);
			for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
				mesh_.property(is_face_located, *freev_it) = false;
			}
			for (int i = 0; i < n; ++ i) { // 循环n个面, 这n个面是假设to->(from, to)之后形成的
				// 判断那所有需要重新定位的free vertices是否在这个triangle(v_from_old, loop[i], loop[i+1])之内
				for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
					if (mesh_.property(is_face_located, *freev_it) == false) {

						OpenMesh::Vec2d pp = mesh_.property(vuv_, *freev_it);//这个free vertex的参数坐标
						OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, v_from_old), mesh_.property(vuv_, loop[i]), mesh_.property(vuv_, loop[(i+1)%n]), pp);
						//std::cout << bc << std::endl;//for test

						if (DGP::is_valid_barycentric_coordinate(bc)) { //要求这个bc的三个分量都是>=0

							mesh_.property(is_face_located, *freev_it) = true;

							// For interior vertex, its' barycentric coordinates is bc according to the triangle(0, i, i+1).
							//std::cout << *freev_it << ": " << v_from_old << "," << loop[i] << "," << loop[i+1] << ": " << bc << "\n"; // for test

							TriMesh::FaceHandle fh = mesh_.face_handle(loop_h[i]);
							if (fh.is_valid() == false) std::cerr << "Error: face invalid." << std::endl;//Only for assert

							// 下面四行就是每一次参数化所应该获得的或是更改的信息.
							mesh_.property(vf_,  *freev_it) = fh;
							mesh_.property(vf0_, *freev_it) = v_from_old; mesh_.property(vf1_, *freev_it) = loop[i]; mesh_.property(vf2_, *freev_it) = loop[(i+1)%n];
							mesh_.property(vbc_, *freev_it) = bc;

							(mesh_.property(fvset_, fh)).push_back(*freev_it);

							//被deleted也就是被parameterized的顶点的3d坐标.
							mesh_collapsed_.set_point(*freev_it, mesh_.point(v_from_old) * bc[0] + mesh_.point(loop[i]) * bc[1] + mesh_.point(loop[(i+1)%n]) * bc[2]);////
						}
					} // 这个free vertex落在这个面上
				} // 一部分free vertex 落在这个面上			
			} // 所有的free vertices 落在这n个面上
			// 当mesh_.property(vp_type_, v_from_old) == DGP::CREASE_VFT时候, free_vertices数组中属于两个特征边的(特征)点的参数化坐标(bc0, 0, bc2)
			for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
				if (mesh_.property(is_face_located, *freev_it) == false) { std::cout << "Error: free_vertices 还没有定位完.\n"; }
			}
			mesh_.remove_property(is_face_located);
			mesh_collapsed_.set_point(v_from_old, mesh_.point(v_from_old));
			free_vertices_size = 0; //重新确认一下是否个数保证
			for (int i = 0; i < n; ++i) {
				free_vertices_size += mesh_.property(fvset_, mesh_.face_handle(loop_h[i])).size();
			}
			if (free_vertices_size != free_vertices.size()) std::cout << "Error: free_vertices's size 在重新定为到新面后就错了.\n";


			// 求新形成的n个边的中点, 这些中点是必须求出来的, 否则没法做下去.
			flatten_vertex_one_ring_resample_edges_midpoint(v_from_old);
			int n_tmp_count = 0;//测试上面按点压平一邻域的方法采样的成功率.
			for (int i = 0; i < n; ++i) {
				if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])))  ++ n_tmp_count;
			}
			std::cout << "Info: flatten vertex one ring resample " << n_tmp_count << "/" << mesh_.valence(v_from_old) << ".\n";

			for (int i = 0; i < n; ++i) {
				// 先判断这个边mesh_.edge_handle(loop_h[i])是否是原子边了, 这是很可能的.
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

				if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) { // 这个边还没有中点采样, 这里的采样是只许成功的.
					flatten_face_resample_3edge_midpoint(mesh_.face_handle(loop_h[i]));
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) { 
						//flatten_face_resample_3edge_midpoint(mesh_.face_handle(loop_h[(i-1)%n]));
					}
					std::cout << "resample loop_h[" << i << "]: " << mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) << ".\n";//输出0表示前面的采样还没有成功.

					TriMesh::FHandle fhi = mesh_.face_handle(loop_h[i]), fhi_1 = mesh_.face_handle(loop_h[(i+n-1)%n]);
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) {
						for (std::vector<TriMesh::VHandle>::iterator it(mesh_.property(fvset_, fhi).begin()), end(mesh_.property(fvset_, fhi).end()); it != end; ++it) {
							TriMesh::VertexHandle vh = *it; //面fhi上的参数化点
							for (TriMesh::VertexOHalfedgeIter voh_it(mesh2_, vh); voh_it; ++voh_it) {
								TriMesh::VHandle v0 = vh, v1 = mesh2_.to_vertex_handle(voh_it);
								if (mesh_.property(vf_, v1) == fhi_1) {
									TriMesh::VHandle v2 = mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(voh_it));
									// 这里其实忽略了v2点可能不是在这个参数域里面的.还可能是没有被deleted的或者是参数化到特征边上的情况.
									if (mesh_.property(vf_, v2) == fhi|| mesh_.property(vf_, v2) == fhi_1) {
										OpenMesh::Vec2d pp = (mesh_.property(vuv_, loop[i]) + mesh_.property(vuv_, v_from_old)) * 0.5;//这个边的中点的参数坐标
										OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, v0), mesh_.property(vuv_, v1), mesh_.property(vuv_, v2), pp);
										//std::cout << bc << std::endl;//for test

										//if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //要求这个bc的三个分量都是>=0
										if (DGP::is_valid_barycentric_coordinate(bc)) { //要求这个bc的三个分量都是>=0
											mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) = true; std::cout << "正常方法搞定了边中点采样.\n";
											mesh_.property(emp_, mesh_.edge_handle(loop_h[i])) = mesh2_.point(v0) * bc[0] + mesh2_.point(v1) * bc[1] + mesh2_.point(v2) * bc[2];
											mesh_.property(emp_closest_vh_, mesh_.edge_handle(loop_h[i])) = v0;
										} 
									} else { // end of if . v2 is inside the valid region
										//std::cout << "Error: free vertex's neighbor v2 isn't in fhi-1.\n"; 
									}
								} else { // end of if.测试
									//std::cout << "Error: free vertex's neighbor isn't in fhi-1.\n"; 
								}
							}
						} // end of for. 循环fhi上的所有参数化点.
					}
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) { std::cout << "Waring: 这边中点采样只是找了最近点.\n"; 
					use_closepoints(mesh_.edge_handle(loop_h[i]));
					}

					// 在实验中发现其实上面的中点采样方法绝大情况下都已经成功了, 下面只是预防而已.
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) {
						// 这里无论用什么办法都要弄出中点来, 暴力方法
						std::cout << "暴力方法 fhi.fvset_.size: " << mesh_.property(fvset_, fhi).size() << ".\n";
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
									} else { // 这里应该正常求出v2的参数化坐标的.
										std::cout << "v2 isn't crease, " << mesh_.property(vp_type_, v2) << ".\t";//v2's type == smooth
										mesh_.property(vuv_, v2) = mesh_.property(vuv_, mesh_.property(vf0_, v2)) * (mesh_.property(vbc_, v2))[0]
										+ mesh_.property(vuv_, mesh_.property(vf1_, v2)) * (mesh_.property(vbc_, v2))[1]
										+ mesh_.property(vuv_, mesh_.property(vf0_, v2)) * (mesh_.property(vbc_, v2))[2];std::cout << "bbb\t";
									}
									std::cout <<  mesh_.property(vuv_, v2) << ".\n";
									OpenMesh::Vec2d pp = (mesh_.property(vuv_, loop[i]) + mesh_.property(vuv_, v_from_old)) * 0.5;//这个边的中点的参数坐标
									OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, v0), mesh_.property(vuv_, v1), mesh_.property(vuv_, v2), pp);
									std::cout << bc << ".\n";//for test

									if (fabs(bc[0]) < 1.0e-10) bc[0] = 0; if (fabs(bc[1]) < 1.0e-10) bc[1] = 0; if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;

									//if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //要求这个bc的三个分量都是>=0
									if ((bc[0] >= 0) && (bc[1] >= 0) && (bc[2] >= 0)) { //要求这个bc的三个分量都是>=0
										mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) = true; std::cout << "True." << mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) << ".\t";
										mesh_.property(emp_, mesh_.edge_handle(loop_h[i])) = mesh2_.point(v0) * bc[0] + mesh2_.point(v1) * bc[1] + mesh2_.point(v2) * bc[2];
										mesh_.property(emp_closest_vh_, mesh_.edge_handle(loop_h[i])) = v0;
									} 

								} 
							} else { // mesh2_.to_vertex_handle(voh_it) 不是被deleted的参数化点,
								if (mesh2_.to_vertex_handle(voh_it) == v_from_old) {}
							} 
						} 
						} // end of for. feature_vers.size()
					} // end of if .这里面使用的是暴力方法.

					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) { std::cout << "Error: 严重, 还有新的边没有完成中点采样.\n"; }
				} // end of if此边还没中点采样
			} //end of for. n个边, 对n个新边中点采样

			// simplified_mesh_ and mesh_分裂新形成的面之间对应上.这在中点细分时候用到.
			mesh_.property(fp_fh_, mesh_.face_handle(loop_h[0])) = simplified_mesh_.face_handle(simplified_mesh_.find_halfedge(v_from, v_to));
			mesh_.property(fp_fh_, mesh_.face_handle(loop_h[n-1])) = simplified_mesh_.face_handle(simplified_mesh_.find_halfedge(v_to, v_from));
			simplified_mesh_.property(fp_fh_sm_, simplified_mesh_.face_handle(simplified_mesh_.find_halfedge(v_from, v_to))) = mesh_.face_handle(loop_h[0]);
			simplified_mesh_.property(fp_fh_sm_, simplified_mesh_.face_handle(simplified_mesh_.find_halfedge(v_to, v_from))) = mesh_.face_handle(loop_h[n-1]);
		} else { /*//  mesh_.is_boundary(v_from_old) == true; // v_to_old -> (v_from_old, v_to_old) 边界分裂
				 // v_l_old或者是v_r_old是is_valid() == false
				 if (mesh_.is_boundary(v_to_old) == false ) std::cout << "Error: v_to_old和v_from_old都是边界点才对的喔.\n";
				 if (!simplified_mesh_.is_boundary(v_to) || !simplified_mesh_.is_boundary(v_from)) std::cout << "Error: v_to和v_from都是边界点才对的喔.\n";
				 std::cout << "分裂边界顶点: " << ".\n";
				 if (v_l_old.is_valid() == false) { // 情况二 
				 TriMesh::HHandle h0 = mesh_.find_halfedge(v_from_old, v_to_old);	
				 if (mesh_.is_boundary(h0) == false) std::cout << "Error: 半边(v_from_old, v_to_old)应该是边界才对.\n"; 
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
				 for (int i=1 ; i<= n-1; ++i) {//求出总的angle
				 OpenMesh::Vec3d d1 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).normalize();
				 OpenMesh::Vec3d d2 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i-1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).normalize();
				 rou_angle += acos(dot(d1, d2));
				 vec_angle.push_back(acos(dot(d1, d2)));
				 } 
				 if (vec_angle.size() != n -1) { std::cout << "Error: vec_angle和loop数组容量不一致.\n"; }
				 double angle_scale_ratio = M_PI / rou_angle; //缩放比率
				 double temp_sum_angle = 0.0, l = 0;
				 for (int i = 1; i <= n-1; ++i) { //设置loop[1...n-1]的参数化坐标
				 temp_sum_angle += (vec_angle[i-1] * angle_scale_ratio);
				 l = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).norm(); 
				 l *= angle_scale_ratio;
				 mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(l*cos(temp_sum_angle), l*sin(temp_sum_angle));
				 }
				 mesh_.property(vuv_, loop[0]) = OpenMesh::Vec2d((mesh_.point(loop[0]) - mesh_.point(v_from_old)).norm() * angle_scale_ratio, 0);
				 mesh_.property(vuv_, v_from_old) = OpenMesh::Vec2d(0, 0); //完成一邻域参数化

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
					free_vertices_size -= 1;//减去那个v_from_old
					if (free_vertices.size() != free_vertices_size) std::cout << "Error: free_vertices's size.\n";
					// 判断free_vertices数组的点落在哪个新面中.
					OpenMesh::VPropHandleT<bool> is_face_located;
					mesh_.add_property(is_face_located);
					for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
					mesh_.property(is_face_located, *freev_it) = false;
					}
					for (int i = 0; i <= n-2; ++ i) { 
					for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
					if (mesh_.property(is_face_located, *freev_it) == false) {
					OpenMesh::Vec2d pp = mesh_.property(vuv_, *freev_it);//这个free vertex的参数坐标
					OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, v_from_old), mesh_.property(vuv_, loop[i]), mesh_.property(vuv_, loop[i+1]), pp);
					if (fabs(bc[0]) < 1.0e-10) bc[0] = 0;  if (fabs(bc[1]) < 1.0e-10) bc[1] = 0;  if (fabs(bc[2]) < 1.0e-10) bc[2] = 0; 

					if ((bc[0] >= 0) && (bc[1] >= 0) && (bc[2] >= 0)) { //要求这个bc的三个分量都是>=0
					mesh_.property(is_face_located, *freev_it) = true;
					TriMesh::FaceHandle fh = mesh_.face_handle(loop_h[i+1]);
					if (fh.is_valid() == false) std::cerr << "Error: face invalid." << std::endl;//Only for assert

					// 下面四行就是每一次参数化所应该获得的或是更改的信息.
					mesh_.property(vf_,  *freev_it) = fh;
					mesh_.property(vf0_, *freev_it) = v_from_old; mesh_.property(vf1_, *freev_it) = loop[i]; mesh_.property(vf2_, *freev_it) = loop[i+1];
					mesh_.property(vbc_, *freev_it) = bc;

					(mesh_.property(fvset_, fh)).push_back(*freev_it);

					TriMesh::Point vfp0 = mesh_.point(v_from_old); TriMesh::Point vfp1 = mesh_.point(loop[i]); TriMesh::Point vfp2 = mesh_.point(loop[i+1]);
					TriMesh::Point vp = vfp0 * bc[0] + vfp1 * bc[1] + vfp2 * bc[2];//被deleted也就是被parameterized的顶点的3d坐标.
					mesh_collapsed_.set_point(*freev_it, vp);////
					}
					} // 这个free vertex落在这个面上
					} // 一部分free vertex 落在这个面上			
					} // 所有的free vertices 落在这n个面上
					for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
					if (mesh_.property(is_face_located, *freev_it) == false) { std::cout << "Error: free_vertices 还没有定位完.\n"; }
					}
					mesh_.remove_property(is_face_located);
					mesh_collapsed_.set_point(v_from_old, mesh_.point(v_from_old));
					free_vertices_size = 0; //重新确认一下是否个数保证
					for (int i = 1; i <= n-1; ++i) {
					free_vertices_size += mesh_.property(fvset_, mesh_.face_handle(loop_h[i])).size();
					}
					if (free_vertices_size != free_vertices.size()) std::cout << "Error: free_vertices's size 在重新定为到新面后就错了.\n";
					// 求边的中点采样
					resample_boundary_edge_midpoint(mesh_, loop_h[0], v_from_old, v_to_old);
					resample_boundary_edge_midpoint(mesh_, mesh_.opposite_halfedge_handle(loop_h[n-1]), loop[n-1], v_from_old);
					for (int i = 1; i <= n-2; ++i) {
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) { // 这个边还没有中点采样, 这里的采样是只许成功的.
					TriMesh::FHandle fhi = mesh_.face_handle(loop_h[i]), fhi_1 = mesh_.face_handle(loop_h[(i+1)%n]);
					for (std::vector<TriMesh::VHandle>::iterator it(mesh_.property(fvset_, fhi).begin()), end(mesh_.property(fvset_, fhi).end()); it != end; ++it) {
					TriMesh::VertexHandle vh = *it; //面fhi上的参数化点
					for (TriMesh::VertexOHalfedgeIter voh_it(mesh2_, vh); voh_it; ++voh_it) {
					TriMesh::VHandle v0 = vh, v1 = mesh2_.to_vertex_handle(voh_it), v2 = mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(voh_it));
					// 这里其实忽略了v2点可能不是在这个参数域里面的.还可能是没有被deleted的或者是参数化到特征边上的情况.
					if (mesh_.property(vf_, v2) == fhi|| mesh_.property(vf_, v2) == fhi_1) {
					OpenMesh::Vec2d pp = (mesh_.property(vuv_, loop[i]) + mesh_.property(vuv_, v_from_old)) * 0.5;//这个边的中点的参数坐标
					OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, v0), mesh_.property(vuv_, v1), mesh_.property(vuv_, v2), pp);
					if (fabs(bc[0]) < 1.0e-10) bc[0] = 0; if (fabs(bc[1]) < 1.0e-10) bc[1] = 0; if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;

					if ((bc[0] >= 0) && (bc[1] >= 0) && (bc[2] >= 0)) { //要求这个bc的三个分量都是>=0
					mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) = true;
					mesh_.property(emp_, mesh_.edge_handle(loop_h[i])) = mesh2_.point(v0) * bc[0] + mesh2_.point(v1) * bc[1] + mesh2_.point(v2) * bc[2];
					mesh_.property(emp_closest_vh_, mesh_.edge_handle(loop_h[i])) = v0;
					} 
					} else { // end of if . v2 is inside the valid region
					std::cout << "Error: free vertex's neighbor v2 isn't in fhi-1.\n"; 
					}
					}
					} // end of for. 循环fhi上的所有参数化点.
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) {
					std::cout << "Error: 边没有中点采样成功.\n";
					}
					}
					}

					} else if (v_r_old.is_valid() == false) { //情况三
					TriMesh::HHandle h0 = mesh_.find_halfedge(v_from_old, v_to_old);	
					if (mesh_.is_boundary(mesh_.opposite_halfedge_handle(h0)) == false) std::cout << "Error: 半边(v_from_old, v_to_old)的对半边应该是边界才对.\n"; 
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
					for (int i=1 ; i<= n-1; ++i) {//求出总的angle
					OpenMesh::Vec3d d1 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).normalize();
					OpenMesh::Vec3d d2 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i-1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).normalize();
					rou_angle += acos(dot(d1, d2));
					vec_angle.push_back(acos(dot(d1, d2)));
					} 
					if (vec_angle.size() != n -1) { std::cout << "Error: vec_angle和loop数组容量不一致.\n"; }
					double angle_scale_ratio = M_PI / rou_angle; //缩放比率
					double temp_sum_angle = 0.0, l = 0;
					for (int i = 1; i <= n-1; ++i) { //设置loop[1...n-1]的参数化坐标
					temp_sum_angle += (vec_angle[i-1] * angle_scale_ratio);
					l = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(v_from_old))).norm(); 
					l *= angle_scale_ratio;
					mesh_.property(vuv_, loop[i]) = OpenMesh::Vec2d(l*cos(temp_sum_angle), l*sin(temp_sum_angle));
					}
					mesh_.property(vuv_, loop[0]) = OpenMesh::Vec2d((mesh_.point(loop[0]) - mesh_.point(v_from_old)).norm() * angle_scale_ratio, 0);
					mesh_.property(vuv_, v_from_old) = OpenMesh::Vec2d(0, 0); //完成一邻域参数化

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
					free_vertices_size -= 1;//减去那个v_from_old
					if (free_vertices.size() != free_vertices_size) std::cout << "Error: free_vertices's size.\n";
					// 判断free_vertices数组的点落在哪个新面中.
					OpenMesh::VPropHandleT<bool> is_face_located;
					mesh_.add_property(is_face_located);
					for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
					mesh_.property(is_face_located, *freev_it) = false;
					}
					for (int i = 0; i <= n-2; ++ i) { 
					for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
					if (mesh_.property(is_face_located, *freev_it) == false) {
					OpenMesh::Vec2d pp = mesh_.property(vuv_, *freev_it);//这个free vertex的参数坐标
					OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, v_from_old), mesh_.property(vuv_, loop[i]), mesh_.property(vuv_, loop[i+1]), pp);
					if (fabs(bc[0]) < 1.0e-10) bc[0] = 0;  if (fabs(bc[1]) < 1.0e-10) bc[1] = 0;  if (fabs(bc[2]) < 1.0e-10) bc[2] = 0; 

					if ((bc[0] >= 0) && (bc[1] >= 0) && (bc[2] >= 0)) { //要求这个bc的三个分量都是>=0
					mesh_.property(is_face_located, *freev_it) = true;
					TriMesh::FaceHandle fh = mesh_.face_handle(loop_h[i]);
					if (fh.is_valid() == false) std::cerr << "Error: face invalid." << std::endl;//Only for assert

					// 下面四行就是每一次参数化所应该获得的或是更改的信息.
					mesh_.property(vf_,  *freev_it) = fh;
					mesh_.property(vf0_, *freev_it) = v_from_old; mesh_.property(vf1_, *freev_it) = loop[i]; mesh_.property(vf2_, *freev_it) = loop[i+1];
					mesh_.property(vbc_, *freev_it) = bc;

					(mesh_.property(fvset_, fh)).push_back(*freev_it);

					TriMesh::Point vfp0 = mesh_.point(v_from_old); TriMesh::Point vfp1 = mesh_.point(loop[i]); TriMesh::Point vfp2 = mesh_.point(loop[i+1]);
					TriMesh::Point vp = vfp0 * bc[0] + vfp1 * bc[1] + vfp2 * bc[2];//被deleted也就是被parameterized的顶点的3d坐标.
					mesh_collapsed_.set_point(*freev_it, vp);////
					}
					} // 这个free vertex落在这个面上
					} // 一部分free vertex 落在这个面上			
					} // 所有的free vertices 落在这n个面上
					for (std::vector<TriMesh::VertexHandle>::iterator freev_it = free_vertices.begin(), freev_end = free_vertices.end(); freev_it != freev_end; ++freev_it) {
					if (mesh_.property(is_face_located, *freev_it) == false) { std::cout << "Error: free_vertices 还没有定位完.\n"; }
					}
					mesh_.remove_property(is_face_located);
					mesh_collapsed_.set_point(v_from_old, mesh_.point(v_from_old));
					free_vertices_size = 0; //重新确认一下是否个数保证
					for (int i = 0; i <= n-2; ++i) {
					free_vertices_size += mesh_.property(fvset_, mesh_.face_handle(loop_h[i])).size();
					}
					if (free_vertices_size != free_vertices.size()) std::cout << "Error: free_vertices's size 在重新定为到新面后就错了.\n";
					// 求边的中点采样
					resample_boundary_edge_midpoint(mesh_, mesh_.opposite_halfedge_handle(loop_h[0]), v_to_old, v_from_old);
					resample_boundary_edge_midpoint(mesh_, loop_h[n-1], v_from_old, loop[n-1]);
					for (int i = 1; i <= n-2; ++i) {
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) { // 这个边还没有中点采样, 这里的采样是只许成功的.
					TriMesh::FHandle fhi = mesh_.face_handle(loop_h[i]), fhi_1 = mesh_.face_handle(loop_h[i-1]);
					for (std::vector<TriMesh::VHandle>::iterator it(mesh_.property(fvset_, fhi).begin()), end(mesh_.property(fvset_, fhi).end()); it != end; ++it) {
					TriMesh::VertexHandle vh = *it; //面fhi上的参数化点
					for (TriMesh::VertexOHalfedgeIter voh_it(mesh2_, vh); voh_it; ++voh_it) {
					TriMesh::VHandle v0 = vh, v1 = mesh2_.to_vertex_handle(voh_it), v2 = mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(voh_it));
					// 这里其实忽略了v2点可能不是在这个参数域里面的.还可能是没有被deleted的或者是参数化到特征边上的情况.
					if (mesh_.property(vf_, v2) == fhi|| mesh_.property(vf_, v2) == fhi_1) {
					OpenMesh::Vec2d pp = (mesh_.property(vuv_, loop[i]) + mesh_.property(vuv_, v_from_old)) * 0.5;//这个边的中点的参数坐标
					OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, v0), mesh_.property(vuv_, v1), mesh_.property(vuv_, v2), pp);
					if (fabs(bc[0]) < 1.0e-10) bc[0] = 0; if (fabs(bc[1]) < 1.0e-10) bc[1] = 0; if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;

					if ((bc[0] >= 0) && (bc[1] >= 0) && (bc[2] >= 0)) { //要求这个bc的三个分量都是>=0
					mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) = true;
					mesh_.property(emp_, mesh_.edge_handle(loop_h[i])) = mesh2_.point(v0) * bc[0] + mesh2_.point(v1) * bc[1] + mesh2_.point(v2) * bc[2];
					mesh_.property(emp_closest_vh_, mesh_.edge_handle(loop_h[i])) = v0;
					} 
					} else { // end of if . v2 is inside the valid region
					std::cout << "Error: free vertex's neighbor v2 isn't in fhi-1.\n"; 
					}
					}
					} // end of for. 循环fhi上的所有参数化点.
					if (mesh_.property(empl_, mesh_.edge_handle(loop_h[i])) == false) {
					std::cout << "Error: 边没有中点采样成功.\n";
					}
					}
					}			 
					}*/
		} // end of if (mesh_.is_boundary(v_from_old).v_to_old -> (v_from_old, v_to_old) 的分内部和边界
		// 边界情况这里好像不能处理, 至少少了:simplified_mesh_ and mesh_分裂新形成的面之间还没有对应上.这在中点细分时候用到.

		midpoint_sample_succeed_ = true;//假设新边都中点采样成功了
		for (TriMesh::VertexEdgeIter ve_it(mesh_, v_from_old); ve_it; ++ve_it) { // 搜索v_from_old的一邻域是否的却都中点采样成功了.
			if (mesh_.property(empl_, ve_it.handle()) == false) midpoint_sample_succeed_ = false;
		}
	} else{ //end of if lc_node_handle != -1 可分裂
		std::cout << "Info: vsplit_refine(v_to): lc_node_handle == -1. 已经上叶子节点了.\n";
	}
	return true;
} // end of function vsplit_refine.
bool DecimationModel::error_driven_selective_refine() {
	std::cout << "进入error_driven_selective_refine().\n";
	// 方法一 从icm入手来找分裂的顶点. icm经过一次loop细分之后, 所有顶点都找在原始模型上的最最近点.
	// 但这方法应该不好. 应该从原始模型入手找哪里误差大，然后对simplified_mesh_做顶点分类.
	find_closestp_icm();

	double largest_err = 0;
	TriMesh::VertexHandle vh_largest_err; 
	int i = 0;
	for (TriMesh::VIter v_it(initial_control_mesh_.vertices_begin()), v_end(initial_control_mesh_.vertices_end()); v_it != v_end; ++v_it) {
		// initial_control_mesh_经过一次loop细分之后和refined_simplified_mesh_的顶点对应.
		if (v_it.handle() != refined_simplified_mesh_.vertex_handle(i)) std::cout << "Error: 这两个mesh的index应该是match的" <<  ".\n"; 
		++ i;
		if (initial_control_mesh_.property(vp_distances_,v_it) > largest_err) {
			largest_err = initial_control_mesh_.property(vp_distances_,v_it);
			vh_largest_err =v_it.handle();
		}
	}
	
	
	std::cout << "largest fitting error: " << largest_err << ".\n";
	std::cout << "simplified mesh's vertices num: " << simplified_mesh_.n_vertices() << ", vh_largest_err: " << vh_largest_err.idx() << ".\n";

	// --------------------------------------------------------
	// 在simplified mesh上找可能被分裂的顶点.
	TriMesh::VertexHandle v0, v1, v2, v3;
	if (vh_largest_err.idx() >= simplified_mesh_.n_vertices()) { //这时vh_largest_err点是在中点细分时候生成的, 我要的是simplified_mesh_上的点
		std::vector<TriMesh::VertexHandle> vec;	//找vh_largest_err的一邻域顶点, 其中有且仅有两个是在simplified_mesh_上的点.
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
	} else { // 保证v0总是在simplified_mesh_上的点
		v0 = vh_largest_err;
	} 
	// -----------------------------------------------------------------------------
	TriMesh::VertexHandle v_to = v0;// v_to顶点对应的v_to_node_handle这个节点Node做分裂.
	if (v1.is_valid()) { // 还要在v0, v1中选择出用哪一个做分裂
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
			TriMesh::VHandle lc_old_vh = voldhandles_[lc_point_handle];    // v0分裂的话就分裂出这个点.看这个点时候在左右两个面上.
			if (mesh_.status(lc_old_vh).deleted() == false) std::cout << "Error: 这时的v0的左节点在mesh_上是被deleted才对的.\n"; 
			if (mesh_.property(vf_, lc_old_vh) == f0_old|| mesh_.property(vf_, lc_old_vh) == f1_old) v_to = v0;
		}
		if (v1_lc != -1) {
			TriMesh::VHandle lc_old_vh = voldhandles_[vhierarchy_.node(v1_lc).point_handle()];
			if (mesh_.status(lc_old_vh).deleted() == false) std::cout << "Error: 这时的vl的左节点在mesh_上是被deleted才对的.\n";
			if (mesh_.property(vf_, lc_old_vh) == f0_old|| mesh_.property(vf_, lc_old_vh) == f1_old) v_to = v1;
		}
		if (v0_lc == -1) v_to = v1; // v0已经是叶子不可分裂了, 那就选别的点吧.
	}

	// ----------------------------
	// 测试已经选择的v_to是否是叶子或者分裂会造成翻转.
	int n_change = 0;
	while (true) {
		if (test_is_leaf_of_vertex_hierarchy(v_to) || test_is_triangle_flipping_after_vsplit(v_to)) {
			// 表示这个v_to要么已经是叶子了, // 表示分裂会造成翻转, 所以要另觅v_to顶点了.
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
	vsplit_refine(v_to); //顶点分裂, 边中点采用, 更新参数化信息.
	// 
	std::cout << "离开error_driven_selective_refine().\n";
	return true;
}
void DecimationModel::error_driven_selective_refine_process() {
	std::cout << "\n";
	std::cout << "========开始误差驱动的refine============\n";//这种提示用40个长度的=,前面第9个位置开始文字,即最先有8个=.
	if (error_driven_selective_refine() == false) std::cout << "Error: error driven refine.\n";
	set_mesh_toberendered(&simplified_mesh_, SIMPLIFIED);
	std::cout << "现在的simplified_mesh_: v " << simplified_mesh_.n_vertices() << ", e " << simplified_mesh_.n_edges() << ", f " << simplified_mesh_.n_faces() << ".\n";
	std::cout << "========结束误差驱动的refine============\n";//这种提示用40个长度的=,前面第9个位置开始文字,即最先有8个=.
}
// ---------------
void DecimationModel::refine_process() {
	std::cout << "\n";
	std::cout << "========开始误差驱动的refine============\n";//这种提示用40个长度的=,前面第9个位置开始文字,即最先有8个=.
	// 2009-05-23, 利用simplified mesh上每一个顶点的误差来确定分裂.
	int n_split = 0; // for test. 每次分裂10个点
	//for (TriMesh::VIter v_it(simplified_mesh_.vertices_begin()), v_end(simplified_mesh_.vertices_end()); v_it != v_end; ++v_it) {
	//	if (simplified_mesh_.property(vp_fitting_error_sm_, v_it) > standard_error_sm_) {
	//		TriMesh::VHandle v_to = v_it.handle();
	//		if (test_is_leaf_of_vertex_hierarchy(v_to) || test_is_triangle_flipping_after_vsplit(v_to)) {
	//			continue;
	//		} else {
	//			vsplit_refine(v_to); //顶点分裂, 边中点采用, 更新参数化信息.
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
			vsplit_refine(vh); //顶点分裂, 边中点采用, 更新参数化信息.
			++n_split;
		}
	}
	std::cout << "n_split: " << n_split << ".\n";

	set_mesh_toberendered(&simplified_mesh_, RESAMPLE_MIDPOINT);
	std::cout << "现在的simplified_mesh_: v " << simplified_mesh_.n_vertices() << ", e " << simplified_mesh_.n_edges() << ", f " << simplified_mesh_.n_faces() << ".\n";
	std::cout << "========结束误差驱动的refine============\n";//这种提示用40个长度的=,前面第9个位置开始文字,即最先有8个=.

}
