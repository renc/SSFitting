
#include "DecimationModel.h" 
//#include "../../../src/CourseExamples/04-Fairing/TaucsSolver.hh"

#include <algorithm> // for max_element algorithm
#include <assert.h>
#include <fstream>


void DecimationModel::resample_original_feature_edge_midpoint()//��ԭʼ�ߺ������߲���
{
	for (TriMesh::EdgeIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end()); e_it != e_end; ++ e_it) 
	{
		if (mesh_.status(e_it).deleted() == false && mesh_.property(empl_, e_it) == false) { // �������ϵı�.// ����ߵ��е㻹û�ж�λ.

			// (1)������ڼ򻯵Ĺ��̵���û�з����仯,�������˵㻹�Ǽ�֮���������.
			TriMesh::VHandle v0 = mesh_.to_vertex_handle(mesh_.halfedge_handle(e_it, 0));
			TriMesh::VHandle v1 = mesh_.to_vertex_handle(mesh_.halfedge_handle(e_it, 1));
			if (((v0 == mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 0))) && (v1 == mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 1))))
			|| ((v1 == mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 0))) && (v0 == mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 1))))) 
			{
					mesh_.property(emp_, e_it.handle()) = (mesh_.point(v0) + mesh_.point(v1)) / 2.0;
					mesh_.property(empl_, e_it) = true;
					mesh_.property(emp_closest_vh_, e_it) = v0;
			}
			// (2)�������ߵ��е�
			if (mesh_.status(e_it).feature() && mesh_.property(empl_, e_it) == false) {
				TriMesh::HHandle h0 = mesh_.halfedge_handle(e_it, 0);
				TriMesh::HHandle h1 = mesh_.halfedge_handle(e_it, 1);
				unsigned int n0 = mesh_.property(hep_heset_, h0).size();
				unsigned int n1 = mesh_.property(hep_heset_, h1).size();
				unsigned int i0 = 0, i1 = 0;
				double sum_len = 0, tmp = 0;
				for (i0 = 0, i1 = 0; i0 < n0; ++i0, ++i1) {
					sum_len += 
						(mesh2_.point(mesh2_.from_vertex_handle(mesh_.property(hep_heset_, h0)[i0]))
						- mesh2_.point(mesh2_.to_vertex_handle(mesh_.property(hep_heset_, h0)[i0]))).norm();
				}
				double half_len = sum_len / 2.0;
				for (i0 = 0; i0 < n0; ++i0) {// �������crease edge���е�
					tmp += 
						(mesh2_.point(mesh2_.from_vertex_handle(mesh_.property(hep_heset_, h0)[i0])) - mesh2_.point(mesh2_.to_vertex_handle(mesh_.property(hep_heset_, h0)[i0]))).norm();
					if (tmp >= half_len) {
						mesh_.property(emp_, e_it.handle()) 
							= mesh2_.point(mesh2_.to_vertex_handle(mesh_.property(hep_heset_, h0)[i0]))
							+ (mesh2_.point(mesh2_.from_vertex_handle(mesh_.property(hep_heset_, h0)[i0])) 
							- mesh2_.point(mesh2_.to_vertex_handle(mesh_.property(hep_heset_, h0)[i0]))).normalize() * (tmp - half_len);
						mesh_.property(empl_, e_it) = true;

						mesh_.property(emp_closest_vh_, e_it) = mesh2_.to_vertex_handle(mesh_.property(hep_heset_, h0)[i0]);
						break;
					}
				}
			} // end of if. feature edge 
		} // end of if edge not deleted
	} // end of for each edge 
	// ����Ƿ����е�crease edge�������ɹ���, ����������������. ��һ�㲻�ᷢ����, ֻ��Ϊ��ȷ��.
	for (TriMesh::EdgeIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end()); e_it != e_end; ++ e_it) {
		if (mesh_.status(e_it).deleted() == false && mesh_.status(e_it).feature()) {
			if (mesh_.property(empl_, e_it) == false) std::cout << "Error: resample_original_feature_edge_midpoint()�л�û����������ߵ��е������.\n";
		}
	} 
}
void DecimationModel::resample_boundary_edge_midpoint(TriMesh &_mesh, TriMesh::HHandle _hh, TriMesh::VHandle _vstart, TriMesh::VHandle _vend) //�Ա߽�߲���
{   // ��ʵ�����������Ǿ���_hh�Ķ԰�ߵ�from��to����
	if (_mesh.is_boundary(_hh) == false) std::cout << "Error: resample_boundary_edge_midpoint ��߷Ǳ߽�.\n";
	if (_mesh.property(empl_, _mesh.edge_handle(_hh)) == false) {
		TriMesh::VHandle vstart = _vstart, vend = _vend; //�߽��ʼĩ����
		TriMesh::HalfedgeHandle boundaryhe;//
		std::vector<TriMesh::VertexHandle> loop;
		double sum_length(0);
		for (TriMesh::VertexOHalfedgeIter vsoh(mesh2_, vstart); vsoh; ++vsoh) {
			if (mesh2_.is_boundary(vsoh)) {
				boundaryhe = vsoh.handle();
				break;//һ������ֻ��һ��������Ǳ߽��, ��ȻҲֻ��һ�������Ǳ߽��.
			}
		} 
		loop.push_back(vstart);
		while (mesh2_.to_vertex_handle(boundaryhe) != vend) {
			loop.push_back(mesh2_.to_vertex_handle(boundaryhe));
			sum_length += (mesh2_.point(mesh2_.from_vertex_handle(boundaryhe)) - mesh2_.point(mesh2_.to_vertex_handle(boundaryhe))).norm();
			boundaryhe = mesh2_.next_halfedge_handle(boundaryhe);
		}
		loop.push_back(vend);//���Ķ���ͱ߳���û�м���
		sum_length += (mesh2_.point(loop[loop.size() - 1] ) - mesh2_.point(loop[loop.size() - 2] )).norm();
		// std::cout << "sum_length: " <<   sum_length << std::endl;//for test
		double half_length = sum_length / 2.0, tmp(0);
		for (int loop_i = 0, loop_n = loop.size(); loop_i <= loop_n-2; loop_i++) {
			tmp += (mesh2_.point(loop[loop_i]) - mesh2_.point(loop[loop_i + 1])).norm();
			if (tmp >= half_length) { 
				_mesh.property(emp_, _mesh.edge_handle(_hh))
					= mesh2_.point(loop[loop_i + 1]) + (mesh2_.point(loop[loop_i]) - mesh2_.point(loop[loop_i + 1])).normalize() * (tmp - half_length);
				_mesh.property(empl_, _mesh.edge_handle(_hh)) = true;
				_mesh.property(emp_closest_vh_, _mesh.edge_handle(_hh)) = loop[loop_i + 1];
				//std::cout << "h0: " <<   _mesh.property(emp_, _mesh.edge_handle(h0)) << std::endl;
				break; 
			}
		}
	} 
}

void DecimationModel::flatten_face_resample_3edge_midpoint(TriMesh::FaceHandle _fh) {
	TriMesh::FHandle f_it = _fh; //��f_it�����,ֻ��Ϊ��ͳһ�ҿ�ʼʱ��ʹ�õ����ֶ���.
	if (mesh_.status(f_it).deleted() == false) { //ѭ��ÿһ��û�б�deleted����,Ҳ����base complex�ϵ���,domain triangle
		// �����Ƕ���f_it����.
		//std::cout << "In face " << f_it << std::endl; // for test
		//////////////////////////////////////////////////////////////////////////
		// (1) ���������� //����һ��������, ��4���ȱ����������.
		TriMesh::FHIter fh_it(mesh_, f_it);
		TriMesh::HalfedgeHandle h0 = fh_it.handle(); ++fh_it;
		TriMesh::HalfedgeHandle h1 = fh_it.handle(); ++fh_it;
		TriMesh::HalfedgeHandle h2 = fh_it.handle();
		TriMesh::VHandle v0(mesh_.to_vertex_handle(h0)), v1(mesh_.to_vertex_handle(h1)), v2(mesh_.to_vertex_handle(h2));
		//std::cout << mesh2_.point(v0) << "; " << mesh2_.point(v1) << "; " << mesh2_.point(v2) << ".\n"; 
		mesh_.property(vuv_, v0) = OpenMesh::Vec2d(1, 0);//v0, v1, v2�����ǲ��ᷢ���ص�����˿������������������.
		mesh_.property(vuv_, v1) = OpenMesh::Vec2d(1.5, 0.8660254);
		mesh_.property(vuv_, v2) = OpenMesh::Vec2d(0.5, 0.8660254);

		TriMesh::FaceHandle f0(-1), f1(-1), f2(-1);
		TriMesh::VHandle    v3, v4, v5;//(-1)
		if (mesh_.is_boundary(mesh_.opposite_halfedge_handle(h0)) == false) {
			f0 = mesh_.face_handle(mesh_.opposite_halfedge_handle(h0));
			v3 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(h0))); 
		} else { //Ӧ�ÿ������h0���ڵı߽�ߵ��е���ԭʼ�����еĲ�����
			resample_boundary_edge_midpoint(mesh_, mesh_.opposite_halfedge_handle(h0), v0, v2); 					
		}
		if (mesh_.is_boundary(mesh_.opposite_halfedge_handle(h1)) == false) {
			f1 = mesh_.face_handle(mesh_.opposite_halfedge_handle(h1));
			v4 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(h1)));	
		} else { //Ӧ�ÿ������h1���ڵı߽�ߵ��е���ԭʼ�����еĲ�����
			resample_boundary_edge_midpoint(mesh_, mesh_.opposite_halfedge_handle(h1), v1, v0);				
		}
		if (mesh_.is_boundary(mesh_.opposite_halfedge_handle(h2)) == false) {//h2�Ķ԰�߲��Ǳ߽�
			f2 = mesh_.face_handle(mesh_.opposite_halfedge_handle(h2));
			v5 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(h2)));
		} else { // h2�Ķ԰���Ǳ߽� //Ӧ�ÿ������h2���ڵı߽�ߵ��е���ԭʼ�����еĲ�����
			resample_boundary_edge_midpoint(mesh_, mesh_.opposite_halfedge_handle(h2), v2, v1);
		} //��֮���f0,f1,f2��Ȼ������.is_valid() == fasle��, ���Ǳ߽���Ľ��.
		TriMesh::FHandle f3, f4, f5, f6, f7, f8;//�������β����������extension
		TriMesh::HHandle h3, h4, h5, h6, h7, h8;
		TriMesh::VHandle v6, v7, v8, v9, v10, v11;
		if (f0.is_valid() == true) {
			TriMesh::HHandle oh0 = mesh_.opposite_halfedge_handle(h0);
			h3 = mesh_.opposite_halfedge_handle(mesh_.next_halfedge_handle(oh0)); h4 = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(oh0));
			if (mesh_.is_boundary(h3) == false) {
				f3 = mesh_.face_handle(h3);		v6 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h3));					
			}
			if (mesh_.is_boundary(h4) == false) {
				f4 = mesh_.face_handle(h4);		v7 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h4));
			}							
		}
		if (f1.is_valid() == true) {
			TriMesh::HHandle oh1 = mesh_.opposite_halfedge_handle(h1);
			h5 = mesh_.opposite_halfedge_handle(mesh_.next_halfedge_handle(oh1)); h6 = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(oh1));
			if (mesh_.is_boundary(h5) == false) {
				f5 = mesh_.face_handle(h5);		v8 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h5));								
			}
			if (mesh_.is_boundary(h6) == false) {
				f6 = mesh_.face_handle(h6);		v9 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h6));	
			} 
		}
		if (f2.is_valid() == true) {
			TriMesh::HHandle oh2 = mesh_.opposite_halfedge_handle(h2);
			h7 = mesh_.opposite_halfedge_handle(mesh_.next_halfedge_handle(oh2)); h8 = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(oh2));
			if (mesh_.is_boundary(h7) == false) {
				f7 = mesh_.face_handle(h7);		v10 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h7));							
			}
			if (mesh_.is_boundary(h8) == false) {
				f8 = mesh_.face_handle(h8);		v11 = mesh_.to_vertex_handle(mesh_.next_halfedge_handle(h8));	
			} 
		}
		// ���(1)�Բ�����Ľ���.
		//////////////////////////////////////////////////////////////////////////
		// (2)
		//������������������ϵ����ж���, �ر����Щ������Ҳ�������������ĳ����.
		std::vector<TriMesh::VertexHandle> fvset(mesh_.property(fvset_, f_it));
		//һЩ��������������������������, ����ȴ��������. �������⴦��. ��Щ�������dart, crease, corner
		if (mesh_.status(mesh_.edge_handle(h0)).feature()) {
			for (int i = 0; i < mesh_.property(hep_heset_, h0).size(); ++i) { 
				fvset.push_back(mesh2_.to_vertex_handle(mesh_.property(hep_heset_, h0)[i]));
			}
		}
		if (mesh_.status(mesh_.edge_handle(h1)).feature()) {
			for (int i = 0; i < mesh_.property(hep_heset_, h1).size(); ++i) {
				fvset.push_back(mesh2_.to_vertex_handle(mesh_.property(hep_heset_, h1)[i]));
			}
		}
		if (mesh_.status(mesh_.edge_handle(h2)).feature()) {
			for (int i = 0; i < mesh_.property(hep_heset_, h2).size(); ++i) {
				fvset.push_back(mesh2_.to_vertex_handle(mesh_.property(hep_heset_, h2)[i]));
			}
		} //ע��: mesh2_.to_vertex_handle(mesh_.property(hep_heset_, h2).size()-1)����ܿ��ܾ��Ǽ������ϵĵ�, û�б�deleted��.
		//
		//std::cout << ".->: " << v3 << " " << v4 << " " << v5 << ".\n";// չ��������for test
		for (std::vector<TriMesh::VertexHandle>::iterator fvset_it(fvset.begin()), it_end(fvset.end()); fvset_it != it_end; ++fvset_it) {
			//ÿ�� ����һ���������������ϵĶ���*fvset_it.
			////std::cout << "��δ���Ķ�����*fvset_it: " << *fvset_it << ".\n"; //for test
			if (mesh_.status(*fvset_it).deleted()) {
				if (mesh_.property(vp_type_, *fvset_it) == DGP::SMOOTH_VFT) { //��crease vertex����smooth��, ����������������, ��������Щcrease v��������������
					if (mesh_.property(vf_, *fvset_it) != f_it) { std::cout << "Error, para check, faces not match.\n";}// only to confirm
					mesh_.property(vuv_, *fvset_it) = mesh_.property(vuv_, mesh_.property(vf0_, *fvset_it)) * (mesh_.property(vbc_, *fvset_it))[0]
					+ mesh_.property(vuv_, mesh_.property(vf1_, *fvset_it)) * (mesh_.property(vbc_, *fvset_it))[1]// 2008.02.22�ŷ������bug,
					+ mesh_.property(vuv_, mesh_.property(vf2_, *fvset_it)) * (mesh_.property(vbc_, *fvset_it))[2];//��ʼʱ��Ȼû�������!

				} else if (mesh_.property(vp_type_, *fvset_it) == DGP::CREASE_VFT) { // Ҫ�ǲ���������h0�ϵĵ�����䵽f0����, ����䵽f1,f2��Ҳ����
					// crease ����������, һ������f_it���ϵ�, ����һ��������f_it����������ϵ�, ������������
					if (mesh_.property(vf_, *fvset_it) == f_it) { // v0, v1, v2������������
					} else if (mesh_.property(vf_, *fvset_it) == f0) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0);
					} else if (mesh_.property(vf_, *fvset_it) == f1) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0);
					} else if (mesh_.property(vf_, *fvset_it) == f2) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508);
					} else { 
						std::cout << "Error: ���ó�����������.\n";std::cout << "In face " << f_it << ", fvset_'s size " << fvset.size(); // for test
						std::cout << ".->: " << f0 << " " << f1 << " " << f2 << ".\n";// չ��������for test
						std::cout << mesh_.property(vf_, *fvset_it) << ".\n";
					}
					mesh_.property(vuv_, *fvset_it) = mesh_.property(vuv_, mesh_.property(vf0_, *fvset_it)) * (mesh_.property(vbc_, *fvset_it))[0]
					+ mesh_.property(vuv_, mesh_.property(vf1_, *fvset_it)) * (mesh_.property(vbc_, *fvset_it))[1]
					+ mesh_.property(vuv_, mesh_.property(vf2_, *fvset_it)) * (mesh_.property(vbc_, *fvset_it))[2];
				} 
			} else {// �����������f_it����������ϵĵ���뵽fvset�����е����, ��������ĩ�Ǹ�������ϼ������ϵĵ�.
				// DGP::DART_VFT, DGP::CREASE_VFT, DGP::CORNER_VFT ��ֻ����v0, v1, v2, ����������Ѿ�������.
			}

			/**/
			// �ж�vj����Ч��
			int kind = 0; 
			bool vj_in_valid_region = true;
			for (TriMesh::VertexOHalfedgeIter voh_it(mesh2_, *fvset_it); voh_it; ++voh_it) {
				TriMesh::VHandle vj = mesh2_.to_vertex_handle(voh_it.handle()); 
				//�ھӵ�vjҪô�Ǳ�deleted�˲�����4������, Ҫô��û�б�deleted����������������, ����������Ч��.
				//std::cout << "vj " << vj << ", ";//for test
				if (mesh_.status(vj).deleted()) { //�ڽӵ�vj�Ǳ�deleted��Ҳ�����б���������base complex��ĳһ�����ϵ�.
					TriMesh::FaceHandle pf = mesh_.property(vf_, vj);//�ڵ�vj��initial parameterize����������������
					if (pf == f_it) {														
					} else if (pf == f0) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0);
					} else if (pf == f1) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0);
					} else if (pf == f2) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); 
					} else if (pf == f3) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0); mesh_.property(vuv_, v6) = OpenMesh::Vec2d(-0.5, 0.8660254);
					} else if (pf == f4) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0); mesh_.property(vuv_, v7) = OpenMesh::Vec2d(0.5, -0.8660254);
					} else if (pf == f5) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0); mesh_.property(vuv_, v8) = OpenMesh::Vec2d(1.5, -0.8660254);
					} else if (pf == f6) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0); mesh_.property(vuv_, v9) = OpenMesh::Vec2d(2.5, 0.8660254);
					} else if (pf == f7) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); mesh_.property(vuv_, v10) = OpenMesh::Vec2d(2, 1.7320508);
					} else if (pf == f8) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); mesh_.property(vuv_, v11) = OpenMesh::Vec2d(0, 1.7320508);
					} else vj_in_valid_region = false; 
					if (vj_in_valid_region == true) {
						mesh_.property(vuv_, vj) = mesh_.property(vuv_, mesh_.property(vf0_, vj)) * (mesh_.property(vbc_, vj))[0]
						+ mesh_.property(vuv_, mesh_.property(vf1_, vj)) * (mesh_.property(vbc_, vj))[1]
						+ mesh_.property(vuv_, mesh_.property(vf2_, vj)) * (mesh_.property(vbc_, vj))[2];
					}

				} else { //�ڽӵ�vjû�б�deleted��Ҳ����������base complex�ϵ�ĳһ������.
					if ((vj == v0) || (vj == v1) || (vj == v2)) { //��vj�͵�������������Ļ���ô�����������Ѿ����������õ���.
					} else if (vj == v3) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0);
					} else if (vj == v4) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0);
					} else if (vj == v5) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); 
					} else if (vj == v6) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0); mesh_.property(vuv_, v6) = OpenMesh::Vec2d(-0.5, 0.8660254); kind = 3; 
					} else if (vj == v7) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0); mesh_.property(vuv_, v7) = OpenMesh::Vec2d(0.5, -0.8660254); kind = 4; 					
					} else if (vj == v8) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0); mesh_.property(vuv_, v8) = OpenMesh::Vec2d(1.5, -0.8660254); kind = 5;
					} else if (vj == v9) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0); mesh_.property(vuv_, v9) = OpenMesh::Vec2d(2.5, 0.8660254);  kind = 6;
					} else if (vj == v10) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); mesh_.property(vuv_, v10) = OpenMesh::Vec2d(2, 1.7320508); kind = 7;
					} else if (vj == v11) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); mesh_.property(vuv_, v11) = OpenMesh::Vec2d(0, 1.7320508); kind = 8;
					} else { vj_in_valid_region = false; }
				}

				TriMesh::VertexHandle vh0 = *fvset_it;//triangle(vh0, vh1, vh2).
				TriMesh::VertexHandle vh1 = vj;
				TriMesh::VertexHandle vh2;
				if (mesh_.is_boundary(mesh2_.next_halfedge_handle(voh_it.handle())) == false) {
					vh2 = mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(voh_it.handle()));// ������������ȷ��ǰ���ǷǱ߽����
					//TriMesh::VertexHandle vh3 = mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(mesh2_.opposite_halfedge_handle(voh_it.handle())));
				} else continue;

				// �ж�������(vi, vj, vh2)��(vh0, vh1, vh2)
				bool vh2_in_valid_region = true;//�ȼ�������������,��Ч
				if (mesh_.status(vh2).deleted()) {
					TriMesh::FHandle vh2f = mesh_.property(vf_, vh2);
					if (vh2f == f_it) { 
					} else if (vh2f == f0) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0);
					} else if (vh2f == f1) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0);
					} else if (vh2f == f2) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); 
					} else if (vh2f == f3) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0); mesh_.property(vuv_, v6) = OpenMesh::Vec2d(-0.5, 0.8660254);
					} else if (vh2f == f4) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0); mesh_.property(vuv_, v7) = OpenMesh::Vec2d(0.5, -0.8660254);
					} else if (vh2f == f5) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0); mesh_.property(vuv_, v8) = OpenMesh::Vec2d(1.5, -0.8660254);
					} else if (vh2f == f6) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0); mesh_.property(vuv_, v9) = OpenMesh::Vec2d(2.5, 0.8660254);
					} else if (vh2f == f7) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); mesh_.property(vuv_, v10) = OpenMesh::Vec2d(2, 1.7320508);
					} else if (vh2f == f8) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); mesh_.property(vuv_, v11) = OpenMesh::Vec2d(0, 1.7320508);
					} else { vh2_in_valid_region = false;
					}
					if (vh2_in_valid_region == true) {
						mesh_.property(vuv_, vh2) = mesh_.property(vuv_, mesh_.property(vf0_, vh2)) * (mesh_.property(vbc_, vh2))[0]
						+ mesh_.property(vuv_, mesh_.property(vf1_, vh2)) * (mesh_.property(vbc_, vh2))[1]
						+ mesh_.property(vuv_, mesh_.property(vf2_, vh2)) * (mesh_.property(vbc_, vh2))[2];	 
					}								
				} else { // vh2 is deleted, 
					if ((vh2 == v0) || (vh2 == v1) || (vh2 == v2)) { 
					} else if (vh2 == v3) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0);
					} else if (vh2 == v4) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0);
					} else if (vh2 == v5) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); 
					} else if (vh2 == v6) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0); mesh_.property(vuv_, v6) = OpenMesh::Vec2d(-0.5, 0.8660254); kind = 3; 
					} else if (vh2 == v7) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0); mesh_.property(vuv_, v7) = OpenMesh::Vec2d(0.5, -0.8660254); kind = 4; 					
					} else if (vh2 == v8) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0); mesh_.property(vuv_, v8) = OpenMesh::Vec2d(1.5, -0.8660254); kind = 5;
					} else if (vh2 == v9) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0); mesh_.property(vuv_, v9) = OpenMesh::Vec2d(2.5, 0.8660254);  kind = 6;
					} else if (vh2 == v10) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); mesh_.property(vuv_, v10) = OpenMesh::Vec2d(2, 1.7320508); kind = 7;
					} else if (vh2 == v11) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); mesh_.property(vuv_, v11) = OpenMesh::Vec2d(0, 1.7320508); kind = 8;
					} else { vh2_in_valid_region = false; }
				} 
				//std::cout << "vh0 " << vh0 << ", " << mesh_.property(vp_type_, vh0) << ", " << "vj " << vj_in_valid_region << ", vh2 " << vh2_in_valid_region << ".\n";
				// ���һ,vj(vh1)���м���ĸ��������� and vh2Ҳ�����м���ĸ���������
				std::vector<TriMesh::HHandle > vec_h;
				vec_h.push_back(h0); vec_h.push_back(h1); vec_h.push_back(h2); 
				vec_h.push_back(h3); vec_h.push_back(h4); vec_h.push_back(h5); vec_h.push_back(h6); vec_h.push_back(h7); vec_h.push_back(h8); 
				std::vector<OpenMesh::Vec2d> vec_mp;
				vec_mp.push_back(OpenMesh::Vec2d(0.75, 0.4330127)); vec_mp.push_back(OpenMesh::Vec2d(1.25, 0.4330127)); vec_mp.push_back(OpenMesh::Vec2d(1, 0.8660254)); 
				vec_mp.push_back(OpenMesh::Vec2d(0.25, 0.4330127)); vec_mp.push_back(OpenMesh::Vec2d(0.5, 0)); 
				vec_mp.push_back(OpenMesh::Vec2d(1.5, 0)); vec_mp.push_back(OpenMesh::Vec2d(1.75, 0.4330127)); 
				vec_mp.push_back(OpenMesh::Vec2d(1.25, 1.2990381)); vec_mp.push_back(OpenMesh::Vec2d(0.75, 1.2990381)); 
				if (vj_in_valid_region) { //vi������һ��vj��Ч//						
					if (vh2_in_valid_region == true) {//vh2������Ч����֮��
						//std::cout << "vh1 valid, vh2 valid. \t";//for test	
						for (int i = 0; i < 9; i ++) { //�����ߵ��е���Ҫ����
							if (i == 0) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0);
							} else if (i == 1) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0);
							} else if (i == 2) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); 
							} else if (i == 3) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0); mesh_.property(vuv_, v6) = OpenMesh::Vec2d(-0.5, 0.8660254);  
							} else if (i == 4) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0); mesh_.property(vuv_, v7) = OpenMesh::Vec2d(0.5, -0.8660254); 
							} else if (i == 5) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0); mesh_.property(vuv_, v8) = OpenMesh::Vec2d(1.5, -0.8660254);
							} else if (i == 6) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0); mesh_.property(vuv_, v9) = OpenMesh::Vec2d(2.5, 0.8660254);
							} else if (i == 7) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); mesh_.property(vuv_, v10) = OpenMesh::Vec2d(2, 1.7320508);				
							} else if (i == 8) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); mesh_.property(vuv_, v11) = OpenMesh::Vec2d(0, 1.7320508); 
							} 
							TriMesh::HHandle hx = vec_h[i]; 
							//if (test_ && i ==0) std::cout << mesh_.property(empl_, mesh_.edge_handle(hx)) << "; " << ;
							if (mesh_.property(empl_, mesh_.edge_handle(hx)) == false) {
								OpenMesh::Vec2d mpx = vec_mp[i];	
								OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, vh0), mesh_.property(vuv_, vh1), mesh_.property(vuv_, vh2), mpx);
								if (fabs(bc[0]) < 1.0e-10) bc[0] = 0; if (fabs(bc[1]) < 1.0e-10) bc[1] = 0; if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;
								if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //Ҫ�����bc��������������>=0
									mesh_.property(emp_, mesh_.edge_handle(hx)) = mesh2_.point(vh0) * bc[0] + mesh2_.point(vh1) * bc[1] + mesh2_.point(vh2) * bc[2];
									mesh_.property(empl_, mesh_.edge_handle(hx)) = true; //std::cout << "i = " << i << ", " << fvset.size() << ".\n";
									mesh_.property(emp_closest_vh_, mesh_.edge_handle(hx)) = vh0;

									//std::cout << ": " <<  v0 << ", " << v1 << ", " << v2 << ".\n";
									//std::cout << mesh2_.point(vh0) << "; " << mesh2_.point(vh1) << "; " << mesh2_.point(vh2) << "; " << bc << ";" << ".\n";
									//std::cout << mesh_.property(vuv_, vh0) << "; " << mesh_.property(vuv_, vh1) << "; " << mesh_.property(vuv_, vh2) << "; " << mpx << ";" << ".\n";//

								}/*
								 if (test_ && i ==0) { 
								 if (mesh_.status(*fvset_it).deleted() && mesh_.property(vf_, *fvset_it) == f_it && (mesh_.property(vf_, vh1) == f0 || mesh_.property(vf_, vh2) == f0)) {
								 std::cout << ": " <<  mesh_.property(vuv_, vh0) << ", " << mesh_.property(vuv_, vh1) << ", " << mesh_.property(vuv_, vh2) << ".\n";
								 std::cout << "" <<  vh0 << ", " << vh1 << ", " << vh2 << "; " << mesh_.property(vuv_, TriMesh::VHandle(1073)) << ": " << mpx << ".\n";
								 test_vertex_ = v0; test_vertex1_ = v1;
									}
									}*/
							} // end of this hx's empl_ == false 
						} // end of for . i = 0~9									
					} else {
						//std::cout << "vh1 valid, vh2 invalid, (" << f_it << ", " << mesh_.property(vf_, vh1) << ", "<< mesh_.property(vf_, vh2) << ").\n";							
					}
				} // end of if vj(vh1)���м���ĸ��������� and vh2Ҳ�����м���ĸ���������	
				if (vh2_in_valid_region) {
					if (vj_in_valid_region == true) {//vh2������Ч����֮��
						//std::cout << "vh1 valid, vh2 valid. \t";//for test	
						for (int i = 0; i < 9; i ++) { //�����ߵ��е���Ҫ����
							if (i == 0) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0);
							} else if (i == 1) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0);
							} else if (i == 2) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); 
							} else if (i == 3) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0); mesh_.property(vuv_, v6) = OpenMesh::Vec2d(-0.5, 0.8660254);  
							} else if (i == 4) { mesh_.property(vuv_, v3) = OpenMesh::Vec2d(0, 0); mesh_.property(vuv_, v7) = OpenMesh::Vec2d(0.5, -0.8660254); 
							} else if (i == 5) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0); mesh_.property(vuv_, v8) = OpenMesh::Vec2d(1.5, -0.8660254);
							} else if (i == 6) { mesh_.property(vuv_, v4) = OpenMesh::Vec2d(2, 0); mesh_.property(vuv_, v9) = OpenMesh::Vec2d(2.5, 0.8660254);
							} else if (i == 7) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); mesh_.property(vuv_, v10) = OpenMesh::Vec2d(2, 1.7320508);				
							} else if (i == 8) { mesh_.property(vuv_, v5) = OpenMesh::Vec2d(1, 1.7320508); mesh_.property(vuv_, v11) = OpenMesh::Vec2d(0, 1.7320508); 
							} 
							TriMesh::HHandle hx = vec_h[i]; 
							//if (test_ && i ==0) std::cout << mesh_.property(empl_, mesh_.edge_handle(hx)) << "; " << ;
							if (mesh_.property(empl_, mesh_.edge_handle(hx)) == false) {
								OpenMesh::Vec2d mpx = vec_mp[i];	
								OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, vh0), mesh_.property(vuv_, vh1), mesh_.property(vuv_, vh2), mpx);
								if (fabs(bc[0]) < 1.0e-10) bc[0] = 0; if (fabs(bc[1]) < 1.0e-10) bc[1] = 0; if (fabs(bc[2]) < 1.0e-10) bc[2] = 0;
								if (!(bc[0]<0) && !(bc[1]<0) && !(bc[2]<0)) { //Ҫ�����bc��������������>=0
									mesh_.property(emp_, mesh_.edge_handle(hx)) = mesh2_.point(vh0) * bc[0] + mesh2_.point(vh1) * bc[1] + mesh2_.point(vh2) * bc[2];
									mesh_.property(empl_, mesh_.edge_handle(hx)) = true; std::cout << "i = " << i << ", " << fvset.size() << ".\n";
									mesh_.property(emp_closest_vh_, mesh_.edge_handle(hx)) = vh0;
 
									//std::cout << ": " <<  v0 << ", " << v1 << ", " << v2 << ".\n";
									//std::cout << mesh2_.point(vh0) << "; " << mesh2_.point(vh1) << "; " << mesh2_.point(vh2) << "; " << bc << ";" << ".\n";
									//std::cout << mesh_.property(vuv_, vh0) << "; " << mesh_.property(vuv_, vh1) << "; " << mesh_.property(vuv_, vh2) << "; " << mpx << ";" << ".\n";//

								}/*
								 if (test_ && i ==0) { 
								 if (mesh_.status(*fvset_it).deleted() && mesh_.property(vf_, *fvset_it) == f_it && (mesh_.property(vf_, vh1) == f0 || mesh_.property(vf_, vh2) == f0)) {
								 std::cout << ": " <<  mesh_.property(vuv_, vh0) << ", " << mesh_.property(vuv_, vh1) << ", " << mesh_.property(vuv_, vh2) << ".\n";
								 std::cout << "" <<  vh0 << ", " << vh1 << ", " << vh2 << "; " << mesh_.property(vuv_, TriMesh::VHandle(1073)) << ": " << mpx << ".\n";
								 test_vertex_ = v0; test_vertex1_ = v1;
									}
									}*/
							} // end of this hx's empl_ == false 
						} // end of for . i = 0~9									
					} else {
						//std::cout << "vh1 valid, vh2 invalid, (" << f_it << ", " << mesh_.property(vf_, vh1) << ", "<< mesh_.property(vf_, vh2) << ").\n";							
					}
				}
			} // end of for.ѭ������˶���*fvset_it��һ���򶥵�vj, ���ж϶���*fvset_it�Ƿ���Ч,Ҳ����˵�ǲ��Ƕ���һ��������.
		} // end of for. ѭ��������������������f_it�ϵ����е� 	

	} // end of if.�ж�f_it�Ƿ�deleted, Ҳ���Ƿ���base complex�ϵ�һ����		
}
void DecimationModel::flatten_face_resample_3edge_midpoint() {
	for (TriMesh::FaceIter f_it(mesh_.faces_begin()), f_end(mesh_.faces_end()); f_it != f_end; ++f_it) {
		if (mesh_.status(f_it).deleted() == false)
			flatten_face_resample_3edge_midpoint(f_it.handle());	
	} // end of for.��������ѭ��һ��, ѭ��ʱ������ϵ������ߵ��е���в���

	int n_empl = 0;
	for (TriMesh::EIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end()); e_it != e_end; ++e_it) {
		if (mesh_.status(e_it).deleted() == false)
			if (mesh_.property(empl_, e_it.handle())) ++ n_empl;
	}
	std::cout << "Info: After flatten_face resampled: " << n_empl << " of " << simplified_mesh_.n_edges() << ".\n";
}
void DecimationModel::flatten_vertex_one_ring_resample_edges_midpoint() {
	for (TriMesh::VIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++ v_it) {
		if (mesh_.status(v_it).deleted() == false) {
			// ��ʱ������߽���һ�����ѹƽ.
			if (mesh_.is_boundary(v_it)) continue;

			flatten_vertex_one_ring_resample_edges_midpoint(v_it.handle());
		}

	} // end of for ��ÿ������.

	int n_empl = 0;
	for (TriMesh::EIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end()); e_it != e_end; ++e_it) {
		if (mesh_.status(e_it).deleted() == false)
			if (mesh_.property(empl_, e_it.handle())) ++ n_empl;
	}
	std::cout << "Info: After flatten_vertex_one_ring resampled: " << n_empl << " of " << simplified_mesh_.n_edges() << ".\n";
}
void DecimationModel::flatten_vertex_one_ring_resample_edges_midpoint(TriMesh::VHandle vh) {
	// ȷ��vh����һ���򶥵�Ĳ�������
	std::vector<TriMesh::VertexHandle> loop; //����ѭ��ؼ�¼vi����(����Ҳ����vh����)��one-ring��������Ķ���
	std::vector<TriMesh::HalfedgeHandle> loop_h;//���ڼ�¼vi����(����Ҳ����vh����)��one-ring��������ĳ����
	std::vector<TriMesh::FHandle> loop_f;
	TriMesh::HHandle hh = (mesh_.voh_iter(vh)).handle();
	TriMesh::HalfedgeHandle h = hh;
	do {	// ���ﲻʹ��VertexOHalfedgeIter��VertexVertexIter��ԭ������Ҫ��֤loop[j]��˳��.loop[0]=to
		loop_h.push_back(h);
		loop_f.push_back(mesh_.face_handle(h));
		loop.push_back(mesh_.to_vertex_handle(h));//std::cout << mesh_.to_vertex_handle(h) << "\t";	 ////	for test	
		h = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_handle(h));
	} while(h != hh);
	if (loop.size() != loop_h.size()) { std::cerr << "Error: collect one-ring elements.\n"; };

	unsigned int n = loop_h.size();//�ж��ٸ��߽綥��
	unsigned int i = 0;
	double length = 0.0, rou_angle = 0.0;//
	std::vector<double> vec_angle;
	for (i=0 ; i<n; ++i) {//����ܵĳ���, �ܳ�, �Լ��Ƕ�֮��. 
		length += (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n]))).norm();

		OpenMesh::Vec3d d1 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[i])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(vh))).normalize();
		OpenMesh::Vec3d d2 = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(vh))).normalize();
		rou_angle += acos(dot(d1, d2));
		vec_angle.push_back(acos(dot(d1, d2)));
	}  //std::cout << "angle " << rou_angle << ".\n";//std::cout << "len: " << length << std::endl; //for test	
	if (vec_angle.size() != n) { std::cout << "Error: ������һ��.\n"; }	

	double angle_scale_ratio = 2 * M_PI / rou_angle; //���ű���
	double temp_sum_angle = 0.0, len = 0;
	for (int i = 0; i < n; ++i) { //һ�����Ĳ���������.
		temp_sum_angle += (vec_angle[i] * angle_scale_ratio);
		len = (OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(loop[(i+1)%n])) - OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_.point(vh))).norm(); 
		len *= angle_scale_ratio; //len = pow(len, angle_scale_ratio);
		mesh_.property(vuv_, loop[(i+1)%n]) = OpenMesh::Vec2d(len*cos(temp_sum_angle), len*sin(temp_sum_angle));
	} 		
	mesh_.property(vuv_, vh) = OpenMesh::Vec2d(0, 0); //���һ��������� 

	// ȷ��vh����һ���������б����������������
	for (std::vector<TriMesh::FHandle>::const_iterator f_it(loop_f.begin()), f_end(loop_f.end()); f_it != f_end; ++ f_it) {
		TriMesh::FVIter fv_it(mesh_, *f_it);
		TriMesh::VHandle v0(fv_it.handle()), v1((++fv_it).handle()), v2((++fv_it).handle());
		std::vector<TriMesh::VertexHandle> fvset = mesh_.property(fvset_, *f_it);
		for (std::vector<TriMesh::VertexHandle>::iterator it(fvset.begin()), it_end(fvset.end()); it != it_end; ++it) {
			if (mesh_.status(*it).deleted() == false) { std::cerr << "Error: \n"; return; } // ��Щ��ǰ�Ѿ����������ĵ�Ӧ�ö���deleted vertices.

			TriMesh::VertexHandle vf0 = mesh_.property(vf0_, *it); //���*it��������ǰ�Ĳ��������������ڵ������������,
			TriMesh::VertexHandle vf1 = mesh_.property(vf1_, *it); //���������������ɵ���, �պ�Ӧ�þ�������ѭ���е��Ǹ���. 
			TriMesh::VertexHandle vf2 = mesh_.property(vf2_, *it);
			if (vf0 != v0 && vf0 != v1 && vf0 != v2) { std::cerr << "Error: parameterziation vf0.���ó���\n"; return; }
			if (vf1 != v0 && vf1 != v1 && vf1 != v2) { std::cerr << "Error: parameterziation vf1.���ó���\n"; return; }
			if (vf2 != v0 && vf2 != v1 && vf2 != v2) { std::cerr << "Error: parameterziation vf2.���ó���\n"; return; }
			mesh_.property(vuv_, *it) = mesh_.property(vuv_, vf0) * (mesh_.property(vbc_, *it))[0]
			+ mesh_.property(vuv_, vf1) * (mesh_.property(vbc_, *it))[1]
			+ mesh_.property(vuv_, vf2) * (mesh_.property(vbc_, *it))[2];				
		}	
	}

	// ��ÿ���ڱ߽����е����.
	for (std::vector<TriMesh::FHandle>::const_iterator f_it(loop_f.begin()), f_end(loop_f.end()); f_it != f_end; ++ f_it) 
	{
		std::vector<TriMesh::VertexHandle> fvset = mesh_.property(fvset_, *f_it);
		for (std::vector<TriMesh::VertexHandle>::iterator it(fvset.begin()), it_end(fvset.end()); it != it_end; ++it) 
		{
			TriMesh::VHandle v0 = *it; //���ڿ��ǵĶ���.
			for (TriMesh::VOHIter voh_it(mesh2_, v0); voh_it; ++voh_it) 
			{
				TriMesh::VHandle v1 = mesh2_.to_vertex_handle(voh_it);
				std::vector<TriMesh::FHandle>::iterator find_result = find(loop_f.begin(), loop_f.end(), mesh_.property(vf_, v1));
				if (find_result != loop_f.end()) { // find succeed.
					TriMesh::VHandle v2 = mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(voh_it)); 
					find_result = find(loop_f.begin(), loop_f.end(), mesh_.property(vf_, v2));
					if (find_result != loop_f.end()) {
						for (std::vector<TriMesh::HHandle>::const_iterator h_it(loop_h.begin()), h_end(loop_h.end()); h_it != h_end; ++ h_it) {
							if (mesh_.property(empl_, mesh_.edge_handle(*h_it)) == false) {
								OpenMesh::Vec2d midpoint((mesh_.property(vuv_, vh) + mesh_.property(vuv_, mesh_.to_vertex_handle(*h_it))) * 0.5);
								OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh_.property(vuv_, v0), mesh_.property(vuv_, v1), mesh_.property(vuv_, v2), midpoint);
								if (DGP::is_valid_barycentric_coordinate(bc)) {
									mesh_.property(empl_, mesh_.edge_handle(*h_it)) = true;
									mesh_.property(emp_, mesh_.edge_handle(*h_it)) = mesh2_.point(v0) * bc[0] + mesh2_.point(v1) * bc[1] + mesh2_.point(v2) * bc[2];
									if (bc[0] > bc[1]) {
										if (bc[0] > bc[2]) mesh_.property(emp_closest_vh_, mesh_.edge_handle(*h_it)) = v0;
										else mesh_.property(emp_closest_vh_, mesh_.edge_handle(*h_it)) = v2;
									} else {
										if (bc[1] > bc[2]) mesh_.property(emp_closest_vh_, mesh_.edge_handle(*h_it)) = v1;
										else mesh_.property(emp_closest_vh_, mesh_.edge_handle(*h_it)) = v2;
									}
								}
							}
						}
					}
				}
			}
		}	
	}
}
void DecimationModel::use_closepoints(TriMesh::EHandle _eh) {
	// Notes: �߽���ǲ�������󲻵��е�������, ���������_eh���Ǳ߽��.
	TriMesh::HHandle h0 = mesh_.halfedge_handle(_eh, 0);
	TriMesh::HHandle h1 = mesh_.halfedge_handle(_eh, 1);
	TriMesh::FHandle f0 = mesh_.face_handle(h0);
	TriMesh::FHandle f1 = mesh_.face_handle(h1);
	TriMesh::Point mp = mesh_.point(mesh_.to_vertex_handle(h1)) + (mesh_.point(mesh_.to_vertex_handle(h0)) - mesh_.point(mesh_.to_vertex_handle(h1))) * 0.5;
	TriMesh::VHandle closest_vertex;
	double dis = (mesh_.point(mesh_.to_vertex_handle(h0)) - mesh_.point(mesh_.to_vertex_handle(h1))).norm();//��ʼ��Ϊ�߳�
	//std::cout <<"use_closepoints: " << dis << ".\n";
	for (TriMesh::VIter v_it(mesh2_.vertices_begin()), v_end(mesh2_.vertices_end()); v_it != v_end; ++v_it) {
		if ((mesh2_.point(v_it) - mp).norm() < dis) {
			closest_vertex = v_it.handle();
			dis = (mesh2_.point(v_it) - mp).norm();
		}
	}
	// �ҵ����Ǹ��������֮��, �ͳ���ͶӰ������һ��������, ���Ƿ���ͶӰ�ɹ�.���������ʦ������.
	// Ҫ��ͶӰ���ɹ�, �Ǿ�ֱ���������ҵ����Ǹ�����Ķ���, ���������ʦԭ��������.
	if (closest_vertex.is_valid()) {
		for (TriMesh::VOHIter voh_it(mesh2_, closest_vertex); voh_it; ++ voh_it) {
			OpenMesh::Vec3d n = mesh2_.calc_face_normal(mesh2_.face_handle(voh_it));
			double d = dot(n, mp - mesh2_.point(closest_vertex));
			OpenMesh::Vec3d point_on_plane = mp - d * n;

			OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(mesh2_.point(closest_vertex), mesh2_.point(mesh2_.to_vertex_handle(voh_it)), 
				mesh2_.point(mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(voh_it))), point_on_plane);
			if (DGP::is_valid_barycentric_coordinate(bc)) {
				mesh_.property(empl_, _eh) = true; 
				mesh_.property(emp_, _eh) = mesh2_.point(closest_vertex) * bc[0] + mesh2_.point(mesh2_.to_vertex_handle(voh_it)) * bc[1]
				+ mesh2_.point(mesh2_.to_vertex_handle(mesh2_.next_halfedge_handle(voh_it))) * bc[2]; 
				mesh_.property(emp_closest_vh_, _eh) = closest_vertex; 
			}
		}
		if (mesh_.property(empl_, _eh) == false) {
			std::cout << "Info: One edge use closetest point as midpoint directly1.\n";
			mesh_.property(empl_, _eh) = true;
			mesh_.property(emp_, _eh) = mesh2_.point(closest_vertex);
			mesh_.property(emp_closest_vh_, _eh) = closest_vertex;
		}
	} else {
		mesh_.property(empl_, _eh) = false; 
		std::cout << "use_closepoints(): vh is not valid, ����û���ҵ�.";//\n
	} 
}
void DecimationModel::project_midpoint_to_triangle() {
	for (TriMesh::EdgeIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end()); e_it != e_end; ++e_it) {
		if (mesh_.status(e_it).deleted() == false && mesh_.property(empl_, e_it) == false) {
			use_closepoints(e_it.handle());
		}
	}
	// ��������ж��ٱ��Ѿ��ɹ�������ص����.
	int n_empl = 0;
	for (TriMesh::EIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end()); e_it != e_end; ++e_it) {
		if (mesh_.status(e_it).deleted() == false)
			if (mesh_.property(empl_, e_it.handle())) ++ n_empl;
	}
	std::cout << "Info: After use closet resampled: " << n_empl << " of " << simplified_mesh_.n_edges() << ".\n"; 
}
void DecimationModel::resample_edge_midpoint_process() {
	std::cout << "\n";
	std::cout << "========��ʼ�е��������================\n";
	// ------------��ʼ�Աߵ��е���в���-------------------
	//std::cout << "��ʼ�Աߵ��е���в���.\n";
	// ��һ�����,����ı�, ����original mesh�оʹ��ڵı�, �򻯹�����û�иı�, �Լ�crease edges.
	resample_original_feature_edge_midpoint();

	// ����չ������һ����ı߽����е����.//
	flatten_vertex_one_ring_resample_edges_midpoint();
	// ����չ��������������߽����е����.
	flatten_face_resample_3edge_midpoint();
	// ��ͶӰ����������
	project_midpoint_to_triangle();
	// ------�����Աߵ��е���в���

	// -------------
	midpoint_sample_succeed_ = true;//�������б߶��ҵ����е�Ĳ�����.
	int count_crosspatches = 0;
	int edge_midpoint_valid(0);//�ж��ٸ���, ����е��Ѿ��ҵ���������?
	for (TriMesh::EdgeIter e_it(mesh_.edges_begin()), e_end(mesh_.edges_end()); e_it != e_end; ++ e_it) { //ѭ��ԭʼ�����е����б�
		if (mesh_.status(e_it).deleted() == true) {
			TriMesh::VHandle v0 = mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 0));
			TriMesh::VHandle v1 = mesh2_.to_vertex_handle(mesh2_.halfedge_handle(e_it, 1));
			if (mesh_.property(vf_, v0) != mesh_.property(vf_, v1)) {
				count_crosspatches++;
			}
		} else { // ����û�б�deleted�ı�, ��Ҫ�Աߵ��е���в���.
			if (mesh_.property(empl_, e_it.handle()) == true) { //�Ѿ�����е����.
				edge_midpoint_valid ++; 
				// std::cout << mesh_.to_vertex_handle(mesh_.halfedge_handle(e_it, 0)) << " " << mesh_.to_vertex_handle(mesh_.halfedge_handle(e_it, 1)) << ".\n";
			} else { // ��û������е��
				midpoint_sample_succeed_ = false; //��һ���ߵ��е�û�в����ɹ��Ļ�Ҳ��ʧ����.
 
				std::cout << "û�в����ɹ��ı�����Ӧ��������: " 
					<< mesh_.face_handle(mesh_.halfedge_handle(e_it, 0)) << " " << mesh_.face_handle(mesh_.halfedge_handle(e_it, 1)) << ".\n";
			}
		}
	}
	std::cout << "para check: E: cross patches: " << count_crosspatches << ", midpoint resampled: " << edge_midpoint_valid << std::endl;

	set_mesh_toberendered(&simplified_mesh_, RESAMPLE_MIDPOINT);
	std::cout << "========�����е��������================\n";
}