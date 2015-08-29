#include "DecimationModel.h"
//#include "../QuarticBezierTriangle/QuadraticBezierTriangles.h"// renc, 20150829
#include "QuadraticBezierTriangles.h" // renc, 20150829

bool DecimationModel::evaluation_get_footpoint_byQAS(const TriMesh& _mesh, TriMesh::FaceHandle _fh_of_initial_control_mesh_, // input
									OpenMesh::Vec3d _Sk,  OpenMesh::Vec3d _bc, // input, 
									OpenMesh::Vec3d &_clp, int _cross_patch_depth) //output 
{
	 
	bool debug = false;
	if (debug) std::cout << "-0";
	TriMesh::FHandle fh = _fh_of_initial_control_mesh_;
	TriMesh::HHandle hh = _mesh.property(fp_hh_icm_, fh);//��Ͱ�ʾ�˲����е�_mesh��ʵ��initial_control_mesh_.
	//TriMesh::VHandle v0 = _mesh.from_vertex_handle(hh);
	//TriMesh::VHandle v1 = _mesh.to_vertex_handle(hh);
	//TriMesh::VHandle v2 = _mesh.to_vertex_handle((_mesh.next_halfedge_handle(hh)));

	if (debug) std::cout << "1"; 
	if (DGP::is_valid_barycentric_coordinate(_bc)) {  
		//std::vector<OpenMesh::Vec3d> f_pijk = qas_->f_pijk(fh, true);
		
		enum return_type { Null, Min_dis, Max_times, Out_vw };
		return_type rt = Null;

		const int Max_Iter_Times = 10000;//���ĵ�������.��Щ�������ǿ��Ե���.
		const double Max_Step_length = 0.5; //��󲽳�

		OpenMesh::Vec3d clp_arr[2];//������, ��������, �Ĵκ���ȡ�Ǵ�.
		for (int i = 0; i < 2; ++ i) {
			double u = _bc[0], v = _bc[1], w = _bc[2]; // initial values for iteration.
			OpenMesh::Vec3d clp;// record the result in this round.
			int iter_times = 0; 
			for (iter_times = 0; iter_times < Max_Iter_Times; ++iter_times) 
			{
				// 0. Jacobian Matrix
				OpenMesh::Vec3d dQAS_dv = qas_->dPos_dv(fh, 1-v-w, v, w, true);
				OpenMesh::Vec3d dQAS_dw = qas_->dPos_dw(fh, 1-v-w, v, w, true); 

				// 1. search direction h = [delta_v, delta_w]T
				double delta_v, delta_w;
				//if (iter_times == 0) std::cout << " : " << qas_->pos(fh, 1-v-w, v, w, true) << "\\ ";
				//���û������,��ӡ�õ��Ľ�����ǵ����⺯��֮ǰ�õ���mesh_.property(vp_eval_pos_, *it),
				//˵��clp����Ҫ=mesh_.property(vp_eval_pos_, *it)��Ϊ��ֵ, 
				solve_func_get_h(dQAS_dv, dQAS_dw, _Sk, qas_->pos(fh, 1-v-w, v, w, true), delta_v, delta_w);
				if (debug) std::cout << "delta_v " << delta_v << ", delta_w " << delta_w << ".\n";
				if (fabs(delta_v) <= 1 || fabs(delta_w) <= 1) {  //��������ķ�������ǶԵ�, ���ǳ���̫����, ����0~1
					// ����ֵ
				} else {
					// �ر���Ԥ������1.#QNAN, ����֪��ô����.2009-07-06���˺�����.
					//double len = sqrt(delta_v * delta_v + delta_w * delta_w);
					delta_v = 0.9;// /= len;
					delta_w = 0.9; // /= len;
				}
				if (delta_v < 0.0001 && delta_w < 0.00001) {
					// ��������Կ����ǵ�����ֹ��, ��ô�̵Ĳ���.
					clp = qas_->pos(fh, 1-v-w, v, w, true);
					rt = Min_dis; break;				
				} 
				else {
					// 2. ��ʼ����.
					double step_len = -0.01; //-0.008;//1;
					if (i == 0) step_len *= 1;
					else step_len *= -1;

					//std::cout << "--\n";
					// 3. Update delta_v and delta_w.// �õ�����step side֮��			
					delta_v *= step_len;
					delta_w *= step_len;

					if ((v + delta_v) < 0) {
						//std::cout << "v < 0. "; //, " << v << ", + " << delta_v << ", " << _v << " ////����û������, ��w<0һ��û������, ����>1ʱ��������.
						if (_cross_patch_depth < Max_Cross_Patch_Depth) {
							TriMesh::HHandle hh_new = initial_control_mesh_.opposite_halfedge_handle(initial_control_mesh_.prev_halfedge_handle(hh));
							TriMesh::FHandle fh_new = initial_control_mesh_.face_handle(hh_new);
							double w_new, v_new;
							update_uvw_when_v_new_is_less_than_zero(v, delta_v, w, delta_w, v_new, w_new);
							OpenMesh::Vec3d bc_new(1-v_new-w_new, v_new, w_new);
							//std::cout << _v << ", " << v << ", " << delta_v << ". " << _w << ", "  w << ", " << delta_w << ".\n";
							if (DGP::is_valid_barycentric_coordinate(bc_new) == false) {
								std::cout << "Error: v < 0, the (u', v', w') is invalid.\n";  return false;
							}
							if (evaluation_get_footpoint_byQAS(initial_control_mesh_, fh_new, _Sk, bc_new, clp, ++_cross_patch_depth)== false) {
								break; //���Ƿ��ش���, �����ҵ���ǰΪֹ, ��ǰ�Ľ��Ҳ����Ч��. //return false;
							} 
						}
						rt = Out_vw; break; 				
					} else if ((w + delta_w) < 0) {
						//std::cout << "w < 0. "; //\n// //����û������, ��v<0һ��û������, ��>1Ҳû������.
						if (_cross_patch_depth < Max_Cross_Patch_Depth) {
							TriMesh::HHandle hh_new = initial_control_mesh_.opposite_halfedge_handle(hh);
							TriMesh::FHandle fh_new = initial_control_mesh_.face_handle(hh_new);
							double w_new, v_new;
							update_uvw_when_w_new_is_less_than_zero(v, delta_v, w, delta_w, v_new, w_new);
							OpenMesh::Vec3d bc_new(1-v_new-w_new, v_new, w_new);
							//std::cout << _v << ", " << v << ", " << delta_v << ". " << _w << ", "  w << ", " << delta_w << ".\n";
							if (DGP::is_valid_barycentric_coordinate(bc_new) == false) {
								std::cout << "Error: v < 0, the (u', v', w') is invalid.\n";  return false;
							}
							if (evaluation_get_footpoint_byQAS(initial_control_mesh_, fh_new, _Sk, bc_new, clp, ++_cross_patch_depth)== false) {
								break; //���Ƿ��ش���, �����ҵ���ǰΪֹ, ��ǰ�Ľ��Ҳ����Ч��. //return false;
							} 
						}				
						rt = Out_vw; break;
					} else if ((v + delta_v + w + delta_w) > 1) {
						//std::cout << "v + w > 1. "; //\n//����û������, ��w<0Ҳû������, ����v<0һ���������
						if (_cross_patch_depth < Max_Cross_Patch_Depth) {
							TriMesh::HHandle hh_new = initial_control_mesh_.opposite_halfedge_handle(initial_control_mesh_.next_halfedge_handle(hh));
							TriMesh::FHandle fh_new = initial_control_mesh_.face_handle(hh_new);
							double v_new, w_new;
							update_uvw_when_v_new_plus_w_new_is_greaterr_than_one(v, delta_v, w, delta_w, v_new, w_new);
							OpenMesh::Vec3d bc_new(1-v_new-w_new, v_new, w_new);
							if (DGP::is_valid_barycentric_coordinate(bc_new) == false) {
								std::cout << "Error: v+w>1, the (u', v', w') is invalid.\n"; return false;
							}
							if (evaluation_get_footpoint_byQAS(initial_control_mesh_, fh_new, _Sk, bc_new, clp, ++_cross_patch_depth)== false) {
								break; //���Ƿ��ش���, �����ҵ���ǰΪֹ, ��ǰ�Ľ��Ҳ����Ч��. //return false;
							} 
						}
						rt = Out_vw; break;
					} else if ((v + delta_v) > 1) { //����û�г��ֹ������.
						//std::cout << "v > 1, " << v << ", + " << delta_v << ".\n"; rt = Out_vw; break;
					} else if ((w + delta_w) > 1) { //����û�г��ֹ������.
						//std::cout << "w > 1.\n"; rt = Out_vw; break;
					} else { 
						//����ͬһ������������
						//std::cout << "v: " << v << ", w: " << w << ". ";
						v += delta_v; 
						w += delta_w;
						if (debug) std::cout << "v': " << v << ", w': " << w << ".\n";

						double distance_before = (clp-_Sk).norm();
						OpenMesh::Vec3d clp_before = clp; //����ǰһ����������� 
						clp = qas_->pos(fh, 1-v-w, v, w, true);//std::cout << "' " << clp; 
						double distance_after = (clp-_Sk).norm(); 
						if (distance_before <= distance_after) {
							//std::cout << "Info: Iteration " << iter_times << ", ends, OK.\n"; //for distance_before < distance_after, 
							if (iter_times == 0 || iter_times == 1) {
								//std::cout << distance_before << ", " << distance_after << ".\n";
								test_array_.push_back(_Sk);
							}
							clp = clp_before; 
							rt = Min_dis; break;
						}
						//if (fabs(distance_after - distance_before) < 1.0e-4) {
						//	clp = clp_before;
						//	rt = Min_dis; break;
						//}
					} 
				} // delta_v and delta_w is meaningful.				
			} // end of for iter_times
			if (iter_times == Max_Iter_Times) {
				//std::cout << "Info: Iteration ends Max_Iter_Times, OK.\n";//���ڳ���������������, 
				std::cout << "Info: Max_Iter_Times\n";
				rt = Max_times;
			} 
			clp_arr[i] = clp; //������һ��closest position. 
		} // end of for i == 0 or 1. 
		if ((_Sk - clp_arr[0]).norm() < (_Sk - clp_arr[1]).norm() ) _clp = clp_arr[0];
		else _clp = clp_arr[1];

		if (rt == Null) {
			return false;
		} else if (rt == Max_times || rt == Min_dis ) {
			return true;
		} else if (rt == Out_vw) {
			return true;
		}
	} else {
		std::cout << "Error: evaluation_get_footpoint_byQAS: _bc is invalid, _bc " << _bc << ".\n";
	}


	//_pp = pp;

	//// ���ڼ��������ϵ�cloesest point position
	//OpenMesh::Vec3d clp = pp; //�����ĳ�ʼ����Ϊevaluation����.
	//if (cloop_stam98_get_closest_points(cN6, v, w, Sk, clp, initial_control_mesh_, h12, cross_patch_depth)) {
	//	_clp = clp;
	//	return true;
	//} else { //û����ȷ��������.
	//	std::cout << "Error: evaluation_get_footpoint(...): non-feature region closest point.\n";
	//	return false;
	//}

	// 


	return false; 
 
}