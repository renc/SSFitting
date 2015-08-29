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
	TriMesh::HHandle hh = _mesh.property(fp_hh_icm_, fh);//这就暗示了参数中的_mesh其实是initial_control_mesh_.
	//TriMesh::VHandle v0 = _mesh.from_vertex_handle(hh);
	//TriMesh::VHandle v1 = _mesh.to_vertex_handle(hh);
	//TriMesh::VHandle v2 = _mesh.to_vertex_handle((_mesh.next_halfedge_handle(hh)));

	if (debug) std::cout << "1"; 
	if (DGP::is_valid_barycentric_coordinate(_bc)) {  
		//std::vector<OpenMesh::Vec3d> f_pijk = qas_->f_pijk(fh, true);
		
		enum return_type { Null, Min_dis, Max_times, Out_vw };
		return_type rt = Null;

		const int Max_Iter_Times = 10000;//最大的迭代次数.这些参数都是可以调的.
		const double Max_Step_length = 0.5; //最大步长

		OpenMesh::Vec3d clp_arr[2];//求两次, 正反方向, 哪次合理取那次.
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
				//这个没有问题,打印得到的结果就是调用这函数之前得到的mesh_.property(vp_eval_pos_, *it),
				//说明clp不需要=mesh_.property(vp_eval_pos_, *it)作为初值, 
				solve_func_get_h(dQAS_dv, dQAS_dw, _Sk, qas_->pos(fh, 1-v-w, v, w, true), delta_v, delta_w);
				if (debug) std::cout << "delta_v " << delta_v << ", delta_w " << delta_w << ".\n";
				if (fabs(delta_v) <= 1 || fabs(delta_w) <= 1) {  //计算出来的方向可能是对的, 但是长度太大了, 超出0~1
					// 正常值
				} else {
					// 特别是预防那种1.#QNAN, 都不知怎么来的.2009-07-06被此害惨了.
					//double len = sqrt(delta_v * delta_v + delta_w * delta_w);
					delta_v = 0.9;// /= len;
					delta_w = 0.9; // /= len;
				}
				if (delta_v < 0.0001 && delta_w < 0.00001) {
					// 这情况可以看做是迭代终止了, 这么短的步长.
					clp = qas_->pos(fh, 1-v-w, v, w, true);
					rt = Min_dis; break;				
				} 
				else {
					// 2. 初始步长.
					double step_len = -0.01; //-0.008;//1;
					if (i == 0) step_len *= 1;
					else step_len *= -1;

					//std::cout << "--\n";
					// 3. Update delta_v and delta_w.// 得到步长step side之后			
					delta_v *= step_len;
					delta_w *= step_len;

					if ((v + delta_v) < 0) {
						//std::cout << "v < 0. "; //, " << v << ", + " << delta_v << ", " << _v << " ////单独没有问题, 和w<0一起没有问题, 但和>1时就有问题.
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
								break; //不是返回错误, 而是找到当前为止, 当前的结果也是有效的. //return false;
							} 
						}
						rt = Out_vw; break; 				
					} else if ((w + delta_w) < 0) {
						//std::cout << "w < 0. "; //\n// //单独没有问题, 和v<0一起没有问题, 和>1也没有问题.
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
								break; //不是返回错误, 而是找到当前为止, 当前的结果也是有效的. //return false;
							} 
						}				
						rt = Out_vw; break;
					} else if ((v + delta_v + w + delta_w) > 1) {
						//std::cout << "v + w > 1. "; //\n//单独没有问题, 和w<0也没有问题, 但和v<0一起就有问题
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
								break; //不是返回错误, 而是找到当前为止, 当前的结果也是有效的. //return false;
							} 
						}
						rt = Out_vw; break;
					} else if ((v + delta_v) > 1) { //好像没有出现过此情况.
						//std::cout << "v > 1, " << v << ", + " << delta_v << ".\n"; rt = Out_vw; break;
					} else if ((w + delta_w) > 1) { //好像没有出现过此情况.
						//std::cout << "w > 1.\n"; rt = Out_vw; break;
					} else { 
						//还在同一个面里面找着
						//std::cout << "v: " << v << ", w: " << w << ". ";
						v += delta_v; 
						w += delta_w;
						if (debug) std::cout << "v': " << v << ", w': " << w << ".\n";

						double distance_before = (clp-_Sk).norm();
						OpenMesh::Vec3d clp_before = clp; //保存前一个最近点坐标 
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
				//std::cout << "Info: Iteration ends Max_Iter_Times, OK.\n";//由于超过给定迭代次数, 
				std::cout << "Info: Max_Iter_Times\n";
				rt = Max_times;
			} 
			clp_arr[i] = clp; //记下这一次closest position. 
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

	//// 求在极限曲面上的cloesest point position
	//OpenMesh::Vec3d clp = pp; //最近点的初始坐标为evaluation坐标.
	//if (cloop_stam98_get_closest_points(cN6, v, w, Sk, clp, initial_control_mesh_, h12, cross_patch_depth)) {
	//	_clp = clp;
	//	return true;
	//} else { //没有正确求出最近点.
	//	std::cout << "Error: evaluation_get_footpoint(...): non-feature region closest point.\n";
	//	return false;
	//}

	// 


	return false; 
 
}