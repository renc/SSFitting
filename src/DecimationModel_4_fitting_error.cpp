#include "DecimationModel.h"
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_min.h>

//#include "../QuarticBezierTriangle/QuadraticBezierTriangles.h"// renc, 20150829
#include "QuadraticBezierTriangles.h" // renc, 20150829

OpenMesh::Vec4d DecimationModel::calc_step_side1_parames_ = OpenMesh::Vec4d(0,0,0,0); 
OpenMesh::Vec3d DecimationModel::calc_step_side1_sk_ = OpenMesh::Vec3d(0,0,0);
std::vector<OpenMesh::Vec3d> DecimationModel::calc_step_side1_c_ = std::vector<OpenMesh::Vec3d>(0, OpenMesh::Vec3d(0,0,0));

double DecimationModel::calc_step_side1(double _h, void *params) {
	double v = currentDV->calc_step_side1_parames_[0], w = currentDV->calc_step_side1_parames_[1];
	double delta_v = currentDV->calc_step_side1_parames_[2], delta_w = currentDV->calc_step_side1_parames_[3];

	OpenMesh::Vec3d pos(0,0,0);
	if (currentDV->get_cloop_evaluation().evaluate_irregular_patch(currentDV->calc_step_side1_c_, v+_h*delta_v, w+_h*delta_w, pos, OpenMesh::Vec3d(0,0,0), OpenMesh::Vec3d(0,0,0)) == true)
	return (pos - currentDV->calc_step_side1_sk_).norm();
	else return 100;
}
void DecimationModel::solve_func_get_h (const OpenMesh::Vec3d dStam_dv, const OpenMesh::Vec3d dStam_dw, 
										const OpenMesh::Vec3d Sk, const OpenMesh::Vec3d clp, double & delta_v, double & delta_w) 
{
	// 下面手动化简并解这个二元方程组(两个方程).JTJh=-JTf
	// J(3*2) = [dStan_dv, dStam_dw]
	// Ah=B
	double B[2];//-JTf = -dot(JT, f)
	B[0] = -dot(dStam_dv, (clp - Sk));//std::cout << dStam_dv << ", " << dStam_dw << ", " << clp << ", " << Sk << ".";
	B[1] = -dot(dStam_dw, (clp - Sk));
	//std::cout << B[0] << ", " << B[1] << ".\n"; // for test //
	double A[2][2];//JTJ
	A[0][0] = dot(dStam_dv, dStam_dv); A[0][1] = dot(dStam_dv, dStam_dw); 
	A[1][0] = dot(dStam_dw, dStam_dv); A[1][1] = dot(dStam_dw, dStam_dw); 
	//std::cout << A[0][0] << ", " << A[0][1] << "; " << A[1][0] << ", " << A[1][1] << ".\n";  // for test////
	if (fabs(A[0][0] - A[1][0]) < 1e-6) {
		delta_v = B[0] * 0.5 / A[0][0];
		delta_w = B[0] * 0.5 / A[0][0];
		return;
	}

	enum IterationType { Steepest_descent, Gauss_newton };
	IterationType itype = Gauss_newton;
	if (itype == Steepest_descent) { //结果是这个即使使用全步长(1or-1)都是收敛极慢的! 迭代50次只是前进了一点点.
		delta_v = B[0];
		delta_w = B[1]; 
	} else {
		delta_v = (B[0]*A[1][1] - B[1]*A[0][1]) / (A[0][0]*A[1][1] - A[0][1]*A[1][0]);	//std::cout << "delta_v: " << delta_v << "; ";////
		delta_w = (B[0] - A[0][0] * delta_v) / A[0][1];									//std::cout << "delta_w: " << delta_w << ".\n";// //
	}		 
}
void DecimationModel::update_uvw_when_v_new_is_less_than_zero(const double v, const double delta_v, const double w, const double delta_w, double &v_new, double &w_new) {
	//正常变化v.
	v_new= fabs(v + delta_v);
	w_new = w + delta_w;
	if (v_new > 1) {
		v_new *= 0.5;
		if (v_new > 1) {
			v_new *= 0.5;
			if (v_new > 1) {
				v_new = 0.7; w_new = 0.2;
			}
		}
	}
	if (w_new > 1) {
		w_new = 1 - v_new - 0.05;
	}
	if (w_new < 0) {
		w_new = 0.001;
	}
	if (v_new + w_new > 1) { //如果正常变化不行的话,
		if (delta_w > 0) w_new = w;	//w先不变,
		if (v_new + w_new > 1) {//还是不行的话, 让v的步长缩小一半.
			v_new *= 0.5;
			if (v_new + w_new > 1) v_new *= 0.5; 
			if (v_new + w_new > 1) v_new *= 0.5; 
			if (v_new + w_new > 1) {				
				v_new = 0.0001; 
				w_new -= 0.0001;
			}
		}
	} //std::cout << v_new << ", " << w_new << ".\t";
}
void DecimationModel::update_uvw_when_w_new_is_less_than_zero(const double v, const double delta_v, const double w, const double delta_w, double &v_new, double &w_new) {
	v_new = v + delta_v;
	w_new = fabs(w + delta_w);
	//std::cout << "B: " << 1-v-w << ", " << v << ", " << w << ". A: " << 1-(w + delta_w+v + delta_v) << ", " << v + delta_v << ", " << w + delta_w 
	//	<< ". Old bc:" << 1-v_new-w_new << ", " << v_new << ", " << w_new << ". ";
	if (v_new + w_new > 1) {
		if (v_new > 0.999) 
			if (delta_v > 0) v_new = v;
		if (v_new + w_new > 1) {
			w_new *= 0.5;
			if (v_new + w_new > 1) w_new *= 0.5;
			if (v_new + w_new > 1) { w_new *= 0.5; }
		}
	}
}
void DecimationModel::update_uvw_when_v_new_plus_w_new_is_greaterr_than_one(const double v, const double delta_v, const double w, const double delta_w, double &v_new, double &w_new) {
	v_new = 1 - (w + delta_w);
	w_new = 1 - (v + delta_v);
	//std::cout << "B: " << 1-v-w << ", " << v << ", " << w << ". delta: " << delta_v << ", " << delta_w << "A: " << 1-( 1- w - delta_w) - (1 - v - delta_v) << ", " << 1 - (w + delta_w) << ", " << 1 - (v + delta_v) << ". ";
	if (v + delta_v > 1) {
		w_new = 0.0001; v_new -= 0.0001; // because: w_new = 1 - (v + delta_v) < 0
	}
	if (w + delta_w > 1) {
		v_new = 0.0001; w_new -= 0.0001;
	}
	//
}////

// Newton Iteration求最近点.
// 求Sk在cN6这邻域的顶点组成的曲面上的最近点_clp, 迭代的初始值为(1 - _v - _w, _v, _w), _clp的初始值为Sk在这曲面片上的精确坐标。
bool DecimationModel::cloop_stam98_get_closest_points (const std::vector<OpenMesh::Vec3d> & cN6, 
													   const double _v, const double _w, OpenMesh::Vec3d Sk, 
													   OpenMesh::Vec3d & _clp, 
													   const TriMesh& _mesh, TriMesh::HHandle h12, int cross_patch_depth) 
{
	//这是irregular的模板. regular也用irregular的方式来求解.
	// (1), search descent direction h(delta_v, delta_w)
	//std::cout << "Sk and initial clp " << Sk << ", " << _clp << ". Initial dis: " << (Sk - _clp).norm() << ".\n";
	std::vector<OpenMesh::Vec3d> c = cN6;
	enum return_type { Null, Min_dis, Max_times, Out_vw };
	return_type rt = Null;

	const int Max_Iter_Times = 2000;//最大的迭代次数.这些参数都是可以调的.
	const double Max_Step_length = 0.5; //最大步长

	OpenMesh::Vec3d clp_arr[2];//求两次, 正反方向, 哪次合理取那次.
	for (int i = 0; i < 2; ++ i) {
		double v = _v, w = _w;
		OpenMesh::Vec3d clp = _clp;
		int iter_times = 0; 
		for (iter_times = 0; iter_times < Max_Iter_Times; ++iter_times) 
		{
			// 0. Jacobian Matrix
			OpenMesh::Vec3d dStam_dv(0, 0, 0), dStam_dw(0, 0, 0); 

			if (cloop_evaluation.evaluate_irregular_patch(c, v, w, OpenMesh::Vec3d(0,0,0), dStam_dv, dStam_dw) == false) {
				std::cout << "Error: iter_times, dStam/dv and dStam/dw: " << iter_times << ", " << dStam_dv << ", " << dStam_dw << ".\n";
			}

			// 1. search direction h = [delta_v, delta_w]T
			double delta_v, delta_w;
			solve_func_get_h(dStam_dv, dStam_dw, Sk, clp, delta_v, delta_w);
			//std::cout << "delta_v " << delta_v << ", delta_w " << delta_w << ".\n";
			if (fabs(delta_v) > 1 || fabs(delta_w) > 1) {  //计算出来的方向可能是对的, 但是长度太大了, 超出0~1
				//double len = sqrt(delta_v * delta_v + delta_w * delta_w);
				delta_v = 0.9;// /= len;
				delta_w = 0.9; // /= len;
			}
			if (delta_v < 0.0001 && delta_w < 0.00001) {
				// 这情况可以看做是迭代终止了, 这么短的步长.
				if (cloop_evaluation.evaluate_irregular_patch(c, v, w, clp, OpenMesh::Vec3d(0,0,0), OpenMesh::Vec3d(0,0,0)) == false) {
					std::cout << "Error: stam98 test regular patch.\n";
				}   
				rt = Min_dis; break;				
			} 
			else {
				// 2. 初始步长.
				double step_len = -0.01; //-0.008;//1;
				if (i == 0) step_len *= 1;
				else step_len *= -1;

				if (delta_v < 0.0001 || delta_w < 0.00001) {
				} else {
					// 下面是使用brent minimization在求迭代的步长, 但从结果来看, 好像蛮力法来求也可以.
					//calc_step_side1_parames_ = OpenMesh::Vec4d(v, w, delta_v, delta_w);
					//calc_step_side1_sk_ = Sk;
					//calc_step_side1_c_ = c;
					//std::cout << currentDV->calc_step_side1_parames_ << ",~~~ " <<currentDV->calc_step_side1_sk_ << ".\n";

					//int status;
					//int iter = 0, max_iter = 100;
					//const gsl_min_fminimizer_type *T;
					//gsl_min_fminimizer *s;
					//double m = 0.01;
					////double a = (-0.9999-v)/delta_v, b = (0.9999-w)/delta_w;
					//double a = -999.999, b = 999.999;
					//gsl_function F;
					//std::cout << "a";
					//F.function = &calc_step_side1;
					//F.params = 0;

					//T = gsl_min_fminimizer_brent;std::cout << "b";
					//s = gsl_min_fminimizer_alloc (T); std::cout << "c";
					//if (v + m*delta_v < 0.000001 || w+m*delta_w < 0.000001 || v + m*delta_v + w+m*delta_w > 0.99999) {
					//	std::cout << "f";
					//} else { 
					//	status = gsl_min_fminimizer_set (s, &F, m, a, b); //设置所求的函数F, 初值m, 范围ab
					//	std::cout << "-status-" << status << " ";
					//	if (status == GSL_SUCCESS) {
					//		std::cout << "d";
					//		printf ("using %s method\n", gsl_min_fminimizer_name (s));std::cout << "e";
					//		printf ("%5s [%9s, %9s] %9s %9s\n",
					//			"iter", "lower", "upper", "min", "err(b-a)");

					//		printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
					//			iter, a, b, m, b - a);

					//		do
					//		{
					//			if (v + m*delta_v < 0.000001 || w+m*delta_w < 0.000001 || v + m*delta_v + w+m*delta_w > 0.99999) break;
					//			iter++;
					//			status = gsl_min_fminimizer_iterate (s);

					//			m = gsl_min_fminimizer_x_minimum (s);
					//			a = gsl_min_fminimizer_x_lower (s);
					//			b = gsl_min_fminimizer_x_upper (s);

					//			status = gsl_min_test_interval (a, b, 0.001, 0.0); // 测试if b-a < 0.001,是的话就status == GSL_SUCCESS

					//			if (status == GSL_SUCCESS) printf ("Converged:\n");

					//			printf ("%5d [%.7f, %.7f] " "%.7f %.7f\n",
					//				iter, a, b, m, b - a);

					//			step_len = m;
					//		}
					//		while (status == GSL_CONTINUE && iter < max_iter);	  
					//	} else {
					//		std::cout << "g";
					//	}

					//}
					//gsl_min_fminimizer_free (s);
				}
				//std::cout << "--\n";
				// 3. Update delta_v and delta_w.// 得到步长step side之后			
				delta_v *= step_len;
				delta_w *= step_len;

				if ((v + delta_v) < 0) {
					std::cout << "v < 0. "; //, " << v << ", + " << delta_v << ", " << _v << " ////单独没有问题, 和w<0一起没有问题, 但和>1时就有问题.
					if (cross_patch_depth < Max_Cross_Patch_Depth) {
						TriMesh::FHandle fh = initial_control_mesh_.face_handle(initial_control_mesh_.opposite_halfedge_handle(initial_control_mesh_.prev_halfedge_handle(h12)));
						double w_new, v_new;
						update_uvw_when_v_new_is_less_than_zero(v, delta_v, w, delta_w, v_new, w_new);
						OpenMesh::Vec3d bc_new(1-v_new-w_new, v_new, w_new);
						//std::cout << _v << ", " << v << ", " << delta_v << ". " << _w << ", "  w << ", " << delta_w << ".\n";
						if (DGP::is_valid_barycentric_coordinate(bc_new) == false) {
							std::cout << "Error: regular v < 0, the (u', v', w') is invalid.\n";  return false;
						}
						OpenMesh::Vec3d p = cN6[0] * (1 - w_new - v_new) + cN6[cN6.size() - 6] * w_new + cN6[cN6.size() - 6 - 1] * v_new;
						OpenMesh::Vec3d pp(0, 0, 0);
						if (evaluation_get_footpoint(initial_control_mesh_, fh, Sk, p,bc_new, pp, clp, ++cross_patch_depth) == false) {
							break; //不是返回错误, 而是找到当前为止, 当前的结果也是有效的. //return false;
						} 
					}
					rt = Out_vw; break; 				
				} else if ((w + delta_w) < 0) {
					std::cout << "w < 0. "; //\n// //单独没有问题, 和v<0一起没有问题, 和>1也没有问题.
					if (cross_patch_depth < Max_Cross_Patch_Depth) {
						TriMesh::FHandle fh = initial_control_mesh_.face_handle(initial_control_mesh_.opposite_halfedge_handle(h12));
						double v_new, w_new;
						update_uvw_when_w_new_is_less_than_zero(v, delta_v, w, delta_w, v_new, w_new);
						OpenMesh::Vec3d bc_new(1-v_new-w_new, v_new, w_new);
						//std::cout << _v << ", " << v << ", " << delta_v << ". " << _w << ", " << w << ", " << delta_w << ".\n";
						if (DGP::is_valid_barycentric_coordinate(bc_new) == false) {
							std::cout << "Error info: " << v << ", " << delta_v << ", " << w << ", " << delta_w << ".\n";
							std::cout << "Error: regular w < 0, the (u', v', w') is invalid, " << 1-v_new-w_new << ", " << v_new << ", " << w_new << ".\n"; 
							return false;
						}
						OpenMesh::Vec3d p = cN6[0] * (1 - w_new - v_new) + cN6[1] * v_new + cN6[2] * w_new;
						
						OpenMesh::Vec3d pp(0, 0, 0);
						if (evaluation_get_footpoint(initial_control_mesh_, fh, Sk, p, bc_new, pp, clp, ++cross_patch_depth) == false) {
							break; //不是返回错误, 而是找到当前为止, 当前的结果也是有效的. //return false;
						} 
					}				
					rt = Out_vw; break;
				} else if ((v + delta_v + w + delta_w) > 1) {
					std::cout << "v + w > 1. "; //\n//单独没有问题, 和w<0也没有问题, 但和v<0一起就有问题
					if (cross_patch_depth < Max_Cross_Patch_Depth) {
						TriMesh::FHandle fh = initial_control_mesh_.face_handle(initial_control_mesh_.opposite_halfedge_handle(initial_control_mesh_.next_halfedge_handle(h12)));
						double v_new, w_new;
						update_uvw_when_v_new_plus_w_new_is_greaterr_than_one(v, delta_v, w, delta_w, v_new, w_new);
						OpenMesh::Vec3d bc_new(1-v_new-w_new, v_new, w_new);
						if (DGP::is_valid_barycentric_coordinate(bc_new) == false) {
							std::cout << "Error: irregular v+w>1, the (u', v', w') is invalid.\n"; return false;
						}
						OpenMesh::Vec3d p = cN6[cN6.size() - 6 + 1] * (1 - w_new - v_new) + cN6[1] * v_new + cN6[cN6.size() - 6] * w_new;
						OpenMesh::Vec3d pp(0, 0, 0); 
						if (evaluation_get_footpoint(initial_control_mesh_, fh, Sk, p, bc_new, pp, clp, ++cross_patch_depth) == false) {
							break; //不是返回错误, 而是找到当前为止, 当前的结果也是有效的. //return false;
						}/**/
					}
					rt = Out_vw; break;
				} else if ((v + delta_v) > 1) { //好像没有出现过此情况.
					std::cout << "v > 1, " << v << ", + " << delta_v << ".\n"; rt = Out_vw; break;
				} else if ((w + delta_w) > 1) { //好像没有出现过此情况.
					std::cout << "w > 1.\n"; rt = Out_vw; break;
				} else { 
					//还在同一个面里面找着
					//std::cout << "v: " << v << ", w: " << w << ". ";
					v += delta_v; 
					w += delta_w;
					//std::cout << "v': " << v << ", w': " << w << ".\n";

					double distance_before = (clp-Sk).norm();
					OpenMesh::Vec3d clp_before = clp; //保存前一个最近点坐标 
					if (cloop_evaluation.evaluate_irregular_patch(c, v, w, clp, OpenMesh::Vec3d(0,0,0), OpenMesh::Vec3d(0,0,0)) == false) {
						std::cout << "Error: stam98 test regular patch.\n";
					}
					double distance_after = (clp-Sk).norm(); 
					if (distance_before <= distance_after) {
						//std::cout << "Info: Iteration ends, OK.\n"; //for distance_before < distance_after, 
						clp = clp_before; 
						rt = Min_dis; break;
					}
				} 
			} // delta_v and delta_w is meaningful.				
		} // end of for iter_times
		if (iter_times == Max_Iter_Times) {
			std::cout << "Info: Iteration ends Max_Iter_Times, OK.\n";//由于超过给定迭代次数, 
			rt = Max_times;
		} 
		clp_arr[i] = clp; //记下这一次closest position.
		if ((Sk - clp).norm() < (Sk - _clp).norm()) { //这就不用找第二次了.
			_clp = clp; return true;
		}
	} // end of for i == 0 or 1. 
	if ((Sk - clp_arr[0]).norm() < (Sk - clp_arr[1]).norm() ) _clp = clp_arr[0];
	else _clp = clp_arr[1];

	if (rt != Null) {
		return true;
	} else {
		return false;
	}	
} 
bool DecimationModel::h94loop_get_closest_points(Hoppe94LoopEval &h94loop_evaluation, bool _is_regular, int feature_type, 
												 const std::vector<OpenMesh::Vec3d> & _c,  const std::vector<bool> &h, 
												 const double _v, const double _w, OpenMesh::Vec3d Sk, OpenMesh::Vec3d & _clp, 
												 const TriMesh& _mesh, TriMesh::HHandle h84orh12, int cross_patch_depth) 
{
	if (feature_type == 0) return false;
	//std::cout << "Sk and initial clp: " << Sk << ", " << clp << ".\n";//for test
	//std::cout << "Initial dis " << (Sk - clp).norm() << ".\n";

	std::vector<OpenMesh::Vec3d> c = _c;
	if (_is_regular) { // feature regular, feature_type 7种(1~7), _c是12个顶点, h是15个特征边.
	} else { // feature irregular, feature type 17种(1~17), _c是N+6个顶点, h是N+9个特征边.

	}
	enum return_type { Null, Min_dis, Max_times, Out_vw, Little_delta };
	return_type rt = Null;

	const int Max_Iter_Times = 4000;//最大的迭代次数.这些参数都是可以调的.
	const double Max_Step_length = 0.5; //最大步长
	double distance_result = 100;
	double distance[2] = {100};
	OpenMesh::Vec3d clp_arr[2] = { OpenMesh::Vec3d(10000, 10000, 10000) }; //在实现的时候那个步长难以定(正负, 大小),所以我分正负求两次, 取效果好的那次.

	for (int i = 0; i < 2; ++i) {
		double v = _v, w = _w;
		OpenMesh::Vec3d clp = _clp;

		int iter_times = 0; 
		for (iter_times = 0; iter_times < Max_Iter_Times; ++iter_times) {		
			// 0. 求出雅克比矩阵的元素
			OpenMesh::Vec3d dStam_dv(0, 0, 0), dStam_dw(0, 0, 0); 
			if (_is_regular) {
				if (h94loop_evaluation.eval_regular_1in7(static_cast<feature_regular_type>(feature_type), c, h, v, w, 
					OpenMesh::Vec3d(0,0,0), dStam_dv, dStam_dw) == false) { 
						std::cerr << "Error: eval_regular_1in7 Der_V Der_W " << feature_type << ". \n"; 
				} 
			} else {
				// irregular
				if (h94loop_evaluation.eval_irregular_1in17(static_cast<feature_irregular_type>(feature_type), c, h, v, w, 
					OpenMesh::Vec3d(0,0,0), dStam_dv, dStam_dw) == false) {
						std::cerr << "Error: eval_irregular_1in17 Der_V Der_W " << feature_type << ".\n";
				} 
			}

			// 1. 解方程JTJh = -JTf得到迭代的方向.
			double delta_v, delta_w;//search direction h = [delta_v, delta_w]T
			solve_func_get_h(dStam_dv, dStam_dw, Sk, clp, delta_v, delta_w);
			//std::cout << "dS/dv " << dStam_dv << ", dS/dw " << dStam_dw << 
			//	"delta_v " << delta_v << ", delta_w " << delta_w << ".\n";//
			if ( (fabs(delta_v) < 1e-6 && fabs(delta_w) < 1e-6) || (delta_v * delta_v + delta_w * delta_w < 1e-12)) {
				rt = Little_delta;
				//std::cout << "Info: Too little delta.\n";
				if (_is_regular) {
					if (h94loop_evaluation.eval_regular_1in7((feature_regular_type)feature_type, c, h, v, w, 
						clp, OpenMesh::Vec3d(0,0,0), OpenMesh::Vec3d(0,0,0)) == false) {
							std::cerr << "Error: eval_regular_1in7 ORI " << feature_type << ".\n";
					} 
				} else { // irregular
					if (h94loop_evaluation.eval_irregular_1in17(static_cast<feature_irregular_type>(feature_type), c, h, v, w, 
						clp, OpenMesh::Vec3d(0,0,0), OpenMesh::Vec3d(0,0,0)) == false) {
							std::cerr << "Error: eval_irregular_1in17 Der_V " << feature_type << ".\n";
					}
				}
				break;
			}
			if (fabs(delta_v) > 1 || fabs(delta_w) > 1) {  //计算出来的方向可能是对的, 但是长度太大了, 超出0~1
				double len = sqrt(delta_v * delta_v + delta_w * delta_w);
				delta_v /= len;
				delta_w /= len;
			}

			// 2. 初始步长.
			double step_len = 0.01;//0.05;//; //-0.002;//1;
			if (i == 0) step_len *= 1;
			else step_len *= -1;

			delta_v *= step_len;
			delta_w *= step_len;
			if ((v + delta_v) < 0) {     
				//std::cout << "feature region v < 0.\n"; //
				if (cross_patch_depth < Max_Cross_Patch_Depth) {
					if (_is_regular) {
						TriMesh::HHandle h84 = h84orh12;
						TriMesh::FHandle fh = initial_control_mesh_.face_handle(initial_control_mesh_.opposite_halfedge_handle(h84));
						double w_new, v_new;
						update_uvw_when_v_new_is_less_than_zero(v, delta_v, w, delta_w, v_new, w_new);
						if (DGP::is_valid_barycentric_coordinate(OpenMesh::Vec3d(1-v_new-w_new, v_new, w_new)) == false) {
							std::cout << "Error: feature regular v < 0, the (u', v', w') is invalid.\n"; 
						}
						OpenMesh::Vec3d bc_new(1-w_new-v_new, v_new, w_new);
						OpenMesh::Vec3d p = c[3] * (1 - w_new - v_new) + c[7] * w_new + c[4] * v_new;
						OpenMesh::Vec3d pp(0, 0, 0);
						evaluation_get_footpoint(initial_control_mesh_, fh, Sk, p, bc_new, pp, clp, ++cross_patch_depth);

					} else {
						TriMesh::HHandle h12 = h84orh12;
						TriMesh::FHandle fh = initial_control_mesh_.face_handle(initial_control_mesh_.opposite_halfedge_handle(initial_control_mesh_.prev_halfedge_handle(h12)));
						double w_new, v_new;
						update_uvw_when_v_new_is_less_than_zero(v, delta_v, w, delta_w, v_new, w_new);
						if (DGP::is_valid_barycentric_coordinate(OpenMesh::Vec3d(1-v_new-w_new, v_new, w_new)) == false) {
							std::cout << "Error: feature regular v < 0, the (u', v', w') is invalid.\n"; 
						}
						OpenMesh::Vec3d bc_new(1-w_new-v_new, v_new, w_new);
						OpenMesh::Vec3d p = c[0] * (1 - w_new - v_new) + c[c.size() - 6] * w_new + c[c.size() - 6 - 1] * v_new;
						OpenMesh::Vec3d pp(0, 0, 0);
						evaluation_get_footpoint(initial_control_mesh_, fh, Sk, p, bc_new, pp, clp, ++cross_patch_depth);

					}/**/ 
				}				
				rt = Out_vw; break;
			} else if ((w + delta_w) < 0) { 
				//std::cout << "feature region w < 0.\n"; //
				if (cross_patch_depth < Max_Cross_Patch_Depth) {
					if (_is_regular) { 
						TriMesh::HHandle h84 = h84orh12;
						TriMesh::FHandle fh = initial_control_mesh_.face_handle(initial_control_mesh_.opposite_halfedge_handle(initial_control_mesh_.next_halfedge_handle(h84)));
						double v_new, w_new;
						update_uvw_when_w_new_is_less_than_zero(v, delta_v, w, delta_w, v_new, w_new);
						if (DGP::is_valid_barycentric_coordinate(OpenMesh::Vec3d(1-v_new-w_new, v_new, w_new)) == false) {
							std::cout << "Error: feature regular w < 0, the (u', v', w') is invalid, " << 1-v_new-w_new << ", " << v_new << ", " << w_new << ".\n"; 
						}
						OpenMesh::Vec3d bc_new(1-w_new-v_new, v_new, w_new);
						OpenMesh::Vec3d p = c[3] * (1 - w_new - v_new) + c[6] * v_new + c[2] * w_new;
						OpenMesh::Vec3d pp(0, 0, 0);
						evaluation_get_footpoint(initial_control_mesh_, fh, Sk, p, bc_new, pp, clp, ++cross_patch_depth); 
					} else {
						TriMesh::HHandle h12 = h84orh12;
						TriMesh::FHandle fh = initial_control_mesh_.face_handle(initial_control_mesh_.opposite_halfedge_handle(h12));
						double w_new, v_new;
						update_uvw_when_w_new_is_less_than_zero(v, delta_v, w, delta_w, v_new, w_new);
						if (DGP::is_valid_barycentric_coordinate(OpenMesh::Vec3d(1-v_new-w_new, v_new, w_new)) == false) {
							std::cout << "Error Info: " << v << ", " << delta_v << ", " << w << ", " << delta_w << ".\n";
							std::cout << "Error: feature irregular w < 0, the (u', v', w') is invalid, " << 1-v_new-w_new << ", " << v_new << ", " << w_new << ".\n"; 
						}
						OpenMesh::Vec3d bc_new(1-w_new-v_new, v_new, w_new);
						OpenMesh::Vec3d p = c[0] * (1 - w_new - v_new) + c[1] * v_new + c[2] * w_new;
						OpenMesh::Vec3d pp(0, 0, 0);
						evaluation_get_footpoint(initial_control_mesh_, fh, Sk, p, bc_new, pp, clp, ++cross_patch_depth);

					}/**/ 
				}				
				rt = Out_vw; break;
			} else if ((v + delta_v + w + delta_w) > 1) { 
				//std::cout << "featurn region v + w > 1." << ".\n"; //
				if (cross_patch_depth < Max_Cross_Patch_Depth) {
					if (_is_regular) {
						TriMesh::HHandle h84 = h84orh12;
						TriMesh::FHandle fh = initial_control_mesh_.face_handle(initial_control_mesh_.opposite_halfedge_handle(initial_control_mesh_.prev_halfedge_handle(h84)));
						double v_new, w_new;
						update_uvw_when_v_new_plus_w_new_is_greaterr_than_one(v, delta_v, w, delta_w, v_new, w_new);
						if (DGP::is_valid_barycentric_coordinate(OpenMesh::Vec3d(1-v_new-w_new, v_new, w_new)) == false) {
							std::cout << "Error: feature irregular v+w>1, the (u', v', w') is invalid.\n"; 
						}
						OpenMesh::Vec3d bc_new(1-w_new-v_new, v_new, w_new);
						OpenMesh::Vec3d p = c[10] * bc_new[0] + c[6] * v_new + c[7] * w_new;
						OpenMesh::Vec3d pp(0, 0, 0);
						evaluation_get_footpoint(initial_control_mesh_, fh, Sk, p, bc_new, pp, clp, ++cross_patch_depth);

					} else {  
						TriMesh::HHandle h12 = h84orh12;
						TriMesh::FHandle fh = initial_control_mesh_.face_handle(initial_control_mesh_.opposite_halfedge_handle(initial_control_mesh_.next_halfedge_handle(h12)));
						double v_new, w_new;
						update_uvw_when_v_new_plus_w_new_is_greaterr_than_one(v, delta_v, w, delta_w, v_new, w_new);
						if (DGP::is_valid_barycentric_coordinate(OpenMesh::Vec3d(1-v_new-w_new, v_new, w_new)) == false) {
							std::cout << "Error: feature irregular v+w>1, the (u', v', w') is invalid.\n"; 
						}
						OpenMesh::Vec3d bc_new(1-w_new-v_new, v_new, w_new);
						OpenMesh::Vec3d p = c[c.size() - 6 + 1] * bc_new[0] + c[1] * v_new + c[c.size() - 6] * w_new;
						OpenMesh::Vec3d pp(0, 0, 0); 
						evaluation_get_footpoint(initial_control_mesh_, fh, Sk, p, bc_new, pp, clp, ++cross_patch_depth);
					}/**/
				}				
				rt = Out_vw; break;
			} else if ((v + delta_v) > 1) { std::cout << "v > 1.\n"; rt = Out_vw; break; //这后两种情况好像没有出现的.
			} else if ((w + delta_w) > 1) { std::cout << "w > 1.\n"; rt = Out_vw; break;
			} else { //还在同一个面里面找着
				//std::cout << "SFG. "; //Same feature region//
				//std::cout << "v: " << v << ", w: " << w << ". ";////
				v += delta_v; 
				w += delta_w; 
				//std::cout << "v': " << v << ", w': " << w << ".\n";////

				double distance_before = (Sk - clp).norm();
				OpenMesh::Vec3d clp_before = clp; //保存前一个最近点坐标

				if (_is_regular) {
					if (h94loop_evaluation.eval_regular_1in7((feature_regular_type)feature_type, c, h, v, w, 
						clp, OpenMesh::Vec3d(0,0,0), OpenMesh::Vec3d(0,0,0)) == false) {
							std::cerr << "Error: eval_regular_1in7 ORI " << feature_type << ".\n";
					} 
				} else { // irregular
					if (h94loop_evaluation.eval_irregular_1in17(static_cast<feature_irregular_type>(feature_type), c, h, v, w, 
						clp, OpenMesh::Vec3d(0,0,0), OpenMesh::Vec3d(0,0,0)) == false) {
							std::cerr << "Error: eval_irregular_1in17 Der_V " << feature_type << ".\n";
					}
				}
				double distance_after = (Sk - clp).norm();//
				//std::cout << "clp': " << clp << ". dis(clp'-Sk): " << distance_after << "; " << distance_before << ". " << iter_times << ".\n";
				//std::cout << distance_after << "; " << distance_before << ". " << iter_times << ", " << v << ", " << w << ".\n";
				if (distance_before <= distance_after) { //
					//std::cout << "Info: Iteration ends, OK." << distance_before << ".\n";//for distance_before < distance_after, 
					clp = clp_before; ////
					rt = Min_dis; break;
				} 
			} // end of if-else v and w's range
		} // end of for iter_times
		// 退出for iter_times的情况:
		// (1) distance after > distance before, 属于正常退出.
		// (2) 达到最大迭代次数, 属于正常退出.
		if (iter_times >= Max_Iter_Times) {////
			std::cout << "Info: Iteration ends Max_Iter_Times, OK." << ".\n";// 由于超过给定迭代次数
			rt = Max_times;//
		} 
		// (3) v, w out of [0, 1] range. 也是正常退出.
		// (4) delta_v, delta_w 太小了也正常退出.

		clp_arr[i] = clp;
		// 提前退出. 但是是不是这一步距离减少了就意味着是最比另一个方向的要小呢?
		if ((Sk - clp).norm() < (Sk - _clp).norm()) {
			_clp = clp; return true;
		}
	} // end of for i == 0, 1
	if ((Sk - clp_arr[0]).norm() < (Sk - clp_arr[1]).norm() ) _clp = clp_arr[0];
	else _clp = clp_arr[1];

	if (rt != Null) return true;
	else return false;	
}
void DecimationModel::evaluation_get_face_type(TriMesh::FaceHandle fh, 
											   bool &_is_feature, bool &_is_regular, int &feature_type, TriMesh::HalfedgeHandle &h) 
{   //这个函数没有放到Hoppe94LoopEval.h中是因为那里的feature_type还分为是否regular, 返回值太烦.
	std::vector<TriMesh::HHandle> h_arr(3, TriMesh::HHandle(-1));
	std::vector<TriMesh::VHandle> v_arr(3, TriMesh::VHandle(-1));
	h_arr[0] = initial_control_mesh_.halfedge_handle(fh); v_arr[0] = initial_control_mesh_.to_vertex_handle(h_arr[0]);
	h_arr[1] = initial_control_mesh_.next_halfedge_handle(h_arr[0]); v_arr[1] = initial_control_mesh_.to_vertex_handle(h_arr[1]);
	h_arr[2] = initial_control_mesh_.next_halfedge_handle(h_arr[1]); v_arr[2] = initial_control_mesh_.to_vertex_handle(h_arr[2]);

	if (initial_control_mesh_.property(vp_type_icm_, v_arr[0]) == DGP::SMOOTH_VFT && initial_control_mesh_.property(vp_type_icm_, v_arr[1]) == DGP::SMOOTH_VFT 
		&& initial_control_mesh_.property(vp_type_icm_, v_arr[2]) == DGP::SMOOTH_VFT) { 
			TriMesh::HalfedgeHandle h12(-1); //没有尖锐特征时候,将regular patch也统一为irregular patch来做.
			for (int i = 0; i < 3; ++i) {
				if (initial_control_mesh_.valence(v_arr[i]) != 6) {
					h12 = initial_control_mesh_.next_halfedge_handle(h_arr[i]); break;//有非规则点的情况
				}
			}
			if (h12.is_valid() == false) h12 = h_arr[1]; //没有非规则点的情况.

			_is_feature = false; _is_regular = false; feature_type = 0; h = h12;
	}  
	else { // feature situation.
		if (initial_control_mesh_.valence(v_arr[0]) == 6 && initial_control_mesh_.valence(v_arr[1]) == 6 && initial_control_mesh_.valence(v_arr[2]) == 6) {  
			// feature regular
			enum feature_regular_type { tNull = 0, t1v0e = 1, t2v0e, t2v1e, t3v0e, t3v1e, t3v2e, t3v3e };//1~7有意义
			feature_regular_type frt = tNull;
			TriMesh::HalfedgeHandle h84(-1);

			int n_feature_vertics = 0;
			if (initial_control_mesh_.property(vp_type_icm_, v_arr[0]) != DGP::SMOOTH_VFT) ++ n_feature_vertics;
			if (initial_control_mesh_.property(vp_type_icm_, v_arr[1]) != DGP::SMOOTH_VFT) ++ n_feature_vertics;
			if (initial_control_mesh_.property(vp_type_icm_, v_arr[2]) != DGP::SMOOTH_VFT) ++ n_feature_vertics;
			if (n_feature_vertics == 0) std::cout << "Error: no feature stituation 应该在上面就已处理的了.\n";
			else if (n_feature_vertics == 1) { // 1种情况
				if (initial_control_mesh_.property(vp_type_icm_, v_arr[0]) != DGP::SMOOTH_VFT) h84 = h_arr[0];
				if (initial_control_mesh_.property(vp_type_icm_, v_arr[1]) != DGP::SMOOTH_VFT) h84 = h_arr[1];
				if (initial_control_mesh_.property(vp_type_icm_, v_arr[2]) != DGP::SMOOTH_VFT) h84 = h_arr[2];

				frt = t1v0e;
			} else if (n_feature_vertics == 2) { // 2种情况
				if (initial_control_mesh_.property(vp_type_icm_, v_arr[0]) != DGP::SMOOTH_VFT && initial_control_mesh_.property(vp_type_icm_, v_arr[1]) != DGP::SMOOTH_VFT) {
					h84 = h_arr[0]; 
					if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[1])).feature()) {
						frt = t2v1e;  
					} else {
						frt = t2v0e;  
					}						
				}
				if (initial_control_mesh_.property(vp_type_icm_, v_arr[1]) != DGP::SMOOTH_VFT && initial_control_mesh_.property(vp_type_icm_, v_arr[2]) != DGP::SMOOTH_VFT) {
					h84 = h_arr[1]; 
					if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[2])).feature()) {
						frt = t2v1e;  
					} else {
						frt = t2v0e;  
					}
				}
				if (initial_control_mesh_.property(vp_type_icm_, v_arr[2]) != DGP::SMOOTH_VFT && initial_control_mesh_.property(vp_type_icm_, v_arr[0]) != DGP::SMOOTH_VFT) {
					h84 = h_arr[2]; 
					if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[0])).feature()) {
						frt = t2v1e;  
					} else {
						frt = t2v0e;  
					}
				}

			} else if (n_feature_vertics == 3) { // 4种情况
				if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[0])).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[1])).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[2])).feature()) {
					h84 = h_arr[0]; frt = t3v3e;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[0])).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[1])).feature() ) { //两个特征边.
					h84 = h_arr[0]; frt = t3v2e;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[1])).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[2])).feature() ) {
					h84 = h_arr[1]; frt = t3v2e;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[2])).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[0])).feature() ) {
					h84 = h_arr[2]; frt = t3v2e;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[0])).feature()) { //一个特征边.
					h84 = h_arr[2]; frt = t3v1e;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[1])).feature()) {
					h84 = h_arr[0]; frt = t3v1e;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[2])).feature()) {
					h84 = h_arr[1]; frt = t3v1e;
				} else { //没有特征边.
					h84 = h_arr[0]; frt = t3v0e;
				} 					
			}  
			if (h84.is_valid() == false) std::cout << "Error: 本不该有的错误, 上面应该已经列举了所有可能情况了.\n";

			_is_feature = true; _is_regular = true; feature_type = (int)frt; h = h84;			
		} else { // irregular
			int extraordinary_vertex_index = -1;
			TriMesh::HHandle h12;
			if (initial_control_mesh_.valence(v_arr[0]) != 6) {
				extraordinary_vertex_index = 0; 
				h12 = h_arr[1];
			}
			if (initial_control_mesh_.valence(v_arr[1]) != 6) {
				extraordinary_vertex_index = 1; 
				h12 = h_arr[2];
			}
			if (initial_control_mesh_.valence(v_arr[2]) != 6) {
				extraordinary_vertex_index = 2; 
				h12 = h_arr[0];
			}
			if (h12.is_valid() == false) {
				std::cout << "Error: feature irregular. 除非3个点都是regular否则不会出现这情况: .\n";
				// << initial_control_mesh_.valence(v_arr[0]) << ", " << initial_control_mesh_.valence(v_arr[1]) << ", " << initial_control_mesh_.valence(v_arr[2]) << "\n";
			}

			enum feature_irregular_type { tNull = 0, 
				t1v0e_up,	t1v0e_left,	t1v0e_right, 
				t2v0e_left, t2v0e_right, t2v0e_bottom,
				t3v0e,
				t2v1e_left, t2v1e_right, t2v1e_bottom,
				t3v1e_left, t3v1e_right, t3v1e_bottom,
				t3v2e_up,	t3v2e_left, t3v2e_right, 
				t3v3e
			};
			feature_irregular_type ft = tNull;
			TriMesh::HHandle h12prev = initial_control_mesh_.prev_halfedge_handle(h12), h12next(initial_control_mesh_.next_halfedge_handle(h12));
			TriMesh::VHandle vhup = initial_control_mesh_.to_vertex_handle(h12prev), vhleft(initial_control_mesh_.to_vertex_handle(h12)), vhright(initial_control_mesh_.to_vertex_handle(h12next));
			if (initial_control_mesh_.property(vp_type_icm_, vhup) == DGP::SMOOTH_VFT) {
				if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12next)).feature()) {
					ft = t2v1e_bottom;
				} else if (initial_control_mesh_.property(vp_type_icm_, vhleft) != DGP::SMOOTH_VFT && initial_control_mesh_.property(vp_type_icm_, vhright) != DGP::SMOOTH_VFT) {
					ft = t2v0e_bottom;
				}  else if (initial_control_mesh_.property(vp_type_icm_, vhleft) != DGP::SMOOTH_VFT ) {
					ft = t1v0e_left;
				}  else if (initial_control_mesh_.property(vp_type_icm_, vhright) != DGP::SMOOTH_VFT ) {
					ft = t1v0e_right;
				} 
			} else {
				if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12prev)).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12)).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12next)).feature()) {
					ft = t3v3e;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12prev)).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12)).feature() ) {
					ft = t3v2e_up; //这个在例子ne....160.off中没有错.
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12)).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12next)).feature()) {
					ft = t3v2e_left;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12next)).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12prev)).feature()) {
					ft = t3v2e_right;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12next)).feature()) {
					ft = t3v1e_bottom;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12)).feature()) {
					if (initial_control_mesh_.property(vp_type_icm_, vhright) != DGP::SMOOTH_VFT) ft = t3v1e_left;
					else ft = t2v1e_left; // nefertiti.160.off来看也没有错.
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12prev)).feature()) {
					if (initial_control_mesh_.property(vp_type_icm_, vhleft) != DGP::SMOOTH_VFT) ft = t3v1e_right;
					else ft = t2v1e_right;// nefertiti.160.off来看也没有错.
				} else if (initial_control_mesh_.property(vp_type_icm_, vhleft) != DGP::SMOOTH_VFT && initial_control_mesh_.property(vp_type_icm_, vhright) != DGP::SMOOTH_VFT) {
					ft = t3v0e;
				}  else if (initial_control_mesh_.property(vp_type_icm_, vhleft) != DGP::SMOOTH_VFT ) {
					ft = t2v0e_left;
				}  else if (initial_control_mesh_.property(vp_type_icm_, vhright) != DGP::SMOOTH_VFT ) {
					ft = t2v0e_right;
				} else if (initial_control_mesh_.property(vp_type_icm_, vhleft) == DGP::SMOOTH_VFT && initial_control_mesh_.property(vp_type_icm_, vhright) == DGP::SMOOTH_VFT) {
					ft = t1v0e_up;
				}  
			}

			_is_feature = true; _is_regular = false; feature_type = (int)ft; h = h12;	
		} //end of if-else regular feature or irregular feature
	} // end of if-else not feature or feature
}
// 求顶点Sk在平面fh上的投影点p对应极限曲面上的精确坐标_pp以及最近点坐标_clp. 注意p是肯定是会投影到fh上的的, 之前我这里犯错了.
bool DecimationModel::evaluation_get_footpoint(const TriMesh& _mesh, TriMesh::FaceHandle _fh_of_initial_control_mesh_, //在此网格的此面上
											   OpenMesh::Vec3d Sk, OpenMesh::Vec3d p, // input: 采样点和迭代搜索的初值
											   OpenMesh::Vec3d _bc, // input: p 相对于_fh_of_initial_control_mesh_的参数坐标
											   OpenMesh::Vec3d &_pp, OpenMesh::Vec3d &_clp, const int cross_patch_depth) // output:前两个.
{ 
	TriMesh::FHandle fh = _fh_of_initial_control_mesh_;

	
	// 先判断面的类型
	bool is_feature = false;// false,没有尖锐特征, true, 有尖锐特征.
	bool is_regular = true; // is regular or irregular.
	int feature_type = 0;   // is_feature == true时候, regular face has 7 cases, irregular 17 cases.
	TriMesh::HalfedgeHandle h; // regular, h84. irregular h12.
	evaluation_get_face_type(fh, is_feature, is_regular, feature_type, h);
	// is_feature, is_regular, feature_type, h
	// false       false       0             h12  //当没有特征的情况下，统一当做非规则面来处理.
	// true        true        1~7           h84   
	// true        false       1~17          h12   
	//if (h.is_valid() == false) std::cout << "Error: h is non-valid.\n";//这不该出现的, 函数应该很简单啊.
	//std::cout << is_feature << ", " << is_regular << ", " << feature_type << ". ";
	//if (is_feature == false || is_regular == false ) continue; 

	if (is_feature == false) { 
		TriMesh::HalfedgeHandle h12 = h; //没有尖锐特征时候,将regular patch也统一为irregular patch来做.
		int N = initial_control_mesh_.valence(initial_control_mesh_.from_vertex_handle(h12));

		std::vector<OpenMesh::Vec3d> cN6(N+6, OpenMesh::Vec3d(0, 0, 0));
		DGP::cloop_stam98_irregular_get_cN6(initial_control_mesh_, h12, cN6);

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(cN6[0], cN6[1], cN6[N], p); // 理应跟参数输入时候的_bc一致!
		if (DGP::is_valid_barycentric_coordinate(bc, true) == false) {
			std::cout << "Error: evaluation_get_footpoint(...): non-feature bc: " << bc << ".\n";
			std::cout << "contu.." << cN6[0] << ", " << cN6[1] << ", " << cN6[N] << ", " << p << ".\n";
			return false;
		}
		double v = bc[1], w = bc[2]; 
		// 求在极限曲面上的evaluation position
		OpenMesh::Vec3d pp(0, 0, 0);
		if (cloop_evaluation.evaluate_irregular_patch(cN6, v, w, pp, OpenMesh::Vec3d(0,0,0), OpenMesh::Vec3d(0,0,0)) == false) {
			std::cout << "Error: 这里错了." << v << ", " << w << " \n";
		}
		_pp = pp;

		// 求在极限曲面上的cloesest point position
		OpenMesh::Vec3d clp = pp; //最近点的初始坐标为evaluation坐标.
		if (cloop_stam98_get_closest_points(cN6, v, w, Sk, clp, initial_control_mesh_, h12, cross_patch_depth)) {
			_clp = clp;
			return true;
		} else { //没有正确求出最近点.
			std::cout << "Error: evaluation_get_footpoint(...): non-feature region closest point.\n";
			return false;
		}

	} else { // feature situation.// is_feature == true;

		if (is_regular == true) {
			TriMesh::HHandle h84 = h;
			std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));// 严重注意: 这12个控制点的排序以stam98 "evaluation of loop subdivision surfaces". Fig1为准.
			std::vector<bool> fe15(15, false);//12个顶点之间的15个边是否是特征边
			DGP::h94loop_eval_regular_get_c12_h15(initial_control_mesh_, h84, c12, fe15);

			OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(c12[3], c12[6], c12[7], p);
			if (DGP::is_valid_barycentric_coordinate(bc, true) == false) {
				std::cout << "Error: evaluation_get_footpoint(...): feature regular bc: " << bc << ", _bc " << _bc << ".\n";
			}
			double v = bc[1], w = bc[2];
			//std::cout << "New bc: " << 1-bc[1]-bc[2] << ", " << bc[1] << ", " << bc[2] << ". ";

			OpenMesh::Vec3d pp(0, 0, 0);
			if (h94loop_evaluation.eval_regular_1in7((feature_regular_type)feature_type, c12, fe15, v, w, 
				pp, OpenMesh::Vec3d(0,0,0), OpenMesh::Vec3d(0,0,0)) == false) {
					std::cout << "Error: feature regular, pp, " <<feature_type << ".\n";
					//return false;
			}
			_pp = pp;
			OpenMesh::Vec3d clp = pp;
			// 根据类型求类型求最近点
			if (fabs(v) < 1e-5 || fabs(w) < 1e-5) { //特征点crease vertices被参数化到特征边上
				//这情况我就不求了, 特殊地将evaluation坐标作为最近点坐标. 
				//原因是这时w=0而v=-1.#IND, 之后求出来的就是无意义的-1.#IND和1.QNAN了.
				//std::cout << "Crease? " << mesh_.property(vp_type_, *it) << ".\n";
				//exit(-1); //for test
				_clp = clp;
				return true;
			} else {
				if (h94loop_get_closest_points(h94loop_evaluation, true, feature_type, c12, fe15, v, w, Sk, clp, initial_control_mesh_, h84, cross_patch_depth) ) {
					_clp = clp; //
					return true; //  
				} else {
					std::cout << "Error: evaluation_get_footpoint(...): feature regular closest point.\n";
					return false;
				}  
			}
		} else { // irregular

			TriMesh::HHandle h12 = h;
			int N = initial_control_mesh_.valence(initial_control_mesh_.from_vertex_handle(h12));  
			if (N == 6) std::cout << "Error: 这里应该是irregular.\n";//其实是不该有此情况的.
			int K = N + 6;
			std::vector<OpenMesh::Vec3d> c(K, OpenMesh::Vec3d(0, 0, 0));// 严重注意: 这K个控制点的排序以stam98 "evaluation of loop subdivision surfaces". Fig2为准.
			std::vector<bool> hN9(N+4+4+1, false);

			DGP::h94loop_eval_irregular_get_cN6_hN9(initial_control_mesh_, h12, c, hN9);

			OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(c[0], c[1], c[N], p);
			if (DGP::is_valid_barycentric_coordinate(bc) == false) {
				std::cout << "Error: evaluation_get_footpoint(...): feature irregular bc: " << bc << ", _bc " << _bc << ".\n";
			}
			double v = bc[1], w = bc[2];

			OpenMesh::Vec3d pp(0, 0, 0);
			h94loop_evaluation.eval_irregular_1in17(feature_irregular_type(feature_type), c, hN9, v, w, 
				pp, OpenMesh::Vec3d(0,0,0), OpenMesh::Vec3d(0,0,0));
			_pp = pp; 
			//std::cout << pp << ", " << (Sk - pp).norm() << "; " << (Sk - _clp).norm() << ".\n";
			OpenMesh::Vec3d clp = pp;
			if (fabs(v) < 1e-6 || fabs(w) < 1e-6) { //特征点crease vertices被参数化到特征边上
				_clp = clp;
				return true;
			} else {
				if (h94loop_get_closest_points(h94loop_evaluation, false, feature_type, c, hN9, v, w, Sk, clp, initial_control_mesh_, h12, cross_patch_depth) ) {
					_clp = clp;
					return true; 
				} else {
					std::cout << "Error: evaluation_get_footpoint(...): feature irregular closest point.\n";
					return false;
				} 
			}						

		} //if-else feature regular-irregular
	} // end of if-else 面的类型

	return false;
}
bool DecimationModel::evaluation_get_all_footpoints() {
	std::cout << "evaluation_get_all_footpoints(), begin.\n";
	// -----------下面用QAS曲面先逼近拟合曲面, 然后在QAS曲面上求最近点
	// 给每一个面指定一个明确的半边, 这样有利于明确一个三角形中三个顶点的顺序.
	initial_control_mesh_.add_property(fp_hh_icm_);
	for (TriMesh::FaceIter f_it(initial_control_mesh_.faces_begin()), f_end(initial_control_mesh_.faces_end()); f_it != f_end; ++f_it) {
		// Note: no boundary here
		TriMesh::FaceHalfedgeIter fh_it(initial_control_mesh_, f_it);
		TriMesh::HHandle hx(fh_it.handle()), hy = initial_control_mesh_.next_halfedge_handle(hx), hz = initial_control_mesh_.next_halfedge_handle(hy);
		int n_count_feature_edges = 0;
		if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(hx)).feature()) ++ n_count_feature_edges;
		if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(hy)).feature()) ++ n_count_feature_edges;
		if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(hz)).feature()) ++ n_count_feature_edges;

		//std::cout << n_count_feature_edges << " " ;
		TriMesh::HHandle h01; //
		if (n_count_feature_edges == 0) h01 = hx; // hy, hz is ok too.
		else if (n_count_feature_edges == 1 || n_count_feature_edges == 3) { // 记录下这个特征边
			if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(hx)).feature()) h01 = hx;
			else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(hy)).feature()) h01 = hy;
			else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(hz)).feature()) h01 = hz; 
		} else if (n_count_feature_edges == 2) { //选择两个特征边中的一个.
			if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(hx)).feature() == false) h01 = hz;
			else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(hy)).feature() == false) h01 = hx;
			else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(hz)).feature() == false) h01 = hy;  
		} if (h01.is_valid() == false) std::cerr << "Error: h01 is supposed to be valid.\n";
		//记录下3个点的枚举顺序, 并且也记录下了要是存在特征边时候那些边就是特征边..
		initial_control_mesh_.property(fp_hh_icm_, f_it) = h01;
	}

	qas_ = new DGP::QuadraticBezierTriangles;
	qas_->create_geometry(initial_control_mesh_, fp_hh_icm_, vp_type_icm_);
	std::cout << "Info: QAS crease geometry end.\n";

	TriMesh limit = initial_control_mesh_;
	OpenMesh::VPropHandleT<DGP::v_feature_type> vp_type_limit;
	limit.add_property(vp_type_limit);
	for (TriMesh::VIter v_it(limit.vertices_begin()), v_end(limit.vertices_end()); v_it != v_end; ++v_it) {
		limit.property(vp_type_limit, v_it) = initial_control_mesh_.property(vp_type_icm_, v_it);
	}
	DGP::Hoppe94LoopSubT<TriMesh> subdi_obj;
	subdi_obj.attach(limit, vp_type_limit);
	subdi_obj.set_limit_position();
	subdi_obj.detach();
	limit.update_normals();

	int try_face_num = 0; 
	int method_1_count =0, method_3_count =0;
	int one_pro =0; int two_pro=0; int three_pro =0; int four_pro = 0;
	int where_get_pro = -1;
	// initial_control_mesh_经过一次hoppe94 loop细分之后结构跟refined_simplified_mesh_是一样的。
	//因为initial_control_mesh_在求出来之后经过了一次细分, 其拓扑与refined_simplified_mesh_是一样,
	//利用这一点, 对于initial_control_mesh_先判断某个面的类型, 再借用refined_simplified_mesh_.property(fvset_rsm_, fh)
	//获得此面的包含的原始顶点, 并利用其重心坐标开始求其在极限曲面上的坐标.
	for (TriMesh::FIter f_it(initial_control_mesh_.faces_begin()), f_end(initial_control_mesh_.faces_end()); f_it != f_end; ++f_it) {
		if (initial_control_mesh_.is_boundary(f_it.handle(), true)) {
			std::cout << "Error: Sorry, 现在还不能处理边界面.\n";
			continue; //暂时不能处理边界面
		}
		//std::cout << "in face " << f_it.handle().idx() << ".\n";

		// 控制网格每一个面包含了哪些原始顶点(被deleted掉的), 这些点的footpoint是什么.
		for (std::vector<TriMesh::VHandle>::const_iterator it = refined_simplified_mesh_.property(fvset_rsm_, f_it).begin(),
			it_end = refined_simplified_mesh_.property(fvset_rsm_, f_it).end(); it != it_end; ++it) { 
				// 这里其实对于crease特征顶点可以直接给出其pp和clp.

				OpenMesh::Vec3d Sk = mesh2_.point(*it);// *it是mesh_上的顶点.作为采样点.
				OpenMesh::Vec3d vp_eval_pos, vp_closest_pos;

				// 方法一.
				//std::cout << "--------Use projection as initial values.\n";
				// (1). 下面尝试用投影的方法求原始顶点Sk在initial_control_mesh_的一个三角面上的投影坐标, 作为初始参数化坐标
				// 这个初始参数化坐标用于求Sk在极限曲面上的精确坐标, 以及最近点坐标.
				OpenMesh::Vec3d p(0, 0, 0);//初始参数化坐标
				OpenMesh::Vec3d bc_proj(0,0,0); //
				bool is_project_ok = false;
				TriMesh::FHandle fh_porject_ok = f_it.handle();
				// 首先*it是否落在initial_control_mesh_的f_it面上, 否则查询f_it的邻域.
				TriMesh::HHandle he_of_face = initial_control_mesh_.property(fp_hh_icm_, f_it);
				TriMesh::VHandle v0 = limit.from_vertex_handle(he_of_face);
				TriMesh::VHandle v1 = limit.to_vertex_handle(he_of_face);
				TriMesh::VHandle v2 = limit.to_vertex_handle(limit.next_halfedge_handle(he_of_face));
				OpenMesh::Vec3d n_of_f_it = limit.calc_face_normal(f_it);//面的法向
				double d = dot(Sk - limit.point(v0), n_of_f_it);//点Sk到面的距离.
				OpenMesh::Vec3d point_on_plane = Sk - (d /*>= 0 ? d: -1.0*d*/) * n_of_f_it;

				OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(limit.point(v0), limit.point(v1), limit.point(v2), point_on_plane);

			
				if (DGP::is_valid_barycentric_coordinate(bc)) {
					is_project_ok = true;
					fh_porject_ok = f_it.handle();
					p = point_on_plane;
					bc_proj = bc;
					
					one_pro ++; where_get_pro = 1;
				} else { 
					// 找f_it面的3个邻域面.
					for (TriMesh::FHIter fh_it(limit, f_it); fh_it && is_project_ok == false; ++fh_it) {
						TriMesh::FHandle of = limit.face_handle(limit.opposite_halfedge_handle(fh_it));
						he_of_face = initial_control_mesh_.property(fp_hh_icm_, of);
						v0 = limit.from_vertex_handle(he_of_face);
						v1 = limit.to_vertex_handle(he_of_face);
						v2 = limit.to_vertex_handle(limit.next_halfedge_handle(he_of_face));

						OpenMesh::Vec3d n = limit.calc_face_normal(of);
						if ((dot(n_of_f_it, n)) < 0.5) continue; // < -0.5 这个邻面与f_it面几乎翻转了, 跳过这个面..
						double d = dot(Sk - limit.point(v0), n);//点Sk都这个面的距离
						OpenMesh::Vec3d point_on_plane = Sk - (d/* >= 0 ? d: -1.0*d*/) * n;

						OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(limit.point(v0), limit.point(v1), limit.point(v2), point_on_plane);
						if (DGP::is_valid_barycentric_coordinate(bc)) {
							is_project_ok = true;
							fh_porject_ok = of;
							p = point_on_plane;
							bc_proj = bc; two_pro ++;where_get_pro = 2;
							
							break;
						} 
					}
					// 上面找不到的话, 再在面的3个点的一邻域
					he_of_face = initial_control_mesh_.property(fp_hh_icm_, f_it);
					v0 = limit.from_vertex_handle(he_of_face);
					v1 = limit.to_vertex_handle(he_of_face);
					v2 = limit.to_vertex_handle(limit.next_halfedge_handle(he_of_face));
					std::vector<TriMesh::VHandle> v_array;
					v_array.push_back(v0); v_array.push_back(v1); v_array.push_back(v2); 
					for (int i = 0; i < 3 && is_project_ok == false; ++i) {
						for (TriMesh::VOHIter voh_it(limit, v_array[i]); voh_it && is_project_ok == false; ++ voh_it) {
							OpenMesh::Vec3d n = limit.calc_face_normal(limit.face_handle(voh_it));
							if ((dot(n_of_f_it, n)) < 0.5) continue; // < -0.5 这个邻面与f_it面几乎翻转了, 跳过这个面..
							double d = dot(Sk - limit.point(v_array[i]), n);//点Sk都这个面的距离
							OpenMesh::Vec3d point_on_plane = Sk - (d/* >= 0 ? d: -1.0*d*/) * n;

							OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(limit.point(v_array[i]), 
								limit.point(limit.to_vertex_handle(voh_it)), 
								limit.point(limit.to_vertex_handle(limit.next_halfedge_handle(voh_it))), point_on_plane);
							if (DGP::is_valid_barycentric_coordinate(bc)) {
								is_project_ok = true;
								fh_porject_ok = limit.face_handle(voh_it);
								p = point_on_plane;
								bc_proj = bc; three_pro ++; where_get_pro = 3 ;

								//std::cout << p << "/" << (limit.point(v_array[i])*bc[0]+limit.point(limit.to_vertex_handle(voh_it))*bc[1]+limit.point(limit.to_vertex_handle(limit.next_halfedge_handle(voh_it)))*bc[2]) << ", ";

								break;
							}
						}
					}				
				}
				//is_project_ok = false;// forse to use this 
				if (is_project_ok == false) {// 实在最后不行的话就用回这个
					// 上面投影得到的结果好.  但经测试, 这里不会造成bc无效, 只会是初值不好，得到最后的结果不够小.
					is_project_ok = true;
					fh_porject_ok = f_it.handle();
					p = limit.point(mesh_.property(vf0_of_rsm_, *it)) * mesh_.property(vbc_of_rsm_, *it)[0]
					+ limit.point(mesh_.property(vf1_of_rsm_, *it)) * mesh_.property(vbc_of_rsm_, *it)[1]
					+ limit.point(mesh_.property(vf2_of_rsm_, *it)) * mesh_.property(vbc_of_rsm_, *it)[2]; 
					bc_proj = mesh_.property(vbc_of_rsm_, *it);
					four_pro ++;
				} 
				// make sure that we have at least projection successed.
				mesh_.property(vp_project_, *it) = true;
				mesh_.property(vp_project_pos_, *it) = p;

				//// (2). Sk在极限曲面上的evaluation position and closest position.
				//// 采样点Sk在initial_control_mesh_的面fh上的参数坐标是p,
				TriMesh::HHandle hh = initial_control_mesh_.property(fp_hh_icm_, fh_porject_ok);//这就暗示了参数中的_mesh其实是initial_control_mesh_.
				v0 = initial_control_mesh_.from_vertex_handle(hh);
				v1 = initial_control_mesh_.to_vertex_handle(hh);
				v2 = initial_control_mesh_.to_vertex_handle((initial_control_mesh_.next_halfedge_handle(hh)));
				bc = DGP::calc_barycentric_coordinates(limit.point(v0), limit.point(v1), limit.point(v2), p);

				//	std::cout << p << "/" << (limit.point(v0)*bc[0]+limit.point(v1)*bc[1]+limit.point(v2)*bc[2]) << ", ";

				if (DGP::is_valid_barycentric_coordinate(bc)) {
					mesh_.property(vp_eval_pos_, *it) = qas_->pos(fh_porject_ok, bc[0], bc[1], bc[2], true);
					mesh_.property(vp_eval_, *it) = true; 
					//clp = mesh_.property(vp_eval_pos_, *it); //最近点的初值.
			 
				//std::cout << mesh_.property(vp_project_pos_, *it) << "/" << (limit.point(v0)*bc[0]+limit.point(v1)*bc[1]+limit.point(v2)*bc[2]) << ". "; //是一致的.
				//	std::cout << mesh_.property(vp_project_pos_, *it) << "/" << mesh_.property(vp_eval_pos_, *it) << ", "; //差别不大的
					//std::cout << (Sk - mesh_.property(vp_project_pos_, *it)).norm() << ":=" << std::cout << mesh_.property(vp_project_pos_, *it) << "~=" << mesh_.property(vp_eval_pos_, *it) << ", ";

					if ((Sk - mesh_.property(vp_project_pos_, *it)).norm() < 0.00001) { // <1.0e-5
						// 2009-05-21, 之前发现Sk和投影点几乎重合了, sk也就是在面上,但是得到的eval point和最近点和Sk却距离很大.
						mesh_.property(vp_closest_pos_, *it) = mesh_.property(vp_project_pos_, *it);//
						mesh_.property(vp_closest_, *it) = true; 
					} else {
						OpenMesh::Vec3d pp(0, 0, 0), clp(0, 0, 0);
						if (evaluation_get_footpoint_byQAS(initial_control_mesh_, fh_porject_ok, Sk, bc, clp, 0)) {
							mesh_.property(vp_closest_pos_, *it) = clp;//
							mesh_.property(vp_closest_, *it) = true; // 
							//std::cout << (Sk - clp).norm() << " ";
				 
						} else {
							
							//处理方法一, 一个原始采样点Sk找不到最近点就直接退出程序.
							//return false;
							//处理方法二, 将这个最近点指定为0, 跳过这个, 继续找其他Sk的最近点.//2009-05-13
							//vp_eval_pos = OpenMesh::Vec3d(0,0,0);
							//mesh_.property(vp_eval_, *it) = true; 
							//vp_closet_pos =  OpenMesh::Vec3d(0,0,0);//
							//mesh_.property(vp_closest_, *it) = true; 
							//处理方法三, 不管了, 直接跳过. 因为下面会将projection position 付给它.
							if (mesh_.property(vp_project_, *it)) {
								mesh_.property(vp_closest_pos_, *it) = mesh_.property(vp_project_pos_, *it); 
								mesh_.property(vp_closest_, *it) = true;
								std::cout << "Info: 用project点做closest点.\n";
							} else {
								std::cout << "Error: evaluation_get_footpoint 找不到最近点.\n";
							}
						 
						} 
					}					
				} else { // 2009-05-18-2126 测试了上面要是单用parameter的结果的话好像没有这情况的..
					std::cout << "Error: 上面才求完的bc现在就不能用了, 不正常. where_get_pro " << where_get_pro << ".\n";
				}

				if (mesh_.property(vp_closest_, *it)) {
					// 2009-05-21, a breakthrough, 这样就毁灭了(Sk和project point很近但是和closest point却很远)这情况.
					// To decide whether to use the calculated closest point or the projection point
					if ((Sk - mesh_.property(vp_project_pos_, *it)).sqrnorm() < (Sk - mesh_.property(vp_closest_pos_, *it)).sqrnorm()) {
						mesh_.property(vp_closest_pos_, *it) = mesh_.property(vp_project_pos_, *it);
					} // 
				} else { // 2009-05-21, this makes sure that every *it has its closest point.
					 mesh_.property(vp_closest_pos_, *it) = mesh_.property(vp_project_pos_, *it);
					 mesh_.property(vp_closest_, *it) = true;
				}
				
		} // end of for each vertex 
	} // end of for each face.
	std::cout << "求投影时候的四种情况的比较: " << one_pro << ", " << two_pro << ", " << three_pro << ", " << four_pro << ".\n";
	//std::cout << "Info: use projection vs parameterization: " << method_1_count << ": " << mesh2_.n_vertices() - method_1_count << "." << method_3_count <<".\n";
	// 2009-05-13, fan_6475 -> 200v, 5891vs5xx, 可见大部分都是用project做初始值好.
	
	std::cout << "evaluation_get_all_footpoints(), end.\n";

	// To calculate the fitting error.
	// (1). 计算原始网格mesh_上每一个顶点的fitting error. 
	int m = Xm_vertex_.size(), n = simplified_mesh_.n_vertices(), Xm_idx= 0;
	for (TriMesh::VIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++ v_it) {
		if (mesh_.status(v_it.handle()).deleted()) {
			if (mesh_.property(vp_closest_, v_it)) {
				mesh_.property(vp_fitting_error_, v_it) 
					= (mesh_.property(vp_closest_pos_, v_it.handle()) - mesh2_.point(v_it)).norm() / bb.diagonal_lenght;
			} else {
				std::cout << "Error: 竟然还有顶点还没有求出最近点.\n";
			} 	  
			
		} else {
			// method 1. using one-ring
			int count = 0;
			double sum_fitting_error = 0;// 一邻域的fitting error
			for (TriMesh::VertexVertexIter vv_it(mesh2_, v_it); vv_it; ++ vv_it) { //这些点有某个可能也是没有被删去的点.
				// 这里fix one bug. 2009-05-21, vv_it和v_it混淆了.
				if (mesh_.property(vp_closest_, vv_it)) {
					sum_fitting_error += (mesh_.property(vp_closest_pos_, vv_it.handle()) - mesh2_.point(vv_it)).norm() / bb.diagonal_lenght;
					++ count;
				}	 
			} // count <= mesh2_.valence(v_it);
			if (count > 0)
			mesh_.property(vp_fitting_error_, v_it) = sum_fitting_error / count;
			else mesh_.property(vp_fitting_error_, v_it) = 0;
			if (mesh_.property(vp_fitting_error_, v_it) < (mesh2_.point(v_it) - refined_simplified_mesh_.point(Xm_vertex_[Xm_idx])).norm() / bb.diagonal_lenght) {
				 // 
			} else {
				mesh_.property(vp_fitting_error_, v_it) = (mesh2_.point(v_it) - refined_simplified_mesh_.point(Xm_vertex_[Xm_idx])).norm() / bb.diagonal_lenght;
			} 

			++Xm_idx;
		}
	} 
	if (Xm_idx != n) std::cout << "Error: "; 
	std::cout << "之前没有被删除的点数是" << Xm_idx << ".\n";
	
	vsplit_array_.clear();
	vsplit_array_sm_.clear();
	// 对所有的fitting error来做排序.
	std::vector<double> fitting_error_array; 
	for (TriMesh::VIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++ v_it) {
		fitting_error_array.push_back(mesh_.property(vp_fitting_error_, v_it));
	}
	std::sort(fitting_error_array.begin(), fitting_error_array.end());
	standard_error_ = fitting_error_array[fitting_error_array.size() * 0.9];
	long double sum_error = std::accumulate(fitting_error_array.begin(), fitting_error_array.end(), 0.0L);
	long double avg_error = sum_error / (fitting_error_array.size());
	std::cout << "Info: " << sum_error << " / " << fitting_error_array.size() << " = " << avg_error << ".\n";
	std::cout << "Info: fitting error " << fitting_error_array[0] << ", " << fitting_error_array[1] << ", ... , " 
		<< fitting_error_array[fitting_error_array.size() * 0.5] << ", " << fitting_error_array[fitting_error_array.size()-3] << ".\n";
	if (fitting_error_array[fitting_error_array.size()-3] < 0.001) {
		fitting_terminate_ = true;
	}
	fitting_error_array.clear();

	std::set<TriMesh::VHandle, FErr_v> v_fe_arr;//对顶点按fittting error来排序, 就可以知道误差大的顶点.
	for (TriMesh::VIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++ v_it) {
		v_fe_arr.insert(v_it.handle());
	}
	 
	std::vector<TriMesh::VHandle> temp_array(v_fe_arr.begin(), v_fe_arr.end()); 
	int n_want_to_split = 20;
	std::vector<TriMesh::VHandle> temp_array2 = std::vector<TriMesh::VertexHandle>(temp_array.end() - n_want_to_split, temp_array.end());//get the biggest last 10.
	std::cout << " vh " << temp_array2[temp_array2.size()-1] << ", size " <<temp_array2.size() << "\n";
	//for (int i = 0; i < temp_array2.size(); ++i ) { 
	//	std::cout << (mesh2_.point(temp_array2[i]) - mesh_.property(vp_project_pos_, temp_array2[i])).norm() << ", ";
	//	std::cout << (mesh_.property(vp_fitting_error_, temp_array2[i])) << ". ";
	//	//std::cout << mesh_.property(vp_project_pos_, temp_array2[i]) << "/" << mesh_.property(vp_eval_pos_,temp_array2[i]) << ", "; //是区别不大,
	//}std::cout << ".\n";
	// 从mesh_的最大误差的顶点数组中找到他们对应的vertex hierarchy中的leaf node,
	// 进而找他们对应simplified mesh上的顶点, 用来做split.
	vsplit_array_.clear();
	for (std::vector<TriMesh::VHandle>::const_iterator it(temp_array2.begin()),  it_end(temp_array2.end()); it != it_end; ++it) {
		TriMesh::VHandle vh = *it;
		int node_handle = mesh_.property(vp_leaf_node_handle_, vh);
		DGP::VHierarchyNode node = vhierarchy_.node(node_handle); // 
		while (node.is_active() == false) {  // 
			node = vhierarchy_.node(node.parent_node_handle());
		}  
		vsplit_array_.push_back(node.vertex_handle());
	} 
	std::cout << "Info: vsplit_array_.size: " << vsplit_array_.size() << ", going to be split.\n";

	//(2). 计算简化网格上每一个顶点的fitting error
	fitting_error_array.clear();
	simplified_mesh_.add_property(vp_fitting_error_sm_);
	for (TriMesh::VIter v_it(simplified_mesh_.vertices_begin()), v_end(simplified_mesh_.vertices_end()); v_it != v_end; ++ v_it) {
		int node_handle = simplified_mesh_.property(vp_node_handle_sm_, v_it.handle());
		if (vhierarchy_.node(node_handle).is_active() == false) std::cout << "Error: 基网格上顶点对应的node应该是active才对.\n";

		int n_leaves = 0;
		double sum_fitting_error = 0; // 这棵树的所有叶子节点的fitting error之和.
		std::vector<int> leaf_node_handle_array;
		traverse_binary_tree_count_leaves(node_handle, n_leaves);
		sum_fitting_error = traverse_binary_tree_add_leaves_fitting_error(node_handle);
		//traverse_binary_tree_collect_leaves(node_handle, leaf_node_handle_array);
		//std::cout << v_it.handle().idx() << ": " << n_leaves << ", " << leaf_node_handle_array.size() << ", " << sum_fitting_error << ", " << sum_fitting_error / n_leaves << ".\n";

		simplified_mesh_.property(vp_fitting_error_sm_, v_it) = sum_fitting_error / n_leaves; 
		fitting_error_array.push_back(simplified_mesh_.property(vp_fitting_error_sm_, v_it));
	} 
	std::sort(fitting_error_array.begin(), fitting_error_array.end());
	standard_error_sm_ = fitting_error_array[fitting_error_array.size() * 0.5];
	std::cout << "Info: fitting error sm " << fitting_error_array[0] << ", " << fitting_error_array[1] << ", ... , " 
		<< fitting_error_array[fitting_error_array.size() * 0.9] << ", " << fitting_error_array[fitting_error_array.size()-1] << ".\n";
	fitting_error_array.clear();
	//std::cout << "a";
	std::set<TriMesh::VHandle, FErr_v_sm> v_fe_arr_sm;
	for (TriMesh::VIter v_it(simplified_mesh_.vertices_begin()), v_end(simplified_mesh_.vertices_end()); v_it != v_end; ++ v_it) {
		v_fe_arr_sm.insert(v_it.handle());
	}//std::cout << "b";

	std::vector<TriMesh::VHandle> temp_array_sm(v_fe_arr_sm.begin(), v_fe_arr_sm.end()); // std::cout << "c";
	vsplit_array_sm_.clear();
	int k_idx = 0; 
	for (std::vector<TriMesh::VHandle>::reverse_iterator it = temp_array_sm.rbegin(); k_idx < n_want_to_split && it != temp_array_sm.rend();
		 ++it) 
	{
		 vsplit_array_sm_.push_back(*it); 
		 ++k_idx;
	}
	//vsplit_array_sm_ = std::vector<TriMesh::VHandle>(temp_array_sm.end()-n_want_to_split, temp_array_sm.end());//取误差最大的10.
	//std::cout << "d";
	return true;
}
// ------------------------------------------------------------------------
void DecimationModel::traverse_binary_tree_count_leaves(int _node_handle, int &_n_leaves) {
	DGP::VHierarchyNode node = vhierarchy_.node(_node_handle);
	if (node.lchild_node_handle() == -1) { //leaf
		++_n_leaves;
	} else {
		traverse_binary_tree_count_leaves(node.lchild_node_handle(), _n_leaves);
		traverse_binary_tree_count_leaves(node.rchild_node_handle(), _n_leaves);
	}
}
void DecimationModel::traverse_binary_tree_collect_leaves(int _node_handle, std::vector<int> &leaf_node_handle_array) {
	DGP::VHierarchyNode node = vhierarchy_.node(_node_handle);
	if (node.lchild_node_handle() == -1) { //leaf
		leaf_node_handle_array.push_back(node.lchild_node_handle());
	} else {
		traverse_binary_tree_collect_leaves(node.lchild_node_handle(), leaf_node_handle_array);
		traverse_binary_tree_collect_leaves(node.rchild_node_handle(), leaf_node_handle_array); 
	}
}

double DecimationModel::traverse_binary_tree_add_leaves_fitting_error(int _node_handle) {
	DGP::VHierarchyNode node = vhierarchy_.node(_node_handle);
	if (node.lchild_node_handle() == -1) {
		TriMesh::VHandle vh = voldhandles_[node.point_handle()]; // vertex of mesh_
		if (mesh_.status(vh).deleted()) {
			return mesh_.property(vp_fitting_error_, vh); 
		} else {
			return 0;
		}		
	} else {
		return traverse_binary_tree_add_leaves_fitting_error(node.lchild_node_handle()) 
			+ traverse_binary_tree_add_leaves_fitting_error(node.rchild_node_handle());
	} 
}
// 计算原始网格Sk点都其对应的最近点的距离作为fitting error
void DecimationModel::calc_fitting_error() {
	
}
// -----------------------------------------------------------------------
void DecimationModel::evaluation_process() { //这一步只能在刚完成了求解limit surface之后进行.
	std::cout << "\n";
	std::cout << "========开始计算原始点在极限曲面上的坐标\n";//这种提示用40个长度的=,前面第9个位置开始文字,即最先有8个=.
	
	if (evaluation_get_all_footpoints() == false) {
		std::cout << "Error: 求解原始顶点到极限曲面上的最近点坐标时候出错.\n";
	} else {
		std::cout << "完成求原始顶点在极限曲面上的最近点坐标.\n"; 
		calc_fitting_error();
	}
	std::cout << "========结束计算原始点在极限曲面上的坐标\n";//这种提示用40个长度的=,前面第9个位置开始文字,即最先有8个=. 
}
void DecimationModel::view_fitting_error_mesh2_process() {

	std::vector<double> fitting_error_array; 
	for (TriMesh::VIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++ v_it) {
		fitting_error_array.push_back(mesh_.property(vp_fitting_error_, v_it));
	}
	std::sort(fitting_error_array.begin(), fitting_error_array.end());

	// discard upper and lower 5% = 5/100 = 1/20, 丢掉了一头一尾
	unsigned int n = fitting_error_array.size()-1;// 例如 n = 100
	unsigned int i = n / 20;// i = 5
	double min_fitting_error = fitting_error_array[i];// 5
	double max_fitting_error = 0.01;//貌似0.02也是挺合适的2009-05-21 // 94 // fitting_error_array[n-1-i]


	// define uniform color intervalls [v0,v1,v2,v3,v4]
	double v0, v1, v2, v3, v4;
	v0 = min_fitting_error + 0.0/4.0 * (max_fitting_error - min_fitting_error);
	v1 = min_fitting_error + 1.0/4.0 * (max_fitting_error - min_fitting_error);
	v2 = min_fitting_error + 2.0/4.0 * (max_fitting_error - min_fitting_error);
	v3 = min_fitting_error + 3.0/4.0 * (max_fitting_error - min_fitting_error);
	v4 = min_fitting_error + 4.0/4.0 * (max_fitting_error - min_fitting_error);


	// map fitting error to colors 误差由小变大,颜色由蓝变红
	// blue~ v0	(0,u,255)	v1	(0, 255, 255-u)	v2	(u, 255, 0)	v3	(255, 255-u, 0)	v4	~red
	double fitting_error;
	TriMesh::Color col;
	mesh2_.add_property(vp_color_fitting_error_mesh2_);
	for (TriMesh::VIter v_it=mesh_.vertices_begin(), v_end(mesh_.vertices_end()); v_it!=v_end; ++v_it)
	{ // 开始循环每一个顶点
		fitting_error = mesh_.property(vp_fitting_error_, v_it);
		col = TriMesh::Color(255,255,255);

		unsigned char u;

		if (fitting_error < v0) 
		{
			col = TriMesh::Color(0, 0, 255);// 蓝色
		}
		else if (fitting_error > v4) 
		{
			col = TriMesh::Color(255, 0, 0);// 红色
		}

		else if (fitting_error <= v2) 
		{
			if (fitting_error <= v1) // [v0, v1]
			{
				u = (unsigned char) (255.0 * (fitting_error - v0) / (v1 - v0));
				col = TriMesh::Color(0, u, 255);
			}      
			else // ]v1, v2]
			{
				u = (unsigned char) (255.0 * (fitting_error - v1) / (v2 - v1));
				col = TriMesh::Color(0, 255, 255-u);
			}
		}
		else 
		{
			if (fitting_error <= v3) // ]v2, v3]
			{
				u = (unsigned char) (255.0 * (fitting_error - v2) / (v3 - v2));
				col = TriMesh::Color(u, 255, 0);
			}
			else // ]v3, v4]
			{
				u = (unsigned char) (255.0 * (fitting_error - v3) / (v4 - v3));
				col = TriMesh::Color(255, 255-u, 0);
			}
		}

		mesh2_.property(vp_color_fitting_error_mesh2_, v_it) = col;
		mesh2_.set_color(v_it, col);
	} //结束循环每一个顶点
	// ------------
	is_view_fitting_error_of_mesh2_ = !is_view_fitting_error_of_mesh2_;
	is_view_fitting_error_of_sm_ = false;
}
void DecimationModel::view_fitting_error_sm_process() {

	std::vector<double> fitting_error_array; 
	for (TriMesh::VIter v_it(simplified_mesh_.vertices_begin()), v_end(simplified_mesh_.vertices_end()); v_it != v_end; ++ v_it) {
		fitting_error_array.push_back(simplified_mesh_.property(vp_fitting_error_sm_, v_it));
	}
	std::sort(fitting_error_array.begin(), fitting_error_array.end());

	// discard upper and lower 5% = 5/100 = 1/20, 丢掉了一头一尾
	unsigned int n = fitting_error_array.size()-1;// 例如 n = 100
	unsigned int i = n / 20;// i = 5
	double min_fitting_error = fitting_error_array[i];// 5
	double max_fitting_error =  0.01;// 94 // fitting_error_array[n-1-i];

	// define uniform color intervalls [v0,v1,v2,v3,v4]
	double v0, v1, v2, v3, v4;
	v0 = min_fitting_error + 0.0/4.0 * (max_fitting_error - min_fitting_error);
	v1 = min_fitting_error + 1.0/4.0 * (max_fitting_error - min_fitting_error);
	v2 = min_fitting_error + 2.0/4.0 * (max_fitting_error - min_fitting_error);
	v3 = min_fitting_error + 3.0/4.0 * (max_fitting_error - min_fitting_error);
	v4 = min_fitting_error + 4.0/4.0 * (max_fitting_error - min_fitting_error);


	// map fitting error to colors 误差由小变大,颜色由蓝变红
	// blue~ v0	(0,u,255)	v1	(0, 255, 255-u)	v2	(u, 255, 0)	v3	(255, 255-u, 0)	v4	~red
	double fitting_error;
	TriMesh::Color col;
	simplified_mesh_.add_property(vp_color_fitting_error_sm_);
	for (TriMesh::VIter v_it=simplified_mesh_.vertices_begin(), v_end(simplified_mesh_.vertices_end()); v_it!=v_end; ++v_it)
	{ // 开始循环每一个顶点
		fitting_error = simplified_mesh_.property(vp_fitting_error_sm_, v_it);
		col = TriMesh::Color(255,255,255);

		unsigned char u;

		if (fitting_error < v0) 
		{
			col = TriMesh::Color(0, 0, 255);// 蓝色
		}
		else if (fitting_error > v4) 
		{
			col = TriMesh::Color(255, 0, 0);// 红色
		}

		else if (fitting_error <= v2) 
		{
			if (fitting_error <= v1) // [v0, v1]
			{
				u = (unsigned char) (255.0 * (fitting_error - v0) / (v1 - v0));
				col = TriMesh::Color(0, u, 255);
			}      
			else // ]v1, v2]
			{
				u = (unsigned char) (255.0 * (fitting_error - v1) / (v2 - v1));
				col = TriMesh::Color(0, 255, 255-u);
			}
		}
		else 
		{
			if (fitting_error <= v3) // ]v2, v3]
			{
				u = (unsigned char) (255.0 * (fitting_error - v2) / (v3 - v2));
				col = TriMesh::Color(u, 255, 0);
			}
			else // ]v3, v4]
			{
				u = (unsigned char) (255.0 * (fitting_error - v3) / (v4 - v3));
				col = TriMesh::Color(255, 255-u, 0);
			}
		}

		simplified_mesh_.property(vp_color_fitting_error_sm_, v_it) = col;
		simplified_mesh_.set_color(v_it, col);
	} //结束循环每一个顶点
	// ------------
	is_view_fitting_error_of_sm_ = !is_view_fitting_error_of_sm_;
	is_view_fitting_error_of_mesh2_ = false;
}
