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
	// �����ֶ����򲢽������Ԫ������(��������).JTJh=-JTf
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
	if (itype == Steepest_descent) { //����������ʹʹ��ȫ����(1or-1)��������������! ����50��ֻ��ǰ����һ���.
		delta_v = B[0];
		delta_w = B[1]; 
	} else {
		delta_v = (B[0]*A[1][1] - B[1]*A[0][1]) / (A[0][0]*A[1][1] - A[0][1]*A[1][0]);	//std::cout << "delta_v: " << delta_v << "; ";////
		delta_w = (B[0] - A[0][0] * delta_v) / A[0][1];									//std::cout << "delta_w: " << delta_w << ".\n";// //
	}		 
}
void DecimationModel::update_uvw_when_v_new_is_less_than_zero(const double v, const double delta_v, const double w, const double delta_w, double &v_new, double &w_new) {
	//�����仯v.
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
	if (v_new + w_new > 1) { //��������仯���еĻ�,
		if (delta_w > 0) w_new = w;	//w�Ȳ���,
		if (v_new + w_new > 1) {//���ǲ��еĻ�, ��v�Ĳ�����Сһ��.
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

// Newton Iteration�������.
// ��Sk��cN6������Ķ�����ɵ������ϵ������_clp, �����ĳ�ʼֵΪ(1 - _v - _w, _v, _w), _clp�ĳ�ʼֵΪSk��������Ƭ�ϵľ�ȷ���ꡣ
bool DecimationModel::cloop_stam98_get_closest_points (const std::vector<OpenMesh::Vec3d> & cN6, 
													   const double _v, const double _w, OpenMesh::Vec3d Sk, 
													   OpenMesh::Vec3d & _clp, 
													   const TriMesh& _mesh, TriMesh::HHandle h12, int cross_patch_depth) 
{
	//����irregular��ģ��. regularҲ��irregular�ķ�ʽ�����.
	// (1), search descent direction h(delta_v, delta_w)
	//std::cout << "Sk and initial clp " << Sk << ", " << _clp << ". Initial dis: " << (Sk - _clp).norm() << ".\n";
	std::vector<OpenMesh::Vec3d> c = cN6;
	enum return_type { Null, Min_dis, Max_times, Out_vw };
	return_type rt = Null;

	const int Max_Iter_Times = 2000;//���ĵ�������.��Щ�������ǿ��Ե���.
	const double Max_Step_length = 0.5; //��󲽳�

	OpenMesh::Vec3d clp_arr[2];//������, ��������, �Ĵκ���ȡ�Ǵ�.
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
			if (fabs(delta_v) > 1 || fabs(delta_w) > 1) {  //��������ķ�������ǶԵ�, ���ǳ���̫����, ����0~1
				//double len = sqrt(delta_v * delta_v + delta_w * delta_w);
				delta_v = 0.9;// /= len;
				delta_w = 0.9; // /= len;
			}
			if (delta_v < 0.0001 && delta_w < 0.00001) {
				// ��������Կ����ǵ�����ֹ��, ��ô�̵Ĳ���.
				if (cloop_evaluation.evaluate_irregular_patch(c, v, w, clp, OpenMesh::Vec3d(0,0,0), OpenMesh::Vec3d(0,0,0)) == false) {
					std::cout << "Error: stam98 test regular patch.\n";
				}   
				rt = Min_dis; break;				
			} 
			else {
				// 2. ��ʼ����.
				double step_len = -0.01; //-0.008;//1;
				if (i == 0) step_len *= 1;
				else step_len *= -1;

				if (delta_v < 0.0001 || delta_w < 0.00001) {
				} else {
					// ������ʹ��brent minimization��������Ĳ���, ���ӽ������, ��������������Ҳ����.
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
					//	status = gsl_min_fminimizer_set (s, &F, m, a, b); //��������ĺ���F, ��ֵm, ��Χab
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

					//			status = gsl_min_test_interval (a, b, 0.001, 0.0); // ����if b-a < 0.001,�ǵĻ���status == GSL_SUCCESS

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
				// 3. Update delta_v and delta_w.// �õ�����step side֮��			
				delta_v *= step_len;
				delta_w *= step_len;

				if ((v + delta_v) < 0) {
					std::cout << "v < 0. "; //, " << v << ", + " << delta_v << ", " << _v << " ////����û������, ��w<0һ��û������, ����>1ʱ��������.
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
							break; //���Ƿ��ش���, �����ҵ���ǰΪֹ, ��ǰ�Ľ��Ҳ����Ч��. //return false;
						} 
					}
					rt = Out_vw; break; 				
				} else if ((w + delta_w) < 0) {
					std::cout << "w < 0. "; //\n// //����û������, ��v<0һ��û������, ��>1Ҳû������.
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
							break; //���Ƿ��ش���, �����ҵ���ǰΪֹ, ��ǰ�Ľ��Ҳ����Ч��. //return false;
						} 
					}				
					rt = Out_vw; break;
				} else if ((v + delta_v + w + delta_w) > 1) {
					std::cout << "v + w > 1. "; //\n//����û������, ��w<0Ҳû������, ����v<0һ���������
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
							break; //���Ƿ��ش���, �����ҵ���ǰΪֹ, ��ǰ�Ľ��Ҳ����Ч��. //return false;
						}/**/
					}
					rt = Out_vw; break;
				} else if ((v + delta_v) > 1) { //����û�г��ֹ������.
					std::cout << "v > 1, " << v << ", + " << delta_v << ".\n"; rt = Out_vw; break;
				} else if ((w + delta_w) > 1) { //����û�г��ֹ������.
					std::cout << "w > 1.\n"; rt = Out_vw; break;
				} else { 
					//����ͬһ������������
					//std::cout << "v: " << v << ", w: " << w << ". ";
					v += delta_v; 
					w += delta_w;
					//std::cout << "v': " << v << ", w': " << w << ".\n";

					double distance_before = (clp-Sk).norm();
					OpenMesh::Vec3d clp_before = clp; //����ǰһ����������� 
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
			std::cout << "Info: Iteration ends Max_Iter_Times, OK.\n";//���ڳ���������������, 
			rt = Max_times;
		} 
		clp_arr[i] = clp; //������һ��closest position.
		if ((Sk - clp).norm() < (Sk - _clp).norm()) { //��Ͳ����ҵڶ�����.
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
	if (_is_regular) { // feature regular, feature_type 7��(1~7), _c��12������, h��15��������.
	} else { // feature irregular, feature type 17��(1~17), _c��N+6������, h��N+9��������.

	}
	enum return_type { Null, Min_dis, Max_times, Out_vw, Little_delta };
	return_type rt = Null;

	const int Max_Iter_Times = 4000;//���ĵ�������.��Щ�������ǿ��Ե���.
	const double Max_Step_length = 0.5; //��󲽳�
	double distance_result = 100;
	double distance[2] = {100};
	OpenMesh::Vec3d clp_arr[2] = { OpenMesh::Vec3d(10000, 10000, 10000) }; //��ʵ�ֵ�ʱ���Ǹ��������Զ�(����, ��С),�����ҷ�����������, ȡЧ���õ��Ǵ�.

	for (int i = 0; i < 2; ++i) {
		double v = _v, w = _w;
		OpenMesh::Vec3d clp = _clp;

		int iter_times = 0; 
		for (iter_times = 0; iter_times < Max_Iter_Times; ++iter_times) {		
			// 0. ����ſ˱Ⱦ����Ԫ��
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

			// 1. �ⷽ��JTJh = -JTf�õ������ķ���.
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
			if (fabs(delta_v) > 1 || fabs(delta_w) > 1) {  //��������ķ�������ǶԵ�, ���ǳ���̫����, ����0~1
				double len = sqrt(delta_v * delta_v + delta_w * delta_w);
				delta_v /= len;
				delta_w /= len;
			}

			// 2. ��ʼ����.
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
			} else if ((v + delta_v) > 1) { std::cout << "v > 1.\n"; rt = Out_vw; break; //��������������û�г��ֵ�.
			} else if ((w + delta_w) > 1) { std::cout << "w > 1.\n"; rt = Out_vw; break;
			} else { //����ͬһ������������
				//std::cout << "SFG. "; //Same feature region//
				//std::cout << "v: " << v << ", w: " << w << ". ";////
				v += delta_v; 
				w += delta_w; 
				//std::cout << "v': " << v << ", w': " << w << ".\n";////

				double distance_before = (Sk - clp).norm();
				OpenMesh::Vec3d clp_before = clp; //����ǰһ�����������

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
		// �˳�for iter_times�����:
		// (1) distance after > distance before, ���������˳�.
		// (2) �ﵽ����������, ���������˳�.
		if (iter_times >= Max_Iter_Times) {////
			std::cout << "Info: Iteration ends Max_Iter_Times, OK." << ".\n";// ���ڳ���������������
			rt = Max_times;//
		} 
		// (3) v, w out of [0, 1] range. Ҳ�������˳�.
		// (4) delta_v, delta_w ̫С��Ҳ�����˳�.

		clp_arr[i] = clp;
		// ��ǰ�˳�. �����ǲ�����һ����������˾���ζ���������һ�������ҪС��?
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
{   //�������û�зŵ�Hoppe94LoopEval.h������Ϊ�����feature_type����Ϊ�Ƿ�regular, ����ֵ̫��.
	std::vector<TriMesh::HHandle> h_arr(3, TriMesh::HHandle(-1));
	std::vector<TriMesh::VHandle> v_arr(3, TriMesh::VHandle(-1));
	h_arr[0] = initial_control_mesh_.halfedge_handle(fh); v_arr[0] = initial_control_mesh_.to_vertex_handle(h_arr[0]);
	h_arr[1] = initial_control_mesh_.next_halfedge_handle(h_arr[0]); v_arr[1] = initial_control_mesh_.to_vertex_handle(h_arr[1]);
	h_arr[2] = initial_control_mesh_.next_halfedge_handle(h_arr[1]); v_arr[2] = initial_control_mesh_.to_vertex_handle(h_arr[2]);

	if (initial_control_mesh_.property(vp_type_icm_, v_arr[0]) == DGP::SMOOTH_VFT && initial_control_mesh_.property(vp_type_icm_, v_arr[1]) == DGP::SMOOTH_VFT 
		&& initial_control_mesh_.property(vp_type_icm_, v_arr[2]) == DGP::SMOOTH_VFT) { 
			TriMesh::HalfedgeHandle h12(-1); //û�м�������ʱ��,��regular patchҲͳһΪirregular patch����.
			for (int i = 0; i < 3; ++i) {
				if (initial_control_mesh_.valence(v_arr[i]) != 6) {
					h12 = initial_control_mesh_.next_halfedge_handle(h_arr[i]); break;//�зǹ��������
				}
			}
			if (h12.is_valid() == false) h12 = h_arr[1]; //û�зǹ��������.

			_is_feature = false; _is_regular = false; feature_type = 0; h = h12;
	}  
	else { // feature situation.
		if (initial_control_mesh_.valence(v_arr[0]) == 6 && initial_control_mesh_.valence(v_arr[1]) == 6 && initial_control_mesh_.valence(v_arr[2]) == 6) {  
			// feature regular
			enum feature_regular_type { tNull = 0, t1v0e = 1, t2v0e, t2v1e, t3v0e, t3v1e, t3v2e, t3v3e };//1~7������
			feature_regular_type frt = tNull;
			TriMesh::HalfedgeHandle h84(-1);

			int n_feature_vertics = 0;
			if (initial_control_mesh_.property(vp_type_icm_, v_arr[0]) != DGP::SMOOTH_VFT) ++ n_feature_vertics;
			if (initial_control_mesh_.property(vp_type_icm_, v_arr[1]) != DGP::SMOOTH_VFT) ++ n_feature_vertics;
			if (initial_control_mesh_.property(vp_type_icm_, v_arr[2]) != DGP::SMOOTH_VFT) ++ n_feature_vertics;
			if (n_feature_vertics == 0) std::cout << "Error: no feature stituation Ӧ����������Ѵ������.\n";
			else if (n_feature_vertics == 1) { // 1�����
				if (initial_control_mesh_.property(vp_type_icm_, v_arr[0]) != DGP::SMOOTH_VFT) h84 = h_arr[0];
				if (initial_control_mesh_.property(vp_type_icm_, v_arr[1]) != DGP::SMOOTH_VFT) h84 = h_arr[1];
				if (initial_control_mesh_.property(vp_type_icm_, v_arr[2]) != DGP::SMOOTH_VFT) h84 = h_arr[2];

				frt = t1v0e;
			} else if (n_feature_vertics == 2) { // 2�����
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

			} else if (n_feature_vertics == 3) { // 4�����
				if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[0])).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[1])).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[2])).feature()) {
					h84 = h_arr[0]; frt = t3v3e;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[0])).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[1])).feature() ) { //����������.
					h84 = h_arr[0]; frt = t3v2e;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[1])).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[2])).feature() ) {
					h84 = h_arr[1]; frt = t3v2e;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[2])).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[0])).feature() ) {
					h84 = h_arr[2]; frt = t3v2e;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[0])).feature()) { //һ��������.
					h84 = h_arr[2]; frt = t3v1e;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[1])).feature()) {
					h84 = h_arr[0]; frt = t3v1e;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h_arr[2])).feature()) {
					h84 = h_arr[1]; frt = t3v1e;
				} else { //û��������.
					h84 = h_arr[0]; frt = t3v0e;
				} 					
			}  
			if (h84.is_valid() == false) std::cout << "Error: �������еĴ���, ����Ӧ���Ѿ��о������п��������.\n";

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
				std::cout << "Error: feature irregular. ����3���㶼��regular���򲻻���������: .\n";
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
					ft = t3v2e_up; //���������ne....160.off��û�д�.
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12)).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12next)).feature()) {
					ft = t3v2e_left;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12next)).feature() && initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12prev)).feature()) {
					ft = t3v2e_right;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12next)).feature()) {
					ft = t3v1e_bottom;
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12)).feature()) {
					if (initial_control_mesh_.property(vp_type_icm_, vhright) != DGP::SMOOTH_VFT) ft = t3v1e_left;
					else ft = t2v1e_left; // nefertiti.160.off����Ҳû�д�.
				} else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(h12prev)).feature()) {
					if (initial_control_mesh_.property(vp_type_icm_, vhleft) != DGP::SMOOTH_VFT) ft = t3v1e_right;
					else ft = t2v1e_right;// nefertiti.160.off����Ҳû�д�.
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
// �󶥵�Sk��ƽ��fh�ϵ�ͶӰ��p��Ӧ���������ϵľ�ȷ����_pp�Լ����������_clp. ע��p�ǿ϶��ǻ�ͶӰ��fh�ϵĵ�, ֮ǰ�����ﷸ����.
bool DecimationModel::evaluation_get_footpoint(const TriMesh& _mesh, TriMesh::FaceHandle _fh_of_initial_control_mesh_, //�ڴ�����Ĵ�����
											   OpenMesh::Vec3d Sk, OpenMesh::Vec3d p, // input: ������͵��������ĳ�ֵ
											   OpenMesh::Vec3d _bc, // input: p �����_fh_of_initial_control_mesh_�Ĳ�������
											   OpenMesh::Vec3d &_pp, OpenMesh::Vec3d &_clp, const int cross_patch_depth) // output:ǰ����.
{ 
	TriMesh::FHandle fh = _fh_of_initial_control_mesh_;

	
	// ���ж��������
	bool is_feature = false;// false,û�м�������, true, �м�������.
	bool is_regular = true; // is regular or irregular.
	int feature_type = 0;   // is_feature == trueʱ��, regular face has 7 cases, irregular 17 cases.
	TriMesh::HalfedgeHandle h; // regular, h84. irregular h12.
	evaluation_get_face_type(fh, is_feature, is_regular, feature_type, h);
	// is_feature, is_regular, feature_type, h
	// false       false       0             h12  //��û������������£�ͳһ�����ǹ�����������.
	// true        true        1~7           h84   
	// true        false       1~17          h12   
	//if (h.is_valid() == false) std::cout << "Error: h is non-valid.\n";//�ⲻ�ó��ֵ�, ����Ӧ�úܼ򵥰�.
	//std::cout << is_feature << ", " << is_regular << ", " << feature_type << ". ";
	//if (is_feature == false || is_regular == false ) continue; 

	if (is_feature == false) { 
		TriMesh::HalfedgeHandle h12 = h; //û�м�������ʱ��,��regular patchҲͳһΪirregular patch����.
		int N = initial_control_mesh_.valence(initial_control_mesh_.from_vertex_handle(h12));

		std::vector<OpenMesh::Vec3d> cN6(N+6, OpenMesh::Vec3d(0, 0, 0));
		DGP::cloop_stam98_irregular_get_cN6(initial_control_mesh_, h12, cN6);

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(cN6[0], cN6[1], cN6[N], p); // ��Ӧ����������ʱ���_bcһ��!
		if (DGP::is_valid_barycentric_coordinate(bc, true) == false) {
			std::cout << "Error: evaluation_get_footpoint(...): non-feature bc: " << bc << ".\n";
			std::cout << "contu.." << cN6[0] << ", " << cN6[1] << ", " << cN6[N] << ", " << p << ".\n";
			return false;
		}
		double v = bc[1], w = bc[2]; 
		// ���ڼ��������ϵ�evaluation position
		OpenMesh::Vec3d pp(0, 0, 0);
		if (cloop_evaluation.evaluate_irregular_patch(cN6, v, w, pp, OpenMesh::Vec3d(0,0,0), OpenMesh::Vec3d(0,0,0)) == false) {
			std::cout << "Error: �������." << v << ", " << w << " \n";
		}
		_pp = pp;

		// ���ڼ��������ϵ�cloesest point position
		OpenMesh::Vec3d clp = pp; //�����ĳ�ʼ����Ϊevaluation����.
		if (cloop_stam98_get_closest_points(cN6, v, w, Sk, clp, initial_control_mesh_, h12, cross_patch_depth)) {
			_clp = clp;
			return true;
		} else { //û����ȷ��������.
			std::cout << "Error: evaluation_get_footpoint(...): non-feature region closest point.\n";
			return false;
		}

	} else { // feature situation.// is_feature == true;

		if (is_regular == true) {
			TriMesh::HHandle h84 = h;
			std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));// ����ע��: ��12�����Ƶ��������stam98 "evaluation of loop subdivision surfaces". Fig1Ϊ׼.
			std::vector<bool> fe15(15, false);//12������֮���15�����Ƿ���������
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
			// ���������������������
			if (fabs(v) < 1e-5 || fabs(w) < 1e-5) { //������crease vertices������������������
				//������ҾͲ�����, ����ؽ�evaluation������Ϊ���������. 
				//ԭ������ʱw=0��v=-1.#IND, ֮��������ľ����������-1.#IND��1.QNAN��.
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
			if (N == 6) std::cout << "Error: ����Ӧ����irregular.\n";//��ʵ�ǲ����д������.
			int K = N + 6;
			std::vector<OpenMesh::Vec3d> c(K, OpenMesh::Vec3d(0, 0, 0));// ����ע��: ��K�����Ƶ��������stam98 "evaluation of loop subdivision surfaces". Fig2Ϊ׼.
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
			if (fabs(v) < 1e-6 || fabs(w) < 1e-6) { //������crease vertices������������������
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
	} // end of if-else �������

	return false;
}
bool DecimationModel::evaluation_get_all_footpoints() {
	std::cout << "evaluation_get_all_footpoints(), begin.\n";
	// -----------������QAS�����ȱƽ��������, Ȼ����QAS�������������
	// ��ÿһ����ָ��һ����ȷ�İ��, ������������ȷһ�������������������˳��.
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
		else if (n_count_feature_edges == 1 || n_count_feature_edges == 3) { // ��¼�����������
			if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(hx)).feature()) h01 = hx;
			else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(hy)).feature()) h01 = hy;
			else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(hz)).feature()) h01 = hz; 
		} else if (n_count_feature_edges == 2) { //ѡ�������������е�һ��.
			if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(hx)).feature() == false) h01 = hz;
			else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(hy)).feature() == false) h01 = hx;
			else if (initial_control_mesh_.status(initial_control_mesh_.edge_handle(hz)).feature() == false) h01 = hy;  
		} if (h01.is_valid() == false) std::cerr << "Error: h01 is supposed to be valid.\n";
		//��¼��3�����ö��˳��, ����Ҳ��¼����Ҫ�Ǵ���������ʱ����Щ�߾���������..
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
	// initial_control_mesh_����һ��hoppe94 loopϸ��֮��ṹ��refined_simplified_mesh_��һ���ġ�
	//��Ϊinitial_control_mesh_�������֮�󾭹���һ��ϸ��, ��������refined_simplified_mesh_��һ��,
	//������һ��, ����initial_control_mesh_���ж�ĳ���������, �ٽ���refined_simplified_mesh_.property(fvset_rsm_, fh)
	//��ô���İ�����ԭʼ����, ���������������꿪ʼ�����ڼ��������ϵ�����.
	for (TriMesh::FIter f_it(initial_control_mesh_.faces_begin()), f_end(initial_control_mesh_.faces_end()); f_it != f_end; ++f_it) {
		if (initial_control_mesh_.is_boundary(f_it.handle(), true)) {
			std::cout << "Error: Sorry, ���ڻ����ܴ���߽���.\n";
			continue; //��ʱ���ܴ���߽���
		}
		//std::cout << "in face " << f_it.handle().idx() << ".\n";

		// ��������ÿһ�����������Щԭʼ����(��deleted����), ��Щ���footpoint��ʲô.
		for (std::vector<TriMesh::VHandle>::const_iterator it = refined_simplified_mesh_.property(fvset_rsm_, f_it).begin(),
			it_end = refined_simplified_mesh_.property(fvset_rsm_, f_it).end(); it != it_end; ++it) { 
				// ������ʵ����crease�����������ֱ�Ӹ�����pp��clp.

				OpenMesh::Vec3d Sk = mesh2_.point(*it);// *it��mesh_�ϵĶ���.��Ϊ������.
				OpenMesh::Vec3d vp_eval_pos, vp_closest_pos;

				// ����һ.
				//std::cout << "--------Use projection as initial values.\n";
				// (1). ���波����ͶӰ�ķ�����ԭʼ����Sk��initial_control_mesh_��һ���������ϵ�ͶӰ����, ��Ϊ��ʼ����������
				// �����ʼ����������������Sk�ڼ��������ϵľ�ȷ����, �Լ����������.
				OpenMesh::Vec3d p(0, 0, 0);//��ʼ����������
				OpenMesh::Vec3d bc_proj(0,0,0); //
				bool is_project_ok = false;
				TriMesh::FHandle fh_porject_ok = f_it.handle();
				// ����*it�Ƿ�����initial_control_mesh_��f_it����, �����ѯf_it������.
				TriMesh::HHandle he_of_face = initial_control_mesh_.property(fp_hh_icm_, f_it);
				TriMesh::VHandle v0 = limit.from_vertex_handle(he_of_face);
				TriMesh::VHandle v1 = limit.to_vertex_handle(he_of_face);
				TriMesh::VHandle v2 = limit.to_vertex_handle(limit.next_halfedge_handle(he_of_face));
				OpenMesh::Vec3d n_of_f_it = limit.calc_face_normal(f_it);//��ķ���
				double d = dot(Sk - limit.point(v0), n_of_f_it);//��Sk����ľ���.
				OpenMesh::Vec3d point_on_plane = Sk - (d /*>= 0 ? d: -1.0*d*/) * n_of_f_it;

				OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(limit.point(v0), limit.point(v1), limit.point(v2), point_on_plane);

			
				if (DGP::is_valid_barycentric_coordinate(bc)) {
					is_project_ok = true;
					fh_porject_ok = f_it.handle();
					p = point_on_plane;
					bc_proj = bc;
					
					one_pro ++; where_get_pro = 1;
				} else { 
					// ��f_it���3��������.
					for (TriMesh::FHIter fh_it(limit, f_it); fh_it && is_project_ok == false; ++fh_it) {
						TriMesh::FHandle of = limit.face_handle(limit.opposite_halfedge_handle(fh_it));
						he_of_face = initial_control_mesh_.property(fp_hh_icm_, of);
						v0 = limit.from_vertex_handle(he_of_face);
						v1 = limit.to_vertex_handle(he_of_face);
						v2 = limit.to_vertex_handle(limit.next_halfedge_handle(he_of_face));

						OpenMesh::Vec3d n = limit.calc_face_normal(of);
						if ((dot(n_of_f_it, n)) < 0.5) continue; // < -0.5 ���������f_it�漸����ת��, ���������..
						double d = dot(Sk - limit.point(v0), n);//��Sk�������ľ���
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
					// �����Ҳ����Ļ�, �������3�����һ����
					he_of_face = initial_control_mesh_.property(fp_hh_icm_, f_it);
					v0 = limit.from_vertex_handle(he_of_face);
					v1 = limit.to_vertex_handle(he_of_face);
					v2 = limit.to_vertex_handle(limit.next_halfedge_handle(he_of_face));
					std::vector<TriMesh::VHandle> v_array;
					v_array.push_back(v0); v_array.push_back(v1); v_array.push_back(v2); 
					for (int i = 0; i < 3 && is_project_ok == false; ++i) {
						for (TriMesh::VOHIter voh_it(limit, v_array[i]); voh_it && is_project_ok == false; ++ voh_it) {
							OpenMesh::Vec3d n = limit.calc_face_normal(limit.face_handle(voh_it));
							if ((dot(n_of_f_it, n)) < 0.5) continue; // < -0.5 ���������f_it�漸����ת��, ���������..
							double d = dot(Sk - limit.point(v_array[i]), n);//��Sk�������ľ���
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
				if (is_project_ok == false) {// ʵ������еĻ����û����
					// ����ͶӰ�õ��Ľ����.  ��������, ���ﲻ�����bc��Ч, ֻ���ǳ�ֵ���ã��õ����Ľ������С.
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

				//// (2). Sk�ڼ��������ϵ�evaluation position and closest position.
				//// ������Sk��initial_control_mesh_����fh�ϵĲ���������p,
				TriMesh::HHandle hh = initial_control_mesh_.property(fp_hh_icm_, fh_porject_ok);//��Ͱ�ʾ�˲����е�_mesh��ʵ��initial_control_mesh_.
				v0 = initial_control_mesh_.from_vertex_handle(hh);
				v1 = initial_control_mesh_.to_vertex_handle(hh);
				v2 = initial_control_mesh_.to_vertex_handle((initial_control_mesh_.next_halfedge_handle(hh)));
				bc = DGP::calc_barycentric_coordinates(limit.point(v0), limit.point(v1), limit.point(v2), p);

				//	std::cout << p << "/" << (limit.point(v0)*bc[0]+limit.point(v1)*bc[1]+limit.point(v2)*bc[2]) << ", ";

				if (DGP::is_valid_barycentric_coordinate(bc)) {
					mesh_.property(vp_eval_pos_, *it) = qas_->pos(fh_porject_ok, bc[0], bc[1], bc[2], true);
					mesh_.property(vp_eval_, *it) = true; 
					//clp = mesh_.property(vp_eval_pos_, *it); //�����ĳ�ֵ.
			 
				//std::cout << mesh_.property(vp_project_pos_, *it) << "/" << (limit.point(v0)*bc[0]+limit.point(v1)*bc[1]+limit.point(v2)*bc[2]) << ". "; //��һ�µ�.
				//	std::cout << mesh_.property(vp_project_pos_, *it) << "/" << mesh_.property(vp_eval_pos_, *it) << ", "; //��𲻴��
					//std::cout << (Sk - mesh_.property(vp_project_pos_, *it)).norm() << ":=" << std::cout << mesh_.property(vp_project_pos_, *it) << "~=" << mesh_.property(vp_eval_pos_, *it) << ", ";

					if ((Sk - mesh_.property(vp_project_pos_, *it)).norm() < 0.00001) { // <1.0e-5
						// 2009-05-21, ֮ǰ����Sk��ͶӰ�㼸���غ���, skҲ����������,���ǵõ���eval point��������Skȴ����ܴ�.
						mesh_.property(vp_closest_pos_, *it) = mesh_.property(vp_project_pos_, *it);//
						mesh_.property(vp_closest_, *it) = true; 
					} else {
						OpenMesh::Vec3d pp(0, 0, 0), clp(0, 0, 0);
						if (evaluation_get_footpoint_byQAS(initial_control_mesh_, fh_porject_ok, Sk, bc, clp, 0)) {
							mesh_.property(vp_closest_pos_, *it) = clp;//
							mesh_.property(vp_closest_, *it) = true; // 
							//std::cout << (Sk - clp).norm() << " ";
				 
						} else {
							
							//������һ, һ��ԭʼ������Sk�Ҳ���������ֱ���˳�����.
							//return false;
							//��������, ����������ָ��Ϊ0, �������, ����������Sk�������.//2009-05-13
							//vp_eval_pos = OpenMesh::Vec3d(0,0,0);
							//mesh_.property(vp_eval_, *it) = true; 
							//vp_closet_pos =  OpenMesh::Vec3d(0,0,0);//
							//mesh_.property(vp_closest_, *it) = true; 
							//��������, ������, ֱ������. ��Ϊ����Ὣprojection position ������.
							if (mesh_.property(vp_project_, *it)) {
								mesh_.property(vp_closest_pos_, *it) = mesh_.property(vp_project_pos_, *it); 
								mesh_.property(vp_closest_, *it) = true;
								std::cout << "Info: ��project����closest��.\n";
							} else {
								std::cout << "Error: evaluation_get_footpoint �Ҳ��������.\n";
							}
						 
						} 
					}					
				} else { // 2009-05-18-2126 ����������Ҫ�ǵ���parameter�Ľ���Ļ�����û���������..
					std::cout << "Error: ����������bc���ھͲ�������, ������. where_get_pro " << where_get_pro << ".\n";
				}

				if (mesh_.property(vp_closest_, *it)) {
					// 2009-05-21, a breakthrough, �����ͻ�����(Sk��project point�ܽ����Ǻ�closest pointȴ��Զ)�����.
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
	std::cout << "��ͶӰʱ�����������ıȽ�: " << one_pro << ", " << two_pro << ", " << three_pro << ", " << four_pro << ".\n";
	//std::cout << "Info: use projection vs parameterization: " << method_1_count << ": " << mesh2_.n_vertices() - method_1_count << "." << method_3_count <<".\n";
	// 2009-05-13, fan_6475 -> 200v, 5891vs5xx, �ɼ��󲿷ֶ�����project����ʼֵ��.
	
	std::cout << "evaluation_get_all_footpoints(), end.\n";

	// To calculate the fitting error.
	// (1). ����ԭʼ����mesh_��ÿһ�������fitting error. 
	int m = Xm_vertex_.size(), n = simplified_mesh_.n_vertices(), Xm_idx= 0;
	for (TriMesh::VIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++ v_it) {
		if (mesh_.status(v_it.handle()).deleted()) {
			if (mesh_.property(vp_closest_, v_it)) {
				mesh_.property(vp_fitting_error_, v_it) 
					= (mesh_.property(vp_closest_pos_, v_it.handle()) - mesh2_.point(v_it)).norm() / bb.diagonal_lenght;
			} else {
				std::cout << "Error: ��Ȼ���ж��㻹û����������.\n";
			} 	  
			
		} else {
			// method 1. using one-ring
			int count = 0;
			double sum_fitting_error = 0;// һ�����fitting error
			for (TriMesh::VertexVertexIter vv_it(mesh2_, v_it); vv_it; ++ vv_it) { //��Щ����ĳ������Ҳ��û�б�ɾȥ�ĵ�.
				// ����fix one bug. 2009-05-21, vv_it��v_it������.
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
	std::cout << "֮ǰû�б�ɾ���ĵ�����" << Xm_idx << ".\n";
	
	vsplit_array_.clear();
	vsplit_array_sm_.clear();
	// �����е�fitting error��������.
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

	std::set<TriMesh::VHandle, FErr_v> v_fe_arr;//�Զ��㰴fittting error������, �Ϳ���֪������Ķ���.
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
	//	//std::cout << mesh_.property(vp_project_pos_, temp_array2[i]) << "/" << mesh_.property(vp_eval_pos_,temp_array2[i]) << ", "; //�����𲻴�,
	//}std::cout << ".\n";
	// ��mesh_��������Ķ����������ҵ����Ƕ�Ӧ��vertex hierarchy�е�leaf node,
	// ���������Ƕ�Ӧsimplified mesh�ϵĶ���, ������split.
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

	//(2). �����������ÿһ�������fitting error
	fitting_error_array.clear();
	simplified_mesh_.add_property(vp_fitting_error_sm_);
	for (TriMesh::VIter v_it(simplified_mesh_.vertices_begin()), v_end(simplified_mesh_.vertices_end()); v_it != v_end; ++ v_it) {
		int node_handle = simplified_mesh_.property(vp_node_handle_sm_, v_it.handle());
		if (vhierarchy_.node(node_handle).is_active() == false) std::cout << "Error: �������϶����Ӧ��nodeӦ����active�Ŷ�.\n";

		int n_leaves = 0;
		double sum_fitting_error = 0; // �����������Ҷ�ӽڵ��fitting error֮��.
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
	//vsplit_array_sm_ = std::vector<TriMesh::VHandle>(temp_array_sm.end()-n_want_to_split, temp_array_sm.end());//ȡ�������10.
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
// ����ԭʼ����Sk�㶼���Ӧ�������ľ�����Ϊfitting error
void DecimationModel::calc_fitting_error() {
	
}
// -----------------------------------------------------------------------
void DecimationModel::evaluation_process() { //��һ��ֻ���ڸ���������limit surface֮�����.
	std::cout << "\n";
	std::cout << "========��ʼ����ԭʼ���ڼ��������ϵ�����\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=.
	
	if (evaluation_get_all_footpoints() == false) {
		std::cout << "Error: ���ԭʼ���㵽���������ϵ����������ʱ�����.\n";
	} else {
		std::cout << "�����ԭʼ�����ڼ��������ϵ����������.\n"; 
		calc_fitting_error();
	}
	std::cout << "========��������ԭʼ���ڼ��������ϵ�����\n";//������ʾ��40�����ȵ�=,ǰ���9��λ�ÿ�ʼ����,��������8��=. 
}
void DecimationModel::view_fitting_error_mesh2_process() {

	std::vector<double> fitting_error_array; 
	for (TriMesh::VIter v_it(mesh_.vertices_begin()), v_end(mesh_.vertices_end()); v_it != v_end; ++ v_it) {
		fitting_error_array.push_back(mesh_.property(vp_fitting_error_, v_it));
	}
	std::sort(fitting_error_array.begin(), fitting_error_array.end());

	// discard upper and lower 5% = 5/100 = 1/20, ������һͷһβ
	unsigned int n = fitting_error_array.size()-1;// ���� n = 100
	unsigned int i = n / 20;// i = 5
	double min_fitting_error = fitting_error_array[i];// 5
	double max_fitting_error = 0.01;//ò��0.02Ҳ��ͦ���ʵ�2009-05-21 // 94 // fitting_error_array[n-1-i]


	// define uniform color intervalls [v0,v1,v2,v3,v4]
	double v0, v1, v2, v3, v4;
	v0 = min_fitting_error + 0.0/4.0 * (max_fitting_error - min_fitting_error);
	v1 = min_fitting_error + 1.0/4.0 * (max_fitting_error - min_fitting_error);
	v2 = min_fitting_error + 2.0/4.0 * (max_fitting_error - min_fitting_error);
	v3 = min_fitting_error + 3.0/4.0 * (max_fitting_error - min_fitting_error);
	v4 = min_fitting_error + 4.0/4.0 * (max_fitting_error - min_fitting_error);


	// map fitting error to colors �����С���,��ɫ�������
	// blue~ v0	(0,u,255)	v1	(0, 255, 255-u)	v2	(u, 255, 0)	v3	(255, 255-u, 0)	v4	~red
	double fitting_error;
	TriMesh::Color col;
	mesh2_.add_property(vp_color_fitting_error_mesh2_);
	for (TriMesh::VIter v_it=mesh_.vertices_begin(), v_end(mesh_.vertices_end()); v_it!=v_end; ++v_it)
	{ // ��ʼѭ��ÿһ������
		fitting_error = mesh_.property(vp_fitting_error_, v_it);
		col = TriMesh::Color(255,255,255);

		unsigned char u;

		if (fitting_error < v0) 
		{
			col = TriMesh::Color(0, 0, 255);// ��ɫ
		}
		else if (fitting_error > v4) 
		{
			col = TriMesh::Color(255, 0, 0);// ��ɫ
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
	} //����ѭ��ÿһ������
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

	// discard upper and lower 5% = 5/100 = 1/20, ������һͷһβ
	unsigned int n = fitting_error_array.size()-1;// ���� n = 100
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


	// map fitting error to colors �����С���,��ɫ�������
	// blue~ v0	(0,u,255)	v1	(0, 255, 255-u)	v2	(u, 255, 0)	v3	(255, 255-u, 0)	v4	~red
	double fitting_error;
	TriMesh::Color col;
	simplified_mesh_.add_property(vp_color_fitting_error_sm_);
	for (TriMesh::VIter v_it=simplified_mesh_.vertices_begin(), v_end(simplified_mesh_.vertices_end()); v_it!=v_end; ++v_it)
	{ // ��ʼѭ��ÿһ������
		fitting_error = simplified_mesh_.property(vp_fitting_error_sm_, v_it);
		col = TriMesh::Color(255,255,255);

		unsigned char u;

		if (fitting_error < v0) 
		{
			col = TriMesh::Color(0, 0, 255);// ��ɫ
		}
		else if (fitting_error > v4) 
		{
			col = TriMesh::Color(255, 0, 0);// ��ɫ
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
	} //����ѭ��ÿһ������
	// ------------
	is_view_fitting_error_of_sm_ = !is_view_fitting_error_of_sm_;
	is_view_fitting_error_of_mesh2_ = false;
}
