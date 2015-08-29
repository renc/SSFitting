#ifndef dgp_loop_evaluator_h
#define dgp_loop_evaluator_h

#include <vector>
#include "OpenMeshAll.h"
#include "OpenMeshMeshType.h"
#include "Stam98ClassicLoopEval.h" 
#include "VertexFeatureType.h"
#include "Hoppe94LoopEval.h"

namespace DGP {

	template <class TriMeshT>
	class LoopEvaluatorT : public OpenMesh::Utils::Noncopyable {
	public:
		// If every vertex's valence is 6, then this face is regular, or irregular。
		static enum RegularTypeFace { IRREGULAR_FACE = 0, REGULAR_FACE = 1 };
		
		// 这些typedef定义只是为了下面使用的方便, 否则每次都要写typename很烦.
		typedef typename TriMeshT::HalfedgeHandle	HalfedgeHandle;
		typedef typename TriMeshT::EdgeHandle		EdgeHandle;
		typedef typename TriMeshT::VertexHandle		VertexHandle;
		typedef typename TriMeshT::FaceHandle		FaceHandle;

		LoopEvaluatorT(TriMeshT &_m) : mesh_(_m) { //这就暗含了一个要求, 一个LoopEvaluatorT实例只绑定了一个模型.
			cloop_evaluation = new Stam98ClassicLoopEval;   
			h94loop_evaluation  = new Hoppe94LoopEval(cloop_evaluation);
		}
		~LoopEvaluatorT() { delete cloop_evaluation; cloop_evaluation = NULL; }
		Stam98ClassicLoopEval * stam98ClassicLoopEvalObject() { return cloop_evaluation; }

		// 利用stam98的工作对classic loop的evaluation.统一了regular and irregular patches.
		// u ----------> w=1  注意参数坐标(v, w)是相对应那3个有序点来说的.
		// |     _hh=h84(regular)
		// |
		// |_hh = h0=h12(irregular)
		// \/ v = 1
		// 对于没有考虑特征情况的classic loop subdivsion surface使用Stam98的方法来求
		// 面mesh_.face_handle(_hh)上重心坐标为(1-_v-_w, _v, _w)的精确坐标, 输出为_Sk and _p.
		// Note: 
		// _hh指向的三角形最多只有一个非规则点, 只能是_hh的from顶点,
		//   也就是说in regular patch _hh is h84's next halfedge, in irregular patch _hh is h12 directly.
		//   这就暗含了_hh已经指定了这个面枚举点的顺序次序.
		// _v, _w重心坐标必须有效范围之内, 并且是对应_hh的三角形来说的.		
		bool cloop_eval_pos(HalfedgeHandle _hh, // in regular patch _hh is h84's prev halfedge, in irregular patch _hh is h12 directly
			float _v, float _w, // input parameters.
			OpenMesh::Vec3d &_Sk, OpenMesh::Vec3d &_p,  OpenMesh::Vec3d &_dv,  OpenMesh::Vec3d &_dw) // output parameters.
		{
			_Sk = OpenMesh::Vec3d(0, 0, 0); _p = OpenMesh::Vec3d(0, 0, 0);
			HalfedgeHandle h0 = _hh;
			HalfedgeHandle h1 = mesh_.next_halfedge_handle(h0);
			HalfedgeHandle h2 = mesh_.next_halfedge_handle(h1);
			VertexHandle v0 = mesh_.from_vertex_handle(h0);
			VertexHandle v1 = mesh_.to_vertex_handle(h0);
			VertexHandle v2 = mesh_.to_vertex_handle(h1);
			if (mesh_.valence(v1) != 6 || mesh_.valence(v2) != 6) {
				//只能是参数_hh的from_vertex_handle可能非规则, 其它两点都要求规则
				std::cerr << "Error: cloop_eval_pos(), v1 and v2 are supposed to be regular.\n";
			}
			
			if (mesh_.valence(v0) == 6) { // regular patch
				HalfedgeHandle h84 = h2; 
				std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
				cloop_stam98_regular_get_c12(mesh_, h84, c12);

				//_Sk = c12[3] * (1-_v-_w) + c12[6] * _v + c12[7] * _w;
				//if (cloop_evaluation->evaluate_regular_patch(c12, _v, _w, _p, _dv, _dw) == false) {
				//	std::cout << "Error: LoopEvaluatorT::cloop_eval_pos, regular.\n";
				//	return false;
				//}
				return cloop_eval_pos(REGULAR_FACE, c12, _v, _w, _Sk, _p, _dv, _dw); 
			} 
			else { // irregular patch
				HalfedgeHandle h12 = _hh; 
				int N = mesh_.valence(mesh_.from_vertex_handle(h12));

				std::vector<OpenMesh::Vec3d> cN6(N+6, OpenMesh::Vec3d(0, 0, 0));
				cloop_stam98_irregular_get_cN6(mesh_, h12, cN6);

				//_Sk = cN6[0] * (1-_v-_w) + cN6[1] * _v + cN6[N] * _w;
				//if (cloop_evaluation->evaluate_irregular_patch(cN6, _v, _w, _p, _dv, _dw) == false) {
				//	std::cout << "Error: LoopEvaluatorT::cloop_eval_pos, irregular.\n";
				//	return false; 
				//}
				return cloop_eval_pos(IRREGULAR_FACE, cN6, _v, _w, _Sk, _p, _dv, _dw);
			}
			return true;
		}
		bool cloop_eval_pos(RegularTypeFace _regular, const std::vector<OpenMesh::Vec3d> &_c, // these parameters can be derived from _hh.
			float _v, float _w, // input parameters.
			OpenMesh::Vec3d &_Sk, OpenMesh::Vec3d &_p,  OpenMesh::Vec3d &_dv,  OpenMesh::Vec3d &_dw) // output parameters.
		{
			_Sk = OpenMesh::Vec3d(0, 0, 0); _p = OpenMesh::Vec3d(0, 0, 0);
			// regular patch
			if (_regular == REGULAR_FACE) {
				std::vector<OpenMesh::Vec3d> c12(_c); 
				_Sk = c12[3] * (1-_v-_w) + c12[6] * _v + c12[7] * _w;

				if (cloop_evaluation->evaluate_regular_patch(c12, _v, _w, _p, _dv, _dw) == false) {
					std::cout << "Error: LoopEvaluatorT::cloop_eval_pos, regular.\n";
					return false;
				}
			} // irregular patch
			else {
				int N = _c.size() -6;// vertex 's valence.
				std::vector<OpenMesh::Vec3d> cN6(_c); 
				_Sk = cN6[0] * (1-_v-_w) + cN6[1] * _v + cN6[N] * _w; 

				if (cloop_evaluation->evaluate_irregular_patch(cN6, _v, _w, _p, _dv, _dw) == false) {
					std::cout << "Error: LoopEvaluatorT::cloop_eval_pos, irregular.\n";
					return false; 
				}				 
			}
			return true;
		}
		// 用Newton iteration方法来求最近点.
		// 对于没有考虑特征情况的classic loop subdivision surface使用Newton Iterator
		// 来求一个给定点Sk在_hh所指向的那一片极限曲面上的最近点. 牛顿迭代的初始值是h0=(_v, _w).
		// Note; 调用前确保_hh指向的三角形最多只有一个非规则点, 只能是_hh的from顶点.
		//       确保_v, _w在重心坐标的有效范围之内.
		//       这只是一个测试版本因为例如搜索到边界时候呢，没有处理的.
		bool cloop_closest_pos(OpenMesh::Vec3d _Sk, HalfedgeHandle _hh, OpenMesh::Vec3d &_clp, float _v = 0.25, float _w = 0.25) 
		{
			std::vector<HalfedgeHandle> h_arr(3, HalfedgeHandle(-1));//面上的3个半边
			std::vector<VertexHandle> v_arr(3, VertexHandle(-1));//面上的3个顶点
			h_arr[0] = _hh; v_arr[0] = mesh_.to_vertex_handle(h_arr[0]);
			h_arr[1] = mesh_.next_halfedge_handle(h_arr[0]); v_arr[1] = mesh_.to_vertex_handle(h_arr[1]);
			h_arr[2] = mesh_.next_halfedge_handle(h_arr[1]); v_arr[2] = mesh_.to_vertex_handle(h_arr[2]); 
			if (mesh_.valence(v_arr[0]) != 6 || mesh_.valence(v_arr[1]) != 6) {
				std::cerr << "Error: cloop_closest_pos, 只能是参数_hh的from_vertex_handle可能非规则, 其它两点都要求规则.\n";
				return false;
			}

			HalfedgeHandle h12 = _hh; //the h12's meaning, check the paper "Evaluation of Loop Subdivision Surfaces". 
			int N = mesh_.valence(mesh_.from_vertex_handle(h12));

			std::vector<OpenMesh::Vec3d> cN6(N+6, OpenMesh::Vec3d(0, 0, 0));
			cloop_stam98_irregular_get_cN6(mesh_, h12, cN6);

			bool is_test = false;
			std::cout << "Search ...\n";
			// (1), search descent direction h(delta_v, delta_w)
			const std::vector<OpenMesh::Vec3d> & c = cN6;
			double v = _v, w = _w;

			const int Max_Iter_Times = 150;//最大的迭代次数.
			const double Max_Step_length = 0.5; //最大步长
			int iter_times = 0; 
			for (iter_times = 0; iter_times < Max_Iter_Times; ++iter_times) 
			{
				OpenMesh::Vec3d dStam_dv(0, 0, 0), dStam_dw(0, 0, 0); 

				if (cloop_evaluation->evaluate_irregular_patch(c, v, w, OpenMesh::Vec3d(0,0,0), dStam_dv, dStam_dw) == false) 
				{ //Stam(v, w)对w求导.It's a vector:[dStam(v, w)1/dw, dStam(v, w)2/dw, dStam(v, w)3/dw]
					std::cout << "Error: dStam/dv, dStam/dw.\n";
				}
				if (is_test) std::cout << dStam_dv << ", " << dStam_dw << ".\n";
				// 下面手动化简并解这个二元方程组(两个方程).JTJh=-JTf
				double B[2];//-JTf
				B[0] = -dot(dStam_dv, (_clp - _Sk));
				B[1] = -dot(dStam_dw, (_clp - _Sk)); //std::cout << B[0] << ", " << B[1] << ".\n"; // for test
				double A[2][2];//JTJ
				A[0][0] = dot(dStam_dv, dStam_dv); A[0][1] = dot(dStam_dv, dStam_dw); if (is_test) std::cout << A[0][0] << ", " << A[0][1] << "; "; // for test//
				A[1][0] = dot(dStam_dw, dStam_dv); A[1][1] = dot(dStam_dw, dStam_dw); if (is_test) std::cout << A[1][0] << ", " << A[1][1] << ".\n"; // for test//
				

				enum IterationType { Steepest_descent, Gauss_newton };
				IterationType itype = Gauss_newton;
				double delta_v, delta_w;//search direction h = [delta_v, delta_w]T
				double step_len = 1;//-0.02;//初始步长.
				if (itype == Steepest_descent) { //结果是这个即使使用全步长(1or-1)都是收敛极慢的! 迭代50次只是前进了一点点.
					delta_v = B[0];
					delta_w = B[1];
					step_len = -1;
				} else {
					if (fabs(A[0][0] - A[0][1]) < 1.e-6) {
						delta_v = B[0] * 0.5 / A[0][0];
						delta_w = B[0] * 0.5 / A[0][0];
					}
					else {
						delta_v = (B[0]*A[1][1] - B[1]*A[0][1]) / (A[0][0]*A[1][1] - A[0][1]*A[1][0]);	if (is_test) std::cout << "delta_v: " << delta_v << "; ";//
						delta_w = (B[0] - A[0][0] * delta_v) / A[0][1];									if (is_test) std::cout << "delta_w: " << delta_w << ".\n";//
					}
					step_len = -0.01;
				}		

				//sd0 = c[0]*(1-v-w) + c[1]*v + c[c.size()-6]*w;
				//sd1 = c[0]*(1-v-w-delta_v-delta_w) + c[1]*(v+delta_v) + c[c.size()-6]*(w+delta_w);

				delta_v *= step_len;
				delta_w *= step_len;
				if ((v + delta_v) < 0) {
					std::cout << "v < 0.\n"; break;
				} else if ((v + delta_v) > 1) {
					std::cout << "v > 1.\n"; break;
				} else if ((w + delta_w) < 0) {
					std::cout << "w < 0.\n"; break;
				} else if ((w + delta_w) > 1) {
					std::cout << "w > 1.\n"; break;
				} else if ((v + delta_v + w + delta_w) > 1) {
					std::cout << "v + w > 1.\n"; break;
				} else { 
					//还在同一个面里面找着
					//std::cout << "v: " << v << ", w: " << w << ". ";
					v += delta_v; 
					w += delta_w;//std::cout << "v': " << v << ", w': " << w << ".\n";

					double distance_before = (_clp-_Sk).norm();
					OpenMesh::Vec3d clp_k_1 = _clp; //保存前一个最近点坐标 
					if (cloop_evaluation->evaluate_irregular_patch(c, v, w, _clp, OpenMesh::Vec3d(0,0,0), OpenMesh::Vec3d(0,0,0)) == false) {
						std::cout << "Error: stam98 test regular patch.\n";
					}
					double distance_after = (_clp-_Sk).norm();
					if (is_test) std::cout << "_clp: " << _clp << ". dis(_clp-_Sk): " << distance_after << ".\n";
					if (distance_before < distance_after) {
						if (is_test) std::cout << "Iteration ends, for distance_before < distance_after, OK.\n";
						_clp = clp_k_1; 
						break;
					}
				}
			}
			if (iter_times == Max_Iter_Times) 
				if (is_test) std::cout << "Info: 由于超过给定迭代次数而停止的.\n";
		}

	
		// v0            v2
		//-------------> w
		//|
		//|_hh=h0
		//|
		//\/ v            If there is 1 irregular vertex, it must be v0.
		// v1             So, make sure the _hh's from vertex is irregular if these is one.
		//bool h94loop_eval_pos(Stam98ClassicLoopEval::BasisFuncType _basisfun_type, // input parameters, ORI, Der_V, Der_W
		//	HalfedgeHandle _hh, // in regular patch _hh is h84's prev halfedge, in irregular patch _hh is h12 directly
		//	float _v, float _w, // input parameters.
		//	OpenMesh::VPropHandleT<v_feature_type> _vp_type, // vertex property of this mesh, make sure this is added before.
		//	OpenMesh::Vec3d &_Sk, OpenMesh::Vec3d &_p, OpenMesh::Vec3d &_dv, OpenMesh::Vec3d &_dw) 
		//{
		//	std::cout << "Error: This function is replaced by the overrided function.\n"; 
		//	return false;
		//} // end of function h94loop_eval_pos(...)
		bool h94loop_eval_pos(
			HalfedgeHandle _hh, // in regular patch _hh is h84's next halfedge, in irregular patch _hh is h12 directly
			float _v, float _w, // input parameters.
			OpenMesh::VPropHandleT<v_feature_type> _vp_type, // vertex property of this mesh, make sure this is added before.
			OpenMesh::Vec3d &_Sk, OpenMesh::Vec3d &_p, OpenMesh::Vec3d &_dv, OpenMesh::Vec3d &_dw )
		{
			_Sk = OpenMesh::Vec3d(0, 0, 0); _p = OpenMesh::Vec3d(0, 0, 0);
			HalfedgeHandle h0 = _hh;
			HalfedgeHandle h1 = mesh_.next_halfedge_handle(h0);
			HalfedgeHandle h2 = mesh_.next_halfedge_handle(h1);
			VertexHandle v0 = mesh_.from_vertex_handle(h0);
			VertexHandle v1 = mesh_.to_vertex_handle(h0);
			VertexHandle v2 = mesh_.to_vertex_handle(h1);
			if (mesh_.valence(v1) != 6 || mesh_.valence(v2) != 6) {
				std::cerr << "Error: h94loop_eval_pos(), v1 and v2 are supposed to be regular.\n";
			}
			_Sk = mesh_.point(v0) * (1-_v-_w) + mesh_.point(v1) * _v + mesh_.point(v2) * _w;

			if (mesh_.valence(v0) == 6) {
				// feature, regular.
				// Make sure to the match one of the 7 kinds of pattern.
				feature_regular_type frt = rtNull;
				HalfedgeHandle h84(-1);
				float v = _v, w = _w; //These need to be corrected according to the h84.

				int n_feature_vertics = 0;
				if (mesh_.property(_vp_type, v0) != DGP::SMOOTH_VFT) ++ n_feature_vertics;
				if (mesh_.property(_vp_type, v1) != DGP::SMOOTH_VFT) ++ n_feature_vertics;
				if (mesh_.property(_vp_type, v2) != DGP::SMOOTH_VFT) ++ n_feature_vertics;
				if (n_feature_vertics == 0) {
					std::cout << "Error: no feature stituation, then should be dealed in cloop_eval_pos() instead.\n";
				} else if (n_feature_vertics == 1) { // 1种情况
					if (mesh_.property(_vp_type, v0) != DGP::SMOOTH_VFT) {
						h84 = h2;
						v = _v; w = _w;
					} else if (mesh_.property(_vp_type, v1) != DGP::SMOOTH_VFT) { 
						h84 = h0;
						v = _w; w = 1 - _v - _w;
					} else if (mesh_.property(_vp_type, v2) != DGP::SMOOTH_VFT) {
						h84 = h1;
						v = 1 - _v - _w; w = _v;
					}

					frt = rt1v0e;
				} else if (n_feature_vertics == 2) { // 2种情况
					if (mesh_.property(_vp_type, v0) != DGP::SMOOTH_VFT && mesh_.property(_vp_type, v1) != DGP::SMOOTH_VFT) {
						h84 = h2;
						v = _v; w = _w;
						if (mesh_.status(mesh_.edge_handle(h0)).feature()) {
							frt = rt2v1e;  
						} else {
							frt = rt2v0e;  
						}						
					}
					if (mesh_.property(_vp_type, v1) != DGP::SMOOTH_VFT && mesh_.property(_vp_type, v2) != DGP::SMOOTH_VFT) {
						h84 = h0; 
						v = _w; w = 1 - _v - _w;
						if (mesh_.status(mesh_.edge_handle(h1)).feature()) {
							frt = rt2v1e;  
						} else {
							frt = rt2v0e;  
						}
					}
					if (mesh_.property(_vp_type, v2) != DGP::SMOOTH_VFT && mesh_.property(_vp_type, v0) != DGP::SMOOTH_VFT) {
						h84 = h1; 
						v = 1 - _v - _w; w = _v;
						if (mesh_.status(mesh_.edge_handle(h2)).feature()) {
							frt = rt2v1e;  
						} else {
							frt = rt2v0e;  
						}
					}

				} else if (n_feature_vertics == 3) { // 4种情况
					if (mesh_.status(mesh_.edge_handle(h0)).feature() && mesh_.status(mesh_.edge_handle(h1)).feature() && mesh_.status(mesh_.edge_handle(h2)).feature()) {
						h84 = h2; v = _v; w = _w; frt = rt3v3e;
					} else if (mesh_.status(mesh_.edge_handle(h0)).feature() && mesh_.status(mesh_.edge_handle(h1)).feature() ) { //两个特征边.
						h84 = h0; v = _w; w = 1 - _v - _w; frt = rt3v2e;
					} else if (mesh_.status(mesh_.edge_handle(h1)).feature() && mesh_.status(mesh_.edge_handle(h2)).feature() ) {
						h84 = h1; v = 1 - _v - _w; w = _v; frt = rt3v2e;
					} else if (mesh_.status(mesh_.edge_handle(h2)).feature() && mesh_.status(mesh_.edge_handle(h0)).feature() ) {
						h84 = h2; v = _v; w = _w; frt = rt3v2e;
					} else if (mesh_.status(mesh_.edge_handle(h0)).feature()) { //一个特征边.
						h84 = h2; v = _v; w = _w; frt = rt3v1e;
					} else if (mesh_.status(mesh_.edge_handle(h1)).feature()) {
						h84 = h0; v = _w; w = 1 - _v - _w; frt = rt3v1e;
					} else if (mesh_.status(mesh_.edge_handle(h2)).feature()) {
						h84 = h1; v = 1 - _v - _w; w = _v; frt = rt3v1e;
					} else { //没有特征边.
						h84 = h2; v = _v; w = _w; frt = rt3v0e;
					} 					
				}  
				if (h84.is_valid() == false) std::cout << "Error: h94loop_eval_pos(): 本不该有的错误.\n";


				std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));// 严重注意: 这12个控制点的排序以stam98 "evaluation of loop subdivision surfaces". Fig1为准.
				std::vector<bool> fe15(15, false);//12个顶点之间的15个边是否是特征边
				DGP::h94loop_eval_regular_get_c12_h15(mesh_, h84, c12, fe15);

				//return h94loop_evaluation->eval_regular_1in7(ORI, frt, c12, fe15, v, w, _p);这个接口只是求出pos, 没有切向量.
				return h94loop_evaluation->eval_regular_1in7(frt, c12, fe15, v, w, _p, _dv, _dw);

			} else {
				// feature, irregular.
				TriMesh::HHandle h12 = h0;//先找到这个边. 

				// 枚举N+6个点.
				int N = mesh_.valence(mesh_.from_vertex_handle(h12));  if (N == 6) std::cout << "Error: 这里应该是irregular.\n";
				int K = N + 6;
				std::vector<OpenMesh::Vec3d> c(K, OpenMesh::Vec3d(0, 0, 0));// 严重注意: 这K个控制点的排序以stam98 "evaluation of loop subdivision surfaces". Fig2为准.
				std::vector<bool> hN9(N+4+4+1, false);

				DGP::h94loop_eval_irregular_get_cN6_hN9(mesh_, h12, c, hN9);

				feature_irregular_type ft = irtNull;
				TriMesh::HHandle h12prev = mesh_.prev_halfedge_handle(h12), h12next(mesh_.next_halfedge_handle(h12));
				TriMesh::VHandle vhup = mesh_.to_vertex_handle(h12prev), vhleft(mesh_.to_vertex_handle(h12)), vhright(mesh_.to_vertex_handle(h12next));
				if (mesh_.property(_vp_type, vhup) == DGP::SMOOTH_VFT) {
					if (mesh_.status(mesh_.edge_handle(h12next)).feature()) {
						ft = irt2v1e_bottom;
					} else if (mesh_.property(_vp_type, vhleft) != DGP::SMOOTH_VFT && mesh_.property(_vp_type, vhright) != DGP::SMOOTH_VFT) {
						ft = irt2v0e_bottom;
					}  else if (mesh_.property(_vp_type, vhleft) != DGP::SMOOTH_VFT ) {
						ft = irt1v0e_left;
					}  else if (mesh_.property(_vp_type, vhright) != DGP::SMOOTH_VFT ) {
						ft = irt1v0e_right;
					} 
				} else {
					if (mesh_.status(mesh_.edge_handle(h12prev)).feature() && mesh_.status(mesh_.edge_handle(h12)).feature() && mesh_.status(mesh_.edge_handle(h12next)).feature()) {
						ft = irt3v3e;
					} else if (mesh_.status(mesh_.edge_handle(h12prev)).feature() && mesh_.status(mesh_.edge_handle(h12)).feature() ) {
						ft = irt3v2e_up; //这个在例子ne....160.off中没有错.
					} else if (mesh_.status(mesh_.edge_handle(h12)).feature() && mesh_.status(mesh_.edge_handle(h12next)).feature()) {
						ft = irt3v2e_left;
					} else if (mesh_.status(mesh_.edge_handle(h12next)).feature() && mesh_.status(mesh_.edge_handle(h12prev)).feature()) {
						ft = irt3v2e_right;
					} else if (mesh_.status(mesh_.edge_handle(h12next)).feature()) {
						ft = irt3v1e_bottom;
					} else if (mesh_.status(mesh_.edge_handle(h12)).feature()) {
						if (mesh_.property(_vp_type, vhright) != DGP::SMOOTH_VFT) ft = irt3v1e_left;
						else ft = irt2v1e_left; // nefertiti.160.off来看也没有错.
					} else if (mesh_.status(mesh_.edge_handle(h12prev)).feature()) {
						if (mesh_.property(_vp_type, vhleft) != DGP::SMOOTH_VFT) ft = irt3v1e_right;
						else ft = irt2v1e_right;// nefertiti.160.off来看也没有错.
					} else if (mesh_.property(_vp_type, vhleft) != DGP::SMOOTH_VFT && mesh_.property(_vp_type, vhright) != DGP::SMOOTH_VFT) {
						ft = irt3v0e;
					}  else if (mesh_.property(_vp_type, vhleft) != DGP::SMOOTH_VFT ) {
						ft = irt2v0e_left;
					}  else if (mesh_.property(_vp_type, vhright) != DGP::SMOOTH_VFT ) {
						ft = irt2v0e_right;
					} else if (mesh_.property(_vp_type, vhleft) == DGP::SMOOTH_VFT && mesh_.property(_vp_type, vhright) == DGP::SMOOTH_VFT) {
						ft = irt1v0e_up;
					}  
				}

				return h94loop_evaluation->eval_irregular_1in17(ft, c, hN9, _v, _w, _p, _dv, _dw);
			}
			return false;
		}
	private:
		TriMeshT &mesh_;
		Stam98ClassicLoopEval *cloop_evaluation;	
		Hoppe94LoopEval *h94loop_evaluation;

	}; // end of class LoopEvaluatorT


	// 下面是几个辅助函数.
	// 半边h84/h12, N+6个顶点的序号都以Evaluation of Loop Subdivision Surfaces上的注释为准.
	// 在
	template <typename TriMeshT> 
	void cloop_stam98_regular_get_c12(const TriMeshT& _mesh, const typename TriMeshT::HHandle h84, std::vector<OpenMesh::Vec3d> &c12) {
		typename TriMeshT::HHandle h47 = _mesh.next_halfedge_handle(h84);
		typename TriMeshT::HHandle h78 = _mesh.next_halfedge_handle(h47);
		typename TriMeshT::HHandle h85 = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(h84));
		typename TriMeshT::HHandle h74 = _mesh.opposite_halfedge_handle(h47);
		typename TriMeshT::HHandle h41 = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(_mesh.next_halfedge_handle(h74)));
		typename TriMeshT::HHandle h42 = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(h41));
		typename TriMeshT::HHandle h73 = _mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(h74));
		typename TriMeshT::HHandle h76 = _mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(h73));
		typename TriMeshT::HHandle h711 = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(h78));
		typename TriMeshT::HHandle h811 = _mesh.opposite_halfedge_handle(_mesh.next_halfedge_handle(h711));

		c12[0] = _mesh.point(_mesh.to_vertex_handle(h41));
		c12[1] = _mesh.point(_mesh.to_vertex_handle(h42));
		c12[2] = _mesh.point(_mesh.to_vertex_handle(h73));
		c12[3] = _mesh.point(_mesh.to_vertex_handle(h84));
		c12[4] = _mesh.point(_mesh.to_vertex_handle(h85));
		c12[5] = _mesh.point(_mesh.to_vertex_handle(h76));
		c12[6] = _mesh.point(_mesh.to_vertex_handle(h47));
		c12[7] = _mesh.point(_mesh.to_vertex_handle(h78));
		c12[8] = _mesh.point(_mesh.to_vertex_handle(_mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(h85))));
		c12[9] = _mesh.point(_mesh.to_vertex_handle(_mesh.next_halfedge_handle(h76)));
		c12[10] = _mesh.point(_mesh.to_vertex_handle(h711));
		c12[11] = _mesh.point(_mesh.to_vertex_handle(_mesh.next_halfedge_handle(h811)));
	}
	template <typename TriMeshT> 
	void cloop_stam98_irregular_get_cN6(const TriMeshT& _mesh, const typename TriMeshT::HHandle h12, std::vector<OpenMesh::Vec3d> &cN6)  {
		int N = cN6.size() - 6;

		cN6[0] = _mesh.point(_mesh.from_vertex_handle(h12));
		typename TriMeshT::HHandle tem_h = h12;
		for (int i = 1; i <= N; ++i) {
			cN6[i] = _mesh.point(_mesh.to_vertex_handle(tem_h));
			tem_h = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(tem_h));
			if (tem_h == h12) break;
		}
		typename TriMeshT::HHandle h2N2 = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(_mesh.next_halfedge_handle(h12)));
		cN6[N+1] = _mesh.point(_mesh.to_vertex_handle(h2N2));
		typename TriMeshT::HHandle h2N3 = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(h2N2));
		cN6[N+2] = _mesh.point(_mesh.to_vertex_handle(h2N3));
		typename TriMeshT::HHandle h2N4 = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(h2N3));
		cN6[N+3] = _mesh.point(_mesh.to_vertex_handle(h2N4));
		typename TriMeshT::HHandle hN1N5 = _mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(_mesh.opposite_halfedge_handle(_mesh.next_halfedge_handle(h2N2))));
		cN6[N+4] = _mesh.point(_mesh.to_vertex_handle(hN1N5));
		cN6[N+5] = _mesh.point(_mesh.to_vertex_handle(_mesh.next_halfedge_handle(hN1N5)));

	}
	template <typename TriMeshT> 
	void h94loop_eval_regular_get_c12_h15(const TriMeshT& _mesh, const typename TriMeshT::HHandle h84, std::vector<OpenMesh::Vec3d>& c12, std::vector<bool>& h15) {
		// 枚举序号严格按照Fig 1. Of Stam98 Evaluation of loop subdivision surfaces.
		typename TriMeshT::HHandle h47 = _mesh.next_halfedge_handle(h84);
		typename TriMeshT::HHandle h78 = _mesh.next_halfedge_handle(h47);
		typename TriMeshT::HHandle h85 = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(h84));
		typename TriMeshT::HHandle h74 = _mesh.opposite_halfedge_handle(h47);
		typename TriMeshT::HHandle h41 = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(_mesh.next_halfedge_handle(h74)));
		typename TriMeshT::HHandle h42 = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(h41));
		typename TriMeshT::HHandle h73 = _mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(h74));
		typename TriMeshT::HHandle h76 = _mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(h73));
		typename TriMeshT::HHandle h711 = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(h78));
		typename TriMeshT::HHandle h811 = _mesh.opposite_halfedge_handle(_mesh.next_halfedge_handle(h711));

		c12[0] = _mesh.point(_mesh.to_vertex_handle(h41));
		c12[1] = _mesh.point(_mesh.to_vertex_handle(h42));
		c12[2] = _mesh.point(_mesh.to_vertex_handle(h73));
		c12[3] = _mesh.point(_mesh.to_vertex_handle(h84));
		c12[4] = _mesh.point(_mesh.to_vertex_handle(h85));
		c12[5] = _mesh.point(_mesh.to_vertex_handle(h76));
		c12[6] = _mesh.point(_mesh.to_vertex_handle(h47));
		c12[7] = _mesh.point(_mesh.to_vertex_handle(h78));
		c12[8] = _mesh.point(_mesh.to_vertex_handle(_mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(h85))));
		c12[9] = _mesh.point(_mesh.to_vertex_handle(_mesh.next_halfedge_handle(h76)));
		c12[10] = _mesh.point(_mesh.to_vertex_handle(h711));
		c12[11] = _mesh.point(_mesh.to_vertex_handle(_mesh.next_halfedge_handle(h811)));


		if (_mesh.status(_mesh.edge_handle(h41)).feature()) h15[0] = true; 
		if (_mesh.status(_mesh.edge_handle(h42)).feature()) h15[1] = true;  
		if (_mesh.status(_mesh.edge_handle(_mesh.next_halfedge_handle(h74))).feature()) h15[2] = true; 
		if (_mesh.status(_mesh.edge_handle(_mesh.opposite_halfedge_handle(_mesh.next_halfedge_handle(h85)))).feature()) h15[3] = true; 
		if (_mesh.status(_mesh.edge_handle(h73)).feature()) h15[4] = true;
		if (_mesh.status(_mesh.edge_handle(h47)).feature()) h15[5] = true;
		if (_mesh.status(_mesh.edge_handle(h84)).feature()) h15[6] = true;
		if (_mesh.status(_mesh.edge_handle(h85)).feature()) h15[7] = true;
		if (_mesh.status(_mesh.edge_handle(h76)).feature()) h15[8] = true;
		if (_mesh.status(_mesh.edge_handle(h78)).feature()) h15[9] = true;
		if (_mesh.status(_mesh.edge_handle(_mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(h85)))).feature()) h15[10] = true;
		if (_mesh.status(_mesh.edge_handle(_mesh.prev_halfedge_handle(h76))).feature()) h15[11] = true;
		if (_mesh.status(_mesh.edge_handle(h711)).feature()) h15[12] = true;
		if (_mesh.status(_mesh.edge_handle(h811)).feature()) h15[13] = true;
		if (_mesh.status(_mesh.edge_handle(_mesh.prev_halfedge_handle(h811))).feature()) h15[14] = true;
	}
	template <typename TriMeshT> 
	void h94loop_eval_irregular_get_cN6_hN9(const TriMeshT &_mesh, const typename TriMeshT::HHandle h12, std::vector<OpenMesh::Vec3d>& cN6, std::vector<bool>& hN9) {
		// 枚举序号严格按照Fig 1. Of Stam98 Evaluation of loop subdivision surfaces.
		int N = cN6.size() - 6;

		cN6[0] = _mesh.point(_mesh.from_vertex_handle(h12));
		typename TriMeshT::HHandle tem_h = h12;
		for (int i = 1; i <= N; ++i) {
			cN6[i] = _mesh.point(_mesh.to_vertex_handle(tem_h));
			tem_h = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(tem_h));
			if (tem_h == h12) break;
		}
		typename TriMeshT::HHandle h2N2 = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(_mesh.next_halfedge_handle(h12)));
		cN6[N+1] = _mesh.point(_mesh.to_vertex_handle(h2N2));
		typename TriMeshT::HHandle h2N3 = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(h2N2));
		cN6[N+2] = _mesh.point(_mesh.to_vertex_handle(h2N3));
		typename TriMeshT::HHandle h2N4 = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(h2N3));
		cN6[N+3] = _mesh.point(_mesh.to_vertex_handle(h2N4));
		typename TriMeshT::HHandle hN1N5 = _mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(_mesh.opposite_halfedge_handle(_mesh.next_halfedge_handle(h2N2))));
		cN6[N+4] = _mesh.point(_mesh.to_vertex_handle(hN1N5));
		cN6[N+5] = _mesh.point(_mesh.to_vertex_handle(_mesh.next_halfedge_handle(hN1N5)));

		//std::cout << "_N: " << N << ".\n"; // for test only. 
		std::vector<typename TriMeshT::HHandle > h(N+4+4+1);
		tem_h = h12;
		for (int i = 0; i < N; ++i) { // h[0], ... h[N-1]
			h[i] = tem_h; 
			tem_h = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(tem_h)); 
		}
		if (tem_h != h[0]) std::cout << "Error: at this time, they should be equal.\n";
		tem_h = _mesh.next_halfedge_handle(h[0]);
		for (int i = 0; i < 5; ++i) { // h[N], ... h[N+4]
			h[N + i] = tem_h;
			tem_h = _mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(tem_h)); 
		}
		if (tem_h != _mesh.opposite_halfedge_handle(h[0])) std::cout << "Error: at this time, they should be equal again.\n";
		tem_h = _mesh.opposite_halfedge_handle(_mesh.next_halfedge_handle(h[N+1]));
		for (int i = 5; i <= 8; ++i) {
			h[N + i] = tem_h;
			tem_h = _mesh.opposite_halfedge_handle(_mesh.prev_halfedge_handle(tem_h));
		}

		if ((h[0] !=h12) || (h[N+1] !=h2N2) || (h[N+3] !=h2N4) || (h[N+6] !=hN1N5)) {
			std::cout << "Error: h94loop_eval_irregular() h[].\n"; return;
		}


		for (int i = 0; i < N+4+4+1; ++i) {
			hN9[i] = _mesh.status(_mesh.edge_handle(h[i])).feature();	 
		} 
	}

	// Identify whether this face is regular or not,
	// if yes, return one halfedge arbitrarily;
	// if no, return the halfedge which from vertex is irregular.
	template <typename TriMeshT>
	bool is_regular(const TriMeshT& _mesh, const typename TriMeshT::FaceHandle _f, typename TriMeshT::HalfedgeHandle &_h) {
		typename TriMeshT::HalfedgeHandle h0 = _mesh.halfedge_handle(_f);
		typename TriMeshT::HalfedgeHandle h1 = _mesh.next_halfedge_handle(h0);
		typename TriMeshT::HalfedgeHandle h2 = _mesh.next_halfedge_handle(h1);

		if (_mesh.valence(_mesh.from_vertex_handle(h0)) != 6) {
			_h = h0; return false;
		}
		if (_mesh.valence(_mesh.from_vertex_handle(h1)) != 6) {
			_h = h1; return false;
		}
		if (_mesh.valence(_mesh.from_vertex_handle(h2)) != 6) {
			_h = h2; return false;
		}
		_h = h0; return true;
	}
} // end of namesplace DGP

#endif //dgp_loop_evaluator_h