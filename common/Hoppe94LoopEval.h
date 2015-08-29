// ***************************************************************
//  
//  Hoppe94LoopEval 
//  Copyright (C) 2008 - by rencanjiang All Rights Reserved
//  rencanjiang@163.com	
//	-------------------------------------------------------------
//  
//  version: 1.0
//  data: 04/27/2008
// ***************************************************************
//  
//  File description:
//		尝试怎么evaluation of Hoppe94's feature perserving Loop subdivision surfaces.
//		注意这工作要结合我之前写的Stam98ClassicLoopEval来做.
//		2008.05.14还没有处理边界情况.
// ***************************************************************

#ifndef Hoppe94LoopEval_H_
#define Hoppe94LoopEval_H_

#include <iostream>
#include <vector>

#include "barycentric.h"
#include "Stam98ClassicLoopEval.h"//这个是基础啊.

enum feature_regular_type { rtNull = 0, rt1v0e = 1, rt2v0e = 2, rt2v1e = 3, rt3v0e = 4, rt3v1e = 5, rt3v2e = 6, rt3v3e = 7 };//1~7有意义
enum feature_irregular_type { irtNull = 0, 
irt1v0e_up,	irt1v0e_left,	irt1v0e_right, 
irt2v0e_left, irt2v0e_right, irt2v0e_bottom,
irt3v0e = 7,
irt2v1e_left, irt2v1e_right, irt2v1e_bottom,
irt3v1e_left, irt3v1e_right, irt3v1e_bottom,
irt3v2e_up,	irt3v2e_left, irt3v2e_right, 
irt3v3e
}; //1 ~ 17有意义.

class Hoppe94LoopEval {

public:
	//第一个default constructor是不要用了,统一用第二个.
	Hoppe94LoopEval() { cloop_evaluator_ = new Stam98ClassicLoopEval; dv_ = OpenMesh::Vec3d(0,0,0); dw_ = OpenMesh::Vec3d(0,0,0); }
	Hoppe94LoopEval(Stam98ClassicLoopEval *_e) { cloop_evaluator_ = _e; dv_ = OpenMesh::Vec3d(0,0,0); dw_ = OpenMesh::Vec3d(0,0,0); }
	~Hoppe94LoopEval() { if (cloop_evaluator_) delete cloop_evaluator_; }

	//在使用下面任何一个eval函数前都必须设置基函数的类型.
	//这两个函数我觉得应该作为private, 不要假设user能设置.
	void set_basis_function_type(Stam98ClassicLoopEval::BasisFuncType _t) { cloop_evaluator_->set_basis_function_type(_t); };
	Stam98ClassicLoopEval::BasisFuncType basis_function_type() { return cloop_evaluator_->basis_func_type(); };

	Stam98ClassicLoopEval *stam98classicloopeval() { return cloop_evaluator_; };//只能用, 不能deleted.

private:
	//regular patch, 1 feature vertex, 0 feature edges.
	bool eval_regular_1v0e_old(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, double _v, double _w, OpenMesh::Vec3d &_p); 
	//regular patch, 1 feature vertex, 0 feature edges.
	bool eval_regular_1v0e(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, double _v, double _w, OpenMesh::Vec3d &_p); 
	//regular patch, 2 feature vertex, 0 feature edges.
	bool eval_regular_2v0e(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, double _v, double _w, OpenMesh::Vec3d &_p); 
	//regular patch, 3 feature vertex, 0 feature edges.
	bool eval_regular_3v0e(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, double _v, double _w, OpenMesh::Vec3d &_p); 
	//regular patch, 2 feature vertex, 1 feature edges. 
	bool eval_regular_2v1e(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, double _v, double _w, OpenMesh::Vec3d &_p); 
	//regular patch, 3 feature vertex, 2 feature edges. 
	bool eval_regular_3v2e(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, double _v, double _w, OpenMesh::Vec3d &_p); 
	//regular patch, 3 feature vertex, 3 feature edges. 
	bool eval_regular_3v3e(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, double _v, double _w, OpenMesh::Vec3d &_p);
	//regular patch, 3 feature vertex, 1 feature edges. 
	bool eval_regular_3v1e(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, double _v, double _w, OpenMesh::Vec3d &_p); 
public:
	// 是上面的统一接口, _type为1到7, _basic_func_type为0/1/2.
	bool eval_regular_1in7(feature_regular_type _type, 
		const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, 
		double _v, double _w, 
		OpenMesh::Vec3d &_p, OpenMesh::Vec3d &_dv, OpenMesh::Vec3d &_dw);

private:
	//irregular patch, 1 feature vertex, 0 feature edges.
	bool eval_irregular_1v0e_up(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p); 
	//irregular patch, 1 feature vertex, 0 feature edges.
	bool eval_irregular_1v0e_left(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p); 
	//irregular patch, 1 feature vertex, 0 feature edges.
	bool eval_irregular_1v0e_right(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p); 
	//irregular patch, 2 feature vertex, 0 feature edges.
	bool eval_irregular_2v0e_left(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p); 
	//irregular patch, 2 feature vertex, 0 feature edges.
	bool eval_irregular_2v0e_right(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p); 
	//irregular patch, 2 feature vertex, 0 feature edges.
	bool eval_irregular_2v0e_bottom(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p); 
	//irregular patch, 3 feature vertex, 0 feature edges.
	bool eval_irregular_3v0e(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p); 
	//irregular patch, 2 feature vertex, 1 feature edges.
	bool eval_irregular_2v1e_left(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w,	OpenMesh::Vec3d &_p); 
	//irregular patch, 2 feature vertex, 1 feature edges.
	bool eval_irregular_2v1e_right(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w,	OpenMesh::Vec3d &_p); 
	//irregular patch, 2 feature vertex, 1 feature edges.
	bool eval_irregular_2v1e_bottom(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w,	OpenMesh::Vec3d &_p); 
	//irregular patch, 3 feature vertex, 1 feature edges.
	bool eval_irregular_3v1e_left(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p); 
	//irregular patch, 3 feature vertex, 1 feature edges.
	bool eval_irregular_3v1e_right(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p); 
	//irregular patch, 3 feature vertex, 1 feature edges.
	bool eval_irregular_3v1e_bottom(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p); 
	//irregular patch, 3 feature vertex, 2 feature edges.
	bool eval_irregular_3v2e_up(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p); 
	//irregular patch, 3 feature vertex, 2 feature edges.
	bool eval_irregular_3v2e_left(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p); 
	//irregular patch, 3 feature vertex, 2 feature edges.
	bool eval_irregular_3v2e_right(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p); 
	//irregular patch, 3 feature vertex, 3 feature edges.
	bool eval_irregular_3v3e(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p); 
public:
	// 这是上面的统一接口, _type是1到17有效.
	bool eval_irregular_1in17(feature_irregular_type _type, 
		const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, 
		OpenMesh::Vec3d &_p, OpenMesh::Vec3d &_dv, OpenMesh::Vec3d &_dw);
private:
	// 核心还是需要Stam98ClassicLoopEval来工作的, 自己维护一个这对象.
	Stam98ClassicLoopEval *cloop_evaluator_; //
	OpenMesh::Vec3d dv_, dw_;//每次invoke Stam98ClassicLoopEval求坐标pos的时候同时求出这两个切向量.
	//它们是值是暂时有效的, 每次调用求值都会立刻被改变。

	//下面两个函数被遗弃了, 只是一开始时候的尝试.
	void eval_regular_1v0e_c33(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15,  int _level, std::vector<OpenMesh::Vec3d> & _c33); 
	int eval_regular_1v0e_tri_idx(double _v, double _w, int _n_stam98_level, double & _new_v, double &_new_w);
	//这个函数是为一个的出口点, 调用Stam98ClassicLoopEval.
	bool eval_regular_1in4_0v0e(const std::vector<OpenMesh::Vec3d> &_c12, double _v, double _w, OpenMesh::Vec3d &_p);
	//下面这7个函数被eval_regular_c18函数代替了.
	void eval_regular_1v0e_c18(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, std::vector<OpenMesh::Vec3d> & _c18); 
	void eval_regular_2v0e_c18(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, std::vector<OpenMesh::Vec3d> & _c18); 
	void eval_regular_3v0e_c18(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, std::vector<OpenMesh::Vec3d> & _c18); 
	void eval_regular_2v1e_c18(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, std::vector<OpenMesh::Vec3d> & _c18); 
	void eval_regular_3v2e_c18(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, std::vector<OpenMesh::Vec3d> & _c18);  
	void eval_regular_3v3e_c18(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, std::vector<OpenMesh::Vec3d> & _c18);  
	void eval_regular_3v1e_c18(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, std::vector<OpenMesh::Vec3d> & _c18);  

	// 从12个顶点获得细分一次之后的18个顶点.
	void eval_regular_c18(const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> & _h15, std::vector<OpenMesh::Vec3d> & _c18);
	void eval_irregular_cN12(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, std::vector<OpenMesh::Vec3d> & _cN12);

	// ---------------------------------------------------
	// 带尖锐特征的细分曲线参数化, private函数,不是用户所需考虑的.
	bool eval_curve_crease_crease(const OpenMesh::Vec3d& P_0_iminus1, const OpenMesh::Vec3d& P_0_i, const OpenMesh::Vec3d& P_0_iplus1, const OpenMesh::Vec3d& P_0_iplus2, 
		double _t, OpenMesh::Vec3d &_p);
	bool eval_curve_corner_crease(const OpenMesh::Vec3d& P_0_i, const OpenMesh::Vec3d& P_0_iplus1, const OpenMesh::Vec3d& P_0_iplus2, 
		double _t, OpenMesh::Vec3d &_p);
	bool eval_curve_corner_corner(const OpenMesh::Vec3d& P_0_i, const OpenMesh::Vec3d& P_0_iplus1, double _t, OpenMesh::Vec3d &_p);
};

#endif // Hoppe94LoopEval_H_

