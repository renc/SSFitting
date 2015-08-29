// ***************************************************************
//  
//  Stam98ClassicLoopEval 
//  Copyright (C) 2008 - by rencanjiang All Rights Reserved
//  rencanjiang@163.com	
//	-------------------------------------------------------------
//  
//  version: 1.0
//  data: 04/12/2008
// ***************************************************************
//  
//  File description: 
//		Check stam98's paper "Evaluation of Loop Subdivision Surfaces" for details.
//      但是最好在LoopEvaluatorT.h中通过统一的接口来调用，可以省去很多麻烦.
// ***************************************************************

#ifndef  stam98classicloop_h
#define  stam98classicloop_h

#include <iostream>
#include <vector>
#include "OpenMeshAll.h"

typedef struct {
	double * val; //one dimension array.
	double * vecI;//one dimension array.
	double ** Phi;//two dimension array.
} EVALSTRUCT;

class Stam98ClassicLoopEval {
public:
	//原来的12个basis functions只要稍微改改就可以求他们对v和w的导数了.
	//Usage:直接是Stam98ClassicLoopEval::ORI, 而不需要Stam98ClassicLoopEval::BasisFuncType::ORI.
	static enum BasisFuncType { ORI = 0, Der_V = 1, Der_W = 2 };
public:
	Stam98ClassicLoopEval();
	~Stam98ClassicLoopEval() {
		if (ev != NULL) {
			for (int i = 0; i < Nmax; ++i) {
				//free ev[i];
			}
		}
		ev = NULL;
	}
	void set_basis_function_type(BasisFuncType _t) {
		basis_func_type_ = _t;
	}
	BasisFuncType basis_func_type() { return basis_func_type_; }

	// 函数evaluate_regular_patch应该是函数evaluate_irregular_patch的特例, 只是这里各自的_c控制点排序不一致而已.
	// 也就是说如果控制点统一用irregular patch的排序方式, 那么可以都用evaluate_irregular_patch就可以了.
	// u -------> w
	// |
	// |
	// \/ v
	// 1. 分别单独求pos, der_v, der_w的话, 需要先设定set_basis_function_type(...). 
	//    这个对useer的假设应该去除, basis_func_type()这个函数应该作为private.
	// 2. 一次来求pos, der_v, 和der_w的话, 就不需要设定basis function的类型了. 建议使用这个函数, 代价多很小, 但是代码更安全.
	bool evaluate_regular_patch(const std::vector<OpenMesh::Vec3d> &_c12, const double _v, const double _w, OpenMesh::Vec3d &_p);//_c.size() == 12, 注意点的排序.
	bool evaluate_regular_patch(const std::vector<OpenMesh::Vec3d> &_c12, const double _v, const double _w, OpenMesh::Vec3d &_p, OpenMesh::Vec3d &_dv, OpenMesh::Vec3d &_dw);//_c.size() == 12, 注意点的排序.
	bool evaluate_irregular_patch(const std::vector<OpenMesh::Vec3d> &_cN6, const double _v, const double _w, OpenMesh::Vec3d &_p);//_c.size() == N+6, 注意点的排序
	bool evaluate_irregular_patch(const std::vector<OpenMesh::Vec3d> &_cN6, const double _v, const double _w, OpenMesh::Vec3d &_p, OpenMesh::Vec3d &_dv, OpenMesh::Vec3d &_dw);//_c.size() == N+6, 注意点的排序
	// 当basis function type是ORI时候求出来的点_p是在(v, w)位置上的Stam(v, w)值,
	// 当basis function type是Dev_V时候求出来的点_p是在(v, w)位置上的Stam(v, w)对v求导.[dStam(v, w)1/dv, dStam(v, w)2/dv, dStam(v, w)3/dv],
	// 当basis function type是Dev_W时候求出来的点_p是在(v, w)位置上的Stam(v, w)对w求导.[dStam(v, w)1/dw, dStam(v, w)2/dw, dStam(v, w)3/dw]
	// 所以一定要清楚准确设置基函数类型.

private:
	//bool evaluate_irregular_patch(const std::vector<OpenMesh::Vec3d> &_c, const double _v, const double _w, OpenMesh::Vec3d &_p, Stam98ClassicLoopEval::BasisFuncType _t);//_c.size() == N+6, 注意点的排序

	bool read_eval(); //使用这个类的对象之前请先调用这个函数.

private:
	
	bool has_read_eval_;
	Stam98ClassicLoopEval::BasisFuncType basis_func_type_;

	int IX(int i, int j, int n) { return ((i)+(n)*(j)); }//第i列, 第j行
	EVALSTRUCT ** ev;
	int Nmax;//lpdata50Nt.dat中支持的最大的顶点的valence
	//ev数组能支持Nmax的阶, 但其实ev数组的size=Nmax-2, index=[0,1,2...Nmax-3], 分别对应的valence=[3,4,5...Nmax].
	EVALSTRUCT ** read_eval ( int * pNmax ); 
	void project_points(std::vector<OpenMesh::Vec3d> &_cp, const std::vector<OpenMesh::Vec3d> &_c, const int _N);
	void eval_surf(OpenMesh::Vec3d& _p, double _v, double _w, std::vector<OpenMesh::Vec3d> &_cp, const int _N);
	double eval_basis(const int _N, const int _k, const int _i, const double _v, const double _w);
	std::vector<double> calc_basis_funs(double _v, double _w);
	
public:
	
private:
	
};

#endif // stam98classicloop_h