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
//      ���������LoopEvaluatorT.h��ͨ��ͳһ�Ľӿ������ã�����ʡȥ�ܶ��鷳.
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
	//ԭ����12��basis functionsֻҪ��΢�ĸľͿ��������Ƕ�v��w�ĵ�����.
	//Usage:ֱ����Stam98ClassicLoopEval::ORI, ������ҪStam98ClassicLoopEval::BasisFuncType::ORI.
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

	// ����evaluate_regular_patchӦ���Ǻ���evaluate_irregular_patch������, ֻ��������Ե�_c���Ƶ�����һ�¶���.
	// Ҳ����˵������Ƶ�ͳһ��irregular patch������ʽ, ��ô���Զ���evaluate_irregular_patch�Ϳ�����.
	// u -------> w
	// |
	// |
	// \/ v
	// 1. �ֱ𵥶���pos, der_v, der_w�Ļ�, ��Ҫ���趨set_basis_function_type(...). 
	//    �����useer�ļ���Ӧ��ȥ��, basis_func_type()�������Ӧ����Ϊprivate.
	// 2. һ������pos, der_v, ��der_w�Ļ�, �Ͳ���Ҫ�趨basis function��������. ����ʹ���������, ���۶��С, ���Ǵ������ȫ.
	bool evaluate_regular_patch(const std::vector<OpenMesh::Vec3d> &_c12, const double _v, const double _w, OpenMesh::Vec3d &_p);//_c.size() == 12, ע��������.
	bool evaluate_regular_patch(const std::vector<OpenMesh::Vec3d> &_c12, const double _v, const double _w, OpenMesh::Vec3d &_p, OpenMesh::Vec3d &_dv, OpenMesh::Vec3d &_dw);//_c.size() == 12, ע��������.
	bool evaluate_irregular_patch(const std::vector<OpenMesh::Vec3d> &_cN6, const double _v, const double _w, OpenMesh::Vec3d &_p);//_c.size() == N+6, ע��������
	bool evaluate_irregular_patch(const std::vector<OpenMesh::Vec3d> &_cN6, const double _v, const double _w, OpenMesh::Vec3d &_p, OpenMesh::Vec3d &_dv, OpenMesh::Vec3d &_dw);//_c.size() == N+6, ע��������
	// ��basis function type��ORIʱ��������ĵ�_p����(v, w)λ���ϵ�Stam(v, w)ֵ,
	// ��basis function type��Dev_Vʱ��������ĵ�_p����(v, w)λ���ϵ�Stam(v, w)��v��.[dStam(v, w)1/dv, dStam(v, w)2/dv, dStam(v, w)3/dv],
	// ��basis function type��Dev_Wʱ��������ĵ�_p����(v, w)λ���ϵ�Stam(v, w)��w��.[dStam(v, w)1/dw, dStam(v, w)2/dw, dStam(v, w)3/dw]
	// ����һ��Ҫ���׼ȷ���û���������.

private:
	//bool evaluate_irregular_patch(const std::vector<OpenMesh::Vec3d> &_c, const double _v, const double _w, OpenMesh::Vec3d &_p, Stam98ClassicLoopEval::BasisFuncType _t);//_c.size() == N+6, ע��������

	bool read_eval(); //ʹ�������Ķ���֮ǰ���ȵ����������.

private:
	
	bool has_read_eval_;
	Stam98ClassicLoopEval::BasisFuncType basis_func_type_;

	int IX(int i, int j, int n) { return ((i)+(n)*(j)); }//��i��, ��j��
	EVALSTRUCT ** ev;
	int Nmax;//lpdata50Nt.dat��֧�ֵ����Ķ����valence
	//ev������֧��Nmax�Ľ�, ����ʵev�����size=Nmax-2, index=[0,1,2...Nmax-3], �ֱ��Ӧ��valence=[3,4,5...Nmax].
	EVALSTRUCT ** read_eval ( int * pNmax ); 
	void project_points(std::vector<OpenMesh::Vec3d> &_cp, const std::vector<OpenMesh::Vec3d> &_c, const int _N);
	void eval_surf(OpenMesh::Vec3d& _p, double _v, double _w, std::vector<OpenMesh::Vec3d> &_cp, const int _N);
	double eval_basis(const int _N, const int _k, const int _i, const double _v, const double _w);
	std::vector<double> calc_basis_funs(double _v, double _w);
	
public:
	
private:
	
};

#endif // stam98classicloop_h