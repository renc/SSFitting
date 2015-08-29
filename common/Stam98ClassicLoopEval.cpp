
#include <cstdio>
#include <cstdlib>
#include "Stam98ClassicLoopEval.h"
#include "lpdata50NT.h" //我用此数据头文件代替了读lpdata50NT.dat, 后者有时候忘了
// 跟次程序放在一起, 也可以保证这个对象不会重复读那文件.
Stam98ClassicLoopEval::Stam98ClassicLoopEval() : has_read_eval_(false), basis_func_type_(ORI) {
	// 之前是使用度lpdata50NT.dat这文件得到数据, 现在我将那些数据都输出并放在数组double lpdata50NT[]中了,
	// 好处是只建立一个对象保证不用反复读那文件。

	int idx = 0;  
	int Nmax, N, K;
	Nmax = 50; 
	ev = (EVALSTRUCT **) malloc ( (Nmax-2)*sizeof(EVALSTRUCT *) );

	for (int i=0 ; i<Nmax-2 ; i++ ) // i [0, 47]
	{
		N = i+3; // [3, 50]
		K = N+6;

		ev[i] = (EVALSTRUCT *) malloc ( sizeof(EVALSTRUCT) );
		ev[i]->val = (double *) malloc ( K*sizeof(double) );
		ev[i]->vecI = (double *) malloc ( K*K*sizeof(double) );
		ev[i]->Phi = (double **) malloc ( 3*sizeof(double) );
		ev[i]->Phi[0] = (double *) malloc ( K*12*sizeof(double) );
		ev[i]->Phi[1] = (double *) malloc ( K*12*sizeof(double) );
		ev[i]->Phi[2] = (double *) malloc ( K*12*sizeof(double) );

		for (int j = 0; j < K; ++j) {			
			ev[i]->val[j] = lpdata50NT[idx++];
		} 
		for (int j = 0; j < K*K; ++j) {	
			ev[i]->vecI[j]  = lpdata50NT[idx++];
		} 
		for (int j = 0; j < K*12; ++j) {	
			ev[i]->Phi[0][j]  = lpdata50NT[idx++];
		} 
		for (int j = 0; j < K*12; ++j) {	
			ev[i]->Phi[1][j]  = lpdata50NT[idx++];
		} 
		for (int j = 0; j < K*12; ++j) {	
			ev[i]->Phi[2][j]  = lpdata50NT[idx++];
		}
	} 
	has_read_eval_ = true;
}
bool Stam98ClassicLoopEval::evaluate_regular_patch(const std::vector<OpenMesh::Vec3d> &_c, const double _v, const double _w, OpenMesh::Vec3d &_p) {
	if (_c.size() != 12) std::cout << "Error: stam98, evaluation of regular patch need 12 control points.\n";
	if (_v == 0 && _w ==0) {
		// 这里要给出极限点，这里要求出并返回.
		//std::cout << "_v and _w is 0, please use the limit position directly.\n";
		//return false;
		if (basis_func_type_ == ORI) { // _p直接就是极限坐标 
			_p =  _c[3] * 0.5 + (_c[0]+ _c[1]+ _c[2]+ _c[4]+ _c[6]+ _c[7]) * 0.0833;
			return true;
		}  // else if (basis_func_type_ == Der_V or Der_W)可以正常地来求. 
	} else {
		_p = OpenMesh::Vec3d(0, 0, 0);
		std::vector<double> b;
		b = calc_basis_funs(_v, _w);  

		// _p = s(v, w) = cT * b(v, w)
		for (int i = 0; i < 12; ++i) { // _p.x 
			_p[0] += _c[i][0] * b[i];
		}
		for (int i = 0; i < 12; ++i) { // _p.y
			_p[1] += _c[i][1] * b[i];
		}
		for (int i = 0; i < 12; ++i) { // _p.z
			_p[2] += _c[i][2] * b[i];
		} 
	}
	
	return true;
}
bool Stam98ClassicLoopEval::evaluate_regular_patch(const std::vector<OpenMesh::Vec3d> &_c12, 
												   const double _v, const double _w, 
		OpenMesh::Vec3d &_p, OpenMesh::Vec3d &_dv, OpenMesh::Vec3d &_dw)//_c.size() == 12, 注意点的排序.
{
	if (_c12.size() != 12) std::cout << "Error: stam98, evaluation of regular patch need 12 control points.\n";
	// 1. get the evaluation _p.
	if (_v == 0 && _w ==0) {
		_p =  _c12[3] * 0.5 + (_c12[0]+ _c12[1]+ _c12[2]+ _c12[4]+ _c12[6]+ _c12[7]) * 0.0833;
		// basis_func_type_ == Der_V or Der_W可以正常地来求. 
	} else {
		_p = OpenMesh::Vec3d(0, 0, 0);
		set_basis_function_type(ORI);
		std::vector<double> b;
		b = calc_basis_funs(_v, _w);
		for (int i = 0; i < 12; ++i) {
			_p += _c12[i] * b[i];
		}
	}

	// 2. get the dStam/dv
	_dv = OpenMesh::Vec3d(0, 0, 0);
	set_basis_function_type(Der_V);
	std::vector<double> b;
	b = calc_basis_funs(_v, _w);
	for (int i = 0; i < 12; ++i) {
		_dv += _c12[i] * b[i];
	}
	// 3. get the dStam/dw
	_dw = OpenMesh::Vec3d(0, 0, 0);
	set_basis_function_type(Der_W); 
	b = calc_basis_funs(_v, _w);
	for (int i = 0; i < 12; ++i) {
		_dw += _c12[i] * b[i];
	}	
	return true;
}

bool Stam98ClassicLoopEval::evaluate_irregular_patch(const std::vector<OpenMesh::Vec3d> &_c, 
													 const double _v, const double _w, OpenMesh::Vec3d &_p) 
{
	if (has_read_eval_ == false) {
		if (read_eval() == false) return false;
	}

	if (_v == 0 && _w ==0) {
		//std::cout << "_v and _w is 0, please use the limit position directly.\n";
		//return false;
		if (basis_func_type_ == ORI) { // _p直接就是极限坐标 
			int k = _c.size() - 6;// k is the valence of the extraordinary vertex.
			double inv_v = 1.0/double(k); double t = (3.0 + 2.0 * cos( 2.0 * M_PI * inv_v) ); // double beta = inv_v * (40.0 - t * t)/64.0;
			double k_alpha = 1.0 / ( 24.0/(40.0-t*t) + 1); double alpha = inv_v * k_alpha;
			_p =  _c[0] * (1 - k_alpha);
			for ( int i = 1; i <= k; ++i)
				_p += _c[i] * alpha;
			return true;
		} // else if (basis_func_type_ == Der_V or Der_W)可以正常地来求. 
	} else {
		_p = OpenMesh::Vec3d(0, 0, 0);

		int K = _c.size();// = N + 6;
		int N = K - 6;// N is the valence of the extraordinary vertex.
		int M = K + 6; 

		std::vector<OpenMesh::Vec3d> _cp;
		_cp.resize(K);
		project_points(_cp, _c, N);   //projection of the K initial control vertices onto the eigenspace

		eval_surf(_p, _v, _w, _cp, N);
	}
	return true;
}


bool Stam98ClassicLoopEval::evaluate_irregular_patch(const std::vector<OpenMesh::Vec3d> &_c, 
													 const double _v, const double _w, 
													 OpenMesh::Vec3d &_p, OpenMesh::Vec3d &_dv, OpenMesh::Vec3d &_dw) 
{
	double v = _v, w = _w;
	if (fabs(_v) < 1.0e-5) v = 0;//原来是1.0e-6,现在改为1.0e-5, 否则像_bc(-1.37927e-006 0.19194 0.808062)
	if (fabs(_w) < 1.0e-5) w = 0;//就被认为是无效的. 而事实上第一项该为0, 后两项和为1.
 	
	if ( v < 0 || w < 0) { std::cout << "Error: evaluate_irregular_patch: v, w < 0.\n"; return false; }
	if ( v + w > 1) { 
		if (v+w < 1.00001) {
			if (v < w) w -= 0.00001; else v -= 0.00001;
		} else {
			std::cout << "Error: evaluate_irregular_patch: " << v << " + " << w << " > 1.\n"; return false; 
		}
	}
	if (has_read_eval_ == false) {
		if (read_eval() == false) return std::cout << "Error: has_read_eval_ false.\n"; false;
	}

	_p = OpenMesh::Vec3d(0, 0, 0);
	_dv = OpenMesh::Vec3d(0, 0, 0);
	_dw = OpenMesh::Vec3d(0, 0, 0);

	int K = _c.size();// = N + 6;
	int N = K - 6;// N is the valence of the extraordinary vertex.
	int M = K + 6; 

	std::vector<OpenMesh::Vec3d> _cp;
	_cp.resize(K);
	project_points(_cp, _c, N);   //projection of the K initial control vertices onto the eigenspace

	if (v == 0 && w ==0) { 
		int k = _c.size() - 6;// k is the valence of the extraordinary vertex.
		double inv_v = 1.0/double(k); double t = (3.0 + 2.0 * cos( 2.0 * M_PI * inv_v) ); // double beta = inv_v * (40.0 - t * t)/64.0;
		double k_alpha = 1.0 / ( 24.0/(40.0-t*t) + 1); double alpha = inv_v * k_alpha;
		_p =  _c[0] * (1 - k_alpha);
		for ( int i = 1; i <= k; ++i) 
			_p += _c[i] * alpha; 
		// (basis_func_type_ == Der_V or Der_W)可以正常地来求. 
	} else {
		set_basis_function_type(ORI);
		eval_surf(_p, v, w, _cp, N);
		
	} 
	
	set_basis_function_type(Der_V);
	eval_surf(_dv, v, w, _cp, N);
	set_basis_function_type(Der_W);
	eval_surf(_dw, v, w, _cp, N);
	
	return true;
}

void Stam98ClassicLoopEval::project_points(std::vector<OpenMesh::Vec3d> &_cp, const std::vector<OpenMesh::Vec3d> &_c, const int _N) {
	// 这个函数的原型在stam98's paper中.
	if (_cp.size() != _N + 6) std::cout << "Error: project_points().\n";
	int K = _N + 6;
	for (int i = 0; i < K; ++i) {
		_cp[i] = OpenMesh::Vec3d(0, 0, 0);
		for (int j = 0; j < K; ++j) {
			_cp[i] += ev[_N-3]->vecI[IX(i,j,K)] * _c[j];//从这里看好像文件lpdata50Nt.dat里面是正序存的.
		}
	}
}
void Stam98ClassicLoopEval::eval_surf(OpenMesh::Vec3d& _p, double _v, double _w, std::vector<OpenMesh::Vec3d> &_cp, const int _N) {
	// 这个函数的原型在stam98's paper中. 
	// determine in which domain omega(n, k) the parameter lies
	int n = (int)floor(1 - log(_v + _w) / log(2.0));
	double pow2 = pow(2.0, n-1);
	_v *= pow2; _w *= pow2;
	int k = 0; //k = 0, 1, 2表示irregular parth一分为四后的三个小的regular patch
	if (_v > 0.5) {
		k = 0; _v = 2 * _v - 1; _w = 2 * _w;
	} else if ( _w > 0.5) {
		k = 2; _v = 2 *_v; _w = 2 * _w -1;
	} else {
		k = 1; _v = 1 - 2 * _v; _w = 1 - 2 * _w;
	}
	// now evaluate the surface 
	_p = OpenMesh::Vec3d(0, 0, 0);
	for (int i = 0; i < _N + 6; ++i) {
		_p += pow(ev[_N-3]->val[i], n-1) * eval_basis(_N, k, i, _v, _w) * _cp[i];
	}
}
double Stam98ClassicLoopEval::eval_basis(const int _N, const int _k, const int _i, const double _v, const double _w) {
	// _N是奇异点valence, _k取012, 
	std::vector<double> b; //12 basis functions
	b = calc_basis_funs(_v, _w); 

	double out = 0;
	for (int j = 0; j < 12; ++j) {
		out += ev[_N - 3]->Phi[_k][IX(_i, j, _N+6)] * b[j]; 
	}
	return out;
}
std::vector<double> Stam98ClassicLoopEval::calc_basis_funs(double _v, double _w) {
	// 严重注意: 见原文的附录A, 这12个基函数是和Fig1的点序号相对应的, 否则就错了.
	std::vector<double> b;//12*1
	b.resize(12);
	
	/*// 计算方法一, 里面的公式和stam98的附录A中的一摸一样.
	double u = 1 - _v - _w, v = _v, w = _w;
    double u2 = u*u, u3 = u2*u, u4 = u3*u;
	double v2 = v*v, v3 = v2*v, v4 = v3*v;
	double w2 = w*w, w3 = w2*w, w4 = w3*w;

	b[0] = u4 + 2*u3*v;
	b[1] = u4 + 2*u3*w;
	b[2] = u4 + 2*u3*w + 6*u3*v + 6*u2*v*w + 12*u2*v2 + 6*u*v2*w + 6*u*v3 + 2*v3*w +v4;
	b[3] = 6*u4 + 24*u3*w + 24*u2*w2 + 8*u*w3 + w4 + 24*u3*v + 60*u2*v*w + 36*u*v*w2 + 6*v*w3 + 24*u2*v2 + 36*u*v2*w + 12*v2*w2 + 8*u*v3 + 6*v3*w + v4;
	b[4] = u4 + 6*u3*w + 12*u2*w2 + 6*u*w3 + w4 + 2*u3*v + 6*u2*v*w + 6*u*v*w2 + 2*v*w3;
	b[5] = 2*u*v3 + v4;
	b[6] = u4 + 6*u3*w + 12*u2*w2 + 6*u*w3 + w4 + 8*u3*v + 36*u2*v*w + 36*u*v*w2 + 8*v*w3 + 24*u2*v2 + 60*u*v2*w + 24*v2*w2 + 24*u*v3 + 24*v3*w + 6*v4;
	b[7] = u4 + 8*u3*w + 24*u2*w2 + 24*u*w3 + 6*w4 + 6*u3*v + 36*u2*v*w + 60*u*v*w2 + 24*v*w3 + 12*u2*v2 + 36*u*v2*w + 24*v2*w2 + 6*u*v3 + 8*v3*w + v4;
	b[8] = 2*u*w3 + w4;
	b[9] = 2*v3*w + v4;
	b[10] = 2*u*w3 + w4 + 6*u*v*w2 + 6*v*w3 + 6*u*v2*w + 12*v2*w2 + 2*u*v3 + 6*v3*w + v4;
	b[11] = w4 + 2*v*w3;*/

	// 计算方法二, 只是将上面的式子化成矩阵形式
	std::vector<double> Mvw;//15*1
	Mvw.resize(15);
	if (basis_func_type_ == ORI) { //用于求stam(v, w)
		Mvw[0] = 1; 
		Mvw[1] = _v;			Mvw[2] = _w; 
		Mvw[3] = _v*_v;			Mvw[4] = _v*_w;			Mvw[5] = _w*_w; 
		Mvw[6] = _v*_v*_v;		Mvw[7] = _v*_v*_w;		Mvw[8] = _v*_w*_w;		Mvw[9] = _w*_w*_w; 
		Mvw[10] = _v*_v*_v*_v;	Mvw[11] = _v*_v*_v*_w;	Mvw[12] = _v*_v*_w*_w;	Mvw[13] = _v*_w*_w*_w;	Mvw[14] = _w*_w*_w*_w; 
	} else if (basis_func_type_ == Der_V) { // 用于求dStam(v,w)/dv
		Mvw[0] = 0; 
		Mvw[1] = 1;				Mvw[2] = 0; 
		Mvw[3] = 2*_v;			Mvw[4] = _w;			Mvw[5] = 0; 
		Mvw[6] = 3*_v*_v;		Mvw[7] = 2*_v*_w;		Mvw[8] = _w*_w;			Mvw[9] = 0; 
		Mvw[10] = 4*_v*_v*_v;	Mvw[11] = 3*_v*_v*_w;	Mvw[12] = 2*_v*_w*_w;	Mvw[13] = _w*_w*_w;		Mvw[14] = 0;  
	} else if (basis_func_type_ == Der_W) { // 用于求dStam(v,w)/dw
		Mvw[0] = 0; 
		Mvw[1] = 0;				Mvw[2] = 1; 
		Mvw[3] = 0;				Mvw[4] = _v;			Mvw[5] = 2*_w; 
		Mvw[6] = 0;				Mvw[7] = _v*_v;			Mvw[8] = _v*_w*2;		Mvw[9] = 3*_w*_w; 
		Mvw[10] = 0;			Mvw[11] = _v*_v*_v;		Mvw[12] = _v*_v*_w*2;	Mvw[13] = _v*_w*_w*3;	Mvw[14] = 4*_w*_w*_w;  
	}
	
	int Mb[12][15] = {
		1, -2, -4, 0, 6,		6, 2, 0, -6, -4,		-1, -2, 0, 2, 1,
		1, -4, -2, 6, 6,		0, -4, -6, 0, 2,		1, 2, 0, -2, -1,
		1, 2, -2, 0, -6,		0, -4, 0, 6, 2,			2, 4, 0, -2, -1,
		6, 0, 0, -12, -12,		-12, 8, 12, 12, 8,		-1, -2, 0, -2, -1,
		1, -2, 2, 0, -6,		0, 2, 6, 0, -4,			-1, -2, 0, 4, 2,
		0, 0, 0, 0, 0,			0, 2, 0, 0, 0,			-1, -2, 0, 0, 0, 
		1, 4, 2, 6, 6,			0, -4, -6, -12, -4,		-1, -2, 0, 4, 2,
		1, 2, 4, 0, 6,			6, -4, -12, -6, -4,		2, 4, 0, -2, -1,
		0, 0, 0, 0, 0,			0, 0, 0, 0, 2,			0, 0, 0, -2, -1,
		0, 0, 0, 0, 0,			0, 0, 0, 0, 0,			1, 2, 0, 0, 0,
		0, 0, 0, 0, 0,			0, 2, 6, 6, 2,			-1, -2, 0, -2, -1,
		0, 0, 0, 0, 0,			0, 0, 0, 0, 0,			0, 0, 0, 2, 1 
	};
	for (int i = 0; i < 12; ++i) {
		for (int j = 0; j < 15; ++j)
			b[i] += Mb[i][j] * Mvw[j]; 
	}

	// --- 上面无论用计算方法一还是二, 这里都需要乘以1/12的.
	double inv_12 = 1.0 / 12.0;
	for (int i = 0; i < 12; ++i) 
		b[i] *= inv_12;

	return b;
} 
bool Stam98ClassicLoopEval::read_eval() {
	ev = read_eval ( &Nmax );
	if (ev == NULL) {
		std::cout << "Error: real_eval return NULL.\n";
		return false;
	}
	//std::cout << "Nmax: " << Nmax << ".\n"; // for test// 50
	has_read_eval_ = true;
	return true;
}
EVALSTRUCT** Stam98ClassicLoopEval::read_eval( int * pNmax )
{
	EVALSTRUCT ** ev;
	FILE * f;
	int Nmax, i, N, K;
 
	if ( !(f=fopen("lpdata50NT.dat", "rb")) ) { //这个文档请放在这个应用程序的同目录下.
		if ( !(f=fopen("D:\\libraries\\MeshCourse07-Code\\MSVC\\release\\lpdata50NT.dat", "rb")) ) { 
				std::cout << "Error: read eval() fail to open file 'lpdata50NT.dat'.\n";	
			return ( NULL ); 
		} else std::cout << "Info: open file 'lpdata50NT.dat' successfully.\n";
	} else {
		std::cout << "Info: open file 'lpdata50NT.dat' successfully.\n";
	}

	fread ( &Nmax, sizeof(int), 1, f );//读进最大的个数, 好像就是50

	ev = (EVALSTRUCT **) malloc ( (Nmax-2)*sizeof(EVALSTRUCT *) );

	for ( i=0 ; i<Nmax-2 ; i++ ) // i [0, 47]
	{
		N = i+3; // [3, 50]
		K = N+6;

		ev[i] = (EVALSTRUCT *) malloc ( sizeof(EVALSTRUCT) );
		ev[i]->val = (double *) malloc ( K*sizeof(double) );
		ev[i]->vecI = (double *) malloc ( K*K*sizeof(double) );
		ev[i]->Phi = (double **) malloc ( 3*sizeof(double) );
		ev[i]->Phi[0] = (double *) malloc ( K*12*sizeof(double) );
		ev[i]->Phi[1] = (double *) malloc ( K*12*sizeof(double) );
		ev[i]->Phi[2] = (double *) malloc ( K*12*sizeof(double) );

		fread ( ev[i]->val, sizeof(double), K, f );
		fread ( ev[i]->vecI, sizeof(double), K*K, f );
		fread ( ev[i]->Phi[0], sizeof(double), K*12, f );
		fread ( ev[i]->Phi[1], sizeof(double), K*12, f );
		fread ( ev[i]->Phi[2], sizeof(double), K*12, f );
	}

	fclose ( f );

	*pNmax = Nmax;

	return ( ev );
}


// -------------------
