
#include "Hoppe94LoopEval.h"



// ================================================
void Hoppe94LoopEval::eval_regular_1v0e_c33(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, int _n_stam98_level, std::vector<OpenMesh::Vec3d> & _c33) {
	//std::cout << "Hoppe94LoopEval::eval_regular_1v0e_c33() begin.\n";
	double d1_16 = 0.0625, d1_8 = 0.125, d3_8 = 0.375, d5_8 = 0.625;
	// first time subdiv, 12 -> 18 
	double A1[18][12] = {
		//0     1       2       3       4       5       6       7       8       9      10       11   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 0  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 1  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 2
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 3   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 4
		{ 0,	0,	d3_8,	d1_8,		0,	d1_8,	d3_8,		0,		0,		0,		0,		0 }, // index 5   
		{ 0,	0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8,		0,		0,		0,		0 }, // index 6
		{ 0,	0,		0,	d3_8,	d1_8,		0,	d1_8,	d3_8,		0,		0,		0,		0 }, // index 7  
		{ 0,	0,		0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8,		0,		0,		0 }, // index 8
		{ 0,	0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8,		0,		0 }, // index 9   
		{ 0,	0,	d1_16,	d1_16,		0,	d1_16,	d5_8,	d1_16,		0,	d1_16,	d1_16,		0 }, // index 10
		{ 0,	0,		0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8,		0 }, // index 11  
		{ 0,	0,		0,	d1_16,	d1_16,		0,	d1_16,	d5_8,	d1_16,		0,	d1_16,	d1_16 }, // index 12
		{ 0,	0,		0,		0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8 }, // index 13  
		{ 0,	0,		0,		0,		0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8,		0 }, // index 14
		{ 0,	0,		0,		0,		0,		0,	d3_8,	d1_8,		0,	d1_8,	d3_8,		0 }, // index 15  
		{ 0,	0,		0,		0,		0,		0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8 }, // index 16 
		{ 0,	0,		0,		0,		0,		0,		0,	d3_8,	d1_8,		0,	d1_8,	d3_8 }, // index 17  
	};
	int n_count_feature_edges = 0;
	std::vector<int> fv_tmp;
	if (_h[0]) { // new control point 0.
		A1[0][0] = A1[0][3] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(0);
	} else {
		A1[0][0] = A1[0][3] = d3_8;  A1[0][1] = A1[0][2] = d1_8;  
	}
	if (_h[1]) { // new control point 1.
		A1[1][1] = A1[1][3] = 0.5;   ++n_count_feature_edges; fv_tmp.push_back(1);
	} else {
		A1[1][1] = A1[1][3] = d3_8;  A1[1][0] = A1[1][4] = d1_8;
	}
	if (_h[2]) { // new control point 2.
		A1[2][2] = A1[2][3] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(2);
	} else {
		A1[2][2] = A1[2][3] = d3_8;  A1[2][0] = A1[2][6] = d1_8;
	}
	if (_h[3]) { // new control point 4.
		A1[4][3] = A1[4][4] = 0.5;   ++n_count_feature_edges; fv_tmp.push_back(4);
	} else {
		A1[4][3] = A1[4][4] = d3_8;  A1[4][1] = A1[4][7] = d1_8;  
	}
	if (n_count_feature_edges <= 1) { // smooth or dart // new control point 3.
		A1[3][3] = d5_8;  A1[3][0] = A1[3][1] = A1[3][2] = A1[3][4] = A1[3][6] = A1[3][7] = d1_16;
	} else if (n_count_feature_edges == 2) { // crease
		A1[3][3] = 0.75;  A1[3][fv_tmp[0]] = A1[3][fv_tmp[1]] = d1_8; 
	} else { // corner
		A1[3][3] = 1;
	} //这个矩阵我真的检查许多遍了.
	// second time subdiv 18 -> 33
	double A2[33][18] = { 0 };
	for (int i = 0; i < 9; ++i) { //前9行和A1是一样的.
		for (int j = 0; j < 12; ++j) {
			A2[i][j] = A1[i][j];
		}
	}	 
	A2[9][5] = A2[9][6] = d3_8; A2[9][2] = A2[9][10] = d1_8; 
	A2[10][6] = d5_8; A2[10][2] = A2[10][3] = A2[10][5] = A2[10][7] = A2[10][10] = A2[10][11] = d1_16;
	A2[11][6] = A2[11][7] = d3_8; A2[11][3] = A2[11][11] = d1_8;
	A2[12][7] = d5_8; A2[12][3] = A2[12][4] = A2[12][6] = A2[12][8] = A2[12][11] = A2[12][12] = d1_16;
	A2[13][7] = A2[13][8] = d3_8; A2[13][4] = A2[13][12] = d1_8;
	A2[14][5] = A2[14][10] = d3_8; A2[14][6] = A2[14][9] = d1_8;
	A2[15][6] = A2[15][10] = d3_8; A2[15][5] = A2[15][11] = d1_8;
	A2[16][6] = A2[16][11] = d3_8; A2[16][7] = A2[16][10] = d1_8;
	A2[17][7] = A2[17][11] = d3_8; A2[17][6] = A2[17][12] = d1_8;
	A2[18][7] = A2[18][12] = d3_8; A2[18][8] = A2[18][11] = d1_8;
	A2[19][8] = A2[19][12] = d3_8; A2[19][7] = A2[19][13] = d1_8;
	A2[20][9] = A2[20][10] = d3_8; A2[20][5] = A2[20][14] = d1_8;
	A2[21][10] = d5_8; A2[21][5] = A2[21][6] = A2[21][9] = A2[21][11] = A2[21][14] = A2[21][15] = d1_16;
	A2[22][10] = A2[22][11] = d3_8; A2[22][6] = A2[22][15] = d1_8;
	A2[23][11] = d5_8; A2[23][6] = A2[23][7] = A2[23][10] = A2[23][12] = A2[23][15] = A2[23][16] = d1_16;
	A2[24][11] = A2[24][12] = d3_8; A2[24][7] = A2[24][16] = d1_8;
	A2[25][12] = d5_8; A2[25][7] = A2[25][8] = A2[25][11] = A2[25][13] = A2[25][16] = A2[25][17] = d1_16;
	A2[26][12] = A2[26][13] = d3_8; A2[26][8] = A2[26][17] = d1_8;
	A2[27][10] = A2[27][14] = d3_8; A2[27][9] = A2[27][15] = d1_8;
	A2[28][10] = A2[28][15] = d3_8; A2[28][11] = A2[28][14] = d1_8;
	A2[29][11] = A2[29][15] = d3_8; A2[29][10] = A2[29][16] = d1_8;
	A2[30][11] = A2[30][16] = d3_8; A2[30][12] = A2[30][15] = d1_8;
	A2[31][12] = A2[31][16] = d3_8; A2[31][11] = A2[31][17] = d1_8;
	A2[32][12] = A2[32][17] = d3_8; A2[32][13] = A2[32][16] = d1_8;//这个矩阵也检查多次了.
	//A33_12 = A33_18 * A18_12 
	double A3[33][12] = { 0};
	for (int i = 0; i < 33; ++i) {
		for (int j = 0; j < 12; ++j) {
			for (int k = 0; k < 18; ++k) {
				A3[i][j] += A2[i][k] * A1[k][j];
			}
		}
	}
	/*
	if (_n_stam98_level == 1) {

		// _c -> c33, from 12 control points to 33 control points.
		for (int i = 0; i < 33; ++i) {
			for (int j = 0; j < 12; ++j) {
				_c33[i] +=  A3[i][j] * _c[j];
			}
		} 
	} else {

	}*/
		// when _n_stam98_level > 1, we need to pick 18 points out of the 33 points, 
		// and use the new 18 points to create new 33 points.
		double A18_33[18][33] = { 0 };
		A18_33[0][0] = 1; A18_33[1][1] = 1; 
		A18_33[2][2] = 1; A18_33[3][3] = 1; A18_33[4][4] = 1; 
		A18_33[5][5] = 1; A18_33[6][6] = 1; A18_33[7][7] = 1;  A18_33[8][8] = 1; 
		A18_33[9][9] = 1; A18_33[10][10] = 1; A18_33[11][11] = 1;  A18_33[12][12] = 1; A18_33[13][13] = 1; 
		A18_33[14][15] = 1; A18_33[15][16] = 1;  A18_33[16][17] = 1; A18_33[17][18] = 1; 

		double A33_33[33][33] = {0}; // A33_18 * A18_33
		for (int i = 0; i < 33; ++i) {
			for (int j = 0; j < 33; ++j) {
				for (int k = 0; k < 18; ++k) {
					A33_33[i][j] += A2[i][k] * A18_33[k][j];
				}
			}
		}
		 
		double AI[33][33] = { 0 };//when _n_stam98_level == 1, this is just the unit matrix.
		for (int i = 0; i < 33; ++i) {
			for (int j = 0; j < 33; ++j) {
				if (i == j) AI[i][j] = 1;
			}
		}

		for (int i_level = 0; i_level < _n_stam98_level -1; ++ i_level) { // when _n_stam98_level=1, not go into for.
			double Atmp[33][33] = { 0};
			for (int i = 0; i < 33; ++i) {
				for (int j = 0; j < 33; ++j) {
					for (int k = 0; k < 33; ++k) {
						Atmp[i][j] += A33_33[k][j] * AI[i][k];
					} 
				}
			}
			for (int i = 0; i < 33; ++i) {
				for (int j = 0; j < 33; ++j) {
					AI[i][j] = Atmp[i][j];
				}
			}
		}

		double A4[33][12] = {0 };
		for (int i = 0; i < 33; ++i) {
			for (int j = 0; j < 12; ++j) {
				for (int k = 0; k < 33; ++k) {
					A4[i][j] += AI[i][k] * A3[k][j];
				}
			}
		}
		// _c -> c33, from 12 control points to 33 control points.
		for (int i = 0; i < 33; ++i) {
			for (int j = 0; j < 12; ++j) {
				_c33[i] +=  A4[i][j] * _c[j];
			}
		} 
	
}
int Hoppe94LoopEval::eval_regular_1v0e_tri_idx(double _v, double _w, int _n_stam98_level, double & _new_v, double &_new_w) {
	double pow2 = pow((double)2, _n_stam98_level -1);
	_v *= pow2; _w *= pow2; // 将_v, _w放大到单位为1.

	int tri_idx = -1; 
	// Method 1. 0.5 < (_v + _w) <=1
	if (_v > 0.5) {
		if (_v + _w <= 0.75) {
			tri_idx = 0;
		} else {
			if (_v > 0.75) tri_idx = 5;
			else if (_w > 0.25) tri_idx = 7;
			else tri_idx = 6;
		}
	} else if (_w > 0.5) {
		if (_v + _w <= 0.75) {
			tri_idx = 4;
		} else {
			if (_v > 0.25) tri_idx = 9;
			else if (_w > 0.75) tri_idx = 11;
			else tri_idx = 10;
		}
	} else {
		if (_v + _w > 0.75) tri_idx = 8;
		else {
			if (_v <= 0.25) tri_idx = 3;
			else if (_w <= 0.25) tri_idx = 1;
			else tri_idx = 2;
		}
	}
	if (tri_idx < 0 || tri_idx > 11) std::cout << "Error: eval_regular_1v0e_tri_idx(), tri_idx.\n";
	
	double delta_v = 0, delta_w =0; 
	switch(tri_idx)
	{
	case 0: delta_v = _v - 0.5; delta_w = _w - 0; break;
	case 1: delta_v = 0.5 - _v; delta_w = 0.25 - _w; break;
	case 2: delta_v = _v - 0.25; delta_w = _w - 0.25; break;
	case 3: delta_v = 0.25 - _v; delta_w = 0.5 -_w; break;
	case 4: delta_v = _v - 0; delta_w = _w - 0.5; break;

	case 5: delta_v = _v - 0.75; delta_w = _w - 0; break;
	case 6: delta_v = 0.75 - _v; delta_w = 0.25 - _w; break;
	case 7: delta_v = _v - 0.5; delta_w = _w - 0.25; break;
	case 8: delta_v = 0.5 -_v; delta_w = 0.5 - _w; break;
	case 9: delta_v = _v - 0.25; delta_w = _w - 0.5; break;
	case 10: delta_v = 0.25 -_v; delta_w = 0.75 -_w; break;
	case 11: delta_v = _v - 0; delta_w = _w - 0.75; break;
	default: std::cout <<"switch error.\n"; break;
	};
	_new_v = delta_v * 4; //经历了2层细分, 这里归一放大2^2倍. 
	_new_w = delta_w * 4; 
	std::cout << "Method 1: tri_idx: " << tri_idx << ". (new_v, new_w)= (" << _new_v << ", " << _new_w << ").\n";

	// Method 2. calculate the barycentric coordiantes to determine
	tri_idx = -1;
	_new_v = _new_w = 0;
	OpenMesh::Vec2d p10(0.5, 0), p11(0.25, 0.25), p12(0, 0.5),
		p15(0.75, 0), p16(0.5, 0.25), p17(0.25, 0.5), p18(0, 0.75),
		p21(1, 0), p22(0.75, 0.25), p23(0.5, 0.5), p24(0.25, 0.75), p25(0, 1);
	OpenMesh::Vec2d pp(_v, _w); OpenMesh::Vec3d bc;
	if (DGP::is_valid_barycentric_coordinate(DGP::calc_barycentric_coordinates(p10, p15, p16, pp))) { tri_idx = 0; bc = DGP::calc_barycentric_coordinates(p10, p15, p16, pp); }
	if (DGP::is_valid_barycentric_coordinate(DGP::calc_barycentric_coordinates(p16, p11, p10, pp))) { tri_idx = 1; bc = DGP::calc_barycentric_coordinates(p16, p11, p10, pp); }
	if (DGP::is_valid_barycentric_coordinate(DGP::calc_barycentric_coordinates(p11, p16, p17, pp))) { tri_idx = 2; bc = DGP::calc_barycentric_coordinates(p11, p16, p17, pp); }
	if (DGP::is_valid_barycentric_coordinate(DGP::calc_barycentric_coordinates(p17, p12, p11, pp))) { tri_idx = 3; bc = DGP::calc_barycentric_coordinates(p17, p12, p11, pp); }
	if (DGP::is_valid_barycentric_coordinate(DGP::calc_barycentric_coordinates(p12, p17, p18, pp))) { tri_idx = 4; bc = DGP::calc_barycentric_coordinates(p12, p17, p18, pp); }
	if (DGP::is_valid_barycentric_coordinate(DGP::calc_barycentric_coordinates(p15, p21, p22, pp))) { tri_idx = 5; bc = DGP::calc_barycentric_coordinates(p15, p21, p22, pp); }
	if (DGP::is_valid_barycentric_coordinate(DGP::calc_barycentric_coordinates(p22, p16, p15, pp))) { tri_idx = 6; bc = DGP::calc_barycentric_coordinates(p22, p16, p15, pp); }
	if (DGP::is_valid_barycentric_coordinate(DGP::calc_barycentric_coordinates(p16, p22, p23, pp))) { tri_idx = 7; bc = DGP::calc_barycentric_coordinates(p16, p22, p23, pp); }
	if (DGP::is_valid_barycentric_coordinate(DGP::calc_barycentric_coordinates(p23, p17, p16, pp))) { tri_idx = 8; bc = DGP::calc_barycentric_coordinates(p23, p17, p16, pp); }
	if (DGP::is_valid_barycentric_coordinate(DGP::calc_barycentric_coordinates(p17, p23, p24, pp))) { tri_idx = 9; bc = DGP::calc_barycentric_coordinates(p17, p23, p24, pp); }
	if (DGP::is_valid_barycentric_coordinate(DGP::calc_barycentric_coordinates(p24, p18, p17, pp))) { tri_idx = 10; bc = DGP::calc_barycentric_coordinates(p24, p18, p17, pp); }
	if (DGP::is_valid_barycentric_coordinate(DGP::calc_barycentric_coordinates(p18, p24, p25, pp))) { tri_idx = 11; bc = DGP::calc_barycentric_coordinates(p18, p24, p25, pp); }
	_new_v = bc[1];
	_new_w = bc[2];
	std::cout << "Method 2: tri_idx: " << tri_idx << ". (new_v, new_w)= (" << _new_v << ", " << _new_w << ").\n";

	return tri_idx;
}
bool Hoppe94LoopEval::eval_regular_1v0e_old(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, double _v, double _w, OpenMesh::Vec3d &_p) {
	std::cout << "\n";
	std::cout << "Hoppe94LoopEval::Begin eval_regular_1v0e_old().\n";
	if (_v + _w < 0 ) std::cout << "Error: eval_regular_1v0e_old(): v+w < 0.\n"; 
	if (_v + _w > 1) { std::cout << "Error: eval_regular_1v0e_old(), 1 < v+w now.\n"; }
	if (_c.size() != 12) std::cout << "Error: Hoppe94LoopEval::eval_regular_1v0e_old(), 12 control points are needed.\n";


	int n_stam98_level = (int)floor(1 - log(_v + _w) / log(2.0));
	int n_level = n_stam98_level + 1;
	std::cout << "n_stam98_level: " << n_stam98_level << ", but one more subdiv (" << n_level << ") is neeeded here.\n";

	// 第一步, 
	std::vector<OpenMesh::Vec3d> c33(33);   // 从12个初始控制点计算出第二层细分下的33个控制点
	eval_regular_1v0e_c33(_c, _h, n_stam98_level, c33);


	// 第二步, 那个函数里面用了两种方法求出的结果一致, 哪里错了呢?
	int tri_idx = -1; // 0 ~ 11 is valid.
	double new_v = 0, new_w = 0;          
	tri_idx = eval_regular_1v0e_tri_idx(_v, _w, n_stam98_level, new_v, new_w); // 看(_v, _w)是在哪个小三角形中, 并且返回归一(单位化)放大的(new_v, new_w).

	// 至此, tri_idx指出是12中哪个小三角形, new_v, new_w又给出了参数坐标, 现在可以从33个点中选出12给作为控制顶点.
	// 作图可以数出下面12个pick模板
	// 第三步
	int pick[12][12] = { //12个选择模板: 前面的12是指tri_idx指向的12个小三角形, 后面的12个序号表示33个控制点中取哪些.
		{ 5, 6, 9, 10, 11, 14, 15, 16, 17, 21, 22, 23 },
		{ 23, 22, 17, 16, 15, 12, 11, 10, 9, 7, 6, 5 },
		{ 6, 7, 10, 11, 12, 15, 16, 17, 18, 22, 23, 24 },
		{ 24, 23, 18, 17, 16, 13, 12, 11, 10, 8, 7, 6 },
		{ 7, 8, 11, 12, 13, 16, 17, 18, 19, 23, 24, 25 },
		{ 9, 10, 14, 15, 16, 20, 21, 22, 23, 27, 28, 29 },
		{ 29, 28, 23, 22, 21, 17, 16, 15, 14, 11, 10, 9 },
		{ 10, 11, 15, 16, 17, 21, 22, 23, 24, 28, 29, 30 },
		{ 30, 29, 24, 23, 22, 18, 17, 16, 15, 12, 11, 10 }, // cell 8
		{ 11, 12, 16, 17, 18, 22, 23, 24, 25, 29, 30, 31 },
		{ 31, 30, 25, 24, 23, 19, 18, 17, 16, 13, 12, 11 },
		{ 12, 13, 17, 18, 19, 23, 24, 25, 26, 30, 31, 32 }
	};

	std::vector<OpenMesh::Vec3d> new_c12(12); //根据tri_idx确定使用的pick mask来从33个控制点中选择12个
	for (int i = 0; i < 12; ++i) {
		new_c12[i] = c33[pick[tri_idx][i]];
	}
	OpenMesh::Vec3d Sk = new_c12[3] * (1-new_v-new_w) + new_c12[6] * new_v + new_c12[7] * new_w;
	std::cout << "Sk " << Sk << ".\n";

	Stam98ClassicLoopEval cloop_evaluation;
	if (cloop_evaluation.evaluate_regular_patch(new_c12, new_v, new_w, _p) == false) {
		std::cout << "Error: stam98 test regular patch.\n";
	}
	std::cout << "_p " << _p << ".\n";

	std::cout << "Hoppe94LoopEval::eval_regular_1v0e_old() End.\n";

	return true;
}
// --------
void Hoppe94LoopEval::eval_regular_1v0e_c18(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, std::vector<OpenMesh::Vec3d> &_c18) {
	// one time subdivision: 12 -> 18 
	double d1_16 = 0.0625, d1_8 = 0.125, d3_8 = 0.375, d5_8 = 0.625;
	// first time subdiv, 12 -> 18 
	double A1[18][12] = {
		//0     1       2       3       4       5       6       7       8       9      10       11   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 0  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 1  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 2
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 3   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 4
		{ 0,	0,	d3_8,	d1_8,		0,	d1_8,	d3_8,		0,		0,		0,		0,		0 }, // index 5   
		{ 0,	0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8,		0,		0,		0,		0 }, // index 6
		{ 0,	0,		0,	d3_8,	d1_8,		0,	d1_8,	d3_8,		0,		0,		0,		0 }, // index 7  
		{ 0,	0,		0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8,		0,		0,		0 }, // index 8
		{ 0,	0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8,		0,		0 }, // index 9   
		{ 0,	0,	d1_16,	d1_16,		0,	d1_16,	d5_8,	d1_16,		0,	d1_16,	d1_16,		0 }, // index 10
		{ 0,	0,		0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8,		0 }, // index 11  
		{ 0,	0,		0,	d1_16,	d1_16,		0,	d1_16,	d5_8,	d1_16,		0,	d1_16,	d1_16 }, // index 12
		{ 0,	0,		0,		0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8 }, // index 13  
		{ 0,	0,		0,		0,		0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8,		0 }, // index 14
		{ 0,	0,		0,		0,		0,		0,	d3_8,	d1_8,		0,	d1_8,	d3_8,		0 }, // index 15  
		{ 0,	0,		0,		0,		0,		0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8 }, // index 16 
		{ 0,	0,		0,		0,		0,		0,		0,	d3_8,	d1_8,		0,	d1_8,	d3_8 }, // index 17  
	};
	// 下面确定第0, 1, 2, 4, 3号点的生成模板
	int n_count_feature_edges = 0;
	std::vector<int> fv_tmp;
	if (_h[0]) { // new control point 0.
		A1[0][0] = A1[0][3] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(0);
	} else {
		A1[0][0] = A1[0][3] = d3_8;  A1[0][1] = A1[0][2] = d1_8;  
	}
	if (_h[1]) { // new control point 1.
		A1[1][1] = A1[1][3] = 0.5;   ++n_count_feature_edges; fv_tmp.push_back(1);
	} else {
		A1[1][1] = A1[1][3] = d3_8;  A1[1][0] = A1[1][4] = d1_8;
	}
	if (_h[2]) { // new control point 2.
		A1[2][2] = A1[2][3] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(2);
	} else {
		A1[2][2] = A1[2][3] = d3_8;  A1[2][0] = A1[2][6] = d1_8;
	}
	if (_h[3]) { // new control point 4.
		A1[4][3] = A1[4][4] = 0.5;   ++n_count_feature_edges; fv_tmp.push_back(4);
	} else {
		A1[4][3] = A1[4][4] = d3_8;  A1[4][1] = A1[4][7] = d1_8;  
	}
	if (n_count_feature_edges <= 1) { // smooth or dart // new control point 3.
		A1[3][3] = d5_8;  A1[3][0] = A1[3][1] = A1[3][2] = A1[3][4] = A1[3][6] = A1[3][7] = d1_16;
	} else if (n_count_feature_edges == 2) { // crease
		A1[3][3] = 0.75;  A1[3][fv_tmp[0]] = A1[3][fv_tmp[1]] = d1_8; 
	} else { // corner
		A1[3][3] = 1;
	}
	
	for (int i = 0; i < 18; ++i) { 
		_c18[i] = OpenMesh::Vec3d(0, 0, 0);
		for (int j = 0; j < 12; ++j) {
			_c18[i] += A1[i][j] * _c[j];
		} 
	}
}
bool Hoppe94LoopEval::eval_regular_1v0e(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, double _v, double _w, OpenMesh::Vec3d &_p) {
	//std::cout << "\n";
	//std::cout << "Hoppe94LoopEval::Begin eval_regular_1v0e().\n";
	if (_v + _w < 0 ) {
		std::cout << "Error: eval_regular_1v0e(): v+w < 0, " << _v << ", " << _w << ".\n"; 
		return false;
	}
	if (_v + _w > 1) { 
		std::cout << "Error: eval_regular_1v0e(), 1 < v+w, " << _v << ", " << _w << "\n"; 
		return false;
	}
	if (_c.size() != 12) {
		std::cout << "Error: Hoppe94LoopEval::eval_regular_1v0e(), 12 control points are needed.\n";
		return false;
	}
	//std::cout << "_v_w" << _v <<"&"<< _w << " ";
	//for (int i = 0; i < 12; ++i) 
		//std::cout << "c " << _c[i] << ".\n";
	//for (int i = 0; i < 15; ++i) 
		//std::cout << "h " << _h[i] << ".\n";

	// ------特殊处理 _v and _w == 0, 根据点的类型来直接求出其值.
	if (_v < 1e-6 && _w < 1e-6) {
		std::vector<int> fv_tmp;
		if (_h[0]) fv_tmp.push_back(0);
		if (_h[1]) fv_tmp.push_back(1);
		if (_h[2]) fv_tmp.push_back(2);
		if (_h[3]) fv_tmp.push_back(4); 
		if (fv_tmp.size() >= 3) {
			_p = _c[3]; return true;
		} else if (fv_tmp.size() == 2) {
			_p = _c[3] * 0.66667 + (_c[fv_tmp[0]] + _c[fv_tmp[1]]) * 0.166667; return true;
		} else if (fv_tmp.size() == 1) { // dart, valence == 6
			_p = _c[3] * 0.5 + (_c[0] + _c[1] + _c[2] + _c[4] + _c[6] + _c[7]) * 0.08333; return true;
		} else {
			std::cout << "Error: Hoppe94LoopEval::eval_regular_1v0e(): ==0 是DGP::DGP::SMOOTH_VFT, 做为非特征考虑的,不在这考虑.\n";
			return false;
		}
	}
	// -----
	std::vector<OpenMesh::Vec3d> c18(18, OpenMesh::Vec3d(0, 0, 0));   // 从12个初始控制点计算出第一层细分下的18个控制点
	//eval_regular_1v0e_c18(_c, _h, c18);
	eval_regular_c18(_c, _h, c18);
	//std::cout << "hEre";

	if (_v + _w <= 0.5) { 
		//方法一, 递归调用. 方法简单但是低效益.
		int pick[12] = {0, 1,    2, 3, 4,    5, 6, 7, 8,    10, 11, 12}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0)); // 
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		v = _v * 2; w = _w * 2;

		std::vector<bool> h(15, false);
		h[0] = _h[0]; h[1] = _h[1]; h[2] = _h[2]; h[3] = _h[3]; 
		
		//std::cout << "in v + w <= 0.5." << _v << ", " << _w << ".\n";
		return eval_regular_1v0e(c12, h, v, w, _p); 
		//方法二, 
	}	
	//std::cout << "Here";
	if (_v > 0.5) { //
		int pick[12] = {2, 3,     5, 6, 7,    9, 10, 11, 12,   14, 15, 16};
		std::vector<OpenMesh::Vec3d> c12(12);
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		v = (_v - 0.5) * 2; w = _w * 2;

		//std::cout << "in v > 0.5." << v << ", " << w << ".\n";
		return eval_regular_1in4_0v0e(c12, v, w, _p);
	}
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) { 
		int pick[12] = { 5, 2,    10, 6, 3,    15, 11, 7, 4,    16, 12, 8};  
		std::vector<OpenMesh::Vec3d> c12; //
		for (int i = 0; i < 12; ++i) {
			c12.push_back(c18[pick[i]]); //std::cout << c12[i] << ".\n";
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];// the same as:
		////v = (_v + _w - 0.5) * 2; w = (0.5 - _v) * 2; // u = (0.5 - _w) * 2;
		//std::cout << "in v + w > 0.5." << _v << ", " << _w << ".\n";
		return eval_regular_1in4_0v0e(c12, v, w, _p);  
	} 

	if (_w > 0.5) { 
		int pick[12] = { 3, 4,    6, 7, 8,    10, 11, 12, 13,    15, 16, 17}; //pick 12 vertices from the 18 vertices, which created before.
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		v = _v * 2; w = (_w - 0.5) * 2;
		//std::cout << "in w > 0.5." << _v << ", " << _w << ".\n";
		return eval_regular_1in4_0v0e(c12, v, w, _p);  
	} 
	std::cout << "Hoppe94LoopEval::eval_regular_1v0e() End." << _v + _w << ".\n";
	return false;
}
void Hoppe94LoopEval::eval_regular_2v0e_c18(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, std::vector<OpenMesh::Vec3d> &_c18) {
	// one time subdivision: 12 -> 18 
	double d1_16 = 0.0625, d1_8 = 0.125, d3_8 = 0.375, d5_8 = 0.625;
	// first time subdiv, 12 -> 18 
	double A1[18][12] = {
		//0     1       2       3       4       5       6       7       8       9      10       11   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 0  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 1  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 2
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 3   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 4
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 5   
		{ 0,	0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8,		0,		0,		0,		0 }, // index 6
		{ 0,	0,		0,	d3_8,	d1_8,		0,	d1_8,	d3_8,		0,		0,		0,		0 }, // index 7  
		{ 0,	0,		0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8,		0,		0,		0 }, // index 8
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 9   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 10
		{ 0,	0,		0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8,		0 }, // index 11  
		{ 0,	0,		0,	d1_16,	d1_16,		0,	d1_16,	d5_8,	d1_16,		0,	d1_16,	d1_16 }, // index 12
		{ 0,	0,		0,		0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8 }, // index 13  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 14
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 15  
		{ 0,	0,		0,		0,		0,		0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8 }, // index 16 
		{ 0,	0,		0,		0,		0,		0,		0,	d3_8,	d1_8,		0,	d1_8,	d3_8 }, // index 17  
	};
	// 下面确定第0, 1, 2, 4, 3号点的生成模板
	int n_count_feature_edges = 0;
	std::vector<int> fv_tmp;
	if (_h[0]) { // new control point 0.
		A1[0][0] = A1[0][3] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(0);
	} else {
		A1[0][0] = A1[0][3] = d3_8;  A1[0][1] = A1[0][2] = d1_8;  
	}
	if (_h[1]) { // new control point 1.
		A1[1][1] = A1[1][3] = 0.5;   ++n_count_feature_edges; fv_tmp.push_back(1);
	} else {
		A1[1][1] = A1[1][3] = d3_8;  A1[1][0] = A1[1][4] = d1_8;
	}
	if (_h[2]) { // new control point 2.
		A1[2][2] = A1[2][3] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(2);
	} else {
		A1[2][2] = A1[2][3] = d3_8;  A1[2][0] = A1[2][6] = d1_8;
	}
	if (_h[3]) { // new control point 4.
		A1[4][3] = A1[4][4] = 0.5;   ++n_count_feature_edges; fv_tmp.push_back(4);
	} else {
		A1[4][3] = A1[4][4] = d3_8;  A1[4][1] = A1[4][7] = d1_8;  
	}
	if (n_count_feature_edges <= 1) { // smooth or dart // new control point 3.
		A1[3][3] = d5_8;  A1[3][0] = A1[3][1] = A1[3][2] = A1[3][4] = A1[3][6] = A1[3][7] = d1_16;
	} else if (n_count_feature_edges == 2) { // crease
		A1[3][3] = 0.75;  A1[3][fv_tmp[0]] = A1[3][fv_tmp[1]] = d1_8; 
	} else { // corner
		A1[3][3] = 1;
	}
	// 下面确定第5, 9, 14, 15, 10号点的生成模板
	n_count_feature_edges = 0; fv_tmp.clear();
	if (_h[4]) { // new control point 5.
		A1[5][2] = A1[5][6] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(2);
	} else {
		A1[5][2] = A1[5][6] = d3_8;  A1[5][3] = A1[5][5] = d1_8;  
	}
	if (_h[8]) { // new control point 9.
		A1[9][5] = A1[9][6] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(5);
	} else {
		A1[9][5] = A1[9][6] = d3_8;  A1[9][2] = A1[9][9] = d1_8;  
	}
	if (_h[11]) { // new control point 14.
		A1[14][6] = A1[14][9] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(9);
	} else {
		A1[14][6] = A1[14][9] = d3_8;  A1[14][5] = A1[14][10] = d1_8;  
	}
	if (_h[12]) { // new control point 15.
		A1[15][6] = A1[15][10] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(10);
	} else {
		A1[15][6] = A1[15][10] = d3_8;  A1[15][7] = A1[15][9] = d1_8;  
	}
	if (n_count_feature_edges <= 1) { // smooth or dart // new control point 10.
		A1[10][6] = d5_8;  A1[10][2] = A1[10][3] = A1[10][5] = A1[10][7] = A1[10][9] = A1[10][10] = d1_16;
	} else if (n_count_feature_edges == 2) { // crease
		A1[10][6] = 0.75;  A1[10][fv_tmp[0]] = A1[10][fv_tmp[1]] = d1_8; 
	} else { // corner
		A1[10][6] = 1;
	}

	for (int i = 0; i < 18; ++i) {
		_c18[i] = OpenMesh::Vec3d(0, 0, 0);
		for (int j = 0; j < 12; ++j) {
			_c18[i] += A1[i][j] * _c[j];
		} 
	}
}
bool Hoppe94LoopEval::eval_regular_2v0e(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, double _v, double _w, OpenMesh::Vec3d &_p) {
	if ( _v + _w < 0 || _v + _w > 1) std::cout << "Error: Hoppe94LoopEval::eval_regular_2v0e(): _v + _w < 0 or > 1.\n";
	if (_c.size() != 12) std::cout << "Error: Hoppe94LoopEval::eval_regular_2v0e(): _c.size() != 12.\n";
	if (_h.size() != 15) std::cout << "Error: Hoppe94LoopEval::eval_regular_2v0e(): _h.size() != 15.\n";
	//std::cout << "Hoppe94LoopEval::eval_regular_2v0e()...\n";

	std::vector<OpenMesh::Vec3d> c18(18, OpenMesh::Vec3d(0, 0, 0));//
	//eval_regular_2v0e_c18(_c, _h, c18);
	eval_regular_c18(_c, _h, c18); 

	if (_w > 0.5) {
		int pick[12] = { 3, 4,    6, 7, 8,    10, 11, 12, 13,    15, 16, 17}; //pick 12 vertices from the 18 vertices, which created before.
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		v = _v * 2; w = (_w - 0.5) * 2;
		return eval_regular_1in4_0v0e(c12, v, w, _p);  
	}
	
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) { // middle 1/4 
		int pick[12] = { 5, 2,    10, 6, 3,    15, 11, 7, 4,    16, 12, 8}; //pick 12 vertices from the 18 vertices, which created before.
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		//OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));
		//v = bc[1]; w = bc[2];// the same as:
		v = (_v + _w - 0.5) * 2; w = (0.5 - _v) * 2; // u = (0.5 - _w) * 2;
		return eval_regular_1in4_0v0e(c12, v, w, _p);  
	}

	// 当点靠近那两个特征点(的其中一个). 
	if (_v + _w <= 0.5) { // 点靠近特征点A// top 1/4		
		int pick[12] = {0, 1,    2, 3, 4,    5, 6, 7, 8,    10, 11, 12}; 
		std::vector<OpenMesh::Vec3d> c12(12); // 
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		v = _v * 2; w = _w * 2;

		std::vector<bool> h(15, false);
		h[0] = _h[0]; h[1] = _h[1]; h[2] = _h[2]; h[3] = _h[3]; 
					
		return eval_regular_1v0e(c12, h, v, w, _p);
	} 
	
	if (_v > 0.5) { //点更靠近特征点B
		int pick[12] = {14, 9,    15, 10, 5,     16, 11, 6, 2,    12, 7, 3};
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(1, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];
		std::vector<bool> h(15, false);
		h[0] = _h[11]; h[1] = _h[8]; h[2] = _h[12]; h[3] = _h[4]; 
		return eval_regular_1v0e(c12, h, v, w, _p);  //  
	}
	return false; 
}
void Hoppe94LoopEval::eval_regular_3v0e_c18(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, std::vector<OpenMesh::Vec3d> &_c18) {
	// one time subdivision: 12 -> 18 
	double d1_16 = 0.0625, d1_8 = 0.125, d3_8 = 0.375, d5_8 = 0.625;
	// first time subdiv, 12 -> 18 
	double A1[18][12] = {
		//0     1       2       3       4       5       6       7       8       9      10       11   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 0  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 1  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 2
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 3   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 4
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 5   
		{ 0,	0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8,		0,		0,		0,		0 }, // index 6
		{ 0,	0,		0,	d3_8,	d1_8,		0,	d1_8,	d3_8,		0,		0,		0,		0 }, // index 7  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 8
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 9   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 10
		{ 0,	0,		0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8,		0 }, // index 11 
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 12
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 13  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 14
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 15  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 16 
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 17  
	};

	// 下面确定第0, 1, 2, 4, 3号点的生成模板
	int n_count_feature_edges = 0;
	std::vector<int> fv_tmp;
	if (_h[0]) { // new control point 0.
		A1[0][0] = A1[0][3] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(0);
	} else {
		A1[0][0] = A1[0][3] = d3_8;  A1[0][1] = A1[0][2] = d1_8;  
	}
	if (_h[1]) { // new control point 1.
		A1[1][1] = A1[1][3] = 0.5;   ++n_count_feature_edges; fv_tmp.push_back(1);
	} else {
		A1[1][1] = A1[1][3] = d3_8;  A1[1][0] = A1[1][4] = d1_8;
	}
	if (_h[2]) { // new control point 2.
		A1[2][2] = A1[2][3] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(2);
	} else {
		A1[2][2] = A1[2][3] = d3_8;  A1[2][0] = A1[2][6] = d1_8;
	}
	if (_h[3]) { // new control point 4.
		A1[4][3] = A1[4][4] = 0.5;   ++n_count_feature_edges; fv_tmp.push_back(4);
	} else {
		A1[4][3] = A1[4][4] = d3_8;  A1[4][1] = A1[4][7] = d1_8;  
	}
	if (n_count_feature_edges <= 1) { // smooth or dart // new control point 3.
		A1[3][3] = d5_8;  A1[3][0] = A1[3][1] = A1[3][2] = A1[3][4] = A1[3][6] = A1[3][7] = d1_16;
	} else if (n_count_feature_edges == 2) { // crease
		A1[3][3] = 0.75;  A1[3][fv_tmp[0]] = A1[3][fv_tmp[1]] = d1_8; 
	} else { // corner
		A1[3][3] = 1;
	}
	// 下面确定第5, 9, 14, 15, 10号点的生成模板
	n_count_feature_edges = 0; fv_tmp.clear();
	if (_h[4]) { // new control point 5.
		A1[5][2] = A1[5][6] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(2);
	} else {
		A1[5][2] = A1[5][6] = d3_8;  A1[5][3] = A1[5][5] = d1_8;  
	}
	if (_h[8]) { // new control point 9.
		A1[9][5] = A1[9][6] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(5);
	} else {
		A1[9][5] = A1[9][6] = d3_8;  A1[9][2] = A1[9][9] = d1_8;  
	}
	if (_h[11]) { // new control point 14.
		A1[14][6] = A1[14][9] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(9);
	} else {
		A1[14][6] = A1[14][9] = d3_8;  A1[14][5] = A1[14][10] = d1_8;  
	}
	if (_h[12]) { // new control point 15.
		A1[15][6] = A1[15][10] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(10);
	} else {
		A1[15][6] = A1[15][10] = d3_8;  A1[15][7] = A1[15][9] = d1_8;  
	}
	if (n_count_feature_edges <= 1) { // smooth or dart // new control point 10.
		A1[10][6] = d5_8;  A1[10][2] = A1[10][3] = A1[10][5] = A1[10][7] = A1[10][9] = A1[10][10] = d1_16;
	} else if (n_count_feature_edges == 2) { // crease
		A1[10][6] = 0.75;  A1[10][fv_tmp[0]] = A1[10][fv_tmp[1]] = d1_8; 
	} else { // corner
		A1[10][6] = 1;
	}
	// 
	n_count_feature_edges = 0; fv_tmp.clear();
	if (_h[7]) { // new control point 8.
		A1[8][4] = A1[8][7] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(4);
	} else {
		A1[8][4] = A1[8][7] = d3_8;  A1[8][3] = A1[8][8] = d1_8;  
	}
	if (_h[10]) { // new control point 13.
		A1[13][7] = A1[13][8] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(8);
	} else {
		A1[13][7] = A1[13][8] = d3_8;  A1[13][4] = A1[13][11] = d1_8;  
	}
	if (_h[13]) { // new control point 16.
		A1[16][7] = A1[16][10] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(10);
	} else {
		A1[16][7] = A1[16][10] = d3_8;  A1[16][6] = A1[16][11] = d1_8;  
	}
	if (_h[14]) { // new control point 17.
		A1[17][7] = A1[17][11] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(11);
	} else {
		A1[17][7] = A1[17][11] = d3_8;  A1[17][8] = A1[17][10] = d1_8;  
	}
	if (n_count_feature_edges <= 1) { // dart // new control point 12.
		A1[12][7] = d5_8;  A1[12][3] = A1[12][4] = A1[12][6] = A1[12][8] = A1[12][10] = A1[12][11] = d1_16;
	} else if (n_count_feature_edges == 2) { // crease
		A1[12][7] = 0.75;  A1[12][fv_tmp[0]] = A1[12][fv_tmp[1]] = d1_8; 
	} else { // corner
		A1[12][7] = 1;
	}

	for (int i = 0; i < 18; ++i) {
		_c18[i] = OpenMesh::Vec3d(0, 0, 0);
		for (int j = 0; j < 12; ++j) {
			_c18[i] += A1[i][j] * _c[j];
		} 
	}
}
bool Hoppe94LoopEval::eval_regular_3v0e(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_v + _w < 0 || _v + _w > 1) std::cout << "Error: Hoppe94LoopEval::eval_regular_3v0e(): _v + _w < 0 or > 1.\n";

	std::vector<OpenMesh::Vec3d> c18(18, OpenMesh::Vec3d(0, 0, 0));
	eval_regular_3v0e_c18(_c, _h, c18);

	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) { // middle 1/4 
		int pick[12] = { 5, 2,    10, 6, 3,    15, 11, 7, 4,    16, 12, 8}; //pick 12 vertices from the 18 vertices, which created before.
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		//OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));
		//v = bc[1]; w = bc[2];// the same as:
		v = (_v + _w - 0.5) * 2; w = (0.5 - _v) * 2; // u = (0.5 - _w) * 2;
		return eval_regular_1in4_0v0e(c12, v, w, _p);  
	}

	if (_v + _w <= 0.5) { // 点靠近特征点A// top 1/4		
		int pick[12] = {0, 1,    2, 3, 4,    5, 6, 7, 8,    10, 11, 12}; 
		std::vector<OpenMesh::Vec3d> c12(12); // 
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		v = _v * 2; w = _w * 2;

		std::vector<bool> h(15, false);
		h[0] = _h[0]; h[1] = _h[1]; h[2] = _h[2]; h[3] = _h[3]; 

		return eval_regular_1v0e(c12, h, v, w, _p);
	} 

	if (_v > 0.5) { //点更靠近特征点B
		int pick[12] = {14, 9,    15, 10, 5,     16, 11, 6, 2,    12, 7, 3};
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(1, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];

		std::vector<bool> h(15, false);
		h[0] = _h[11]; h[1] = _h[8]; h[2] = _h[12]; h[3] = _h[4]; 
		return eval_regular_1v0e(c12, h, v, w, _p);  //  
	}

	if (_w > 0.5) { // low right side (1/4) can be evaluation
		int pick[12] = { 13, 17,    8, 12, 16,    4, 7, 11, 15,    3, 6, 10}; //pick 12 vertices from the 18 vertices, which created before.
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];

		std::vector<bool> h(15, false);
		h[0] = _h[10]; h[1] = _h[14];
		h[2] = _h[7]; h[3] = _h[13]; 
		return eval_regular_1v0e(c12, h, v, w, _p); 
	}

	return false; 
}
void Hoppe94LoopEval::eval_regular_2v1e_c18(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, std::vector<OpenMesh::Vec3d> &_c18) {
	// one time subdivision: 12 -> 18 
	double d1_16 = 0.0625, d1_8 = 0.125, d3_8 = 0.375, d5_8 = 0.625;
	// first time subdiv, 12 -> 18 
	double A1[18][12] = {
		//0     1       2       3       4       5       6       7       8       9      10       11   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 0  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 1  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 2
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 3   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 4
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 5   
		{ 0,	0,		0,		0.5,	0,		0,		0.5,	0,		0,		0,		0,		0 }, // index 6
		{ 0,	0,		0,	d3_8,	d1_8,		0,	d1_8,	d3_8,		0,		0,		0,		0 }, // index 7  
		{ 0,	0,		0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8,		0,		0,		0 }, // index 8
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 9   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 10
		{ 0,	0,		0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8,		0 }, // index 11  
		{ 0,	0,		0,	d1_16,	d1_16,		0,	d1_16,	d5_8,	d1_16,		0,	d1_16,	d1_16 }, // index 12
		{ 0,	0,		0,		0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8 }, // index 13  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 14
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 15  
		{ 0,	0,		0,		0,		0,		0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8 }, // index 16 
		{ 0,	0,		0,		0,		0,		0,		0,	d3_8,	d1_8,		0,	d1_8,	d3_8 }, // index 17  
	};
	// 下面确定第0, 1, 2, 4, 3号点的生成模板
	int n_count_feature_edges = 0, fv_tmp = -1;
	if (_h[0]) { // new control point 0.
		A1[0][0] = A1[0][3] = 0.5;  ++n_count_feature_edges; fv_tmp = 0;
	} else {
		A1[0][0] = A1[0][3] = d3_8;  A1[0][1] = A1[0][2] = d1_8;  
	}
	if (_h[1]) { // new control point 1.
		A1[1][1] = A1[1][3] = 0.5;   ++n_count_feature_edges; fv_tmp = 1;
	} else {
		A1[1][1] = A1[1][3] = d3_8;  A1[1][0] = A1[1][4] = d1_8;
	}
	if (_h[2]) { // new control point 2.
		A1[2][2] = A1[2][3] = 0.5;  ++n_count_feature_edges; fv_tmp = 2;
	} else {
		A1[2][2] = A1[2][3] = d3_8;  A1[2][0] = A1[2][6] = d1_8;
	}
	if (_h[3]) { // new control point 4.
		A1[4][3] = A1[4][4] = 0.5;   ++n_count_feature_edges; fv_tmp = 4;
	} else {
		A1[4][3] = A1[4][4] = d3_8;  A1[4][1] = A1[4][7] = d1_8;  
	}
	if (n_count_feature_edges == 0) { // smooth or dart // new control point 3.
		A1[3][3] = d5_8;  A1[3][0] = A1[3][1] = A1[3][2] = A1[3][4] = A1[3][6] = A1[3][7] = d1_16;
	} else if (n_count_feature_edges == 1) { // crease
		A1[3][3] = 0.75;  A1[3][fv_tmp] = A1[3][6] = d1_8; 
	} else { // corner
		A1[3][3] = 1;
	}
	// 下面确定第5, 9, 14, 15, 10号点的生成模板
	n_count_feature_edges = 0; fv_tmp = -1;
	if (_h[4]) { // new control point 5.
		A1[5][2] = A1[5][6] = 0.5;  ++n_count_feature_edges; fv_tmp = 2;
	} else {
		A1[5][2] = A1[5][6] = d3_8;  A1[5][3] = A1[5][5] = d1_8;  
	}
	if (_h[8]) { // new control point 9.
		A1[9][5] = A1[9][6] = 0.5;  ++n_count_feature_edges; fv_tmp = 5;
	} else {
		A1[9][5] = A1[9][6] = d3_8;  A1[9][2] = A1[9][9] = d1_8;  
	}
	if (_h[11]) { // new control point 14.
		A1[14][6] = A1[14][9] = 0.5;  ++n_count_feature_edges; fv_tmp = 9;
	} else {
		A1[14][6] = A1[14][9] = d3_8;  A1[14][5] = A1[14][10] = d1_8;  
	}
	if (_h[12]) { // new control point 15.
		A1[15][6] = A1[15][10] = 0.5;  ++n_count_feature_edges; fv_tmp = 10;
	} else {
		A1[15][6] = A1[15][10] = d3_8;  A1[15][7] = A1[15][9] = d1_8;  
	}
	if (n_count_feature_edges == 0) { // smooth or dart // new control point 10.
		A1[10][6] = d5_8;  A1[10][2] = A1[10][3] = A1[10][5] = A1[10][7] = A1[10][9] = A1[10][10] = d1_16;
	} else if (n_count_feature_edges == 1) { // crease
		A1[10][6] = 0.75;  A1[10][fv_tmp] = A1[10][3] = d1_8; 
	} else { // corner
		A1[10][6] = 1;
	}
	
	for (int i = 0; i < 18; ++i) {
		_c18[i] = OpenMesh::Vec3d(0, 0, 0);
		for (int j = 0; j < 12; ++j) {
			_c18[i] += A1[i][j] * _c[j];
		} 
	}
}
bool Hoppe94LoopEval::eval_regular_1in4_0v0e(const std::vector<OpenMesh::Vec3d> &_c, double _v, double _w, OpenMesh::Vec3d &_p) {
	//std::cout << "Here: Hoppe94LoopEval::eval_regular_1in4_0v0e().\n"; // for test
	// subdiv one time, 12 -> 18.
	double d1_16 = 0.0625, d1_8 = 0.125, d3_8 = 0.375, d5_8 = 0.625;
	double A1[18][12] = {
		//0     1       2       3       4       5       6       7       8       9      10       11   
		{ d3_8,	d1_8,d1_8,	d3_8,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 0  
		{ d1_8,	d3_8,	0,	d3_8,	d1_8,		0,		0,		0,		0,		0,		0,		0 }, // index 1  
		{ d1_8,	0,	d3_8,	d3_8,		0,		0,	d1_8,		0,		0,		0,		0,		0 }, // index 2
		{ d1_16,d1_16,d1_16,d5_8,	d1_16,		0,	d1_16,	d1_16,		0,		0,		0,		0 }, // index 3   
		{ 0,d1_8,		0,	d3_8,	d3_8,		0,		0,	d1_8,		0,		0,		0,		0 }, // index 4
		{ 0,	0,	d3_8,	d1_8,		0,	d1_8,	d3_8,		0,		0,		0,		0,		0 }, // index 5   
		{ 0,	0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8,		0,		0,		0,		0 }, // index 6
		{ 0,	0,		0,	d3_8,	d1_8,		0,	d1_8,	d3_8,		0,		0,		0,		0 }, // index 7  
		{ 0,	0,		0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8,		0,		0,		0 }, // index 8
		{ 0,	0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8,		0,		0 }, // index 9   
		{ 0,	0,	d1_16,	d1_16,		0,	d1_16,	d5_8,	d1_16,		0,	d1_16,	d1_16,		0 }, // index 10
		{ 0,	0,		0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8,		0 }, // index 11  
		{ 0,	0,		0,	d1_16,	d1_16,		0,	d1_16,	d5_8,	d1_16,		0,	d1_16,	d1_16 }, // index 12
		{ 0,	0,		0,		0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8 }, // index 13  
		{ 0,	0,		0,		0,		0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8,		0 }, // index 14
		{ 0,	0,		0,		0,		0,		0,	d3_8,	d1_8,		0,	d1_8,	d3_8,		0 }, // index 15  
		{ 0,	0,		0,		0,		0,		0,	d1_8,	d3_8,		0,		0,	d3_8,	d1_8 }, // index 16 
		{ 0,	0,		0,		0,		0,		0,		0,	d3_8,	d1_8,		0,	d1_8,	d3_8 }, // index 17  
	};

	std::vector<OpenMesh::Vec3d> c18(18, OpenMesh::Vec3d(0, 0, 0));
	for (int i = 0; i < 18; ++i) {
		for (int j = 0; j < 12; ++j) {
			c18[i] += A1[i][j] * _c[j];
		}
	}

	int pick[4][12] = {
		{0, 1,     2, 3, 4,    5, 6, 7, 8,    10, 11, 12},  
		{2, 3,     5, 6, 7,    9, 10, 11, 12,    14, 15, 16},
		{16, 15,   12, 11, 10,    8, 7, 6, 5,    4, 3, 2},
		{3, 4,     6, 7, 8,    10, 11, 12, 13,    15, 16, 17},
	};

	//std::cout << "1in4: " << _v << ", " << _w << ".\n";
	int tri_idx = -1;
	double v_new = 0, w_new = 0;
	if (_v + _w <= 0.5) {
		tri_idx = 0;
		v_new = _v * 2; w_new = _w * 2;
	} else if (_v > 0.5) {
		tri_idx = 1;
		v_new = (_v - 0.5) * 2; w_new = _w * 2;
	} else if (_w > 0.5) { //****
		tri_idx = 3;
		v_new = _v * 2; w_new = (_w - 0.5) * 2;
	} else {
		tri_idx = 2;
		v_new = (0.5 - _v) * 2; w_new = (0.5 - _w) * 2;
	}

	std::vector<OpenMesh::Vec3d> c12_new(12, OpenMesh::Vec3d(0, 0, 0));
	for (int i = 0; i < 12; ++i) {
		c12_new[i] = c18[pick[tri_idx][i]];
	}

	//OpenMesh::Vec3d Sk = c12_new[3] * (1- v_new - w_new) + c12_new[6] * v_new + c12_new[7] * w_new;
	//std::cout << "1in4 Sk " << Sk << ".\n";

	//if (cloop_evaluator_->evaluate_regular_patch(c12_new, v_new, w_new, _p) == false) {
	//	std::cout << "Error: stam98 test regular patch.\n";
	//}
	if (cloop_evaluator_->evaluate_regular_patch(c12_new, v_new, w_new, _p, dv_, dw_) == false) {
		std::cout << "Error: stam98 test regular patch.\n";
	}
	//std::cout << "1in4 _p " << _p << ".\n";
	return true;
}
bool Hoppe94LoopEval::eval_regular_2v1e(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_v + _w < 0) std::cout << "Error: Hoppe94LoopEval::eval_regular_2v1e(): _v + _w < 0, is " << _v << ", " << _w << ".\n";
	if (_v + _w > 1.00001) std::cout << "Error: Hoppe94LoopEval::eval_regular_2v1e(): _v + _w > 1, is " << _v << ", " << _w << ".\n";
	if (_c.size() != 12) std::cout << "Error: Hoppe94LoopEval::eval_regular_2v1e(): _c.size() != 12.\n";
	if (_h.size() != 15) std::cout << "Error: Hoppe94LoopEval::eval_regular_2v1e(): _h.size() != 15.\n";
	// 这情况下_c[3 and 6]是特征点(dart, crease, and corner), _h[5]是特征边crease edge.

	// 特殊处理
	if (_w < 1e-6) { //对_w = 0这个边的曲线evaluation
		enum vtype {SMOOTH = -1, DART, CREASE, CORNER};
		vtype vt3 = SMOOTH, vt6 = SMOOTH; //_c[3], _c[6]的类型.

		std::vector<int> fv_tmp3;
		if (_h[0]) fv_tmp3.push_back(0);
		if (_h[1]) fv_tmp3.push_back(1);
		if (_h[2]) fv_tmp3.push_back(2);
		if (_h[3]) fv_tmp3.push_back(4); 
		if (fv_tmp3.size() >= 2) { //_c[3] is corner, 
			vt3 = CORNER;
		} else if (fv_tmp3.size() == 1) { // _c[3] is crease, 
			vt3 = CREASE;
		} else if (fv_tmp3.size() == 0) { // _c[3] is dart
			vt3 = DART;
		} else {
			std::cout << "Error: Hoppe94LoopEval::eval_regular_2v1e(): ==0 _c[3]是SMOOTH, 做为非特征考虑的,不在这考虑.\n";
			return false;
		} 
		
		std::vector<int> fv_tmp6;
		if (_h[4]) fv_tmp6.push_back(2);
		if (_h[8]) fv_tmp6.push_back(5);
		if (_h[11]) fv_tmp6.push_back(9);
		if (_h[12]) fv_tmp6.push_back(10); 
		if (fv_tmp6.size() >= 2) { //_c[6] is corner, 
			vt6 = CORNER;
		} else if (fv_tmp6.size() == 1) { // _c[6] is crease, 
			vt6 = CREASE;
		} else if (fv_tmp6.size() == 0) { // _c[6] is dart
			vt6 = DART;
		} else {
			std::cout << "Error: Hoppe94LoopEval::eval_regular_2v1e(): ==0 _c[6]是SMOOTH, 做为非特征考虑的,不在这考虑.\n";
			return false;
		} 
		
		if (vt3 == CREASE && vt6 == CREASE) 
			return eval_curve_crease_crease(_c[fv_tmp3[0]], _c[3], _c[6], _c[fv_tmp6[0]], _v, _p);

		if (vt3 == CORNER && vt6 == CREASE) 
			return eval_curve_corner_crease(_c[3], _c[6], _c[fv_tmp6[0]], _v, _p);
		else if (vt3 == CREASE && vt6 == CORNER) 
			return eval_curve_corner_crease(_c[6], _c[3], _c[fv_tmp3[0]], 1 - _v, _p);

		if (vt3 == CORNER && vt6 == CORNER) 
			return eval_curve_corner_corner(_c[3], _c[6], _v, _p);
		
		if (_v < 1e-6 && _w < 1e-6) { //太接近_c[3]
			if (vt3 == CORNER) _p = _c[3];
			else if (vt3 == CREASE) _p = _c[3] * 0.66667 + (_c[6] + _c[fv_tmp3[0]]) * 0.166667;
			else if (vt3 == DART) _p = _c[3] * 0.5 + (_c[0] + _c[1] + _c[2] + _c[4] + _c[6] + _c[7]) * 1.0/12.0;
			else { std::cout << "Error: Hoppe94LoopEval::eval_regular_2v1e(): 类型出错.\n"; return false; }
			return true;
		}
		if (_v > 0.999999) { //太接近_c[6]
			if (vt6 == CORNER) _p = _c[6];
			else if (vt6 == CREASE) _p = _c[6] * 0.66667 + (_c[3] + _c[fv_tmp6[0]]) * 0.166667;
			else if (vt6 == DART) _p = _c[6] * 0.5 + (_c[2] + _c[3] + _c[5] + _c[7] + _c[9] + _c[10]) * 1.0/12.0;
			else { std::cout << "Error: Hoppe94LoopEval::eval_regular_2v1e(): 类型出错.\n"; return false; }
			return true; 
		}		
	}
	// 实在没办法了就出绝招 
	//DGP::change_slight_vw(_v, _w);
	
	std::vector<OpenMesh::Vec3d> c18(18, OpenMesh::Vec3d(0, 0, 0));//
	//eval_regular_2v1e_c18(_c, _h, c18);
	eval_regular_c18(_c, _h, c18);
	 
	// 下面分4中情况来考虑.
	if (_w > 0.5) { // low right side (1/4) can be evaluation
		int pick[12] = { 3, 4,    6, 7, 8,    10, 11, 12, 13,    15, 16, 17}; //pick 12 vertices from the 18 vertices, which created before.
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		v = _v * 2; w = (_w - 0.5) * 2; 
		return eval_regular_1in4_0v0e(c12, v, w, _p); 
	}

	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) { // low middle 1/4 
		int pick[12] = { 5, 2,    10, 6, 3,    15, 11, 7, 4,    16, 12, 8}; //pick 12 vertices from the 18 vertices, which created before.
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0)); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];// the same as:
		////v = (_v + _w - 0.5) * 2; w = (0.5 - _v) * 2; // u = (0.5 - _w) * 2;
		 
		std::vector<bool> h(15, false);
		h[2] = true; h[3] = true; 
		return eval_regular_1v0e(c12, h, v, w, _p); 
	}

	if (_v + _w <= 0.5) { // top 1/4
		int pick[12] = {0, 1,    2, 3, 4,    5, 6, 7, 8,    10, 11, 12}; 
		std::vector<OpenMesh::Vec3d> c12(12); // 
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		v = _v * 2; w = _w * 2;

		std::vector<bool> h(15, false);
		h[0] = _h[0]; h[1] = _h[1]; h[2] = _h[2]; h[3] = _h[3]; 
		h[5] = true; h[11] = true; 
		return eval_regular_2v1e(c12, h, v, w, _p);  // recursion
	}

	if (_v > 0.5) { // low left 1/4
		int pick[12] = {2, 3,    5, 6, 7,    9, 10, 11, 12,    14, 15, 16};
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		v = (_v - 0.5 ) * 2; w = _w * 2; 

		std::vector<bool> h(15, false);
		h[1] = true; h[5] = true;
		h[4] = _h[4]; h[8] = _h[8]; h[11] = _h[11]; h[12] = _h[12];  
		return eval_regular_2v1e(c12, h, v, w, _p);  // recursion
	} 
	
	return false;
}

void Hoppe94LoopEval::eval_regular_3v2e_c18(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, std::vector<OpenMesh::Vec3d> &_c18) {
	// one time subdivision: 12 -> 18 
	double d1_16 = 0.0625, d1_8 = 0.125, d3_8 = 0.375, d5_8 = 0.625;
	// first time subdiv, 12 -> 18 
	double A1[18][12] = {
		//0     1       2       3       4       5       6       7       8       9      10       11   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 0  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 1  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 2
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 3   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 4
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 5   
		{ 0,	0,		0,		0.5,	0,		0,		0.5,	0,		0,		0,		0,		0 }, // index 6
		{ 0,	0,		0,		0.5,	0,		0,		0,		0.5,	0,		0,		0,		0 }, // index 7  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 8
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 9   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 10
		{ 0,	0,		0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8,		0 }, // index 11  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 12
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 13  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 14
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 15  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 16 
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 17  
	};
	
	int n_count_feature_edges = 0, fv_tmp = -1;
	if (_h[0]) { // new control point 0.
		A1[0][0] = A1[0][3] = 0.5;  ++n_count_feature_edges; fv_tmp = 0;
	} else {
		A1[0][0] = A1[0][3] = d3_8;  A1[0][1] = A1[0][2] = d1_8;  
	}
	if (_h[1]) { // new control point 1.
		A1[1][1] = A1[1][3] = 0.5;   ++n_count_feature_edges; fv_tmp = 1;
	} else {
		A1[1][1] = A1[1][3] = d3_8;  A1[1][0] = A1[1][4] = d1_8;
	}
	if (_h[2]) { // new control point 2.
		A1[2][2] = A1[2][3] = 0.5;  ++n_count_feature_edges; fv_tmp = 2;
	} else {
		A1[2][2] = A1[2][3] = d3_8;  A1[2][0] = A1[2][6] = d1_8;
	}
	if (_h[3]) { // new control point 4.
		A1[4][3] = A1[4][4] = 0.5;   ++n_count_feature_edges; fv_tmp = 4;
	} else {
		A1[4][3] = A1[4][4] = d3_8;  A1[4][1] = A1[4][7] = d1_8;  
	}
	if (n_count_feature_edges == 0) { // crease // new control point 3.
		A1[3][3] = 0.75;  A1[3][6] = A1[3][7] = d1_8; 
	} else if (n_count_feature_edges >= 1) { // corner
		A1[3][3] = 1;
	} 
	// 
	n_count_feature_edges = 0; fv_tmp = -1;
	if (_h[4]) { // new control point 5.
		A1[5][2] = A1[5][6] = 0.5;  ++n_count_feature_edges; fv_tmp = 2;
	} else {
		A1[5][2] = A1[5][6] = d3_8;  A1[5][3] = A1[5][5] = d1_8;  
	}
	if (_h[8]) { // new control point 9.
		A1[9][5] = A1[9][6] = 0.5;  ++n_count_feature_edges; fv_tmp = 5;
	} else {
		A1[9][5] = A1[9][6] = d3_8;  A1[9][2] = A1[9][9] = d1_8;  
	}
	if (_h[11]) { // new control point 14.
		A1[14][6] = A1[14][9] = 0.5;  ++n_count_feature_edges; fv_tmp = 9;
	} else {
		A1[14][6] = A1[14][9] = d3_8;  A1[14][5] = A1[14][10] = d1_8;  
	}
	if (_h[12]) { // new control point 15.
		A1[15][6] = A1[15][10] = 0.5;  ++n_count_feature_edges; fv_tmp = 10;
	} else {
		A1[15][6] = A1[15][10] = d3_8;  A1[15][7] = A1[15][9] = d1_8;  
	}
	if (n_count_feature_edges == 0) { // dart // new control point 10.
		A1[10][6] = d5_8;  A1[10][2] = A1[10][3] = A1[10][5] = A1[10][7] = A1[10][9] = A1[10][10] = d1_16;
	} else if (n_count_feature_edges == 1) { // crease
		A1[10][6] = 0.75;  A1[10][fv_tmp] = A1[10][3] = d1_8; 
	} else { // corner
		A1[10][6] = 1;
	}
	// 
	n_count_feature_edges = 0; fv_tmp = -1;
	if (_h[7]) { // new control point 8.
		A1[8][4] = A1[8][7] = 0.5;  ++n_count_feature_edges; fv_tmp = 4;
	} else {
		A1[8][4] = A1[8][7] = d3_8;  A1[8][3] = A1[8][8] = d1_8;  
	}
	if (_h[10]) { // new control point 13.
		A1[13][7] = A1[13][8] = 0.5;  ++n_count_feature_edges; fv_tmp = 8;
	} else {
		A1[13][7] = A1[13][8] = d3_8;  A1[13][4] = A1[13][11] = d1_8;  
	}
	if (_h[13]) { // new control point 16.
		A1[16][7] = A1[16][10] = 0.5;  ++n_count_feature_edges; fv_tmp = 10;
	} else {
		A1[16][7] = A1[16][10] = d3_8;  A1[16][6] = A1[16][11] = d1_8;  
	}
	if (_h[14]) { // new control point 17.
		A1[17][7] = A1[17][11] = 0.5;  ++n_count_feature_edges; fv_tmp = 11;
	} else {
		A1[17][7] = A1[17][11] = d3_8;  A1[17][8] = A1[17][10] = d1_8;  
	}
	if (n_count_feature_edges == 0) { // dart // new control point 12.
		A1[12][7] = d5_8;  A1[12][3] = A1[12][4] = A1[12][6] = A1[12][8] = A1[12][10] = A1[12][11] = d1_16;
	} else if (n_count_feature_edges == 1) { // crease
		A1[12][7] = 0.75;  A1[12][fv_tmp] = A1[12][3] = d1_8; 
	} else { // corner
		A1[12][7] = 1;
	}

	for (int i = 0; i < 18; ++i) {
		_c18[i] = OpenMesh::Vec3d(0, 0, 0);
		for (int j = 0; j < 12; ++j) {
			_c18[i] += A1[i][j] * _c[j];
		} 
	}
}
bool Hoppe94LoopEval::eval_regular_3v2e(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_v + _w < 0 || _v + _w > 1) std::cout << "Error: Hoppe94LoopEval::eval_regular_3v2e(): _v + _w < 0 or > 1.\n";
	if (_c.size() != 12) std::cout << "Error: Hoppe94LoopEval::eval_regular_3v2e(): _c.size() != 12.\n";
	if (_h.size() != 15) std::cout << "Error: Hoppe94LoopEval::eval_regular_3v2e(): _h.size() != 15.\n";
	//std::cout << "Hoppe94LoopEval::eval_regular_3v2e()...\n";

	if (_v < 1e-6 && _w < 1e-6) { //太接近_c[3]
		std::vector<int > fv_tmp;
		if (_h[0]) fv_tmp.push_back(0);
		if (_h[1]) fv_tmp.push_back(1);
		if (_h[2]) fv_tmp.push_back(2);
		if (_h[3]) fv_tmp.push_back(4); 
		if (fv_tmp.size() == 1) { 
			_p = _c[3]; return true;
		} else if (fv_tmp.size() == 0) {
			_p = _c[3] * 0.66667 + (_c[6] + _c[7]) * 0.166667; return true;
		} else { 
			std::cout << "Error: Hoppe94LoopEval::eval_regular_3v2e(): 类型错误.\n"; return false; 
		}
	}
	// 特殊处理	
	//DGP::change_slight_vw(_v, _w);
	// -----
	std::vector<OpenMesh::Vec3d> c18(18, OpenMesh::Vec3d(0, 0, 0));//
	//eval_regular_3v2e_c18(_c, _h, c18);
	eval_regular_c18(_c, _h, c18);

	// 下面也是分四种情况来做.
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) { // low middle 1/4 
		int pick[12] = { 4, 8,    3, 7, 12,    2, 6, 11, 16,    5, 10, 15}; //pick 12 vertices from the 18 vertices, which created before.
		std::vector<OpenMesh::Vec3d> c12; //
		for (int i = 0; i < 12; ++i) {
			c12.push_back(c18[pick[i]]);
			//std::cout << c12[i] << ".\n";
		}
		double v = 0, w = 0;
		//OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));
		//v = bc[1]; w = bc[2];// the same as:
		v = (0.5 - _w) * 2; w = (_v + _w - 0.5) * 2; // u = (0.5 - _v) * 2;		
		
		std::vector<bool> h(15, false);
		h[2] = true; h[3] = true; h[4] = true; h[12] = true;
		return eval_regular_2v0e(c12, h, v, w, _p); 
	} 

	if (_v > 0.5) { // low left 1/4
		int pick[12] = {2, 3,    5, 6, 7,    9, 10, 11, 12,    14, 15, 16};
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		v = (_v - 0.5 ) * 2; w = _w * 2; 

		std::vector<bool> h(15, false);
		h[1] = true; 
		h[4] = _h[4]; h[5] = true; h[8] = _h[8]; h[11] = _h[11]; h[12] = _h[12]; 
		return eval_regular_2v1e(c12, h, v, w, _p);  // recursion
	} 

	if (_w > 0.5) { // low right side (1/4) can be evaluation
		int pick[12] = { 13, 17,    8, 12, 16,    4, 7, 11, 15,    3, 6, 10}; //pick 12 vertices from the 18 vertices, which created before.
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];

		std::vector<bool> h(15, false);
		h[0] = _h[10]; h[1] = _h[14];
		h[2] = _h[7]; h[3] = _h[13];
		h[5] = true; h[11] = true;
		return eval_regular_2v1e(c12, h, v, w, _p); 
	}

	if (_v + _w <= 0.5) { // top 1/4
		int pick[12] = {0, 1,    2, 3, 4,    5, 6, 7, 8,    10, 11, 12}; 
		std::vector<OpenMesh::Vec3d> c12(12); // 
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		v = _v * 2; w = _w * 2;

		std::vector<bool> h(15, false);
		h[0] = _h[0]; h[1] = _h[1]; h[2] = _h[2]; h[3] = _h[3]; 
		h[5] = true; h[6] = true; h[11] = true; h[14] = true;
		return eval_regular_3v2e(c12, h, v, w, _p);  // recursion
	}
	return false;
}
void Hoppe94LoopEval::eval_regular_3v3e_c18(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, std::vector<OpenMesh::Vec3d> &_c18) {
	// one time subdivision: 12 -> 18 
	double d1_16 = 0.0625, d1_8 = 0.125, d3_8 = 0.375, d5_8 = 0.625;
	// first time subdiv, 12 -> 18 
	double A1[18][12] = {
		//0     1       2       3       4       5       6       7       8       9      10       11   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 0  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 1  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 2
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 3   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 4
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 5   
		{ 0,	0,		0,		0.5,	0,		0,		0.5,	0,		0,		0,		0,		0 }, // index 6
		{ 0,	0,		0,		0.5,	0,		0,		0,		0.5,	0,		0,		0,		0 }, // index 7  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 8
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 9   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 10
		{ 0,	0,		0,		0,		0,		0,		0.5,	0.5,	0,		0,		0,		0 }, // index 11  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 12
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 13  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 14
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 15  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 16 
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 17  
	};

	int n_count_feature_edges = 0, fv_tmp = -1;
	if (_h[0]) { // new control point 0.
		A1[0][0] = A1[0][3] = 0.5;  ++n_count_feature_edges; fv_tmp = 0;
	} else {
		A1[0][0] = A1[0][3] = d3_8;  A1[0][1] = A1[0][2] = d1_8;  
	}
	if (_h[1]) { // new control point 1.
		A1[1][1] = A1[1][3] = 0.5;   ++n_count_feature_edges; fv_tmp = 1;
	} else {
		A1[1][1] = A1[1][3] = d3_8;  A1[1][0] = A1[1][4] = d1_8;
	}
	if (_h[2]) { // new control point 2.
		A1[2][2] = A1[2][3] = 0.5;  ++n_count_feature_edges; fv_tmp = 2;
	} else {
		A1[2][2] = A1[2][3] = d3_8;  A1[2][0] = A1[2][6] = d1_8;
	}
	if (_h[3]) { // new control point 4.
		A1[4][3] = A1[4][4] = 0.5;   ++n_count_feature_edges; fv_tmp = 4;
	} else {
		A1[4][3] = A1[4][4] = d3_8;  A1[4][1] = A1[4][7] = d1_8;  
	}
	if (n_count_feature_edges == 0) { // crease // new control point 3.
		A1[3][3] = 0.75;  A1[3][6] = A1[3][7] = d1_8; 
	} else if (n_count_feature_edges >= 1) { // corner
		A1[3][3] = 1;
	} 
	// 
	n_count_feature_edges = 0; fv_tmp = -1;
	if (_h[4]) { // new control point 5.
		A1[5][2] = A1[5][6] = 0.5;  ++n_count_feature_edges; fv_tmp = 2;
	} else {
		A1[5][2] = A1[5][6] = d3_8;  A1[5][3] = A1[5][5] = d1_8;  
	}
	if (_h[8]) { // new control point 9.
		A1[9][5] = A1[9][6] = 0.5;  ++n_count_feature_edges; fv_tmp = 5;
	} else {
		A1[9][5] = A1[9][6] = d3_8;  A1[9][2] = A1[9][9] = d1_8;  
	}
	if (_h[11]) { // new control point 14.
		A1[14][6] = A1[14][9] = 0.5;  ++n_count_feature_edges; fv_tmp = 9;
	} else {
		A1[14][6] = A1[14][9] = d3_8;  A1[14][5] = A1[14][10] = d1_8;  
	}
	if (_h[12]) { // new control point 15.
		A1[15][6] = A1[15][10] = 0.5;  ++n_count_feature_edges; fv_tmp = 10;
	} else {
		A1[15][6] = A1[15][10] = d3_8;  A1[15][7] = A1[15][9] = d1_8;  
	}
	if (n_count_feature_edges == 0) { // crease // new control point 10.
		A1[10][6] = 0.75;  A1[10][3] = A1[10][7] = d1_8; 
	} else if (n_count_feature_edges >= 1) {  // corner
		A1[10][6] = 1;
	} 
	// 
	n_count_feature_edges = 0; fv_tmp = -1;
	if (_h[7]) { // new control point 8.
		A1[8][4] = A1[8][7] = 0.5;  ++n_count_feature_edges; fv_tmp = 4;
	} else {
		A1[8][4] = A1[8][7] = d3_8;  A1[8][3] = A1[8][8] = d1_8;  
	}
	if (_h[10]) { // new control point 13.
		A1[13][7] = A1[13][8] = 0.5;  ++n_count_feature_edges; fv_tmp = 8;
	} else {
		A1[13][7] = A1[13][8] = d3_8;  A1[13][4] = A1[13][11] = d1_8;  
	}
	if (_h[13]) { // new control point 16.
		A1[16][7] = A1[16][10] = 0.5;  ++n_count_feature_edges; fv_tmp = 10;
	} else {
		A1[16][7] = A1[16][10] = d3_8;  A1[16][6] = A1[16][11] = d1_8;  
	}
	if (_h[14]) { // new control point 17.
		A1[17][7] = A1[17][11] = 0.5;  ++n_count_feature_edges; fv_tmp = 11;
	} else {
		A1[17][7] = A1[17][11] = d3_8;  A1[17][8] = A1[17][10] = d1_8;  
	}
	if (n_count_feature_edges == 0) { // crease // new control point 12.
		A1[12][7] = 0.75;  A1[12][6] = A1[12][3] = d1_8;   
	} else if (n_count_feature_edges >= 1) { // corner
		A1[12][7] = 1;
	}

	for (int i = 0; i < 18; ++i) {
		_c18[i] = OpenMesh::Vec3d(0, 0, 0);
		for (int j = 0; j < 12; ++j) {
			_c18[i] += A1[i][j] * _c[j];
		} 
	}
}
bool Hoppe94LoopEval::eval_regular_3v3e(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_v + _w < 0 || _v + _w > 1) std::cout << "Error: Hoppe94LoopEval::eval_regular_3v3e(): _v + _w < 0 or > 1.\n";
	if (_c.size() != 12) std::cout << "Error: Hoppe94LoopEval::eval_regular_3v3e(): _c.size() != 12.\n";
	if (_h.size() != 15) std::cout << "Error: Hoppe94LoopEval::eval_regular_3v3e(): _h.size() != 15.\n";

	std::vector<OpenMesh::Vec3d> c18(18, OpenMesh::Vec3d(0, 0, 0));//
	eval_regular_3v3e_c18(_c, _h, c18);

	// 下面也是分四种情况来做.
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) { // low middle 1/4 
		int pick[12] = { 4, 8,    3, 7, 12,    2, 6, 11, 16,    5, 10, 15}; //pick 12 vertices from the 18 vertices, which created before.
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];// the same as:
		////v = (0.5 - _w) * 2; w = (_v + _w - 0.5) * 2; // u = (0.5 - _v) * 2;		

		std::vector<bool> h(15, false);
		h[2] = true; h[3] = true; h[4] = true; h[7] = true; h[12] = true; h[13] = true;
		return eval_regular_3v0e(c12, h, v, w, _p); 
	} 

	if (_v + _w <= 0.5) { // top 1/4
		int pick[12] = {0, 1,    2, 3, 4,    5, 6, 7, 8,    10, 11, 12}; 
		std::vector<OpenMesh::Vec3d> c12(12); // 
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		v = _v * 2; w = _w * 2;

		std::vector<bool> h(15, false);
		h[0] = _h[0]; h[1] = _h[1]; h[2] = _h[2]; h[3] = _h[3]; 
		h[5] = true; h[6] = true; h[11] = true; h[14] = true;
		return eval_regular_3v2e(c12, h, v, w, _p);  //  
	}

	if (_v > 0.5) { // low left 1/4
		int pick[12] = {14, 9,    15, 10, 5,     16, 11, 6, 2,    12, 7, 3};
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(1, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];

		std::vector<bool> h(15, false);
		h[0] = _h[11]; h[1] = _h[8]; h[2] = _h[12]; h[3] = _h[4]; 
		h[5] = true; h[6] = true; h[11] = true; h[14] = true;
		return eval_regular_3v2e(c12, h, v, w, _p);  //  
	} 

	if (_w > 0.5) { // low right side (1/4) can be evaluation
		int pick[12] = { 13, 17,    8, 12, 16,    4, 7, 11, 15,    3, 6, 10}; //pick 12 vertices from the 18 vertices, which created before.
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];

		std::vector<bool> h(15, false);
		h[0] = _h[10]; h[1] = _h[14]; h[2] = _h[7]; h[3] = _h[13]; 
		h[5] = true; h[6] = true; h[11] = true; h[14] = true;
		return eval_regular_3v2e(c12, h, v, w, _p);  //  
	}

	
	return false;
}

void Hoppe94LoopEval::eval_regular_3v1e_c18(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, std::vector<OpenMesh::Vec3d> &_c18) {
	// one time subdivision: 12 -> 18 
	double d1_16 = 0.0625, d1_8 = 0.125, d3_8 = 0.375, d5_8 = 0.625;
	// first time subdiv, 12 -> 18 
	double A1[18][12] = {
		//0     1       2       3       4       5       6       7       8       9      10       11   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 0  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 1  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 2
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 3   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 4
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 5   
		{ 0,	0,		0,		0.5,	0,		0,		0.5,	0,		0,		0,		0,		0 }, // index 6
		{ 0,	0,		0,	d3_8,	d1_8,		0,	d1_8,	d3_8,		0,		0,		0,		0 }, // index 7  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 8
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 9   
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 10
		{ 0,	0,		0,	d1_8,		0,		0,	d3_8,	d3_8,		0,		0,	d1_8,		0 }, // index 11  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 12
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 13  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 14
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 15  
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 16 
		{ 0,	0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0 }, // index 17  
	};

	int n_count_feature_edges = 0, fv_tmp = -1;
	if (_h[0]) { // new control point 0.
		A1[0][0] = A1[0][3] = 0.5;  ++n_count_feature_edges; fv_tmp = 0;
	} else {
		A1[0][0] = A1[0][3] = d3_8;  A1[0][1] = A1[0][2] = d1_8;  
	}
	if (_h[1]) { // new control point 1.
		A1[1][1] = A1[1][3] = 0.5;   ++n_count_feature_edges; fv_tmp = 1;
	} else {
		A1[1][1] = A1[1][3] = d3_8;  A1[1][0] = A1[1][4] = d1_8;
	}
	if (_h[2]) { // new control point 2.
		A1[2][2] = A1[2][3] = 0.5;  ++n_count_feature_edges; fv_tmp = 2;
	} else {
		A1[2][2] = A1[2][3] = d3_8;  A1[2][0] = A1[2][6] = d1_8;
	}
	if (_h[3]) { // new control point 4.
		A1[4][3] = A1[4][4] = 0.5;   ++n_count_feature_edges; fv_tmp = 4;
	} else {
		A1[4][3] = A1[4][4] = d3_8;  A1[4][1] = A1[4][7] = d1_8;  
	}
	if (n_count_feature_edges == 0) { // new control point 3.
		A1[3][3] = d5_8;   A1[3][0] =  A1[3][1] =  A1[3][2] =  A1[3][4] =  A1[3][6] =  A1[3][7] = d1_16;   
	} else if (n_count_feature_edges == 1) { // crease
		A1[3][3] = 0.75;  A1[3][6] = A1[3][fv_tmp] = d1_8; 
	} else { // corner
		A1[3][3] = 1;
	}
	// 
	n_count_feature_edges = 0; fv_tmp = -1;
	if (_h[4]) { // new control point 5.
		A1[5][2] = A1[5][6] = 0.5;  ++n_count_feature_edges; fv_tmp = 2;
	} else {
		A1[5][2] = A1[5][6] = d3_8;  A1[5][3] = A1[5][5] = d1_8;  
	}
	if (_h[8]) { // new control point 9.
		A1[9][5] = A1[9][6] = 0.5;  ++n_count_feature_edges; fv_tmp = 5;
	} else {
		A1[9][5] = A1[9][6] = d3_8;  A1[9][2] = A1[9][9] = d1_8;  
	}
	if (_h[11]) { // new control point 14.
		A1[14][6] = A1[14][9] = 0.5;  ++n_count_feature_edges; fv_tmp = 9;
	} else {
		A1[14][6] = A1[14][9] = d3_8;  A1[14][5] = A1[14][10] = d1_8;  
	}
	if (_h[12]) { // new control point 15.
		A1[15][6] = A1[15][10] = 0.5;  ++n_count_feature_edges; fv_tmp = 10;
	} else {
		A1[15][6] = A1[15][10] = d3_8;  A1[15][7] = A1[15][9] = d1_8;  
	}
	if (n_count_feature_edges == 0) { // dart // new control point 10.
		A1[10][6] = d5_8;  A1[10][2] = A1[10][3] = A1[10][5] = A1[10][7] = A1[10][9] = A1[10][10] = d1_16;
	} else if (n_count_feature_edges == 1) { // crease
		A1[10][6] = 0.75;  A1[10][fv_tmp] = A1[10][3] = d1_8; 
	} else { // corner
		A1[10][6] = 1;
	}
	// 
	n_count_feature_edges = 0; 
	std::vector<int > fv_tmp_array;
	if (_h[7]) { // new control point 8.
		A1[8][4] = A1[8][7] = 0.5;  ++n_count_feature_edges; fv_tmp_array.push_back(4);
	} else {
		A1[8][4] = A1[8][7] = d3_8;  A1[8][3] = A1[8][8] = d1_8;  
	}
	if (_h[10]) { // new control point 13.
		A1[13][7] = A1[13][8] = 0.5;  ++n_count_feature_edges; fv_tmp_array.push_back(8);
	} else {
		A1[13][7] = A1[13][8] = d3_8;  A1[13][4] = A1[13][11] = d1_8;  
	}
	if (_h[13]) { // new control point 16.
		A1[16][7] = A1[16][10] = 0.5;  ++n_count_feature_edges; fv_tmp_array.push_back(10);
	} else {
		A1[16][7] = A1[16][10] = d3_8;  A1[16][6] = A1[16][11] = d1_8;  
	}
	if (_h[14]) { // new control point 17.
		A1[17][7] = A1[17][11] = 0.5;  ++n_count_feature_edges; fv_tmp_array.push_back(11);
	} else {
		A1[17][7] = A1[17][11] = d3_8;  A1[17][8] = A1[17][10] = d1_8;  
	}
	if (n_count_feature_edges <= 1) { // dart // new control point 12.
		A1[12][7] = d5_8;  A1[12][3] = A1[12][4] = A1[12][6] = A1[12][8] = A1[12][10] = A1[12][11] = d1_16;
	} else if (n_count_feature_edges == 2) { // crease
		A1[12][7] = 0.75;  A1[12][fv_tmp_array[0]] = A1[12][fv_tmp_array[1]] = d1_8; 
	} else { // corner
		A1[12][7] = 1;
	}

	for (int i = 0; i < 18; ++i) {
		_c18[i] = OpenMesh::Vec3d(0, 0, 0);
		for (int j = 0; j < 12; ++j) {
			_c18[i] += A1[i][j] * _c[j];
		} 
	}
}
bool Hoppe94LoopEval::eval_regular_3v1e(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_v + _w < 0 || _v + _w > 1) std::cout << "Error: Hoppe94LoopEval::eval_regular_3v1e(): _v + _w < 0 or > 1.\n";
	if (_c.size() != 12) std::cout << "Error: Hoppe94LoopEval::eval_regular_3v1e(): _c.size() != 12.\n";
	if (_h.size() != 15) std::cout << "Error: Hoppe94LoopEval::eval_regular_3v1e(): _h.size() != 15.\n";
	//std::cout << "Hoppe94LoopEval::eval_regular_3v1e()...\n";

	std::vector<OpenMesh::Vec3d> c18(18, OpenMesh::Vec3d(0, 0, 0));//
	eval_regular_3v1e_c18(_c, _h, c18);

	// 下面也是分四种情况来做.
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) { // low middle 1/4 
		int pick[12] = { 5, 2,    10, 6, 3,     15, 11, 7, 4,    16, 12, 8}; //pick 12 vertices from the 18 vertices, which created before.
		std::vector<OpenMesh::Vec3d> c12; //
		for (int i = 0; i < 12; ++i) {
			c12.push_back(c18[pick[i]]);
			//std::cout << c12[i] << ".\n";
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];// the same as:	

		std::vector<bool> h(15, false);
		h[2] = true; h[3] = true;
		return eval_regular_1v0e(c12, h, v, w, _p); 
	} 

	if (_v > 0.5) { // low left 1/4
		int pick[12] = {2, 3,    5, 6, 7,    9, 10, 11, 12,    14, 15, 16};
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		v = (_v - 0.5 ) * 2; w = _w * 2; 

		std::vector<bool> h(15, false);
		h[1] = true; 
		h[4] = _h[4]; h[5] = true; h[8] = _h[8]; h[11] = _h[11]; h[12] = _h[12]; 
		return eval_regular_2v1e(c12, h, v, w, _p);  // recursion
	} 

	if (_w > 0.5) { // low right side (1/4) can be evaluation
		int pick[12] = { 13, 17,    8, 12, 16,    4, 7, 11, 15,    3, 6, 10}; //pick 12 vertices from the 18 vertices, which created before.
		std::vector<OpenMesh::Vec3d> c12(12); //
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];

		std::vector<bool> h(15, false);
		h[0] = _h[10]; h[1] = _h[14];
		h[2] = _h[7]; h[3] = _h[13];
		return eval_regular_1v0e(c12, h, v, w, _p); 
	}

	if (_v + _w <= 0.5) { // top 1/4
		int pick[12] = {0, 1,    2, 3, 4,    5, 6, 7, 8,    10, 11, 12}; 
		std::vector<OpenMesh::Vec3d> c12(12); // 
		for (int i = 0; i < 12; ++i) {
			c12[i] = c18[pick[i]];
		}
		double v = 0, w = 0;
		v = _v * 2; w = _w * 2;

		std::vector<bool> h(15, false);
		h[0] = _h[0]; h[1] = _h[1]; h[2] = _h[2]; h[3] = _h[3]; 
		h[5] = true;  h[11] = true;  
		return eval_regular_2v1e(c12, h, v, w, _p);  // recursion
	}
	return false;
}

void Hoppe94LoopEval::eval_regular_c18(const std::vector<OpenMesh::Vec3d> &_c, const std::vector<bool> &_h, std::vector<OpenMesh::Vec3d> &_c18) {
	// one time subdivision: 12 -> 18 
	double d1_16 = 0.0625, d1_8 = 0.125, d3_8 = 0.375, d5_8 = 0.625;
	// first time subdiv, 12 -> 18 
	double A[18][12] = { 0 }; // 每一行对应一个新生成的顶点.
	
	int n_count_feature_edges = 0;
	std::vector<int >fv_tmp;
	if (_h[0]) { // new control point 0.
		A[0][0] = A[0][3] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(0);
	} else {
		A[0][0] = A[0][3] = d3_8;  A[0][1] = A[0][2] = d1_8;  
	}
	if (_h[1]) { // new control point 1.
		A[1][1] = A[1][3] = 0.5;   ++n_count_feature_edges; fv_tmp.push_back(1);
	} else {
		A[1][1] = A[1][3] = d3_8;  A[1][0] = A[1][4] = d1_8;
	}
	if (_h[2]) { // new control point 2.
		A[2][2] = A[2][3] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(2);
	} else {
		A[2][2] = A[2][3] = d3_8;  A[2][0] = A[2][6] = d1_8;
	}
	if (_h[3]) { // new control point 4.
		A[4][3] = A[4][4] = 0.5;   ++n_count_feature_edges; fv_tmp.push_back(4);
	} else {
		A[4][3] = A[4][4] = d3_8;  A[4][1] = A[4][7] = d1_8;  
	}
	if (_h[5]) {
		++n_count_feature_edges; fv_tmp.push_back(6);
	}
	if (_h[6]) {
		++n_count_feature_edges; fv_tmp.push_back(7);
	}
	if (n_count_feature_edges <= 1) { // new control point 3.
		A[3][3] = d5_8;   A[3][0] =  A[3][1] =  A[3][2] =  A[3][4] =  A[3][6] =  A[3][7] = d1_16;   
	} else if (n_count_feature_edges == 2) { // crease
		A[3][3] = 0.75;  A[3][fv_tmp[0]] = A[3][fv_tmp[1]] = d1_8; 
	} else { // corner
		A[3][3] = 1;
	} 
	// 
	n_count_feature_edges = 0; fv_tmp.clear();
	if (_h[4]) { // new control point 5.
		A[5][2] = A[5][6] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(2);
	} else {
		A[5][2] = A[5][6] = d3_8;  A[5][3] = A[5][5] = d1_8;  
	}
	if (_h[5]) {
		A[6][3] = A[6][6] = 0.5; ++n_count_feature_edges; fv_tmp.push_back(3);
	} else {
		A[6][3] = A[6][6] = d3_8; A[6][2] = A[6][7] = d1_8; 
	}
	if (_h[8]) { // new control point 9.
		A[9][5] = A[9][6] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(5);
	} else {
		A[9][5] = A[9][6] = d3_8;  A[9][2] = A[9][9] = d1_8;  
	}
	if (_h[9]) {
		A[11][6] = A[11][7] = 0.5; ++n_count_feature_edges; fv_tmp.push_back(7);
	} else {
		A[11][6] = A[11][7] = d3_8; A[11][3] = A[11][10] = d1_8;   
	}
	if (_h[11]) { // new control point 14.
		A[14][6] = A[14][9] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(9);
	} else {
		A[14][6] = A[14][9] = d3_8;  A[14][5] = A[14][10] = d1_8;  
	}
	if (_h[12]) { // new control point 15.
		A[15][6] = A[15][10] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(10);
	} else {
		A[15][6] = A[15][10] = d3_8;  A[15][7] = A[15][9] = d1_8;  
	}
	if (n_count_feature_edges <= 1) { // dart // new control point 10.
		A[10][6] = d5_8;  A[10][2] = A[10][3] = A[10][5] = A[10][7] = A[10][9] = A[10][10] = d1_16;
	} else if (n_count_feature_edges == 2) { // crease
		A[10][6] = 0.75;  A[10][fv_tmp[0]] = A[10][fv_tmp[1]] = d1_8; 
	} else { // corner
		A[10][6] = 1;
	}
	// 
	n_count_feature_edges = 0; 
	fv_tmp.clear();
	if (_h[6]) {
		A[7][3] = A[7][7] = 0.5; ++n_count_feature_edges; fv_tmp.push_back(3);
	} else {
		A[7][3] = A[7][7] = d3_8; A[7][4] = A[7][6] = d1_8;  
	}
	if (_h[7]) { // new control point 8.
		A[8][4] = A[8][7] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(4);
	} else {
		A[8][4] = A[8][7] = d3_8;  A[8][3] = A[8][8] = d1_8;  
	}
	if (_h[9]) {
		++ n_count_feature_edges; fv_tmp.push_back(6);
	}
	if (_h[10]) { // new control point 13.
		A[13][7] = A[13][8] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(8);
	} else {
		A[13][7] = A[13][8] = d3_8;  A[13][4] = A[13][11] = d1_8;  
	}
	if (_h[13]) { // new control point 16.
		A[16][7] = A[16][10] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(10);
	} else {
		A[16][7] = A[16][10] = d3_8;  A[16][6] = A[16][11] = d1_8;  
	}
	if (_h[14]) { // new control point 17.
		A[17][7] = A[17][11] = 0.5;  ++n_count_feature_edges; fv_tmp.push_back(11);
	} else {
		A[17][7] = A[17][11] = d3_8;  A[17][8] = A[17][10] = d1_8;  
	}
	if (n_count_feature_edges <= 1) { // dart // new control point 12.
		A[12][7] = d5_8;  A[12][3] = A[12][4] = A[12][6] = A[12][8] = A[12][10] = A[12][11] = d1_16;
	} else if (n_count_feature_edges == 2) { // crease
		A[12][7] = 0.75;  A[12][fv_tmp[0]] = A[12][fv_tmp[1]] = d1_8; 
	} else { // corner
		A[12][7] = 1;
	}

	for (int i = 0; i < 18; ++i) {
		_c18[i] = OpenMesh::Vec3d(0, 0, 0);
		for (int j = 0; j < 12; ++j) {
			_c18[i] += A[i][j] * _c[j];
		} 
	}
}

// =================================================

// -----------
void Hoppe94LoopEval::eval_irregular_cN12(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, std::vector<OpenMesh::Vec3d> & _cN12) {
	int N = _cN6.size() - 6;
	if (N != _hN9.size() - 9 || N != _cN12.size() - 12) std::cout << "Error: Hoppe94LoopEval::eval_irregular_cN12(), size wrong.\n";

	std::vector<std::vector<double> > A(_cN12.size(), std::vector<double>(_cN6.size(), 0));
	double d1_8 = 0.125, d3_8 = 0.375, d5_8 = 0.625, d1_16 = 0.0625;

	if (_hN9[0]) {
		A[1][0] = A[1][1] = 0.5;
	} else {
		A[1][0] = A[1][1] = d3_8;  A[1][N] = A[1][2] = d1_8;  
	}
	if (_hN9[N-1]) {
		A[N][0] = A[N][N] = 0.5;
	} else {
		A[N][0] = A[N][N] = d3_8;  A[N][N-1] = A[N][1] = d1_8;  
	}
	for (int i = 1; i <= N -2; ++i) {
		if (_hN9[i]) {
			A[i+1][0] =  A[i+1][i+1] = 0.5;  
		} else {
			A[i+1][0] =  A[i+1][i+1] = d3_8;  A[i+1][i] =  A[i+1][i+2] = d1_8;    
		}
	}
	std::vector<int > vec;
	for (int i = 0; i < N; ++i) {
		if (_hN9[i]) vec.push_back(i+1); 
	}
	if (vec.size() >= 3) {
		A[0][0] = 1;
	} else if (vec.size() == 2) {
		A[0][0] = 0.75; A[0][vec[0]] = A[0][vec[1]] = d1_8; 
	} else { // smooth & dart
		double beta = 0;
		if (N == 3) {
			beta = 3.0 / 16.0;
		} else {
			beta = (1.0/(N*64.0)) * (40 - (3 + 2 * cos(2*M_PI/N)) * (3 + 2 * cos(2*M_PI/N)));  
		}
		A[0][0] = 1 - N * beta;
		for (int i = 1; i <= N; ++i) {
			A[0][i] = beta; 
		} 
	}

	int count = 0;
	vec.clear();
	if (_hN9[N]) {
		++ count; vec.push_back(N);
		A[N+1][1] = A[N+1][N] = 0.5; 
	} else {
		A[N+1][1] = A[N+1][N] = d3_8;    A[N+1][0] = A[N+1][N+1] = d1_8;  
	} 
	if (_hN9[N+1]) {
		++ count; vec.push_back(N+1);
		A[N+6][1] = A[N+6][N+1] = 0.5;
	} else {
		A[N+6][1] = A[N+6][N+1] = d3_8;  A[N+6][N] = A[N+6][N+2] = d1_8;  
	}
	if (_hN9[N+2]) {
		++ count; vec.push_back(N+2);
		A[N+7][1] = A[N+7][N+2] = 0.5;
	} else {
		A[N+7][1] = A[N+7][N+2] = d3_8;  A[N+7][N+1] = A[N+7][N+3] = d1_8;  
	}
	if (_hN9[N+3]) {
		++ count; vec.push_back(N+3);
		A[N+8][1] = A[N+8][N+3] = 0.5;
	} else {
		A[N+8][1] = A[N+8][N+3] = d3_8;  A[N+8][2] = A[N+8][N+2] = d1_8;  
	}
	if (_hN9[N+4]) {
		++ count; vec.push_back(2);
		A[N+3][1] = A[N+3][2] = 0.5;
	} else {
		A[N+3][1] = A[N+3][2] = d3_8;  A[N+3][0] = A[N+3][N+3] = d1_8;  
	}
	if (_hN9[0]) {
		++ count; vec.push_back(0);
	} 
	if (count >= 3) {
		A[N+2][1] = 1;
	} else if (count == 2) {
		A[N+2][1] = 0.75; A[N+2][vec[0]] = A[N+2][vec[1]] = d1_8;
	} else {
		A[N+2][1] = d5_8;  A[N+2][0] = A[N+2][2] = A[N+2][N] = A[N+2][N+1] = A[N+2][N+2] = A[N+2][N+3] = d1_16;
	}

	count = 0; vec.clear();
	if (_hN9[N+5]) {
		++count; vec.push_back(N+1);
		A[N+9][N] = A[N+9][N+1] = 0.5; 
	} else {
		A[N+9][N] = A[N+9][N+1] = d3_8;  A[N+9][1] = A[N+9][N+4] = d1_8;   
	}
	if (_hN9[N+6]) {
		++count; vec.push_back(N+4);
		A[N+10][N] = A[N+10][N+4] = 0.5; 
	} else {
		A[N+10][N] = A[N+10][N+4] = d3_8;  A[N+10][N+1] = A[N+10][N+5] = d1_8;   
	}
	if (_hN9[N+7]) {
		++count; vec.push_back(N+5);
		A[N+11][N] = A[N+11][N+5] = 0.5; 
	} else {
		A[N+11][N] = A[N+11][N+5] = d3_8;  A[N+11][N-1] = A[N+11][N+4] = d1_8;   
	}
	if (_hN9[N+8]) {
		++count; vec.push_back(N-1);
		A[N+5][N] = A[N+5][N-1] = 0.5; 
	} else {
		A[N+5][N] = A[N+5][N-1] = d3_8;  A[N+5][0] = A[N+5][N+5] = d1_8;   
	}
	if (_hN9[N-1]) {
		++count; vec.push_back(0);
	}
	if (_hN9[N]) {
		++count; vec.push_back(1);
	}

	if (count >= 3) {
		A[N+4][N] = 1;
	} else if (count == 2) {
		A[N+4][N] = 0.75; A[N+4][vec[0]] = A[N+4][vec[1]] = d1_8;
	} else {
		A[N+4][N] = d5_8;  A[N+4][0] = A[N+4][1] = A[N+4][N+1] = A[N+4][N+4] = A[N+4][N+5] = A[N+4][N-1] = d1_16;
	}

	for (int i = 0; i < N+12; ++i) {
		for (int j = 0; j < N+6; ++j) {
			_cN12[i] += A[i][j] * _cN6[j];			
		}		//std::cout << _cN12[i] << ".\n";
	}
}
bool Hoppe94LoopEval::eval_irregular_1v0e_up(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_1v0e_up(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_1v0e_up(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}

	int N = _cN6.size() - 6;// valence of the extraordinary vertex.
	//特殊处理这个特征点.
	if (_v < 1e-6 && _w < 1e-6) {
		std::vector<int> fv_tmp;
		for (int i = 1; i <= N-2; ++i) {
			if (_hN9[i]) fv_tmp.push_back(i+1);
		}
		if (fv_tmp.size() >= 3) {
			_p = _cN6[0]; return true;
		} else if (fv_tmp.size() == 2) {
			_p = _cN6[0] * 0.66667 + (_cN6[fv_tmp[0]] + _cN6[fv_tmp[1]]) * 0.166667; return true;
		} else if (fv_tmp.size() == 1) { //dart
			int k = N;
			double inv_v = 1.0/double(k); double t = (3.0 + 2.0 * cos( 2.0 * M_PI * inv_v) ); // double beta = inv_v * (40.0 - t * t)/64.0;
			double k_alpha = 1.0 / ( 24.0/(40.0-t*t) + 1); double alpha = inv_v * k_alpha;

			OpenMesh::Vec3d limit_pos = OpenMesh::Vec3d(0, 0, 0);
			for (int i = 1; i <= N; ++i) {
				limit_pos += _cN6[i] * alpha; 
			}
			limit_pos += _cN6[0] * (1.0 - k_alpha);
			_p = limit_pos; return true;
		} else if (fv_tmp.size() == 0) {
			std::cout << "Error: Hoppe94LoopEval::eval_irregular_1v0e_up(), 应该在eval_irregular_1v0e_left/right()中处理吧.\n"; return false;
		} 
	}
	
	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------

	if (_v + _w  <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6);
		std::vector<bool> hN9_pick(N+9, false);
		for (int i = 0; i <= N-1; ++i) {
			hN9_pick[i] = _hN9[i];
		}
		hN9_pick[N+2] = _hN9[0];
		hN9_pick[N+6] = _hN9[N-1];
		return eval_irregular_1v0e_up(cN6_pick, hN9_pick, _v *2, _w *2, _p);
	}  
	if (_v > 0.5) {
		int pick[12] = {2, 0,   N+3, 1, N,   N+8, N+2, N+1, N+4,    N+7, N+6, N+9};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		return eval_regular_1in4_0v0e(c12, (_v - 0.5) * 2, _w * 2, _p); 
	} 
	if (_w > 0.5) {
		int pick[12] = {0, N-1,   1, N, N+5,   N+2, N+1, N+4, N+11,   N+6, N+9, N+10};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		return eval_regular_1in4_0v0e(c12, _v  * 2, (_w - 0.5)* 2, _p); 
	}
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N+3, 2,   N+2, 1, 0,   N+6, N+1, N, N-1,     N+9, N+4, N+5}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));
		return eval_regular_1in4_0v0e(c12, bc[1], bc[2], _p); 
	}
	return false;
}
bool Hoppe94LoopEval::eval_irregular_1v0e_left(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_1v0e_left(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		if (_v + _w > 1) {
			if (_v + _w < 1 + 1e-6) {
				_v -= 0.5e-6; _w -= 0.5e-6; 
			} else {
				std::cout << "Error: Hoppe94LoopEval::eval_irregular_1v0e_left(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; 
				return false; 
			}
		}		
	}

	int N = _cN6.size() - 6;// valence of the extraordinary vertex.

	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------
	//if (_v + _w <= 0.5) {
	//	std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6); 
	//	return cloop_evaluator_->evaluate_irregular_patch(cN6_pick, _v * 2, _w * 2, _p);
	//}
	if (_v + _w <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6); 
		return cloop_evaluator_->evaluate_irregular_patch(cN6_pick, _v * 2, _w * 2, _p, dv_, dw_);
	}
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N+3, 2,   N+2, 1, 0,   N+6, N+1, N, N-1,     N+9, N+4, N+5}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));
		return eval_regular_1in4_0v0e(c12, bc[1], bc[2], _p); 
	}
	if (_v > 0.5) {
		int pick[12] = {N+7, N+8,   N+6, N+2, N+3,   N+9, N+1, 1, 2,    N+4, N, 0};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(1, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+2]; h[1] = _hN9[N+3]; h[2] = _hN9[N+1]; h[3] = _hN9[N+4]; 
		return eval_regular_1v0e(c12, h, v, w, _p);  //  
	} 
	if (_w > 0.5) {
		int pick[12] = {0, N-1,   1, N, N+5,   N+2, N+1, N+4, N+11,   N+6, N+9, N+10};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		return eval_regular_1in4_0v0e(c12, _v  * 2, (_w - 0.5)* 2, _p); 
	}
	return false;
}
bool Hoppe94LoopEval::eval_irregular_1v0e_right(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p){
	//irregular patch, 1 feature vertex, 0 feature edges. 特征点在右下角
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_1v0e_right(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_1v0e_right(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}

	int N = _cN6.size() - 6;// valence of the extraordinary vertex.

	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------
	//if (_v + _w <= 0.5) {
	//	std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6); 
	//	return cloop_evaluator_->evaluate_irregular_patch(cN6_pick, _v * 2, _w * 2, _p);
	//}
	if (_v + _w <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6); 
		return cloop_evaluator_->evaluate_irregular_patch(cN6_pick, _v * 2, _w * 2, _p, dv_, dw_);
	}
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N+3, 2,   N+2, 1, 0,   N+6, N+1, N, N-1,     N+9, N+4, N+5}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));
		return eval_regular_1in4_0v0e(c12, bc[1], bc[2], _p); 
	}
	if (_v > 0.5) {
		int pick[12] = {2, 0,   N+3, 1, N,   N+8, N+2, N+1, N+4,    N+7, N+6, N+9};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		return eval_regular_1in4_0v0e(c12, (_v - 0.5) * 2, _w * 2, _p); 
	} 
	if (_w > 0.5) {
		int pick[12] = {N+11, N+10,   N+5, N+4, N+9,   N-1, N, N+1, N+6,   0, 1, N+2};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+7]; h[1] = _hN9[N+6]; h[2] = _hN9[N+8]; h[3] = _hN9[N+5]; 
		return eval_regular_1v0e(c12, h, bc[1], bc[2], _p);  
	}
	return false;
}
bool Hoppe94LoopEval::eval_irregular_2v0e_left(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_2v0e_left(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_2v0e_left(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}
	DGP::change_slight_vw(_v, _w);

	int N = _cN6.size() - 6;// valence of the extraordinary vertex.
	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
    // -----------------------
	 
	if (_v + _w  <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6);
		std::vector<bool> hN9_pick(N+9, false);
		for (int i = 0; i <= N-1; ++i) {
			hN9_pick[i] = _hN9[i];
		}
		hN9_pick[N+2] = _hN9[0];
		hN9_pick[N+6] = _hN9[N-1];
		return eval_irregular_1v0e_up(cN6_pick, hN9_pick, _v *2, _w *2, _p);
	}  
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N+3, 2,   N+2, 1, 0,   N+6, N+1, N, N-1,     N+9, N+4, N+5}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));
		return eval_regular_1in4_0v0e(c12, bc[1], bc[2], _p); 
	}
	if (_v > 0.5) {
		int pick[12] = {N+7, N+8,   N+6, N+2, N+3,   N+9, N+1, 1, 2,    N+4, N, 0};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(1, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+2]; h[1] = _hN9[N+3]; h[2] = _hN9[N+1]; h[3] = _hN9[N+4]; 
		return eval_regular_1v0e(c12, h, v, w, _p);  //  
	} 
	if (_w > 0.5) {
		int pick[12] = {0, N-1,   1, N, N+5,   N+2, N+1, N+4, N+11,   N+6, N+9, N+10};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		return eval_regular_1in4_0v0e(c12, _v  * 2, (_w - 0.5)* 2, _p); 
	}
	return false;
}
bool Hoppe94LoopEval::eval_irregular_2v0e_right(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p) {
	//irregular patch, 2 feature vertex, 0 feature edges.
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_2v0e_right(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_2v0e_right(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}
	DGP::change_slight_vw(_v, _w);

	int N = _cN6.size() - 6;// valence of the extraordinary vertex.
	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------

	if (_v + _w  <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6);
		std::vector<bool> hN9_pick(N+9, false);
		for (int i = 0; i <= N-1; ++i) {
			hN9_pick[i] = _hN9[i];
		}
		hN9_pick[N+2] = _hN9[0];
		hN9_pick[N+6] = _hN9[N-1];
		return eval_irregular_1v0e_up(cN6_pick, hN9_pick, _v *2, _w *2, _p);
	}  
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N+3, 2,   N+2, 1, 0,   N+6, N+1, N, N-1,     N+9, N+4, N+5}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));
		return eval_regular_1in4_0v0e(c12, bc[1], bc[2], _p); 
	}
	if (_v > 0.5) {
		int pick[12] = {2, 0,   N+3, 1, N,   N+8, N+2, N+1, N+4,    N+7, N+6, N+9};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		return eval_regular_1in4_0v0e(c12, (_v - 0.5) * 2, _w * 2, _p); 
	} 
	if (_w > 0.5) {
		int pick[12] = {N+11, N+10,   N+5, N+4, N+9,   N-1, N, N+1, N+6,   0, 1, N+2};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+7]; h[1] = _hN9[N+6]; h[2] = _hN9[N+8]; h[3] = _hN9[N+5]; 
		return eval_regular_1v0e(c12, h, bc[1], bc[2], _p);  
	}

	return false;
}
bool Hoppe94LoopEval::eval_irregular_2v0e_bottom(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_2v0e_right(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_2v0e_right(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}
	int N = _cN6.size() - 6;// valence of the extraordinary vertex.
	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------
	//if (_v + _w <= 0.5) {
	//	std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6); 
	//	return cloop_evaluator_->evaluate_irregular_patch(cN6_pick, _v * 2, _w * 2, _p);
	//}
	if (_v + _w <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6); 
		return cloop_evaluator_->evaluate_irregular_patch(cN6_pick, _v * 2, _w * 2, _p, dv_, dw_);
	}
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N+3, 2,   N+2, 1, 0,   N+6, N+1, N, N-1,     N+9, N+4, N+5}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));
		return eval_regular_1in4_0v0e(c12, bc[1], bc[2], _p); 
	}
	if (_v > 0.5) {
		int pick[12] = {N+7, N+8,   N+6, N+2, N+3,   N+9, N+1, 1, 2,    N+4, N, 0};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(1, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+2]; h[1] = _hN9[N+3]; h[2] = _hN9[N+1]; h[3] = _hN9[N+4]; 
		return eval_regular_1v0e(c12, h, v, w, _p);  //  
	} 
	if (_w > 0.5) {
		int pick[12] = {N+11, N+10,   N+5, N+4, N+9,   N-1, N, N+1, N+6,   0, 1, N+2};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+7]; h[1] = _hN9[N+6]; h[2] = _hN9[N+8]; h[3] = _hN9[N+5]; 
		return eval_regular_1v0e(c12, h, bc[1], bc[2], _p);  
	}

	return false;
}
bool Hoppe94LoopEval::eval_irregular_3v0e(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v0e(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v0e(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}
	int N = _cN6.size() - 6;// valence of the extraordinary vertex.
	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------

	if (_v + _w  <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6);
		std::vector<bool> hN9_pick(N+9, false);
		for (int i = 0; i <= N-1; ++i) {
			hN9_pick[i] = _hN9[i];
		}
		hN9_pick[N+2] = _hN9[0];
		hN9_pick[N+6] = _hN9[N-1];
		return eval_irregular_1v0e_up(cN6_pick, hN9_pick, _v *2, _w *2, _p);
	}  
	if (_v > 0.5) {
		int pick[12] = {N+7, N+8,   N+6, N+2, N+3,   N+9, N+1, 1, 2,    N+4, N, 0};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(1, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));
		
		std::vector<bool> h(15, false);
		h[0] = _hN9[N+2]; h[1] = _hN9[N+3]; h[2] = _hN9[N+1]; h[3] = _hN9[N+4]; 
		return eval_regular_1v0e(c12, h, bc[1], bc[2], _p);  //  
	} 
	if (_w > 0.5) {
		int pick[12] = {N+11, N+10,   N+5, N+4, N+9,   N-1, N, N+1, N+6,   0, 1, N+2};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+7]; h[1] = _hN9[N+6]; h[2] = _hN9[N+8]; h[3] = _hN9[N+5]; 
		return eval_regular_1v0e(c12, h, bc[1], bc[2], _p);  
	}
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N+3, 2,   N+2, 1, 0,   N+6, N+1, N, N-1,     N+9, N+4, N+5}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));
		return eval_regular_1in4_0v0e(c12, bc[1], bc[2], _p); 
	}
	return false;
}
bool Hoppe94LoopEval::eval_irregular_2v1e_left(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_2v1e_left(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (fabs(_v + _w - 1) < 0.0001) { //有些情况是u=0, 这两个值==1, 但是有一点点偏差, 将那误差忽略掉.
		double len = fabs(_v + _w - 1);
		_v -= len * 0.5;
		_w -= len * 0.5;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_2v1e_left(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}
	int N = _cN6.size() - 6;// valence of the extraordinary vertex.
	
	//特殊处理这个特征点.
	if (_v < 1e-6 && _w < 1e-6) {
		std::vector<int> fv_tmp;
		for (int i = 1; i <= N-2; ++i) {
			if (_hN9[i]) fv_tmp.push_back(i+1);
		}
		if (fv_tmp.size() >= 2) {
			_p = _cN6[0]; return true;
		} else if (fv_tmp.size() == 1) {
			_p = _cN6[0] * 0.66667 + (_cN6[fv_tmp[0]] + _cN6[1]) * 0.166667; return true;
		} else if (fv_tmp.size() == 0) { //dart
			int k = N;
			double inv_v = 1.0/double(k); double t = (3.0 + 2.0 * cos( 2.0 * M_PI * inv_v) ); // double beta = inv_v * (40.0 - t * t)/64.0;
			double k_alpha = 1.0 / ( 24.0/(40.0-t*t) + 1); double alpha = inv_v * k_alpha;

			OpenMesh::Vec3d limit_pos = OpenMesh::Vec3d(0, 0, 0);
			for (int i = 1; i <= N; ++i) {
				limit_pos += _cN6[i] * alpha; 
			}
			limit_pos += _cN6[0] * (1.0 - k_alpha);
			_p = limit_pos; return true;
		} else if (fv_tmp.size() == 0) {
			std::cout << "Error: Hoppe94LoopEval::eval_irregular_1v0e_up(), 应该在eval_irregular_1v0e_left/right()中处理吧.\n"; return false;
		} 
	}

	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------

	if (_v + _w  <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6);
		std::vector<bool> hN9_pick(N+9, false);
		for (int i = 0; i <= N-1; ++i) {
			hN9_pick[i] = _hN9[i];
		}
		hN9_pick[N+2] = _hN9[0];
		hN9_pick[N+6] = _hN9[N-1];
		return eval_irregular_2v1e_left(cN6_pick, hN9_pick, _v * 2, _w * 2, _p);
	}  
	if (_v > 0.5) { 
		int pick[12] = {2, 0,   N+3, 1, N,   N+8, N+2, N+1, N+4,    N+7, N+6, N+9};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		std::vector<bool> h(15, false);
		h[1] = true; h[4] = _hN9[N+4]; h[5] = true; h[8] = _hN9[N+3]; h[11] = _hN9[N+2]; h[12] = _hN9[N+1]; 
		return eval_regular_2v1e(c12, h, (_v-0.5)*2, _w * 2, _p);  //  
	} 
	if (_w > 0.5) {
		int pick[12] = {N+11, N+10,   N+5, N+4, N+9,   N-1, N, N+1, N+6,   0, 1, N+2};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));
		return eval_regular_1in4_0v0e(c12, bc[1], bc[2], _p);
	}
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N+3, 2,   N+2, 1, 0,   N+6, N+1, N, N-1,     N+9, N+4, N+5}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));
		
		std::vector<bool> h(15, false);
		h[2] = true; h[3] = true; 
		return eval_regular_1v0e(c12, h, bc[1], bc[2], _p); 
	}
	return false;
}
bool Hoppe94LoopEval::eval_irregular_2v1e_right(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w,	OpenMesh::Vec3d &_p){
	//irregular patch, 2 feature vertex, 1 feature edges.
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_2v1e_right(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_2v1e_right(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}
	int N = _cN6.size() - 6;// valence of the extraordinary vertex.

	//特殊处理这个特征点.
	if (_v < 1e-6 && _w < 1e-6) {
		std::vector<int> fv_tmp;
		for (int i = 1; i <= N-2; ++i) {
			if (_hN9[i]) fv_tmp.push_back(i+1);
		}
		if (fv_tmp.size() >= 2) {
			_p = _cN6[0]; return true;
		} else if (fv_tmp.size() == 1) {
			_p = _cN6[0] * 0.66667 + (_cN6[fv_tmp[0]] + _cN6[N]) * 0.166667; return true;
		} else if (fv_tmp.size() == 0) { //dart
			int k = N;
			double inv_v = 1.0/double(k); double t = (3.0 + 2.0 * cos( 2.0 * M_PI * inv_v) ); // double beta = inv_v * (40.0 - t * t)/64.0;
			double k_alpha = 1.0 / ( 24.0/(40.0-t*t) + 1); double alpha = inv_v * k_alpha;

			OpenMesh::Vec3d limit_pos = OpenMesh::Vec3d(0, 0, 0);
			for (int i = 1; i <= N; ++i) {
				limit_pos += _cN6[i] * alpha; 
			}
			limit_pos += _cN6[0] * (1.0 - k_alpha);
			_p = limit_pos; return true;
		} else if (fv_tmp.size() == 0) {
			std::cout << "Error: Hoppe94LoopEval::eval_irregular_1v0e_up(), 应该在eval_irregular_1v0e_left/right()中处理吧.\n"; return false;
		} 
	}

	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------

	if (_v + _w  <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6);
		std::vector<bool> hN9_pick(N+9, false);
		for (int i = 0; i <= N-1; ++i) {
			hN9_pick[i] = _hN9[i];
		}
		hN9_pick[N+2] = _hN9[0];
		hN9_pick[N+6] = _hN9[N-1];
		return eval_irregular_2v1e_right(cN6_pick, hN9_pick, _v * 2, _w * 2, _p);
	}  
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N-1, N+5,   0, N, N+4,   2, 1, N+1, N+9,     N+3, N+2, N+6}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[2] = true; h[3] = true; 
		return eval_regular_1v0e(c12, h, bc[1], bc[2], _p); 
	}
	if (_v > 0.5) {
		int pick[12] = {2, 0,   N+3, 1, N,   N+8, N+2, N+1, N+4,    N+7, N+6, N+9};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		return eval_regular_1in4_0v0e(c12, (_v - 0.5) * 2, _w * 2, _p); 
	} 
	if (_w > 0.5) {
		int pick[12] = {N+11, N+10,   N+5, N+4, N+9,   N-1, N, N+1, N+6,   0, 1, N+2};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+7]; h[1] = _hN9[N+6]; h[2] = _hN9[N+8]; h[3] = _hN9[N+5]; 
		h[5] = true; h[11] = true;
		return eval_regular_2v1e(c12, h, bc[1], bc[2], _p);  
	}
	return false;
}
bool Hoppe94LoopEval::eval_irregular_2v1e_bottom(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w,	OpenMesh::Vec3d &_p) {
	//irregular patch, 2 feature vertex, 1 feature edges. 
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_2v1e_bottom(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_2v1e_bottom(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}
	int N = _cN6.size() - 6;// valence of the extraordinary vertex.
	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------

	//if (_v + _w <= 0.5) {
	//	std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6); 
	//	return cloop_evaluator_->evaluate_irregular_patch(cN6_pick, _v * 2, _w * 2, _p);
	//}
	if (_v + _w <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6); 
		return cloop_evaluator_->evaluate_irregular_patch(cN6_pick, _v * 2, _w * 2, _p, dv_, dw_);
	}
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N+9, N+6,     N+4, N+1, N+2,   N+5, N, 1, N+3,   N-1, 0, 2}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[2] = true; h[3] = true; 
		return eval_regular_1v0e(c12, h, bc[1], bc[2], _p); 
	}
	if (_v > 0.5) {
		int pick[12] = {N+7, N+8,   N+6, N+2, N+3,   N+9, N+1, 1, 2,    N+4, N, 0};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(1, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+2]; h[1] = _hN9[N+3]; h[2] = _hN9[N+1]; h[3] = _hN9[N+4]; 
		h[5] = true; h[11] = true;
		return eval_regular_2v1e(c12, h, v, w, _p);  //  
	} 
	if (_w > 0.5) {
		int pick[12] = {N+6, N+2,   N+9, N+1, 1,   N+10, N+4, N, 0,   N+11, N+5, N-1};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[1] = true; h[4] = _hN9[N+5]; h[5] = true; 
		h[8] = _hN9[N+6]; h[11] = _hN9[N+7]; h[12] = _hN9[N+8];
		return eval_regular_2v1e(c12, h, bc[1], bc[2], _p);  
	}
	return false;
}
bool Hoppe94LoopEval::eval_irregular_3v1e_left(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v1e_left(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v1e_left(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}
	int N = _cN6.size() - 6;// valence of the extraordinary vertex.
	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------

	if (_v + _w  <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6);
		std::vector<bool> hN9_pick(N+9, false);
		for (int i = 0; i <= N-1; ++i) {
			hN9_pick[i] = _hN9[i];
		}
		hN9_pick[N+2] = _hN9[0];
		hN9_pick[N+6] = _hN9[N-1];
		return eval_irregular_2v1e_left(cN6_pick, hN9_pick, _v *2, _w *2, _p);
	}  
	if (_v > 0.5) {
		int pick[12] = {2, 0,   N+3, 1, N,   N+8, N+2, N+1, N+4,    N+7, N+6, N+9};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		std::vector<bool> h(15, false);
		h[1] = true; h[4] = _hN9[N+4]; h[5] = true; h[8] = _hN9[N+3]; h[11] = _hN9[N+2]; h[12] = _hN9[N+1]; 
		return eval_regular_2v1e(c12, h, (_v-0.5)*2, _w*2, _p);  //  
	} 
	if (_w > 0.5) {
		int pick[12] = {N+11, N+10,   N+5, N+4, N+9,   N-1, N, N+1, N+6,   0, 1, N+2};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+7]; h[1] = _hN9[N+6]; h[2] = _hN9[N+8]; h[3] = _hN9[N+5]; 
		return eval_regular_1v0e(c12, h, bc[1], bc[2], _p);  
	}
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N+3, 2,   N+2, 1, 0,   N+6, N+1, N, N-1,     N+9, N+4, N+5}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[2] = true; h[3] = true; 
		return eval_regular_1v0e(c12, h, bc[1], bc[2], _p); 
	}
	return false;
}
bool Hoppe94LoopEval::eval_irregular_3v1e_right(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p) {
	//irregular patch, 3 feature vertex, 1 feature edges.
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v1e_right(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v1e_right(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}
	int N = _cN6.size() - 6;// valence of the extraordinary vertex.
	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------
	if (_v + _w  <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6);
		std::vector<bool> hN9_pick(N+9, false);
		for (int i = 0; i <= N-1; ++i) {
			hN9_pick[i] = _hN9[i];
		}
		hN9_pick[N+2] = _hN9[0];
		hN9_pick[N+6] = _hN9[N-1];
		return eval_irregular_2v1e_right(cN6_pick, hN9_pick, _v *2, _w *2, _p);
	}  
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N-1, N+5,   0, N, N+4,   2, 1, N+1, N+9,     N+3, N+2, N+6}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[2] = true; h[3] = true; 
		return eval_regular_1v0e(c12, h, bc[1], bc[2], _p); 
	}
	if (_v > 0.5) {
		int pick[12] = {N+7, N+8,   N+6, N+2, N+3,   N+9, N+1, 1, 2,    N+4, N, 0};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(1, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+2]; h[1] = _hN9[N+3]; h[2] = _hN9[N+1]; h[3] = _hN9[N+4]; 
		return eval_regular_1v0e(c12, h, v, w, _p);  //  
	} 
	if (_w > 0.5) {
		int pick[12] = {N+11, N+10,   N+5, N+4, N+9,   N-1, N, N+1, N+6,   0, 1, N+2};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+7]; h[1] = _hN9[N+6]; h[2] = _hN9[N+8]; h[3] = _hN9[N+5]; 
		h[5] = true; h[11] = true;
		return eval_regular_2v1e(c12, h, bc[1], bc[2], _p);  
	}
	return false;
}
bool Hoppe94LoopEval::eval_irregular_3v1e_bottom(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p) {
	//irregular patch, 3 feature vertex, 1 feature edges.
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v1e_bottom(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v1e_bottom(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}
	int N = _cN6.size() - 6;// valence of the extraordinary vertex.
	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------
	if (_v + _w  <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6);
		std::vector<bool> hN9_pick(N+9, false);
		for (int i = 0; i <= N-1; ++i) {
			hN9_pick[i] = _hN9[i];
		}
		hN9_pick[N+2] = _hN9[0];
		hN9_pick[N+6] = _hN9[N-1];
		return eval_irregular_1v0e_up(cN6_pick, hN9_pick, _v *2, _w *2, _p);
	}  
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N+9, N+6,     N+4, N+1, N+2,   N+5, N, 1, N+3,   N-1, 0, 2}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[2] = true; h[3] = true; 
		return eval_regular_1v0e(c12, h, bc[1], bc[2], _p); 
	}
	if (_v > 0.5) {
		int pick[12] = {N+7, N+8,   N+6, N+2, N+3,   N+9, N+1, 1, 2,    N+4, N, 0};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(1, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+2]; h[1] = _hN9[N+3]; h[2] = _hN9[N+1]; h[3] = _hN9[N+4]; 
		h[5] = true; h[11] = true;
		return eval_regular_2v1e(c12, h, v, w, _p);  //  
	} 
	if (_w > 0.5) {
		int pick[12] = {N+6, N+2,   N+9, N+1, 1,   N+10, N+4, N, 0,   N+11, N+5, N-1};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[1] = true; h[4] = _hN9[N+5]; h[5] = true; 
		h[8] = _hN9[N+6]; h[11] = _hN9[N+7]; h[12] = _hN9[N+8];
		return eval_regular_2v1e(c12, h, bc[1], bc[2], _p);  
	}
	return false;
}
bool Hoppe94LoopEval::eval_irregular_3v3e(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v3e(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v3e(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}
	int N = _cN6.size() - 6;// valence of the extraordinary vertex.
	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------

	if (_v + _w  <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6);
		std::vector<bool> hN9_pick(N+9, false);
		for (int i = 0; i <= N-1; ++i) {
			hN9_pick[i] = _hN9[i];
		}
		hN9_pick[N+2] = _hN9[0];//true
		hN9_pick[N+6] = _hN9[N-1];//true
		return eval_irregular_3v2e_up(cN6_pick, hN9_pick, _v *2, _w *2, _p);
	}  
	if (_v > 0.5) {
		int pick[12] = {N+7, N+8,   N+6, N+2, N+3,   N+9, N+1, 1, 2,    N+4, N, 0};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(1, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+2]; h[1] = _hN9[N+3]; h[2] = _hN9[N+1]; h[3] = _hN9[N+4]; 
		h[5] = true; h[6] = true; h[11] = true; h[14] = true; 
		return eval_regular_3v2e(c12, h, bc[1], bc[2], _p);  //  
	} 
	if (_w > 0.5) {
		int pick[12] = {N+11, N+10,   N+5, N+4, N+9,   N-1, N, N+1, N+6,   0, 1, N+2};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+7]; h[1] = _hN9[N+6]; h[2] = _hN9[N+8]; h[3] = _hN9[N+5]; 
		h[5] = true; h[6] = true; h[11] = true;h[14] = true;
		return eval_regular_3v2e(c12, h, bc[1], bc[2], _p);  
	}
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N-1, N+5,   0, N, N+4,   2, 1, N+1, N+9,     N+3, N+2, N+6}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[2] = true; h[3] = true; h[4] = true; h[7] = true; h[12] = true; h[13] = true;
		return eval_regular_3v0e(c12, h, bc[1], bc[2], _p); 
	}
	return false;
}
bool Hoppe94LoopEval::eval_irregular_3v2e_up(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v2e_up(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v2e_up(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}
	DGP::change_slight_vw(_v, _w);

	int N = _cN6.size() - 6;// valence of the extraordinary vertex.

	//特殊处理这个特征点.
	if (_v < 1e-6 && _w < 1e-6) {
		std::vector<int> fv_tmp;
		for (int i = 1; i <= N-2; ++i) {
			if (_hN9[i]) fv_tmp.push_back(i+1);
		}
		if (fv_tmp.size() >= 1) {
			_p = _cN6[0]; return true;
		} else if (fv_tmp.size() == 0) {
			_p = _cN6[0] * 0.66667 + (_cN6[1] + _cN6[N]) * 0.166667; return true;
		} 
	}

	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------

	if (_v + _w  <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6);
		std::vector<bool> hN9_pick(N+9, false);
		for (int i = 0; i <= N-1; ++i) {
			hN9_pick[i] = _hN9[i];
		}
		hN9_pick[N+2] = _hN9[0];//true
		hN9_pick[N+6] = _hN9[N-1];//true
		return eval_irregular_3v2e_up(cN6_pick, hN9_pick, _v *2, _w *2, _p);
	}  
	if (_v > 0.5) {
		int pick[12] = {2, 0,   N+3, 1, N,   N+8, N+2, N+1, N+4,    N+7, N+6, N+9};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		std::vector<bool> h(15, false);
		h[1] = true; h[4] = _hN9[N+4]; h[5] = true; h[8] = _hN9[N+3]; h[11] = _hN9[N+2]; h[12] = _hN9[N+1]; 
		return eval_regular_2v1e(c12, h, (_v-0.5)*2, _w*2, _p);  //  
	} 
	if (_w > 0.5) {
		int pick[12] = {N+11, N+10,   N+5, N+4, N+9,   N-1, N, N+1, N+6,   0, 1, N+2};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+7]; h[1] = _hN9[N+6]; h[2] = _hN9[N+8]; h[3] = _hN9[N+5]; 
		h[5] = true; h[11] = true;
		return eval_regular_2v1e(c12, h, bc[1], bc[2], _p);  
	}
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N-1, N+5,   0, N, N+4,   2, 1, N+1, N+9,     N+3, N+2, N+6}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[2] = true; h[3] = true; h[4] = true; h[12] = true;
		return eval_regular_2v0e(c12, h, bc[1], bc[2], _p); 
	}
	return false;
}


bool Hoppe94LoopEval::eval_irregular_3v2e_left(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v2e_left(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v2e_left(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}
	int N = _cN6.size() - 6;// valence of the extraordinary vertex.
	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------
	if (_v + _w  <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6);
		std::vector<bool> hN9_pick(N+9, false);
		for (int i = 0; i <= N-1; ++i) {
			hN9_pick[i] = _hN9[i];
		}
		hN9_pick[N+2] = _hN9[0];
		hN9_pick[N+6] = _hN9[N-1];
		return eval_irregular_2v1e_left(cN6_pick, hN9_pick, _v * 2, _w * 2, _p);
	}  
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N+3, 2,   N+2, 1, 0,   N+6, N+1, N, N-1,     N+9, N+4, N+5}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[2] = true; h[3] = true; h[4] = true; h[12] = true; 
		return eval_regular_2v0e(c12, h, bc[1], bc[2], _p); 
	}
	if (_v > 0.5) {
		int pick[12] = {N+7, N+8,   N+6, N+2, N+3,   N+9, N+1, 1, 2,    N+4, N, 0};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(1, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+2]; h[1] = _hN9[N+3]; h[2] = _hN9[N+1]; h[3] = _hN9[N+4]; 
		h[5] = true; h[6] = true; h[11] = true; h[14] = true; 
		return eval_regular_3v2e(c12, h, bc[1], bc[2], _p);  //  
	} 
	if (_w > 0.5) {
		int pick[12] = {N+6, N+2,   N+9, N+1, 1,   N+10, N+4, N, 0,   N+11, N+5, N-1};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[1] = true; h[4] = _hN9[N+5]; h[5] = true; 
		h[8] = _hN9[N+6]; h[11] = _hN9[N+7]; h[12] = _hN9[N+8];
		return eval_regular_2v1e(c12, h, bc[1], bc[2], _p);  
	}
	return false;
}
bool Hoppe94LoopEval::eval_irregular_3v2e_right(const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, double _v, double _w, OpenMesh::Vec3d &_p) {
	if (_cN6.size() - 6 != _hN9.size() - 9) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v2e_right(), _cN6 or _hN9 's size is wrong.\n"; return false;
	}
	if (_v + _w < 0 || _v + _w > 1) {
		std::cout << "Error: Hoppe94LoopEval::eval_irregular_3v2e_right(), _v or _w is invalid, " << _v << ", " << _w << ".\n"; return false;
	}
	int N = _cN6.size() - 6;// valence of the extraordinary vertex.
	std::vector<OpenMesh::Vec3d> cN12(N + 12, OpenMesh::Vec3d(0, 0, 0));
	eval_irregular_cN12(_cN6, _hN9, cN12);
	// -----------------------
	if (_v + _w  <= 0.5) {
		std::vector<OpenMesh::Vec3d> cN6_pick(cN12.begin(), cN12.begin()+N+6);
		std::vector<bool> hN9_pick(N+9, false);
		for (int i = 0; i <= N-1; ++i) {
			hN9_pick[i] = _hN9[i];
		}
		hN9_pick[N+2] = _hN9[0];
		hN9_pick[N+6] = _hN9[N-1];
		return eval_irregular_2v1e_right(cN6_pick, hN9_pick, _v * 2, _w * 2, _p);
	}  
	if (_v + _w > 0.5 && _v <= 0.5 && _w <= 0.5) {
		int pick[12] = {N+9, N+6,     N+4, N+1, N+2,   N+5, N, 1, N+3,   N-1, 0, 2}; 
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}

		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[2] = true; h[3] = true; h[4] = true; h[12] = true; 
		return eval_regular_2v0e(c12, h, bc[1], bc[2], _p); 
	}
	if (_v > 0.5) {
		int pick[12] = {N+7, N+8,   N+6, N+2, N+3,   N+9, N+1, 1, 2,    N+4, N, 0};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		double v = 0, w = 0;
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(1, 0), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(0.5, 0), OpenMesh::Vec2d(_v, _w));
		v = bc[1]; w = bc[2];

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+2]; h[1] = _hN9[N+3]; h[2] = _hN9[N+1]; h[3] = _hN9[N+4]; 
		h[5] = true; h[11] = true;
		return eval_regular_2v1e(c12, h, v, w, _p);  //  
	} 
	if (_w > 0.5) {
		int pick[12] = {N+11, N+10,   N+5, N+4, N+9,   N-1, N, N+1, N+6,   0, 1, N+2};
		std::vector<OpenMesh::Vec3d> c12(12, OpenMesh::Vec3d(0, 0, 0));
		for (int i = 0; i < 12; ++i) {
			c12[i] = cN12[pick[i]]; 
		}
		OpenMesh::Vec3d bc = DGP::calc_barycentric_coordinates(OpenMesh::Vec2d(0, 1), OpenMesh::Vec2d(0, 0.5), OpenMesh::Vec2d(0.5, 0.5), OpenMesh::Vec2d(_v, _w));

		std::vector<bool> h(15, false);
		h[0] = _hN9[N+7]; h[1] = _hN9[N+6]; h[2] = _hN9[N+8]; h[3] = _hN9[N+5]; 
		h[5] = true; h[6] = true; h[11] = true;h[14] = true;
		return eval_regular_3v2e(c12, h, bc[1], bc[2], _p);  
	}
	return false;
}
// ===================================================
// -----------
// -- 为什么的7个regular和17个irregular的函数提供统一的接口
bool Hoppe94LoopEval::eval_regular_1in7(feature_regular_type _type, 
					   const std::vector<OpenMesh::Vec3d> &_c12, const std::vector<bool> &_h15, 
					   double _v, double _w, 
					   OpenMesh::Vec3d &_p, OpenMesh::Vec3d &_dv, OpenMesh::Vec3d &_dw)
{
	bool flag = false; 
	DGP::change_slight_vw(_v, _w);
	switch(_type)
	{ 
	case rt1v0e: flag = eval_regular_1v0e(_c12, _h15, _v, _w, _p);  break; 
	case rt2v0e: flag = eval_regular_2v0e(_c12, _h15, _v, _w, _p);  break; 
	case rt2v1e: flag = eval_regular_2v1e(_c12, _h15, _v, _w, _p);  break; 
	case rt3v0e: flag = eval_regular_3v0e(_c12, _h15, _v, _w, _p);  break; 
	case rt3v1e: flag = eval_regular_3v1e(_c12, _h15, _v, _w, _p);  break; 
	case rt3v2e: flag = eval_regular_3v2e(_c12, _h15, _v, _w, _p);  break; 
	case rt3v3e: flag = eval_regular_3v3e(_c12, _h15, _v, _w, _p);  break; 
	default: std::cout << "Error: eval_regular_1in7, No such feature_regular_type.\n"; flag = false;break; 
	}
	_dv = dv_; _dw = dw_;
	return flag;
}
bool Hoppe94LoopEval::eval_irregular_1in17(feature_irregular_type _type, 
	const std::vector<OpenMesh::Vec3d> &_cN6, const std::vector<bool> &_hN9, 
	double _v, double _w, 
	OpenMesh::Vec3d &_p, OpenMesh::Vec3d &_dv, OpenMesh::Vec3d &_dw)
{
	bool flag = false;
	DGP::change_slight_vw(_v, _w);
	switch(_type) 
	{  						
	case irt1v0e_up:		flag = eval_irregular_1v0e_up(_cN6, _hN9, _v, _w, _p);  break; 
	case irt1v0e_left:		flag = eval_irregular_1v0e_left(_cN6, _hN9, _v, _w, _p);   break; 
	case irt1v0e_right:		flag = eval_irregular_1v0e_right(_cN6, _hN9, _v, _w, _p);  break; 
	case irt2v0e_left:		flag = eval_irregular_2v0e_left(_cN6, _hN9, _v, _w, _p);  break; 
	case irt2v0e_right:		flag = eval_irregular_2v0e_right(_cN6, _hN9, _v, _w, _p);  break; 
	case irt2v0e_bottom:	flag = eval_irregular_2v0e_bottom(_cN6, _hN9, _v, _w, _p);  break; 
	case irt3v0e:			flag = eval_irregular_3v0e(_cN6, _hN9, _v, _w, _p);   break; 
	case irt2v1e_left:		flag = eval_irregular_2v1e_left(_cN6, _hN9, _v, _w, _p);   break; 
	case irt2v1e_right:		flag = eval_irregular_2v1e_right(_cN6, _hN9, _v, _w, _p);  break; 
	case irt2v1e_bottom:	flag = eval_irregular_2v1e_bottom(_cN6, _hN9, _v, _w, _p);   break; 
	case irt3v1e_left:		flag = eval_irregular_3v1e_left(_cN6, _hN9, _v, _w, _p);   break; 
	case irt3v1e_right:		flag = eval_irregular_3v1e_right(_cN6, _hN9, _v, _w, _p);   break; 
	case irt3v1e_bottom:	flag = eval_irregular_3v1e_bottom(_cN6, _hN9, _v, _w, _p);  break; 
	case irt3v2e_up:		flag = eval_irregular_3v2e_up(_cN6, _hN9, _v, _w, _p);  break; 
	case irt3v2e_left:		flag = eval_irregular_3v2e_left(_cN6, _hN9, _v, _w, _p);   break; 
	case irt3v2e_right:		flag = eval_irregular_3v2e_right(_cN6, _hN9, _v, _w, _p);  break; 
	case irt3v3e:			flag = eval_irregular_3v3e(_cN6, _hN9, _v, _w, _p);   break; 
	default: flag = false;break; 
	}	
	_dv = dv_; _dw = dw_;
	return flag;
}
// ===================================================
// -- 带尖锐特征的细分曲线参数化
bool Hoppe94LoopEval::eval_curve_crease_crease(const OpenMesh::Vec3d& P_0_iminus1, const OpenMesh::Vec3d& P_0_i, const OpenMesh::Vec3d& P_0_iplus1, const OpenMesh::Vec3d& P_0_iplus2, 
	double _t, OpenMesh::Vec3d &_p) 
{
	//std::cout << P_0_iminus1 << ", " << P_0_i << ", " << P_0_iplus1 << ", " << P_0_iplus2 << ", " << _t << ", " << _p << ".\n"; // for test.
	// 规则的三次B样条曲线. 能直接evaluation了.
	OpenMesh::Vec3d tmp[4] = {
		(double)-1 * P_0_iminus1 + (double)3 * P_0_i - (double)3 * P_0_iplus1 + P_0_iplus2,
		(double) 3 * P_0_iminus1 - (double)6 * P_0_i + (double)3 * P_0_iplus1,
		(double)-3 * P_0_iminus1                     + (double)3 * P_0_iplus1,
		(double) 1 * P_0_iminus1 + (double)4 * P_0_i + (double)1 * P_0_iplus1
	};
	
	double t = _t, t2 = t*t, t3 = t2*t;
	_p = t3 * tmp[0] + t2 * tmp[1] + t * tmp[2] + tmp[3];
	_p *= 1.0/6.0;
	return true;
}
bool Hoppe94LoopEval::eval_curve_corner_crease(const OpenMesh::Vec3d& P_0_i, const OpenMesh::Vec3d& P_0_iplus1, const OpenMesh::Vec3d& P_0_iplus2, 
	double _t, OpenMesh::Vec3d &_p)
{
	if (_t < 1e-6) { // 太接近角点了,直接给出其值.
		_p = P_0_i; return true;
	}
	if (_t > 0.999999) { // 太接近crease点了,直接给出其值.
		_p = P_0_iplus1 * 0.66667 + (P_0_i + P_0_iplus2) * 0.166667; return true;
	}


	// 
	int k = (int)floor(1 - log(_t) / log(2.0)); // k >= 1.
	double t_new = _t * pow(2.0, k) - 1;

	//std::cout << _t << ", " << k << ", " << t_new << "; "  << ", " << P_0_i << ", " << P_0_iplus1 << ", " << P_0_iplus2 << ", " << _p << ".\n";

	double a = pow(1.0, k-1), b = pow(1.0/2.0, k-1), c = pow(1.0/8.0, k-1);	 
	double A1[4][3] = { {1, 0, 0},
						{0.5, 0.5, 0},
						{0.125, 0.75, 0.125},
						{0, 0.5, 0.5}};  
	double A2[3][3] = { {a, 0, 0},
						{a-b, b, 0},
						{a-2*b+c, 2*b-2*c, c}}; 
	double A[4][3] = { 0 };
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				A[i][j] += A1[i][k] * A2[k][j];
			}
		} 
	}

	OpenMesh::Vec3d p[4] = { OpenMesh::Vec3d(0, 0, 0) };
	for (int i = 0; i < 4; ++i) {
		p[i] = A[i][0] * P_0_i + A[i][1] * P_0_iplus1 + A[i][2] * P_0_iplus2; 
	}

	return eval_curve_crease_crease(p[0], p[1], p[2], p[3], t_new, _p);
}
bool Hoppe94LoopEval::eval_curve_corner_corner(const OpenMesh::Vec3d& P_0_i, const OpenMesh::Vec3d& P_0_iplus1, double _t, OpenMesh::Vec3d &_p) {
	// 退化为直线, 
	_p = (1 - _t) * P_0_i + _t * P_0_iplus1;
	return true;
}