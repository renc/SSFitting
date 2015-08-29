// ***************************************************************
//  
//  barycentric 
//  Copyright (C) 2007 - by rencanjiang All Rights Reserved
//  rencanjiang@163.com	
//	-------------------------------------------------------------
//  
//  version: 1.0
//  data: 12/10/2007
// ***************************************************************
//  
//  File description:
// ***************************************************************

#ifndef dgpstudio_barycentric_h
#define dgpstudio_barycentric_h

#include "OpenMeshAll.h"

namespace DGP {
inline void change_slight_vw(double & _u, double& _v, double& _w) {
	// (u,v,w) 理论上不该做此改动, 但是对于浮点数据需要整理一下, 构成了一个barycentric coordinate是让其有效.
	if (fabs(_u) < 1.0e-5) _u = 0;
	if (fabs(_v) < 1.0e-5) _v = 0;//原来是1.0e-6,现在改为1.0e-5, 否则像_bc(-1.37927e-006 0.19194 0.808062)
	if (fabs(_w) < 1.0e-5) _w = 0;//就被认为是无效的. 而事实上第一项该为0, 后两项和为1.

	if (_u >=1 && _u < 1.00001) _u = 1;
	if (_v >=1 && _v < 1.00001) _v = 1;
	if (_w >=1 && _w < 1.00001) _w = 1;

	if ( _v + _w > 1) { 
		if (_v+_w < 1.00001) {
			if (_v < _w) _w -= 0.00001; else _v -= 0.00001;
		} else {
			//std::cout << "Error: change_slight_vw: " << _v << " + " << _w << " > 1.\n"; 
		}
	}
	if (_u + _v + _w > 1) { 
		if (_u + _v+_w < 1.00001) {
			if (_u < _v && _w < _v) _v -= 0.00001; 
			if (_u < _w && _v < _w) _w -= 0.00001; 
			if (_v < _u && _w < _u) _u -= 0.00001; 
		} else {
			//std::cout << "Error: change_slight_vw: " << _v << " + " << _w << " > 1.\n"; 
		}
	}
	//if (_v < 1e-6) {
	//	_v += 0.00001;
	//	if (_v + _w > 1) _w -= 0.00001;
	//}
	//if (_w < 1e-6) {
	//	_w += 0.00001;
	//	if (_v + _w > 1) _v -= 0.00001;
	//}
	//if (_v + _w > 1 && _v + _w < 1 + 1e-5) {
	//	_v -= 0.5e-5; _w -= 0.5e-5;
	//}
}
inline void change_slight_vw(double& _v, double& _w) {
	double u = 1-_v-_w;
	change_slight_vw(u, _v, _w);
}

// 给出(1)三角形p0p1p2, 注意顶点的循序,(2)三角形面上的一个顶点pp
// 求出这个顶点对于这个三角形的重心坐标
// 返回的时候注意检查重心左边不能是(-1, -1, -1)表示出错了.
inline OpenMesh::Vec3d calc_barycentric_coordinates(const OpenMesh::Vec3d& p0, const OpenMesh::Vec3d& p1, const OpenMesh::Vec3d& p2,
	const OpenMesh::Vec3d& pp,
	bool _debug = false) {
			using namespace OpenMesh;
			const Vec3d v1(p1 - p0);
			const Vec3d v2(p2 - p0);
 
			// I don't need to normalize normal
			const Vec3d norm(cross(v2, v1));
			if (_debug) std::cout << ">" << v1 << ", " << v2 << ", "<< norm << " ";
			// now we need to convert on_plane to barycentric coords.
			// we need edge inward-perpendiculars tv1, tv2, tw12
			const Vec3d tv1((cross(v1, norm)).normalize());
			const Vec3d tv2((cross(norm, v2)).normalize());

			const Vec3d tw12((cross(p2 - p1, norm)).normalize());
			if (_debug) std::cout << ": " << p1 << "; " << (cross(p2 - p1, norm)) << "; " << tw12 << "; ";
			const Vec3d ppp0(pp - p0);

			double b0 = dot(pp - p1, tw12) / dot(-v1, tw12);
			double b1 = dot(ppp0, tv2) / dot(v1, tv2);
			double b2 = dot(ppp0, tv1) / dot(v2, tv1);// b2 = 1 - b0 - b1;

			if ((b0 < 0) && (b1 < 0) && (b2 < 0)) {
				std::cerr << "Error: barycentric coordinates cannot all be negative.\n";
				return Vec3d(-1, -1, -1);
			}
			//std::cout << Vec3d(b0, b1, b2) << std::endl;
			change_slight_vw(b0, b1, b2);
			return Vec3d(b0, b1, b2);		
}
inline OpenMesh::Vec3d calc_barycentric_coordinates(const OpenMesh::Vec2d& _p0, const OpenMesh::Vec2d& _p1, const OpenMesh::Vec2d& _p2,
													const OpenMesh::Vec2d& _pp) {
														using namespace OpenMesh;
														return calc_barycentric_coordinates(Vec3d(_p0[0], _p0[1], 0), Vec3d(_p1[0], _p1[1], 0), 
															Vec3d(_p2[0], _p2[1], 0), Vec3d(_pp[0], _pp[1], 0));
}


inline bool is_valid_barycentric_coordinate(OpenMesh::Vec3d _bc, bool _debug = false) {
	if (fabs(_bc[0]) < 1.0e-5) _bc[0] = 0;//原来是1.0e-6,现在改为1.0e-5, 否则像_bc(-1.37927e-006 0.19194 0.808062)
	if (fabs(_bc[1]) < 1.0e-5) _bc[1] = 0;//就被认为是无效的. 而事实上第一项该为0, 后两项和为1.
	if (fabs(_bc[2]) < 1.0e-5) _bc[2] = 0;//还有(0.303709 -0.000569765 0.696861)
	
	if (_bc[0] > 1 || _bc[1] > 1 || _bc[2] > 1) {
		if (_debug) std::cerr << "Error: bc[i]>1 "; return false; 
	}
	if (_bc[0] + _bc[1] > 1.00001) {
		if (_debug) std::cerr << "Error: bc[0+1]>1 "; return false;
	}
	if (_bc[0] + _bc[2] > 1)  {
		if (_debug) std::cerr << "Error: bc[0+2]>1 "; return false;
	}
	if (_bc[1] + _bc[2] > 1)  {
		if (_debug) std::cerr << "Error: bc[1+2]>1 "; return false;
	}
	if (_bc[0] + _bc[1] + _bc[2] > 1) {
		if (_bc[0] + _bc[1] + _bc[2] < 1.00001) {
			if (_bc[0] <= _bc[1] && _bc[1] <= _bc[2] ) _bc[2] -= 0.00001;
			if (_bc[0] <= _bc[1] && _bc[2] <= _bc[1] ) _bc[1] -= 0.00001; 
			if (_bc[0] >= _bc[1] && _bc[0] >= _bc[2] ) _bc[0] -= 0.00001; 
			if (_bc[0] >= _bc[1] && _bc[0] <= _bc[2] ) _bc[2] -= 0.00001; 
		} else {
			if (_debug) std::cerr << "Error: bc[0+1+2]>1," << _bc[0] + _bc[1] + _bc[2] - 1 << " "; return false; 
		}		
	}
	
	//if (!(_bc[0]<0) && !(_bc[1]<0) && !(_bc[2]<0)) { //要求这个_bc的三个分量都是>=0
	if ((_bc[0] >= 0) && (_bc[1] >= 0) && (_bc[2] >= 0)) { //要求这个_bc的三个分量都是>=0
		return true;
	}
	return false;
}
} // end of DGP
#endif