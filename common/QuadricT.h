//=============================================================================
//                                                                            
//   Example code for the full-day course
//
//   M. Botsch, M. Pauly, L. Kobbelt, P. Alliez, B. Levy,
//   "Geometric Modeling Based on Polygonal Meshes"
//   held at SIGGRAPH 2007, San Diego, USA.
//
//   Copyright (C) 2007 by  Computer Graphics Laboratory, ETH Zurich, 
//                      and Computer Graphics Group,      RWTH Aachen
//
//                                                                            
//-----------------------------------------------------------------------------
//                                                                            
//                                License                                     
//                                                                            
//   This program is free software; you can redistribute it and/or
//   modify it under the terms of the GNU General Public License
//   as published by the Free Software Foundation; either version 2
//   of the License, or (at your option) any later version.
//   
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//   
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 51 Franklin Street, Fifth Floor, 
//   Boston, MA  02110-1301, USA.
//
//   changedg by rencanjiang@163.com
//                                                                            
//=============================================================================
//=============================================================================
//
//  CLASS QuadricT
//  是参考OpenMesh中class QuadricT来修改实现的
//=============================================================================
// ***************************************************************
//  
//  QuadricT 
//  Copyright (C) 2007 - changed by rencanjiang All Rights Reserved
//  rencanjiang@163.com	
//	-------------------------------------------------------------
//  
//  version: 1.0
//  data: 12/10/2007
// ***************************************************************
//  
//  File description:
// ***************************************************************
#ifndef QUADRIC_HH
#define QUADRIC_HH

//== INCLUDES =================================================================

#include "OpenMeshAll.h"
#include <iostream> //

//== CLASS DEFINITION =========================================================

namespace DGP {
	 

/** /class QuadricT 应该是很容易实现的

Stores a quadric as a 4x4 symmetrix matrix. Used by the
error quadric based mesh decimation algorithms.
**/
//这个Scalar是float或者是double类型
template <class Scalar>  
class QuadricT
{
public:


	/// construct with upper triangle of symmetrix 4x4 matrix
	QuadricT(Scalar _a, Scalar _b, Scalar _c, Scalar _d,
						Scalar _e, Scalar _f, Scalar _g,
								   Scalar _h, Scalar _i,
								              Scalar _j)
		     : a(_a), b(_b), c(_c), d(_d), 
		              e(_e), f(_f), g(_g),    
		                     h(_h), i(_i), 
		                            j(_j) {}


	/// constructor from given plane equation: ax+by+cz+d=0, n=[a, b, c] unit normal
	QuadricT( Scalar _a=0.0, Scalar _b=0.0, Scalar _c=0.0, Scalar _d=0.0 )
		     :  a(_a*_a), b(_a*_b),  c(_a*_c),  d(_a*_d),
		                  e(_b*_b),  f(_b*_c),  g(_b*_d),
		                             h(_c*_c),  i(_c*_d),
		                                        j(_d*_d) {}


	/// set all entries to zero
	void clear()  { a = b = c = d = e = f = g = h = i = j = 0.0; }


	/// add quadrics
	QuadricT<Scalar>& operator+=( const QuadricT<Scalar>& _q )
	{
		a += _q.a;  b += _q.b;  c += _q.c;  d += _q.d;
		            e += _q.e;  f += _q.f;  g += _q.g;
		                        h += _q.h;  i += _q.i;
		                                    j += _q.j;
		return *this;
	}
	QuadricT<Scalar>& operator-=( const QuadricT<Scalar>& _q )
	{
		a -= _q.a;  b -= _q.b;  c -= _q.c;  d -= _q.d;
		            e -= _q.e;  f -= _q.f;  g -= _q.g;
		                        h -= _q.h;  i -= _q.i;
		                                    j -= _q.j;
		return *this;
	}


	/// multiply by scalar
	QuadricT<Scalar>& operator*=( Scalar _s)
	{
		a *= _s;  b *= _s;  c *= _s;  d *= _s;
		          e *= _s;  f *= _s;  g *= _s;
				            h *= _s;  i *= _s;
									  j *= _s;
		return *this;
	}


	/// evaluate quadric Q at vector v: v*Q*v, 求出error cost
	template <typename T>
	Scalar operator()(const OpenMesh::VectorT<T,3> _v) const
	{
		Scalar x(_v[0]), y(_v[1]), z(_v[2]);
		return a*x*x + 2.0*b*x*y + 2.0*c*x*z + 2.0*d*x
			         +     e*y*y + 2.0*f*y*z + 2.0*g*y
			                     +     h*z*z + 2.0*i*z
			                                 +     j;
	}
	
	// 加入 Vertex Placement Policies, A可逆的话求最优的新点newV使得errorQ(newV)最小
	// Q = [a, b, c] d
	//     [b, e, f] g
	//     [c, f, h] i
	//     [d, g, i, j    A 是3阶方阵
	/* bool is_invertible_A() // 判断矩阵A是否可逆
	{
		Scalar detA = a*(e*h-f*f) - b*(b*h-f*c) + c*(b*f-e*c);
		if (std::abs(detA) <= 1e-12) return false; // detA = 0 所以A不可逆
		return true;
	}*/
	bool optimal_placement(Scalar &_newvp0, Scalar &_newvp1, Scalar &_newvp2)
	{	//A可逆返回真,否则为假
		Scalar rkInverse[3][3] = { {0.0, 0.0, 0.0},
									{0.0, 0.0, 0.0},
									{0.0, 0.0, 0.0}};   // inverst(A)
		Scalar elt[3][3] = {a, b, c,
							b, e, f,
							c, f, h};
		//std::cout << _newvp0 << " " << _newvp1 << " " << _newvp2 << std::endl;
		rkInverse[0][0] = elt[1][1] * elt[2][2] -
                          elt[1][2] * elt[2][1];
		rkInverse[0][1] = elt[0][2] * elt[2][1] -
						  elt[0][1] * elt[2][2];
		rkInverse[0][2] = elt[0][1] * elt[1][2] -
						  elt[0][2] * elt[1][1];
		rkInverse[1][0] = elt[1][2] * elt[2][0] -
						  elt[1][0] * elt[2][2];
		rkInverse[1][1] = elt[0][0] * elt[2][2] -
						  elt[0][2] * elt[2][0];
		rkInverse[1][2] = elt[0][2] * elt[1][0] -
						  elt[0][0] * elt[1][2];
		rkInverse[2][0] = elt[1][0] * elt[2][1] -
						  elt[1][1] * elt[2][0];
		rkInverse[2][1] = elt[0][1] * elt[2][0] -
						  elt[0][0] * elt[2][1];
		rkInverse[2][2] = elt[0][0] * elt[1][1] -
						  elt[0][1] * elt[1][0];

		Scalar fDet =
			elt[0][0] * rkInverse[0][0] +
			elt[0][1] * rkInverse[1][0] +
			elt[0][2] * rkInverse[2][0];
		//std::cerr << "Q: optimal_placement(): " << fDet << std::endl;//for test
		if (std::abs(fDet) <= 1e-12 )
			return false;//不可逆

		Scalar fInvDet = (float)1.0 / fDet;

		for (int iRow = 0; iRow < 3; iRow++) {
			for (int iCol = 0; iCol < 3; iCol++)
				rkInverse[iRow][iCol] *= fInvDet;
		}
		// (_newvp0, _newvp1,_newvp2 ) = - inverst(A) * [d, g, i]T
		_newvp0 = -(rkInverse[0][0]*d +rkInverse[0][1]*g + rkInverse[0][2]*i);
		_newvp1 = -(rkInverse[1][0]*d +rkInverse[1][1]*g + rkInverse[1][2]*i);
		_newvp2 = -(rkInverse[2][0]*d +rkInverse[2][1]*g + rkInverse[2][2]*i);
		//std::cout << _newvp0 << " " << _newvp1 << " " << _newvp2 << std::endl;
		return true;
	}
private:
	// 一个Quadric就是用10个标量来表示
	Scalar a, b, c, d, 
		      e, f, g, 
		         h, i, 
		            j;
};


/// Quadric using floats
typedef QuadricT<float> Quadricf;

/// Quadric using double
typedef QuadricT<double> Quadricd;

} // end of namespace DGP
//=============================================================================
#endif // QUADRIC_HH defined
//=============================================================================
