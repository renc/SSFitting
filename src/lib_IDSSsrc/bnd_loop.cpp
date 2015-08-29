//////////////////////////////////////////////////////////////////////
// File: bnd_loop.cpp  
//    : a sample routine to compute one component bounding a loop patch, 
//      say x^+, x^- of x.  , using looprange.h and loopdomain.h .

#include <assert.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#define REAL double // �ҽ����е�double��ֱ�Ӹ�Ϊ�˶�double
//typedef double Real;//������ð�

#define MIN_VALENCE 3
#define MAX_VALENCE 16
//const int MIN_VALENCE = 3;
//const int MAX_VALENCE = 16;

#include "bnd_loop.h"
#include "loopdomain.h"
#include "looprange.h"

////////////////////////////////////////////////////////////////
// Return triangle that bounds one (x,y or z) component of the loop patch
//
//  ( For a 3D patch, call this routine three times with 
//    the appropriate input array.
//    for example:
//      for(m=0;m<3;m++) { //m:0,1,2�ֱ��Ӧxyz��������
//        boundLoop(valence, &coeff[m], 3, &upper[m], &lower[m], 3);
//      }
//  )
//
//
// input:
//    valence: of v_0             (the valence of v_1 and v_2 has to be 6)
//    coeff  : the coefficient of each vertex arranged at the following order
//              ... - 9
//             / \ / \
//           n+5 - 0 - 8           where n = valence
//           / \ /.\ / \
//          3 - 1 - 2 - 7 
//           \ / \ / \ /
//            4 - 5 - 6
//
//    stride: the distance of each coefficient in the coeff array
//            for example: 
//            stride = 3 , if the array is packed with xyz positions, 
//            stride = 1 , if the array is packed with only x positions.
// output:
//    upper: upper triangle
//    lower: lower triangle
//    stride_slefe: the distance between each result
//

int BoundLoop(int valence, double *coeff, int stride, 
					   double* upper, double* lower, int stride_slefe)
{
	double top_end, left_end, right_end; 
	
	
	int bas = valence +3; // number of bases, //basisһ����valence+6��, ǰ����Ϊ��, ֻ�����valence+3��
	double *D, *P, *M;
	double D2b[MAX_VALENCE+3];

//    assert(valence>=MIN_VALENCE && valence <= MAX_VALENCE);

	P = Loop_Upper[valence-3];   // upper bound of basis���������bi+
	M = Loop_Lower[valence-3];   // lower bound of basis���������bi-
	D = Loop_Domain[valence-3];  // domain layout//D = &Loop_Domain[valence-3][0];
	
	///////////////////////////////////////////////////////////////
	top_end   = coeff[0*stride]; // save three end points, of the central triangle.
	left_end  = coeff[1*stride]; // �������ʾ��double coeff[]����ֵ��˳�����:
	right_end = coeff[2*stride]; // x0 y0 z0 x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4 ... xn+5 yn+5 zn+5

	// compute the differences of bi to the linear function intepolating
	// three end points. (D[2b], D[2b+1], 1-D[2b]-D[2b+1])���������������������������.
	for (int b=0;b<bas; b++) {
		//std::cout << D[2*b] << ", " << D[2*b +1] << ".\n";
		D2b[b] = //����paper�е�di. 
			coeff[(b+3)*stride]- (D[2*b] * top_end + D[2*b+1] * left_end + (1-D[2*b]-D[2*b+1]) * right_end);
	}

	// compute each break point
	for(int loc=0;loc<3;loc++) {
		double* upper_ptr = &upper[loc*stride_slefe];  // result storage, upper_pter = upper[0], upper[3], upper[6],
		double* lower_ptr = &lower[loc*stride_slefe];  //							  lower[0], lower[3], lower[6].�ֱ��Ӧv0, v1, v3������bounds��
		// 1. initialize to the linear function 
		(*upper_ptr) = (*lower_ptr) =  loc==0? top_end: (loc==1?left_end:right_end);

		// 2. add contribution from every coefficients
		for(int b=0;b<bas;b++) {
			if(D2b[b]>0) {
				(*upper_ptr) += P[b*3+loc]* D2b[b];
				(*lower_ptr) += M[b*3+loc]* D2b[b];
			} else {
				(*upper_ptr) += M[b*3+loc]* D2b[b];
				(*lower_ptr) += P[b*3+loc]* D2b[b];
			}
		}
	}

	return 0;
}

