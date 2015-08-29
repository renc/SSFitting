// EquationSolver.h: interface for the CEquationSolver class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_EQUATIONSOLVER_H__7E4DB6B8_3CEF_4E6E_92BE_D3820B013614__INCLUDED_)
#define AFX_EQUATIONSOLVER_H__7E4DB6B8_3CEF_4E6E_92BE_D3820B013614__INCLUDED_

#include "xVec3.h"	// Added by ClassView
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CEquationSolver  
{
public:
	CEquationSolver(int maxIterationTimes,double epsilon);
	CEquationSolver();
	virtual ~CEquationSolver();

	void OneIteration(double *A, xVec3d *x, xVec3d *b, int n);
	void OneIteration(double *A, double *x, double *b, int n);
	void ConjugateGradient(double *A,xVec3d *x,xVec3d *b,int dimension);
	void nDCopyVector(double *s,double *t,int n);
	void nDScalarDotVector(double alpha,double * a,int n,double *c);
	void nDVectorSubstract(double *a,double *b,int n,double *c);
	void nDVectorAdd(double *a,double *b,int n,double *c);
	void ConjugateGradient(double *A, double *x,double *q,int dimension);

private:
	double CalculateError(double *x, double *y, int n);
	void CopyVector(double *a,double *b,int n);
	double m_epsilon;
	int m_MaxIterationTimes;
	void   nDMatrixDotVector(double *m,double * v,int n,double *result);
	double nDInnerProduct(double *a,double *b,int n);
};

#endif // !defined(AFX_EQUATIONSOLVER_H__7E4DB6B8_3CEF_4E6E_92BE_D3820B013614__INCLUDED_)
