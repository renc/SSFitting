// EquationSolver.cpp: implementation of the CEquationSolver class.
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"
//#include "MeshOpti.h"
#include "EquationSolver.h"

//#ifdef _DEBUG
//#undef THIS_FILE
//static char THIS_FILE[]=__FILE__;
//#define new DEBUG_NEW
//#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CEquationSolver::CEquationSolver()
{
	m_MaxIterationTimes = 50;
	m_epsilon = 1e-10;
}
CEquationSolver::CEquationSolver(int maxIterationTimes,double epsilon)
{
	m_MaxIterationTimes = maxIterationTimes;
	m_epsilon = epsilon;
}

CEquationSolver::~CEquationSolver()
{

}
void CEquationSolver::ConjugateGradient(double *A, xVec3d *x, xVec3d *b, int dimension)
{
	double *tempx, *tempb;
	int i,j;
	tempx =new double[dimension];
	tempb =new double[dimension];

	for (j = 0; j < 3; j++)
	{
		for (i = 0 ; i < dimension; i++)
		{
			tempx[i] = x[i].m_e[j];
			tempb[i] = b[i].m_e[j];
		}
		ConjugateGradient(A,tempx,tempb,dimension);
		for (i = 0 ; i < dimension; i++)
		{
			x[i].m_e[j] = tempx[i];
		}	
	}
	
	delete [] tempx;
	delete [] tempb;
}

void CEquationSolver::ConjugateGradient(double *A, double *x, double *b, int dimension)
{
	int n = dimension;
	double *tempx;

	tempx =    new double[n];
	
	int count = 0;
	while(count++ < m_MaxIterationTimes)
	{
		CopyVector(x,tempx,n);
		OneIteration(A, x, b, n);
		if ( CalculateError(tempx,x,n) < m_epsilon )
			break;
	}
	
	delete [] tempx;
}

double CEquationSolver::nDInnerProduct(double *a, double *b, int n)
{
	int i;
	double x = 0;
	for (i = 0; i < n ; i++)
		x += a[i] * b[i];
	return (x);
}

void CEquationSolver::nDMatrixDotVector(double *m, double *v, int n, double *result)
{   // m是n*n的方阵, v和result都是n维向量.
	int i,j;

	double *x = new double[n];

	for (i = 0; i < n ; i++)
	{
		x[i] = 0;
		for (j = 0; j < n ;j++)
			x[i] += m[ i * n + j] * v[j];
	}
	for (i = 0; i < n ; i++)
	{
		result[i] = x[i];
	}
}

void CEquationSolver::nDVectorAdd(double *a, double *b, int n, double *c)
{
	for ( int i = 0 ; i < n; i++)
		c[i] = a[i] + b[i];
}

void CEquationSolver::nDVectorSubstract(double *a, double *b, int n, double *c)
{
	for ( int i = 0 ; i < n; i++)
		c[i] = a[i] - b[i];
}

void CEquationSolver::nDScalarDotVector(double alpha, double *a, int n,double *c)
{
	for ( int i = 0 ; i < n; i++)
		c[i] = alpha * a[i];
}

void CEquationSolver::nDCopyVector(double *s, double *t, int n)
{
	for ( int i = 0 ; i < n; i++)
		t[i] = s[i];
}

void CEquationSolver::OneIteration(double *A, double *x, double *b, int n)
{

	double *r,*p,*Ap,*temp;
	double alpha,beta,rr,app;

	r=    new double[n];
	p=    new double[n];
	Ap=   new double[n];
	temp= new double[n];
	
	// temp = A * x
	nDMatrixDotVector(A, x, n, temp);
	nDVectorSubstract(b,temp, n, r);
	nDCopyVector(r,p,n);

	rr = nDInnerProduct(r,r,n);             // rr=(r,r)

	nDMatrixDotVector(A,p,n,Ap);	        // Ap = A * p; 	
    app = nDInnerProduct(Ap,p,n);           // App = (Ap,p);
	if ( fabs (app) < 1e-20 )
	{
		//MessageBox("divided by zero.\n",NULL,MB_OK);
		return;
	}
	alpha = rr / app;                   // alpha
	nDScalarDotVector(alpha,p,n,temp);  // temp = alpha * p;
	nDVectorAdd(x,temp,n,x);            // x(j+1) = x(j) + temp;
	nDScalarDotVector(alpha,Ap,n,temp); // temp = alpha * Ap;
	nDVectorSubstract(r,temp,n,r);      // r(j+1) = r(j) - temp;
	if (fabs(rr) < 1e-20)
	{
		//MessageBox("divided by zero.\n",NULL,MB_OK);
		return;
	}
	beta = nDInnerProduct(r,r,n) / rr;   // beta = (r(j+1),r(j+1))/(r(j),r(j));
	nDScalarDotVector(beta,p,n,temp);    // temp = belta * p;
	nDVectorAdd(r,temp,n,p);             // p(j+1) = r(j+1) + temp;

	delete [] r;
	delete [] p;
	delete [] Ap;
	delete [] temp;
}

void CEquationSolver::OneIteration(double *A, xVec3d *x, xVec3d *b, int n)
{
	double *tempx, *tempb;
	int i,j;
	tempx =new double[n];
	tempb =new double[n];

	for (j = 0; j < 3; j++)
	{
		for (i = 0 ; i < n; i++)
		{
			tempx[i] = x[i].m_e[j];
			tempb[i] = b[i].m_e[j];
		}
		OneIteration(A,tempx,tempb,n);
		for (i = 0 ; i < n; i++)
		{
			x[i].m_e[j] = tempx[i];
		}	
	}
	delete [] tempx;
	delete [] tempb;

}

double CEquationSolver::CalculateError(double *x, double *y, int n)
{   
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		sum += fabs(x[i] - y[i]);
	}
	return (sum);
}

void CEquationSolver::CopyVector(double *a, double *b, int n)
{
	for (int i = 0; i < n; i++)
		b[i] = a[i];
}
