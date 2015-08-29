#ifndef __MY_TAUCS_ADDON
#define __MY_TAUCS_ADDON

/****************************************************************/
// Sivan's library!!!
#ifndef WIN32
#include <unistd.h>
#include <pthread.h>
#endif
#define TAUCS_CORE_DOUBLE
extern "C" { // TAUCS is a C library 
	// 这样的目的是包含了taucsaddon.h这个头文件的.cpp中都是以C的形式调用taucs库里面的函数的.
#include "taucs.h" //cpp代码中都这么使用.
}

#include <vector>
#include <map>



typedef double taucsType; //其实在taucs.h中就有typedef double taucs_double;这里定义taucsType只是为了默认使用double类型.
/****************************************************************/

/// Multiplies two symmetric lower matrices. If dimensions don't
/// match or if the arguments are not good in any other way
// returns NULL. Otherwise returns the new matrix.
taucs_ccs_matrix *Mul2SymmetricMatrices(const taucs_ccs_matrix *mat0,
										const taucs_ccs_matrix *mat1);




taucs_ccs_matrix *Mul2NonSymmetricMatrices(const taucs_ccs_matrix *matA,
										   const taucs_ccs_matrix *matB);

// for usage when it's known that the result is symmetric,
// like A^T * A
taucs_ccs_matrix *Mul2NonSymmMatSymmResult(const taucs_ccs_matrix *matA,
                                           const taucs_ccs_matrix *matB);

/// Computes the transpose of a matrix.
taucs_ccs_matrix *MatrixTranspose(const taucs_ccs_matrix *mat);

taucs_ccs_matrix * CreateTaucsMatrixFromColumns(const std::vector< std::map<int,taucsType> > & cols, 
												int nRows,
												int flags);
void CreateColumnsFromTaucsMatrix(const taucs_ccs_matrix *mat, // input
								  std::vector<std::map<int, taucsType> > & cols, int nRows); //output.

// Multiplies matA by x and stores the result
// in b. Assumes all memory has been allocated
// and the sizes match; assumes matA is not
// symmetric!!
void MulNonSymmMatrixVector(const taucs_ccs_matrix *matA,
					        const taucsType * x,
							taucsType * b);

taucs_ccs_matrix * MatrixCopy(const taucs_ccs_matrix *mat);

//这是我扩展的, 目的是为了按列主序输出矩阵的内容看看.
void PrintTaucsCCSMatrix(const taucs_ccs_matrix *mat);
//this is extended by me, 比较两个ccs矩阵是否相等, double的比较较烦
bool EqualTest(const taucs_ccs_matrix *matA, const taucs_ccs_matrix *matB);

#endif // __MY_TAUCS_ADDON


