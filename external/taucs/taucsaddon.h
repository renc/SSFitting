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
	// ������Ŀ���ǰ�����taucsaddon.h���ͷ�ļ���.cpp�ж�����C����ʽ����taucs������ĺ�����.
#include "taucs.h" //cpp�����ж���ôʹ��.
}

#include <vector>
#include <map>



typedef double taucsType; //��ʵ��taucs.h�о���typedef double taucs_double;���ﶨ��taucsTypeֻ��Ϊ��Ĭ��ʹ��double����.
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

//��������չ��, Ŀ����Ϊ�˰������������������ݿ���.
void PrintTaucsCCSMatrix(const taucs_ccs_matrix *mat);
//this is extended by me, �Ƚ�����ccs�����Ƿ����, double�ıȽϽϷ�
bool EqualTest(const taucs_ccs_matrix *matA, const taucs_ccs_matrix *matB);

#endif // __MY_TAUCS_ADDON


