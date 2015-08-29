#pragma warning( disable : 4786 )
#pragma warning( disable : 4503 )

#include <cmath>
#include <map>
#include <algorithm>
#include <vector>
#include <iostream>
using namespace std;

#include "taucsaddon.h"

typedef std::map<std::pair<int,int>, taucsType> Pos2ValueMap;
inline std::pair<int,int> GetPos(int i, int j);

inline std::pair<int,int> GetPos(int i, int j) {
	return (i > j)? std::pair<int,int>(j, i)  :  std::pair<int,int>(i, j);
}

// 下面两个函数注意: matA在taucs_ccs_create函数建立都都不可以使用除TAUCS_DOUBLE之外的TAUCS_SYMMETRIC和TAUCS_LOWER.
// 也就是说matA, matB都不是对称的,也不是因为是对称矩阵所以已经只存储了下三角的数. 也就是说matA和matB只是一般的矩阵.
// 如果matA是matB的转置的话, 其结果是对称正定的, 只存储其下三角, 那么就用下面第二个函数, 
// taucs_ccs_create函数建立其结果矩阵是用到的flags是TAUCS_DOUBLE|TAUCS_SYMMETRIC|TAUCS_LOWER. 
// 这个结果的ccs矩阵可以用在taucs_linsolve()的第一参数中用于求解了. 
// 而下面第一个函数其结果矩阵在taucs_ccs_create函数建立是flags只是TAUCS_DOUBLE, 还不能用在taucs_linsolve()函数中.
// 其结果即使也是对称的话, 也不会只存下三角, 也就是还不满足taucs_linsolve()函数的要求.
// Assuming nothing about the result (the result is not stored symmetric)
taucs_ccs_matrix *Mul2NonSymmetricMatrices(const taucs_ccs_matrix *matA,
                                           const taucs_ccs_matrix *matB) {
	// Compatibility of dimensions        
	//if (matA->m != matB->n) { std::cerr << "Error: Mul2NonSymmetricMatrices, matA->m != matB->n.\n"; return NULL; }
    //上面那句是本来的测试语句, 用来确定matA和matB矩阵可以相乘, 但是其测试好像错了, 下面是我加入的测试语句.
	if (matA->n != matB->m) { std::cerr << "Error: Mul2NonSymmetricMatrices, matA->n != matB->m.\n"; return NULL; }//A的列数应该等于B的行数.

	if ((matA->flags & TAUCS_SYMMETRIC) || //这个如果是真的话, 代表matA的标志包含了TAUCS_SYMMETRIC
		(matB->flags & TAUCS_LOWER) ) 
	{	std::cerr << "Error: Mul2NonSymmetricMatrices, flags errors.\n"; return NULL; }
		
	
	// (m x n)*(n x k) = (m x k)
	int m=matA->m;
	int n=matA->n;
	int k=matB->n;
	
	taucsType biv, valA;
	int rowInd, rowA;
	std::vector<std::map<int, taucsType> > rowsC(k);
	for (int i=0; i<k; ++i) {
		// creating column i of C
		std::map<int, taucsType> & mapRow2Val = rowsC[i];
		// travel on bi
		for (int rowptrBi = matB->colptr[i];rowptrBi < matB->colptr[i+1];++rowptrBi) {
			rowInd = matB->rowind[rowptrBi];
			biv = matB->taucs_values[rowptrBi];
			// make biv*a_{rowInd} and insert into mapRow2Val
			for (int rowptrA=matA->colptr[rowInd];rowptrA<matA->colptr[rowInd+1];++rowptrA) {
				rowA=matA->rowind[rowptrA];
				valA=matA->taucs_values[rowptrA];
				// insert valA*biv into map
				std::map<int, taucsType>::iterator it = mapRow2Val.find(rowA);
				if (it == mapRow2Val.end()) {
					// first time
					mapRow2Val[rowA] = valA*biv;
				}
				else {
					it->second = it->second + valA*biv;
				}
			}
		}
		// now column i is created
	}
	//std::cout << "Here?.\n";
	return CreateTaucsMatrixFromColumns(rowsC,m,TAUCS_DOUBLE);
}

// for usage when it's known that the result is symmetric,
// like A^T * A  // multiply two non-symmetric matrix(当然也可以是symmetric的), the result 必须已知is symmetric. 
taucs_ccs_matrix *Mul2NonSymmMatSymmResult(const taucs_ccs_matrix *matA,
                                           const taucs_ccs_matrix *matB) {
	// Compatibility of dimensions        
	if ((matA->m != matB->n) || (matA->n != matB->m))
		return NULL;
	
	if ((matA->flags & TAUCS_SYMMETRIC) || //真的话, 表示matA是一个对称矩阵, 但是为什么是对称矩阵就错误呢
		(matB->flags & TAUCS_LOWER) ) {
			std::cerr << "Error: Mul2NonSymmMatSymmResult, 1.\n";
		return NULL;
	}
	
	// (m x n)*(n x m) = (m x m)
	int m=matA->m;
	int n=matA->n;
	
	taucsType biv, valA;
	int rowInd, rowA;
	std::vector<std::map<int, taucsType> > rowsC(m);//这m表示列数, 这个结构是以列为主序的,
	for (int i=0; i<m; ++i) {
		// creating column i of C
		std::map<int, taucsType> & mapRow2Val = rowsC[i];
		// travel on bi
		for (int rowptrBi = matB->colptr[i];rowptrBi < matB->colptr[i+1];++rowptrBi) {
			rowInd = matB->rowind[rowptrBi];
			biv = matB->taucs_values[rowptrBi];
			// make biv*a_{rowInd} and insert into mapRow2Val
			// Ignore anything above the diagonal!!
			for (int rowptrA=matA->colptr[rowInd];rowptrA<matA->colptr[rowInd+1];++rowptrA) {
				rowA=matA->rowind[rowptrA];
				if (rowA >= i) {
					valA=matA->taucs_values[rowptrA];
					// insert valA*biv into map
					std::map<int, taucsType>::iterator it = mapRow2Val.find(rowA);
					if (it == mapRow2Val.end()) {
						// first time
						mapRow2Val[rowA] = valA*biv;
					}
					else {
						it->second = it->second + valA*biv;
					}
				}
			}
		}
		// now column i is created
	}
	
	return CreateTaucsMatrixFromColumns(rowsC,m,TAUCS_DOUBLE|TAUCS_SYMMETRIC|TAUCS_LOWER);
}


/// Computes the transpose of a matrix.应该没有错
taucs_ccs_matrix *MatrixTranspose(const taucs_ccs_matrix *mat) {
	taucs_ccs_matrix* ret;
	ret = taucs_ccs_create(mat->n, mat->m, mat->colptr[mat->n], mat->flags);
	if (! ret) { std::cerr << "Error: MatrixTranspose, create ret.\n"; // 我加入的测试语句
		return NULL;
	}

	
	// symmetric matrix -> just copy the matrix, 可见对称矩阵的求转置是很快的.
	if (mat->flags & TAUCS_SYMMETRIC) { //真的话, 表示这是一个对称矩阵. 这与下面的MatrixCopy函数一样.
		if (mat->n != mat->m) std::cerr << "Error: MatrixTranspose, 对称的话应该行列数一样.\n"; // 我加入的测试语句

		memcpy(ret->colptr, mat->colptr, sizeof(int) * (mat->n + 1)); // colptr数组的长度是(mat->n + 1).
		memcpy(ret->rowind, mat->rowind, sizeof(int) * (mat->colptr[mat->n]));//非零元nnz的个数 = mat->colptr[mat->n]
		memcpy(ret->taucs_values, mat->taucs_values, sizeof(taucsType) * (mat->colptr[mat->n]));

		return ret;
	}

	// non-symmetric matrix -> need to build data structure. // 显示怎么从ccs个数转换成正常的行列结构
	// we'll go over the columns and build the rows
	std::vector< std::vector<int> >       rows(mat->m);  //m行, 每一行对应那些列有非零数值
	std::vector< std::vector<taucsType> > values(mat->m);//m行中每一行的非零数值
	for (int c = 0; c < mat->n; ++c) {
		for (int rowi = mat->colptr[c]; rowi < mat->colptr[c+1]; ++rowi) {
			rows[mat->rowind[rowi]].push_back(c);                        //第mat->rowind[rowi]行中的第c列的数值为mat->taucs_values[rowi]
			values[mat->rowind[rowi]].push_back(mat->taucs_values[rowi]);
		}
	}

	// copying the rows as columns in ret
	int cind = 0; //cind中非零元的个数, 同时也是ret->colptr数组中的值
	for (int r = 0; r < mat->m; ++r) { //mat中的第r行, 对应于ret中的第r列
		ret->colptr[r] = cind;

		for (int j = 0; j < (int)rows[r].size(); ++j) {
			ret->rowind[cind]		= rows[r][j];
			ret->taucs_values[cind]	= values[r][j];
			cind++;
		}
	}
	ret->colptr[mat->m] = cind; //nnz
	if (cind != mat->colptr[mat->n]) std::cerr << "Error: MatrixTranspose, 两个矩阵的非零元不一样.";
//	assert(cind == mat->colptr[mat->n]);

	return ret;
}

// 获得一个taucs's compressed columns format matrix
// Note: 输入参数cols本来就是以列为主序的结构而不是常用的行列结构.应该没有错
taucs_ccs_matrix * CreateTaucsMatrixFromColumns(const std::vector< std::map<int,taucsType> > & cols, 
												int nRows,
												int flags) { // flags数据类型
	// count nnz:
	int nCols = (int)cols.size();//多少列

	int nnz = 0;//数组中非零元的个数 number of non-zero
	for (int counter=0; counter < nCols; ++counter) {
		nnz += (int)cols[counter].size();//这一列cols[counter]中有cols[counter].size()这么多个非零元
	}
	
	taucs_ccs_matrix *matC = taucs_ccs_create(nRows,nCols,nnz,flags);//nRows行数, nCols列数, nnz矩阵中非零元的个数, flags数据类型
	if (! matC) //只是根据类型和大小分配空间, 里面还没有有效数值的.
		return NULL;
	
	// copy cols into matC
	std::map<int,taucsType>::const_iterator rit;
	int rowptrC = 0;
	for (int c=0;c<nCols;++c) {   //循环第c列
		matC->colptr[c] = rowptrC;//新的第c列从values数组和rowind数组的第rowptrC位(0-based)开始

		for (rit = cols[c].begin();rit!= cols[c].end();++rit) {      //存放第cols[c]列中的每一个元素rit,这些元素都是非零元.  
			matC->rowind[rowptrC]=rit->first;                        //每一个元素rit是第几行的,
			//int ind = rit->first; //和下面那句一样,我觉得没有用的
			matC->taucs_values[rowptrC]=rit->second;                 //每一个元素rit的内容是什么
			//double val = rit->second;
			++rowptrC;
		}
	}
	matC->colptr[nCols]=nnz;
	return matC;
}
// 是上面函数的逆.
void CreateColumnsFromTaucsMatrix(const taucs_ccs_matrix *mat, // input
								  std::vector<std::map<int, taucsType> > & cols, int nRows) //output.
{
	if (mat->n != cols.size()) {
		std::cerr << "Error: taucs, the culumns size is no the same.\n"; return;
	}
	for (int i = 0; i < mat->n; ++i) { // the i column
		int rowind_begin = mat->colptr[i];
		int rowind_end = mat->colptr[i+1];// this column is from rowind[begin] to rowind[end-1];
		for (int j = rowind_begin; j < rowind_end; ++j) {
			cols[i][mat->rowind[j]] = mat->taucs_values[j];
		}
	}
	int nnz = 0; 
	for (int i = 0; i < cols.size(); ++i) {
		nnz += (int)cols[i].size();
	}
	if (nnz != mat->colptr[mat->n]) {
		std::cerr << "Error: taucs, the nnz is not the same.\n"; return;
	}
}
// Multiplies matA by x and stores the result in b. 
// Assumes all memory has been allocated and the sizes match; 
// assumes matA is not symmetric!!
// 求解matA * x, 结果存于b中. 必须是matA的行数=b的行数, matA的列数=x的行数. 这里只是计算一维的x和b. 应该没有错
void MulNonSymmMatrixVector(const taucs_ccs_matrix *matA,
					        const taucsType * x,
							taucsType * b) {
// make b all zero
	memset(b, 0, matA->m * sizeof(taucsType));

	for (int col = 0; col < matA->n; ++col) {
		// going over column col of matA, multiplying
		// it by x[col] and setting the appropriate values
		// of vector b
		for (int p = matA->colptr[col]; p < matA->colptr[col+1]; ++p) {
			b[matA->rowind[p]] += x[col]*matA->taucs_values[p];
		}
	}
}

taucs_ccs_matrix * MatrixCopy(const taucs_ccs_matrix *mat) {
	taucs_ccs_matrix* ret;
	ret = taucs_ccs_create(mat->m, mat->n, mat->colptr[mat->n], mat->flags);
	if (! ret)
		return NULL;

	
	memcpy(ret->colptr, mat->colptr, sizeof(int) * (mat->n + 1));
	memcpy(ret->rowind, mat->rowind, sizeof(int) * (mat->colptr[mat->n]));
	memcpy(ret->taucs_values, mat->taucs_values, sizeof(taucsType) * (mat->colptr[mat->n]));

	return ret;
}


void PrintTaucsCCSMatrix(const taucs_ccs_matrix *mat)
{
	std::cout << "PrintTaucsCCSMatrix: ";
	int n = mat->n;//列数
	int nnz = mat->colptr[n];//非零元个数.
	std::cout << "m: " << mat->m << ", n: " << n << ", nnz: " << nnz << ";\n";
	for (int i = 0; i < nnz; ++i) {
		std::cout << mat->taucs_values[i] << ", ";
	}
	std::cout << " End.\n";
}

bool EqualTest(const taucs_ccs_matrix *matA, const taucs_ccs_matrix *matB)
{	// 这个函数是很耗时的, 非必要不要用.
	std::cout << "EqualTest: ";
	if (matA->m != matB->m) { std::cout << "numbers of row not equal.\t"; return false; }
	if (matA->n != matB->n) { std::cout << "numbers of column not equal.\t"; return false; }
	if (matA->colptr[matA->n] != matB->colptr[matB->n]) { std::cout << "numbers of nnz not equal.\t"; return false; }
	int nnz = matA->colptr[matA->n];
	for (int i = 0; i < nnz; ++nnz) {
		if (fabs(matA->taucs_values[i] - matB->taucs_values[i]) > 1e-8) {
			std::cout << "values not equal.\n";
			return false;
		}
	}
	std::cout << "End.\n";
	return true;
}