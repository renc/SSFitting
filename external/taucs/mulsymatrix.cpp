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

// ������������ע��: matA��taucs_ccs_create������������������ʹ�ó�TAUCS_DOUBLE֮���TAUCS_SYMMETRIC��TAUCS_LOWER.
// Ҳ����˵matA, matB�����ǶԳƵ�,Ҳ������Ϊ�ǶԳƾ��������Ѿ�ֻ�洢�������ǵ���. Ҳ����˵matA��matBֻ��һ��ľ���.
// ���matA��matB��ת�õĻ�, �����ǶԳ�������, ֻ�洢��������, ��ô��������ڶ�������, 
// taucs_ccs_create�������������������õ���flags��TAUCS_DOUBLE|TAUCS_SYMMETRIC|TAUCS_LOWER. 
// ��������ccs�����������taucs_linsolve()�ĵ�һ���������������. 
// �������һ����������������taucs_ccs_create����������flagsֻ��TAUCS_DOUBLE, ����������taucs_linsolve()������.
// ������ʹҲ�ǶԳƵĻ�, Ҳ����ֻ��������, Ҳ���ǻ�������taucs_linsolve()������Ҫ��.
// Assuming nothing about the result (the result is not stored symmetric)
taucs_ccs_matrix *Mul2NonSymmetricMatrices(const taucs_ccs_matrix *matA,
                                           const taucs_ccs_matrix *matB) {
	// Compatibility of dimensions        
	//if (matA->m != matB->n) { std::cerr << "Error: Mul2NonSymmetricMatrices, matA->m != matB->n.\n"; return NULL; }
    //�����Ǿ��Ǳ����Ĳ������, ����ȷ��matA��matB����������, ��������Ժ������, �������Ҽ���Ĳ������.
	if (matA->n != matB->m) { std::cerr << "Error: Mul2NonSymmetricMatrices, matA->n != matB->m.\n"; return NULL; }//A������Ӧ�õ���B������.

	if ((matA->flags & TAUCS_SYMMETRIC) || //����������Ļ�, ����matA�ı�־������TAUCS_SYMMETRIC
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
// like A^T * A  // multiply two non-symmetric matrix(��ȻҲ������symmetric��), the result ������֪is symmetric. 
taucs_ccs_matrix *Mul2NonSymmMatSymmResult(const taucs_ccs_matrix *matA,
                                           const taucs_ccs_matrix *matB) {
	// Compatibility of dimensions        
	if ((matA->m != matB->n) || (matA->n != matB->m))
		return NULL;
	
	if ((matA->flags & TAUCS_SYMMETRIC) || //��Ļ�, ��ʾmatA��һ���Գƾ���, ����Ϊʲô�ǶԳƾ���ʹ�����
		(matB->flags & TAUCS_LOWER) ) {
			std::cerr << "Error: Mul2NonSymmMatSymmResult, 1.\n";
		return NULL;
	}
	
	// (m x n)*(n x m) = (m x m)
	int m=matA->m;
	int n=matA->n;
	
	taucsType biv, valA;
	int rowInd, rowA;
	std::vector<std::map<int, taucsType> > rowsC(m);//��m��ʾ����, ����ṹ������Ϊ�����,
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


/// Computes the transpose of a matrix.Ӧ��û�д�
taucs_ccs_matrix *MatrixTranspose(const taucs_ccs_matrix *mat) {
	taucs_ccs_matrix* ret;
	ret = taucs_ccs_create(mat->n, mat->m, mat->colptr[mat->n], mat->flags);
	if (! ret) { std::cerr << "Error: MatrixTranspose, create ret.\n"; // �Ҽ���Ĳ������
		return NULL;
	}

	
	// symmetric matrix -> just copy the matrix, �ɼ��Գƾ������ת���Ǻܿ��.
	if (mat->flags & TAUCS_SYMMETRIC) { //��Ļ�, ��ʾ����һ���Գƾ���. ���������MatrixCopy����һ��.
		if (mat->n != mat->m) std::cerr << "Error: MatrixTranspose, �ԳƵĻ�Ӧ��������һ��.\n"; // �Ҽ���Ĳ������

		memcpy(ret->colptr, mat->colptr, sizeof(int) * (mat->n + 1)); // colptr����ĳ�����(mat->n + 1).
		memcpy(ret->rowind, mat->rowind, sizeof(int) * (mat->colptr[mat->n]));//����Ԫnnz�ĸ��� = mat->colptr[mat->n]
		memcpy(ret->taucs_values, mat->taucs_values, sizeof(taucsType) * (mat->colptr[mat->n]));

		return ret;
	}

	// non-symmetric matrix -> need to build data structure. // ��ʾ��ô��ccs����ת�������������нṹ
	// we'll go over the columns and build the rows
	std::vector< std::vector<int> >       rows(mat->m);  //m��, ÿһ�ж�Ӧ��Щ���з�����ֵ
	std::vector< std::vector<taucsType> > values(mat->m);//m����ÿһ�еķ�����ֵ
	for (int c = 0; c < mat->n; ++c) {
		for (int rowi = mat->colptr[c]; rowi < mat->colptr[c+1]; ++rowi) {
			rows[mat->rowind[rowi]].push_back(c);                        //��mat->rowind[rowi]���еĵ�c�е���ֵΪmat->taucs_values[rowi]
			values[mat->rowind[rowi]].push_back(mat->taucs_values[rowi]);
		}
	}

	// copying the rows as columns in ret
	int cind = 0; //cind�з���Ԫ�ĸ���, ͬʱҲ��ret->colptr�����е�ֵ
	for (int r = 0; r < mat->m; ++r) { //mat�еĵ�r��, ��Ӧ��ret�еĵ�r��
		ret->colptr[r] = cind;

		for (int j = 0; j < (int)rows[r].size(); ++j) {
			ret->rowind[cind]		= rows[r][j];
			ret->taucs_values[cind]	= values[r][j];
			cind++;
		}
	}
	ret->colptr[mat->m] = cind; //nnz
	if (cind != mat->colptr[mat->n]) std::cerr << "Error: MatrixTranspose, ��������ķ���Ԫ��һ��.";
//	assert(cind == mat->colptr[mat->n]);

	return ret;
}

// ���һ��taucs's compressed columns format matrix
// Note: �������cols������������Ϊ����Ľṹ�����ǳ��õ����нṹ.Ӧ��û�д�
taucs_ccs_matrix * CreateTaucsMatrixFromColumns(const std::vector< std::map<int,taucsType> > & cols, 
												int nRows,
												int flags) { // flags��������
	// count nnz:
	int nCols = (int)cols.size();//������

	int nnz = 0;//�����з���Ԫ�ĸ��� number of non-zero
	for (int counter=0; counter < nCols; ++counter) {
		nnz += (int)cols[counter].size();//��һ��cols[counter]����cols[counter].size()��ô�������Ԫ
	}
	
	taucs_ccs_matrix *matC = taucs_ccs_create(nRows,nCols,nnz,flags);//nRows����, nCols����, nnz�����з���Ԫ�ĸ���, flags��������
	if (! matC) //ֻ�Ǹ������ͺʹ�С����ռ�, ���滹û����Ч��ֵ��.
		return NULL;
	
	// copy cols into matC
	std::map<int,taucsType>::const_iterator rit;
	int rowptrC = 0;
	for (int c=0;c<nCols;++c) {   //ѭ����c��
		matC->colptr[c] = rowptrC;//�µĵ�c�д�values�����rowind����ĵ�rowptrCλ(0-based)��ʼ

		for (rit = cols[c].begin();rit!= cols[c].end();++rit) {      //��ŵ�cols[c]���е�ÿһ��Ԫ��rit,��ЩԪ�ض��Ƿ���Ԫ.  
			matC->rowind[rowptrC]=rit->first;                        //ÿһ��Ԫ��rit�ǵڼ��е�,
			//int ind = rit->first; //�������Ǿ�һ��,�Ҿ���û���õ�
			matC->taucs_values[rowptrC]=rit->second;                 //ÿһ��Ԫ��rit��������ʲô
			//double val = rit->second;
			++rowptrC;
		}
	}
	matC->colptr[nCols]=nnz;
	return matC;
}
// �����溯������.
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
// ���matA * x, �������b��. ������matA������=b������, matA������=x������. ����ֻ�Ǽ���һά��x��b. Ӧ��û�д�
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
	int n = mat->n;//����
	int nnz = mat->colptr[n];//����Ԫ����.
	std::cout << "m: " << mat->m << ", n: " << n << ", nnz: " << nnz << ";\n";
	for (int i = 0; i < nnz; ++i) {
		std::cout << mat->taucs_values[i] << ", ";
	}
	std::cout << " End.\n";
}

bool EqualTest(const taucs_ccs_matrix *matA, const taucs_ccs_matrix *matB)
{	// ��������Ǻܺ�ʱ��, �Ǳ�Ҫ��Ҫ��.
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