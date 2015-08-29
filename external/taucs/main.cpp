#include <iostream>
//#include "../taucs/taucsaddon.h"
#include "taucsaddon.h"
#pragma comment (lib, ".\\release\\taucs.lib")

int main(int argc, char *argv)
{ 
	// this program solves a linear system Ax=b using TAUCS 
	// 虽然常用是taucs中使用对称正定矩阵, 这里的A应该可以不对称的, 因为下面自己求了At, 使用了法方程
	// 也就是说当A是正定时候, 就可以用taucs了.

	// dimensions
	int num_of_rows = 10;
	int num_of_columns = 10;

	// data structure for non-zero elements of the matrix A 这个数组只放非零元
	std::vector<std::map<int, taucsType> > data_A(num_of_columns);

	// fill in values (in this case, an upper triangular matrix这里以对称矩阵为例, 并且我觉得是下三角阵)
	int i = 0;
	for (i = 0; i < num_of_rows; ++i)
		for (int j = 0; j < num_of_columns; ++j)
		{
			if (i >= j)
				data_A[j][i] = 1; //等同于 data_A[j]这一列中, 第i行的非零元为1 //注意这里只填写非零元
			// notice inversion of indices while accessing data_A
		}

	// create a TAUCS matrix, representing A
	taucs_ccs_matrix *A = CreateTaucsMatrixFromColumns(data_A, num_of_rows, TAUCS_DOUBLE);

	// TAUCS requires A to be symmetric positive definite, so we will solve A^T * A * x = A^T * b instead
    
	taucs_ccs_matrix *A_t = MatrixTranspose(A);
	taucs_ccs_matrix *A_tA = Mul2NonSymmMatSymmResult(A_t, A);//std::cout << "a";
	//taucs_ccs_matrix *A_tA = Mul2NonSymmetricMatrices(A_t, A);//std::cout << "a"; //这两句是我加入的测试, 看是否效果一样, 
	//A_tA->flags += (TAUCS_SYMMETRIC + TAUCS_LOWER);                               //但结果不是, 这导致出错.

	// build 2 different constraints vectors b    //这里给出的例子中b和x都是m行(num_of_rows行), 并且是2维(2列)的
	taucsType *b = new taucsType[2 * num_of_rows];
	for (i = 0; i < num_of_rows; ++i)
		b[i] = i + 1;

	for (i = 0; i < num_of_rows; ++i)
		b[i + num_of_rows] = num_of_rows - i;

	// multiply to get A^T * b
	taucsType *A_tb = new taucsType[2 * num_of_rows];
	MulNonSymmMatrixVector(A_t, b, A_tb); //必须是A_t的行数=A_tb的行数, A_t的列数=b的行数. 这里只是计算一维的b.
	MulNonSymmMatrixVector(A_t, b + num_of_rows, A_tb + num_of_rows);//这里计算b的第二维.

	// result vector
	taucsType *x = new taucsType[2 * num_of_rows];

	// call TAUCS to solve the system
	char* solver_options[] = { "taucs.factor.LLT=true", NULL };
	int error_code = taucs_linsolve(A_tA, NULL, 2, x, A_tb, solver_options, NULL);//这个2是表示x和A_tb都是二维的, 关键是求出第一个参数A_tA
	if (error_code != TAUCS_SUCCESS)
	{
		std::cerr << "Solver failed (" << error_code << ")." << std::endl;
	}
	else
	{
		std::cout << "Equations solved:" << std::endl << "x1 = ";
		for (i = 0; i < num_of_rows; ++i)
			std::cout << x[i] << " ";
		std::cout << std::endl << "x2 = ";
		for (i = 0; i < num_of_rows; ++i)
			std::cout << x[i + num_of_rows] << " ";
		std::cout << std::endl;
	}
	
	// free up the memory
	taucs_ccs_free(A);
	taucs_ccs_free(A_t);
	taucs_ccs_free(A_tA);
	delete [] b;
	delete [] x;
	delete [] A_tb;

	// -------------
	std::cout << "\n"; //下面是一个例子为了测试Mul2NonSymmetricMatrices这个函数.
	//    M1 * A1 = m1a1;
	//[ 1  0  0]  [1 2 2]  [1  2  2]
	//[-4  1  0] *[4 4 2] =[0 -4 -6] //注意了, 结果中的这两个零是计算出来的,
	//[-4  0  1]  [4 6 4]  [0 -2 -4] //是不应该存到结果的ccs个数矩阵里面的.
	std::vector<std::map<int, taucsType> > A1(3);
	A1[0][0] = 1; A1[0][1] = 4; A1[0][2] = 4; //第0列的数据
	A1[1][0] = 2; A1[1][1] = 4; A1[1][2] = 6; //第1列的数据
	A1[2][0] = 2; A1[2][1] = 2; A1[2][2] = 4; //第2列的数据
	std::vector<std::map<int, taucsType> > M1(3);
	M1[0][0] = 1; M1[0][1] =-4; M1[0][2] =-4;
	M1[1][1] = 1;
	M1[2][2] = 1;
	taucs_ccs_matrix *a1 = CreateTaucsMatrixFromColumns(A1, 3, TAUCS_DOUBLE);
	taucs_ccs_matrix *m1 = CreateTaucsMatrixFromColumns(M1, 3, TAUCS_DOUBLE);
	taucs_ccs_matrix *m1a1 = Mul2NonSymmetricMatrices(m1, a1); //这个函数求出的结果好像是有错, 结果中好像有非零元,这应该是不对的啊.
	for (int i = 0; i < m1a1->colptr[m1a1->n]; ++i)
	{
		std::cout << m1a1->taucs_values[i] << ", ";
	} std::cout << "\n";
	// 对函数Mul2NonSymmetricMatrices的测试结果是: 计算是没有问题的, 但是如果在计算中(某一行乘以某一列再相加)产生的结果刚好是零,
	// 这个零也被加入到结果的ccs矩阵中, 但是ccs形式的矩阵里面是不该有零的.
	double *b1 = new double[3]; b1[0] = 3; b1[1] = 6; b1[2] = 10;
	double *m1b1 = new double[3]; 
	MulNonSymmMatrixVector(m1, b1, m1b1);//这个函数就没有问题
	std::cout << m1b1[0] << ", " << m1b1[1] << ", " << m1b1[2] << ".\n "; 

	std::cout << "\n"; //下面是一个例子为了测试Mul2NonSymmMatSymmResult这个函数.
	//[ 1  1]  [1  1]   [2  0]
	//[ 1 -1] *[1 -1] = [0  2] //注意了, 结果中的这两个零是计算出来的, 是不应该存到结果的ccs个数矩阵里面的.
	std::vector<std::map<int, taucsType> > A2(2);
	A2[0][0] = 1; A2[0][1] = 1; //第0列的数据
	A2[1][0] = 1; A2[1][1] =-1; //第1列的数据
	std::vector<std::map<int, taucsType> > M2(2);
	M2[0][0] = 1; M2[0][1] = 1;
	M2[1][0] = 1; M2[1][1] =-1; 
	taucs_ccs_matrix *a2 = CreateTaucsMatrixFromColumns(A2, 2, TAUCS_DOUBLE);
	taucs_ccs_matrix *m2 = CreateTaucsMatrixFromColumns(M2, 2, TAUCS_DOUBLE);
	taucs_ccs_matrix *m2a2 = Mul2NonSymmMatSymmResult(m2, a2); //这个函数求出的结果好像是有错, 结果中好像有非零元,这应该是不对的啊.
	for (int i = 0; i < m2a2->colptr[m2a2->n]; ++i)
	{
		std::cout << m2a2->taucs_values[i] << ", ";
	} std::cout << "\n";

	getchar();
	return 0;
}
/* 猜想, Mul2NonSymmMatSymmResult函数是不是也会遇到Mul2NonSymmetricMatrices这个函数的特殊情况呢?
[a   c]  [a   b]   [a2+c2  ab+cd]
[b   d] *[c   d] = [ab+cd  b2+d2]
那么假如a, b, c, d之前都不是零, 因此存在ccs矩阵里面, 但是刚好结果ab+cd=0呢? 是否也存于结果的矩阵里面呢? 这不是错了吗?
结论: 无论是Mul2NonSymmMatSymmResult函数还是Mul2NonSymmetricMatrices函数, 只要计算的结果不会刚好是零的话, 那就没有问题,
否则这个零值也会被存在ccs矩阵里面, 这是不和要求的.
*/