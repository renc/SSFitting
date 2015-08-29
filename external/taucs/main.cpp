#include <iostream>
//#include "../taucs/taucsaddon.h"
#include "taucsaddon.h"
#pragma comment (lib, ".\\release\\taucs.lib")

int main(int argc, char *argv)
{ 
	// this program solves a linear system Ax=b using TAUCS 
	// ��Ȼ������taucs��ʹ�öԳ���������, �����AӦ�ÿ��Բ��ԳƵ�, ��Ϊ�����Լ�����At, ʹ���˷�����
	// Ҳ����˵��A������ʱ��, �Ϳ�����taucs��.

	// dimensions
	int num_of_rows = 10;
	int num_of_columns = 10;

	// data structure for non-zero elements of the matrix A �������ֻ�ŷ���Ԫ
	std::vector<std::map<int, taucsType> > data_A(num_of_columns);

	// fill in values (in this case, an upper triangular matrix�����ԶԳƾ���Ϊ��, �����Ҿ�������������)
	int i = 0;
	for (i = 0; i < num_of_rows; ++i)
		for (int j = 0; j < num_of_columns; ++j)
		{
			if (i >= j)
				data_A[j][i] = 1; //��ͬ�� data_A[j]��һ����, ��i�еķ���ԪΪ1 //ע������ֻ��д����Ԫ
			// notice inversion of indices while accessing data_A
		}

	// create a TAUCS matrix, representing A
	taucs_ccs_matrix *A = CreateTaucsMatrixFromColumns(data_A, num_of_rows, TAUCS_DOUBLE);

	// TAUCS requires A to be symmetric positive definite, so we will solve A^T * A * x = A^T * b instead
    
	taucs_ccs_matrix *A_t = MatrixTranspose(A);
	taucs_ccs_matrix *A_tA = Mul2NonSymmMatSymmResult(A_t, A);//std::cout << "a";
	//taucs_ccs_matrix *A_tA = Mul2NonSymmetricMatrices(A_t, A);//std::cout << "a"; //���������Ҽ���Ĳ���, ���Ƿ�Ч��һ��, 
	//A_tA->flags += (TAUCS_SYMMETRIC + TAUCS_LOWER);                               //���������, �⵼�³���.

	// build 2 different constraints vectors b    //���������������b��x����m��(num_of_rows��), ������2ά(2��)��
	taucsType *b = new taucsType[2 * num_of_rows];
	for (i = 0; i < num_of_rows; ++i)
		b[i] = i + 1;

	for (i = 0; i < num_of_rows; ++i)
		b[i + num_of_rows] = num_of_rows - i;

	// multiply to get A^T * b
	taucsType *A_tb = new taucsType[2 * num_of_rows];
	MulNonSymmMatrixVector(A_t, b, A_tb); //������A_t������=A_tb������, A_t������=b������. ����ֻ�Ǽ���һά��b.
	MulNonSymmMatrixVector(A_t, b + num_of_rows, A_tb + num_of_rows);//�������b�ĵڶ�ά.

	// result vector
	taucsType *x = new taucsType[2 * num_of_rows];

	// call TAUCS to solve the system
	char* solver_options[] = { "taucs.factor.LLT=true", NULL };
	int error_code = taucs_linsolve(A_tA, NULL, 2, x, A_tb, solver_options, NULL);//���2�Ǳ�ʾx��A_tb���Ƕ�ά��, �ؼ��������һ������A_tA
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
	std::cout << "\n"; //������һ������Ϊ�˲���Mul2NonSymmetricMatrices�������.
	//    M1 * A1 = m1a1;
	//[ 1  0  0]  [1 2 2]  [1  2  2]
	//[-4  1  0] *[4 4 2] =[0 -4 -6] //ע����, ����е����������Ǽ��������,
	//[-4  0  1]  [4 6 4]  [0 -2 -4] //�ǲ�Ӧ�ô浽�����ccs�������������.
	std::vector<std::map<int, taucsType> > A1(3);
	A1[0][0] = 1; A1[0][1] = 4; A1[0][2] = 4; //��0�е�����
	A1[1][0] = 2; A1[1][1] = 4; A1[1][2] = 6; //��1�е�����
	A1[2][0] = 2; A1[2][1] = 2; A1[2][2] = 4; //��2�е�����
	std::vector<std::map<int, taucsType> > M1(3);
	M1[0][0] = 1; M1[0][1] =-4; M1[0][2] =-4;
	M1[1][1] = 1;
	M1[2][2] = 1;
	taucs_ccs_matrix *a1 = CreateTaucsMatrixFromColumns(A1, 3, TAUCS_DOUBLE);
	taucs_ccs_matrix *m1 = CreateTaucsMatrixFromColumns(M1, 3, TAUCS_DOUBLE);
	taucs_ccs_matrix *m1a1 = Mul2NonSymmetricMatrices(m1, a1); //�����������Ľ���������д�, ����к����з���Ԫ,��Ӧ���ǲ��Եİ�.
	for (int i = 0; i < m1a1->colptr[m1a1->n]; ++i)
	{
		std::cout << m1a1->taucs_values[i] << ", ";
	} std::cout << "\n";
	// �Ժ���Mul2NonSymmetricMatrices�Ĳ��Խ����: ������û�������, ��������ڼ�����(ĳһ�г���ĳһ�������)�����Ľ���պ�����,
	// �����Ҳ�����뵽�����ccs������, ����ccs��ʽ�ľ��������ǲ��������.
	double *b1 = new double[3]; b1[0] = 3; b1[1] = 6; b1[2] = 10;
	double *m1b1 = new double[3]; 
	MulNonSymmMatrixVector(m1, b1, m1b1);//���������û������
	std::cout << m1b1[0] << ", " << m1b1[1] << ", " << m1b1[2] << ".\n "; 

	std::cout << "\n"; //������һ������Ϊ�˲���Mul2NonSymmMatSymmResult�������.
	//[ 1  1]  [1  1]   [2  0]
	//[ 1 -1] *[1 -1] = [0  2] //ע����, ����е����������Ǽ��������, �ǲ�Ӧ�ô浽�����ccs�������������.
	std::vector<std::map<int, taucsType> > A2(2);
	A2[0][0] = 1; A2[0][1] = 1; //��0�е�����
	A2[1][0] = 1; A2[1][1] =-1; //��1�е�����
	std::vector<std::map<int, taucsType> > M2(2);
	M2[0][0] = 1; M2[0][1] = 1;
	M2[1][0] = 1; M2[1][1] =-1; 
	taucs_ccs_matrix *a2 = CreateTaucsMatrixFromColumns(A2, 2, TAUCS_DOUBLE);
	taucs_ccs_matrix *m2 = CreateTaucsMatrixFromColumns(M2, 2, TAUCS_DOUBLE);
	taucs_ccs_matrix *m2a2 = Mul2NonSymmMatSymmResult(m2, a2); //�����������Ľ���������д�, ����к����з���Ԫ,��Ӧ���ǲ��Եİ�.
	for (int i = 0; i < m2a2->colptr[m2a2->n]; ++i)
	{
		std::cout << m2a2->taucs_values[i] << ", ";
	} std::cout << "\n";

	getchar();
	return 0;
}
/* ����, Mul2NonSymmMatSymmResult�����ǲ���Ҳ������Mul2NonSymmetricMatrices������������������?
[a   c]  [a   b]   [a2+c2  ab+cd]
[b   d] *[c   d] = [ab+cd  b2+d2]
��ô����a, b, c, d֮ǰ��������, ��˴���ccs��������, ���Ǹպý��ab+cd=0��? �Ƿ�Ҳ���ڽ���ľ���������? �ⲻ�Ǵ�����?
����: ������Mul2NonSymmMatSymmResult��������Mul2NonSymmetricMatrices����, ֻҪ����Ľ������պ�����Ļ�, �Ǿ�û������,
���������ֵҲ�ᱻ����ccs��������, ���ǲ���Ҫ���.
*/