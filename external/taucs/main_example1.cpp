

#include <cstdio>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include <assert.h>
#define TAUCS_CORE_DOUBLE
// #include "taucs.h"
/* taucs.hԭ��ֻ�ǿ�����C����ʹ�õ�, ������Cplusplus������ʹ�õĻ���ҪҪextern "C" { #include "taucs.h" } �������������Ӵ���
���ڱ���...
1>main_example1.cpp
1>��������...
1>main_example1.obj : error LNK2001: �޷��������ⲿ���� "int __cdecl taucs_linsolve(struct taucs_ccs_matrix *,void * *,int,void *,void *,char * * const,void * * const)" (?taucs_linsolve@@YAHPAUtaucs_ccs_matrix@@PAPAXHPAX2QAPADQAPAX@Z)
1>main_example1.obj : error LNK2001: �޷��������ⲿ���� "void __cdecl taucs_logfile(char *)" (?taucs_logfile@@YAXPAD@Z)
1>main_example1.obj : error LNK2001: �޷��������ⲿ���� "void __cdecl taucs_ccs_free(struct taucs_ccs_matrix *)" (?taucs_ccs_free@@YAXPAUtaucs_ccs_matrix@@@Z)
1>main_example1.obj : error LNK2001: �޷��������ⲿ���� "int __cdecl taucs_ccs_write_ijv(struct taucs_ccs_matrix *,char *)" (?taucs_ccs_write_ijv@@YAHPAUtaucs_ccs_matrix@@PAD@Z)
1>main_example1.obj : error LNK2001: �޷��������ⲿ���� "struct taucs_ccs_matrix * __cdecl taucs_ccs_create(int,int,int,int)" (?taucs_ccs_create@@YAPAUtaucs_ccs_matrix@@HHHH@Z)
1>D:\libraries\taucs\Release\example1.exe : fatal error LNK1120: 5 ���޷��������ⲿ����
*/
extern "C" { 
#include "taucs.h" 
}

/* ����Ҳʹ��taucsaddon.h
#include "taucsaddon.h"
#pragma comment (lib, ".\\release\\taucs.lib")
1>���ڱ���...
1>main_example1.cpp
1>C:\Program Files\Microsoft Visual Studio 8\VC\include\malloc.h(240) : error C3861: 'taucs_must_not_call_free_directly': identifier not found
*/

/* Solve a trivial symmetric Ax=b system :

[ 1.0 0.5 0.0 0.0 ]       [ 1.0 ]
[ 0.5 1.0 0.5 0.0 ]       [ 2.0 ]
[ 0.0 0.5 1.0 0.5 ] X   = [ 3.0 ]
[ 0.0 0.0 0.5 1.0 ]       [ 4.0 ]         // answer X = [0 2 0 4]

values:	1.0 0.5 1.0 0.5 1.0 0.5 1.0       // nnz = 7 ������Ԫ
rowind: 0   1   1   2   2   3   3
colptr: 0   2   4   6   7                 // n = 4 ��, length(colptr) = n+1, ��colptr[0]=0, colptr[n]=nnz.
*/
int main(int argc, char* argv[])
{
	int m = 4,n = 4,nnz = 7, i, used;
	taucs_double x[4];
	taucs_double b[4];
	taucs_ccs_matrix * pMatrix;
	taucs_logfile ("stdout");//ֻ��Ϊ����taucs_linsolve()���������һЩ�����Ϣ, �Ǳ�Ҫ��.
	pMatrix = taucs_ccs_create( m, n, nnz, TAUCS_DOUBLE ); //������ʾ����ô�ֶ���ʾ��ָ����������ʼ��һ���Գƾ���
	//flagsһ��ʹ��TAUCS_DOUBLE, ��Ȼ�ǶԳƵ�lower���ں������TAUCS_SYMMETRIC��TAUCS_LOWER.
	pMatrix->colptr[0] = 0;
	pMatrix->colptr[1] = 2;
	pMatrix->colptr[2] = 4;
	pMatrix->colptr[3] = 6;
	pMatrix->colptr[4] = 7;
	used = 0;
	for (i = 0; i < 4; i++ )
	{
		pMatrix->rowind[used] = i;
		pMatrix->taucs_values[used] = 1.0;
		used++;
		if (i+1 <= 3)
		{
			pMatrix->rowind[used] = i+1;
			pMatrix->taucs_values[used] = 0.5;
			used++;
		}
	}
	pMatrix->flags += TAUCS_SYMMETRIC;
	pMatrix->flags += TAUCS_LOWER;
	i = taucs_ccs_write_ijv(pMatrix, "test.txt" );

	b[0] = 1.0; b[1] = 2.0; b[2] = 3.0; b[3] = 4.0;

	char* options[] = { "taucs.factor.LLT=true","taucs.factor.ordering=metis",NULL };
	i = taucs_linsolve (pMatrix, NULL, 1, x, b, options, NULL);
	if (i != TAUCS_SUCCESS)
	{
		printf ("Solution error.\n");
		if (i==TAUCS_ERROR)
			printf ("Generic error.");
		if (i==TAUCS_ERROR_NOMEM)
			printf ("NOMEM error.");
		if (i==TAUCS_ERROR_BADARGS)
			printf ("BADARGS error.");
		if (i==TAUCS_ERROR_MAXDEPTH)
			printf ("MAXDEPTH error.");
		if (i==TAUCS_ERROR_INDEFINITE)
			printf ("NOT POSITIVE DEFINITE error.");
	}
	else
	{
		printf ("Solution success.\n");
		for (i = 0; i < 4; i++)
			printf ("%f\n",x[i]);
	}
	//taucs_dccs_free(pMatrix);
	taucs_ccs_free(pMatrix);

	getchar();
	return 0;
}

