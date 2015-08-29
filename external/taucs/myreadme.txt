
> 2008.05.05
If the sparse matrix is positive definite, Taucs,if not umfpack.
Taucs����������direct solve the spd matrix linear system.
taucsû��sparse matrix��manipulation(inverse, transform, multiplication, etc.).

taucs_linsolve��������ʹ�����ⷽ����(taucsֻ����linear������)�ĺ���, ��ϸ���doc.

taucs��vector���ͺ���, ����Ҳ��doc.

CGAL_wrapperĿ¼���Ǵ�cgal�����Ĵ���.

> 2008.02.27 
�÷�: ǰ�ᵱȻ��������taucs library.
(1) ���Լ�����Ŀ��#include "taucsaddon.h", link��Ŀtaucs���ɵ�taucs.lib
(2) ������÷���, �������ɵĴ���.exe�������ϴ�. ������һ���÷�������
���Լ�����Ŀ��#include "taucsaddon.h", �����mulsymatrix.cpp���Լ�����Ŀ,
����linkʱ��ֻ��Ҫlink��7��lib�Ϳ�����, ������link taucs.lib.

> 2008.02.23 ����example�е�main.cpp��
��Ax=b
	// data structure for non-zero elements of the matrix A
	std::vector<std::map<int, taucsType> > data_A(num_of_columns);
����data_Aֻ���A�еķ���Ԫ

		data_A[0]	data_A[1]	data_A[2]	...	data_A[j]	data_A[num_of_columns -1]	һ��num_of_columns��, ��i��data_A[j]�д��A�е�j�еķ���Ԫ
0		@    
1		@			@
2		@			@			@
...		...
i		@			@			@				8
m-1		@			@			@				@			@

Aһ��m��, ����A�ĵ�j�е�i�е���ֵ��8, ��ô����data_A[j][i] = 8.
��ʵstd::map<int, double> data_A[j]; Ҳ����˵data_A[j]��������һ��map(���Կ�����һά����)


> 2008.02.23
taucs.lib,   Release MT
libf77blas.lib libcblas.lib libatlas.lib liblapack.lib libmetis.lib vcf2c.lib libtaucs.lib(��7��lib����˵ԭ������, Ҳ������meshcourse07��MT-versionĿ¼�µ�)
Note: ԭ���Դ���7��lib�ҷ����ļ���lib.original����, �Ҿ�����7��������meshcourse07��MT-versionĿ¼�µ�Ӧ����һ����.
example.exe, Release MT, ������tmp.cpp, Ҳû��ignore libcmt.lib
����, MT�Ǳ���libcmt.lib��, ���Ǵ�error��ʾ�ɼ�tmp.cpp��ĳ������libcmt.lib�ظ�������,
����ֻ����Ҫlibcmt.lib��������tmp.cpp

taucs.lib,   Release MT
libf77blas.lib libcblas.lib libatlas.lib liblapack.lib libmetis-vc71-mt.lib vcf2c-vc71-mt.lib libtaucs-vc71-mt.lib ��������cgal�н������Ե�.
example.exe, Release MT, �Ӳ�����tmp.cpp������, ��Ҫignore libcmt.lib 

taucs.lib,   Release MD
libf77blas.lib libcblas.lib libatlas.lib liblapack.lib libmetis.lib vcf2c.lib libtaucs.lib(��7��lib��meshcourse07��MDd-MD-versionĿ¼�µ�)
example.exe, Release MD,  �Ӳ�����tmp.cpp������, Ҳû��ignore libcmt.lib

taucs.lib,   Debug MDd
libf77blas.lib libcblas.lib libatlas.lib liblapack.lib libmetis.lib vcf2c.lib libtaucs.lib(��7��lib��meshcourse07��MDd-MD-versionĿ¼�µ�)
example.exe, Debug MDd,  �Ӳ�����tmp.cpp������, Ҳû��ignore libcmt.lib, �����ɳɹ���, ��������ʱ����"û���ҵ�MSVCP80D.dll"�Ĵ���.


> 2007.12.23
��������汾����ԭ��վ�ǹٷ��汾. 
����汾��������ԭ���������˵�����extension, ����������taucsaddon.h��mulsymatrix.cpp, �����˹��ھ����ת�ú�����������˵�����.
����ԭ����
(libf77blas.lib vcf2c.lib libcblas.lib libatlas.lib libtaucs.lib liblapack.lib libmetis.lib)
һ��7��.lib�������µ�taucs.lib.
���������ϵ��÷�, 
To use TAUCS in your Visual Studio .NET you should:

1.      Add the TAUCS static library project to your solution.

2.      In your own project's properties, set the Runtime Library field (found under C/C++ -> Code Generation) to Multi-threaded Debug DLL
        (Ҳ����Debug Multithreading DLL, MDd).

3.      In your own project's properties, add LIBCMT to the Ignore Specific Library field (found under Linker -> Input).

4.      In your own project's properties, add taucs.lib to Linker -> Input -> Additional Dependencies.

5.      In your own project's properties, add the directory where taucs.lib is found (typically ../taucs/Debug) to Linker -> General -> Additional Library Directories.

6.      Add the file tmp.cpp to your project.

