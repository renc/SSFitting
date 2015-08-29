
> 2008.05.05
If the sparse matrix is positive definite, Taucs,if not umfpack.
Taucs就是用来解direct solve the spd matrix linear system.
taucs没有sparse matrix的manipulation(inverse, transform, multiplication, etc.).

taucs_linsolve就是我所使用来解方程组(taucs只能是linear方程组)的函数, 详细情况doc.

taucs的vector类型很弱, 具体也看doc.

CGAL_wrapper目录下是从cgal拿来的代码.

> 2008.02.27 
用法: 前提当然是已有了taucs library.
(1) 在自己的项目中#include "taucsaddon.h", link项目taucs生成的taucs.lib
(2) 上面的用法简单, 但是生成的代码.exe体积好像较大. 所以另一个用法可以是
在自己的项目中#include "taucsaddon.h", 并添加mulsymatrix.cpp进自己的项目,
这样link时候只需要link那7个lib就可以了, 而不是link taucs.lib.

> 2008.02.23 例子example中的main.cpp中
解Ax=b
	// data structure for non-zero elements of the matrix A
	std::vector<std::map<int, taucsType> > data_A(num_of_columns);
数组data_A只存放A中的非零元

		data_A[0]	data_A[1]	data_A[2]	...	data_A[j]	data_A[num_of_columns -1]	一共num_of_columns列, 第i列data_A[j]中存放A中第j列的非零元
0		@    
1		@			@
2		@			@			@
...		...
i		@			@			@				8
m-1		@			@			@				@			@

A一共m行, 假如A的第j列第i行的数值是8, 那么设置data_A[j][i] = 8.
其实std::map<int, double> data_A[j]; 也就是说data_A[j]本来就是一个map(可以看做是一维数组)


> 2008.02.23
taucs.lib,   Release MT
libf77blas.lib libcblas.lib libatlas.lib liblapack.lib libmetis.lib vcf2c.lib libtaucs.lib(这7个lib可以说原来带的, 也可以是meshcourse07中MT-version目录下的)
Note: 原来自带的7个lib我放在文件夹lib.original里面, 我觉得这7个东西和meshcourse07中MT-version目录下的应该是一样的.
example.exe, Release MT, 不加入tmp.cpp, 也没有ignore libcmt.lib
首先, MT是必需libcmt.lib的, 但是从error提示可见tmp.cpp中某东西与libcmt.lib重复定义了,
所以只能是要libcmt.lib而不加入tmp.cpp

taucs.lib,   Release MT
libf77blas.lib libcblas.lib libatlas.lib liblapack.lib libmetis-vc71-mt.lib vcf2c-vc71-mt.lib libtaucs-vc71-mt.lib 后三个从cgal中借来试试的.
example.exe, Release MT, 加不加入tmp.cpp都可以, 需要ignore libcmt.lib 

taucs.lib,   Release MD
libf77blas.lib libcblas.lib libatlas.lib liblapack.lib libmetis.lib vcf2c.lib libtaucs.lib(这7个lib是meshcourse07中MDd-MD-version目录下的)
example.exe, Release MD,  加不加入tmp.cpp都可以, 也没有ignore libcmt.lib

taucs.lib,   Debug MDd
libf77blas.lib libcblas.lib libatlas.lib liblapack.lib libmetis.lib vcf2c.lib libtaucs.lib(这7个lib是meshcourse07中MDd-MD-version目录下的)
example.exe, Debug MDd,  加不加入tmp.cpp都可以, 也没有ignore libcmt.lib, 是生成成功了, 但是运行时后有"没有找到MSVCP80D.dll"的错误.


> 2007.12.23
这里这个版本不是原网站那官方版本. 
这个版本好像是在原基础上做了点扩充extension, 例如增加了taucsaddon.h和mulsymatrix.cpp, 扩充了关于矩阵就转置和两个矩阵相乘的运算.
利用原来的
(libf77blas.lib vcf2c.lib libcblas.lib libatlas.lib libtaucs.lib liblapack.lib libmetis.lib)
一共7个.lib生成了新的taucs.lib.
以下是网上的用法, 
To use TAUCS in your Visual Studio .NET you should:

1.      Add the TAUCS static library project to your solution.

2.      In your own project's properties, set the Runtime Library field (found under C/C++ -> Code Generation) to Multi-threaded Debug DLL
        (也就是Debug Multithreading DLL, MDd).

3.      In your own project's properties, add LIBCMT to the Ignore Specific Library field (found under Linker -> Input).

4.      In your own project's properties, add taucs.lib to Linker -> Input -> Additional Dependencies.

5.      In your own project's properties, add the directory where taucs.lib is found (typically ../taucs/Debug) to Linker -> General -> Additional Library Directories.

6.      Add the file tmp.cpp to your project.

