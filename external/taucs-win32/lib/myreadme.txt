>2007.09.10
一共七个库libtaucs.lib libmetis.lib vcf2c.lib liblapack.lib libcblas.lib libatlas.lib libf77blas.lib
但是其中            
libmetis.lib libtaucs.lib vcf2c.lib 这三个好像是MT生成的, 和MDd生成的OpenMesh库有冲突, 
因此Course Example的Fairing, Parameterization，Repair和Deformation都不能运行，
所以我从cgal-3.2.1中将MDd版本的拷贝过来。
但是OpenMesh库需要用MD生成，而其他所有10项目用MDd生成就可以了。
