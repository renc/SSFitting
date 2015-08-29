> 2008.03.06 从http://www.cise.ufl.edu/~xwu/下载的.
求解loop细分曲面的bounding.

________________________________________________________
This is code is being distributed without any guarantee.

any commercial use of this code need authorization from
the authors:

Xiaobin Wu  (xwu@cise.ufl.edu)
Jorg Peters (jorg@cise.ufl.edu)
SurfLab, University of Florida
Gainesville, FL, U.S.A.
________________________________________________________


---------------------------
1. Source code
---------------------------

looprange.h   : upper and lower bound of each basis function
loopdomain.h  : (u,v) positions of each domain mesh 
bnd_loop.cpp  : the routine to compute bounds for a given loop patch, 
                say x^+, x^- of x. 
                this routine uses looprange.h and loopdomain.h .

__________________________________________________________________

For newest updates, refer to the SurfLab webpage
http://www.cise.ufl.edu/research/SurfLab/

Aug. 30. 2004


