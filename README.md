This program is the implementation of the method proposed in the paper Point-Based Manifold Harmonics.

Acknowledgement:

This program is inspired greatly by the source code and help provided by Dr. Jian Sun (http://www.geomtop.org/sunjian/).
Dependency:

-BOOST C++ Library http://www.boost.org/
-LAPACK/BLAS http://www.netlib.org/lapack/
-CGAL http://www.cgal.org/
    --GMP http://gmplib.org/
    --MPFR http://www.mpfr.org/

Successfully compiled using MinGW GCC 4.5.2 with CGAL 3.8, Boost 1.46.1, GMP 5.0.1, MPFR 3.0.1.
Also known to be working with MinGW-W64.
It should work on most windows system.

Usage:

pcdqb.exe MModel(input) QFile(output) BFile(output)

It takes 1 input of 3D triangular model as M file and output 2 files for Q matrix and B matrix as specified in the paper.
Q file list all elements of the sparse matrix Q in the Matlab manner. In each row, the first number is the row-index starting
from 1; the second number is the column-index starting from 1; the third number is the non-zero element value.
B file lists the diagonal elements of the B matrix as specified in the paper.

M file format is a text-based 3D model format. A python script is provided to convert OBJ model into M format:

obj2mvonly.py objname mfilename

It works with Python 3.2.
