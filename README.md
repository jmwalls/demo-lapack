##LAPACK examples
-----

Occasionally we wish to compute the eigenvalues of sparse PSD matrices.
Unfortunately, this is often a prohibatively expensive operation. General
C/C++ libraries do not support this operation either (e.g., Eigen, GSL).

Fortunately, the problem I'm interested in has additional structure---the
matrices are block tridiagonal (or approximately band diagonal). An eigenvalue
solver for banded matrices is available in LAPACK. So let's test it out.

###To build/execute/clean
-----

$ make  
$ ./lapack-demo  
$ make clean

The python script `demo.py` is provided for verification.

###Authors
-----
Jeff Walls \<jmwalls@umich.edu\>
