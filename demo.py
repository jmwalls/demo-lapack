#!/usr/bin/env python
import sys
import numpy as np

if __name__ == '__main__':
    A = np.array ([[ 5.,  2., -3., -1.,  0.,  0.],
                   [ 2.,  7., -1., -3.,  0.,  0.],
                   [-3., -1.,  6.,  2., -3., -1.],
                   [-1., -3.,  2.,  6., -1., -3.],
                   [ 0.,  0., -3., -1.,  3.,  1.],
                   [ 0.,  0., -1., -3.,  1.,  3.]])
    print 'eigenvalues of A:'
    evals = np.linalg.eigvals (A)
    evals.sort ()
    print '\t',
    for ev in evals:
        print ev,
    print

    print 'determinant: ', np.linalg.det (A)

    sys.exit (0)
