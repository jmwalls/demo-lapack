#ifndef __LINALG_H__
#define __LINALG_H__

#include <vector>
#include <Eigen/Dense>

// returns evals in ascending order
void sym_band_eigvals (const Eigen::MatrixXd& A, size_t k, 
        std::vector<double>* evals);

// numpy's slogdet function: the determinant is computed via LU factorization
// using the LAPACK routine z/dgetrf.
//
// dpbtrf computes the Cholesky factorization of a real symmetric positive
// definite band matrix---determinant is then the product of the diagonal
// elements of the L factor (squared).

double sym_band_det (const Eigen::MatrixXd& A, size_t k);

#endif // __LINALG_H__
