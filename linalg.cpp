#include <cmath>
#include "linalg.h"

using namespace std;

extern "C" {
/* DSBGV computes all the eigenvalues, and optionally, the eigenvectors of a
 * real generalized symmetric-definite banded eigenproblem, of the form
 *      A*x=(lambda)*B*x. 
 * Here A and B are assumed to be symmetric and banded, and B is also positive
 * definite.   
 *
 * JOBZ    (input) CHARACTER*1   
 *          = 'N':  Compute eigenvalues only;   
 *          = 'V':  Compute eigenvalues and eigenvectors.   
 * UPLO    (input) CHARACTER*1   
 *          = 'U':  Upper triangles of A and B are stored;   
 *          = 'L':  Lower triangles of A and B are stored.   
 * 
 * N       (input) INTEGER   
 *          The order of the matrices A and B.  N >= 0.   
 * 
 * KA      (input) INTEGER   
 *          The number of superdiagonals of the matrix A if UPLO = 'U',   
 *          or the number of subdiagonals if UPLO = 'L'. KA >= 0.   
 *          
 * KB      (input) INTEGER   
 *          The number of superdiagonals of the matrix B if UPLO = 'U',   
 *          or the number of subdiagonals if UPLO = 'L'. KB >= 0.   
 *
 * AB      (input/output) DOUBLE PRECISION array, dimension (LDAB, N)   
 *          On entry, the upper or lower triangle of the symmetric band   
 *          matrix A, stored in the first ka+1 rows of the array.  The   
 *          j-th column of A is stored in the j-th column of the array AB   
 *          as follows:   
 *              if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;   
 *              if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).   
 *          On exit, the contents of AB are destroyed.   
 * 
 * LDAB    (input) INTEGER   
 *          The leading dimension of the array AB.  LDAB >= KA+1.   
 * 
 * BB      (input/output) DOUBLE PRECISION array, dimension (LDBB, N)   
 *          On entry, the upper or lower triangle of the symmetric band
 *          matrix B, stored in the first kb+1 rows of the array.  The   j-th
 *          column of B is stored in the j-th column of the array BB   as
 *          follows:   
 *              if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;   
 *              if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).   
 *          On exit, the factor S from the split Cholesky factorization   
 *          B = S**T*S, as returned by DPBSTF.   
 * 
 * LDBB    (input) INTEGER   
 *          The leading dimension of the array BB.  LDBB >= KB+1.   
 * 
 * W       (output) DOUBLE PRECISION array, dimension (N)   
 *          If INFO = 0, the eigenvalues in ascending order.   
 *
 * Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)   
 *          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
 *          eigenvectors, with the i-th column of Z holding the  eigenvector
 *          associated with W(i). The eigenvectors are   normalized so that
 *          Z**T*B*Z = I. 
 *          If JOBZ = 'N', then Z is not referenced.   
 * 
 * LDZ     (input) INTEGER   
 *          The leading dimension of the array Z.  LDZ >= 1, and if   JOBZ =
 *          'V', LDZ >= N.   
 * 
 * WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)   
 *
 * INFO    (output) INTEGER   
 *          = 0:  successful exit   
 *          < 0:  if INFO = -i, the i-th argument had an illegal value   
 *          > 0:  if INFO = i, and i is:   
 *          <= N:  the algorithm failed to converge:   
 *             i off-diagonal elements of an intermediate   
 *             tridiagonal form did not converge to zero;   
 *          > N:   if INFO = N + i, for 1 <= i <= N, then DPBSTF   
 *             returned INFO = i: B is not positive definite.   
 *             The factorization of B could not be completed and   
 *             no eigenvalues or eigenvectors were computed.
 */
void dsbgv_ (const char *jobz, const char *uplo, const int *n, 
        const int *ka, const int *kb, 
        double *AB, const int *ldab, 
        double *BB, const int *ldbb, 
        double *w, 
        double *z__, const int *ldz, 
        double *work, int *info);

/* DPBTRF computes the Cholesky factorization of a real symmetric positive
 * definite band matrix A.   
 * 
 * The factorization has the form   
 *      A = U**T * U,  if UPLO = 'U', or   
 *      A = L  * L**T,  if UPLO = 'L',   
 * where U is an upper triangular matrix and L is lower triangular.   
 * 
 * UPLO    (input) CHARACTER*1   
 *          = 'U':  Upper triangle of A is stored;   
 *          = 'L':  Lower triangle of A is stored.   
 *
 * N       (input) INTEGER   
 *          The order of the matrix A.  N >= 0.   
 *
 * KD      (input) INTEGER   
 *          The number of superdiagonals of the matrix A if UPLO = 'U',   
 *          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.   
 *
 * AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)   
 *          On entry, the upper or lower triangle of the symmetric band   
 *          matrix A, stored in the first KD+1 rows of the array.  The   
 *          j-th column of A is stored in the j-th column of the array AB   
 *          as follows:   
 *              if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;   
 *              if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).   
 *          On exit, if INFO = 0, the triangular factor U or L from the   
 *          Cholesky factorization A = U**T*U or A = L*L**T of the band   
 *          matrix A, in the same storage format as A.   
 * LDAB    (input) INTEGER   
 *          The leading dimension of the array AB.  LDAB >= KD+1.   
 *
 * INFO    (output) INTEGER   
 *          = 0:  successful exit   
 *          < 0:  if INFO = -i, the i-th argument had an illegal value   
 *          > 0:  if INFO = i, the leading minor of order i is not   
 *                positive definite, and the factorization could not be   
 *                completed.
 */
int dpbtrf_ (const char *uplo, const int *n, const int *kd, 
        double *AB, const int *ldab, int *info);
}

static const char jobz = 'N', uplo = 'U';

void sym_band_eigvals (const Eigen::MatrixXd& A, size_t k,
        vector<double>* evals)
{
    eigen_assert (A.rows ()==A.cols ());
    int N = A.rows (); 
    const int ka=k, lda=k+1, kb=0, ldb=1;

    // if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;   
    double **Ab = (double**) malloc (N*sizeof *Ab + N*(lda*sizeof **Ab)),
           **Bb = (double**) malloc (N*sizeof *Bb + N*(ldb*sizeof **Bb));
    double *da = (double*) (Ab + N), *db = (double*) (Bb + N);
    for (int i=0;i<N;++i) {
        Ab[i] = da + i*lda;
        int nels = min (i+1,lda);
        for (int j=0;j<nels;++j) Ab[i][lda-j-1] = A (i-j,i);

        Bb[i] = db + i*ldb;
        Bb[i][0] = 1;
    }

    double *w = (double*) malloc (N*sizeof *w);
    double *work = (double*) malloc (3*N*sizeof *work);
    int info;

    // call lapack
    dsbgv_ (&jobz, &uplo, &N,
            &ka, &kb,
            &Ab[0][0], &lda, &Bb[0][0], &ldb,
            w, 0, &N, work, &info);
    if (info==0) {
        evals->resize (N);
        std::copy (w, w+N, evals->begin ());
    }

    free (work);
    free (w);
    free (Bb);
    free (Ab);
}

double sym_band_det (const Eigen::MatrixXd& A, size_t k)
{
    eigen_assert (A.rows ()==A.cols ());
    int N = A.rows (); 
    const int ka=k, lda=k+1;

    double **Ab = (double**) malloc (N*sizeof *Ab + N*(lda*sizeof **Ab));
    double *da = (double*) (Ab + N);
    for (int i=0;i<N;++i) {
        Ab[i] = da + i*lda;
        int nels = min (i+1,lda);
        for (int j=0;j<nels;++j) Ab[i][lda-j-1] = A (i-j,i);
    }

    int info;
    dpbtrf_ (&uplo, &N, &ka, 
        &Ab[0][0], &lda, &info);
    if (info==0) {
        double det = 1;
        for (int i=0;i<N;++i) {
            det *= Ab[i][k];
        }
        free (Ab);
        return det*det;
    }
    else {
        free (Ab);
        return -1;
    }
}
