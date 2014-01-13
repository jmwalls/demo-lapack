#include <cstdlib>
#include <iostream>

using namespace std;

/*
 * DSBGV computes all the eigenvalues, and optionally, the eigenvectors of a
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
*              tridiagonal form did not converge to zero;   
*           > N:   if INFO = N + i, for 1 <= i <= N, then DPBSTF   
*              returned INFO = i: B is not positive definite.   
*              The factorization of B could not be completed and   
*              no eigenvalues or eigenvectors were computed.
*/
extern "C" {
void dsbgv_ (const char *jobz, const char *uplo, const int *n, 
        const int *ka, const int *kb, 
        double *AB, const int *ldab, 
        double *BB, const int *ldbb, 
        double *w, 
        double *z__, const int *ldz, 
        double *work, int *info);
}

int main(int argc, char *argv[])
{
    cout << "LAPACK demo test..." << endl;

    // A = 6x6 matrix
    //       5  2 -3 -1  0  0
    //       2  7 -1 -3  0  0
    //      -3 -1  6  2 -3 -1
    //      -1 -3  2  6 -1 -3
    //       0  0 -3 -1  3  1
    //       0  0 -1 -3  1  3
    // B = 6x6 identity matrix

    const char jobz = 'N', uplo = 'U';

    const int N = 6; 
    const int ka = 3, lda = 4, kb = 0, ldb = 1;

    double x=0;
    double A[N][lda] = {{ x,  x,  x, 5}, 
                        { x,  x,  2, 7},
                        { x, -3, -1, 6}, 
                        {-1, -3,  2, 6},
                        { 0, -3, -1, 3}, 
                        {-1, -3,  1, 3}};
    double B[N][ldb] = {{1},{1},{1},{1},{1},{1}};

    double w[N] = {0};

    const int ldz = N;
    double z[ldz][N] = {{0}};

    double work[3*N] = {0};
    int info;

    dsbgv_ (&jobz, &uplo, &N,
            &ka, &kb,
            &A[0][0], &lda, &B[0][0], &ldb,
            w, &z[0][0], &ldz, work, &info);

    if (info==0) {
        cout << "eigenvalues of A: " << endl;
        for (int i=0;i<N;++i)
            cout << w[i] << " ";
        cout << endl;
    } else {
        cerr << "dsbgv failed: " << info << endl;
    }

    return EXIT_SUCCESS;
}
