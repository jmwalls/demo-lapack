#include <cstdlib>
#include <iostream>

#include "linalg.h"

using namespace std;

int main(int argc, char *argv[])
{
    cout << "LAPACK demo test..." << endl;

    Eigen::MatrixXd A (6,6);
    A  << 5.,  2., -3., -1.,  0.,  0.,
          2.,  7., -1., -3.,  0.,  0.,
         -3., -1.,  6.,  2., -3., -1.,
         -1., -3.,  2.,  6., -1., -3.,
          0.,  0., -3., -1.,  3.,  1.,
          0.,  0., -1., -3.,  1.,  3.;

    vector<double> evals;
    sym_band_eigvals (A, 3, &evals);
    if (evals.size ()) {
        cout << "found eigenvalues:" << endl << "\t";
        for (auto e : evals) cout << e << " ";
        cout << endl;
    }

    double det = sym_band_det (A, 3);
    cout << "determinant: " << det << endl;

    return EXIT_SUCCESS;
}
