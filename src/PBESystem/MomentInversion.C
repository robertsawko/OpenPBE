#include "MomentInversion.H"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using Eigen::VectorXd;
using Eigen::MatrixXd;

pair wheeler_inversion(const VectorXd &moments)
{
    int N = moments.size() / 2;  // number of nodes of the quadrature approximation

    MatrixXd sigma = MatrixXd::Zero(N + 1, 2 * N);
    sigma.row(1) = moments.head(2 * N);

    VectorXd a = VectorXd::Zero(N);
    VectorXd b = VectorXd::Zero(N);
    a(0) = moments(1) / moments(0);
    b(0) = moments(0);  // This value is insignificant as it's not being used.

    for(int k = 1; k < N; ++k){
        int Ncol = 2 * (N - k);
        sigma.block(k + 1, k, 1, Ncol)
            = sigma.block(k, k + 1, 1, Ncol)
            - a(k - 1) * sigma.block(k, k, 1, Ncol)
            - b(k - 1) * sigma.block(k - 1, k, 1, Ncol);

        a(k) = -sigma(k, k) / sigma(k, k - 1)
            +
            sigma(k+1, k+1) / sigma(k+1, k);

        b(k) = sigma(k+1, k) / sigma(k, k-1);
    }

    VectorXd b_diag = - b.array().abs().sqrt(); //-Eigen::sqrt(Eigen::abs(b));
    //TODO: find a tridiagonal representation
    MatrixXd jacobi = a.asDiagonal();
    for(int i = 0; i < b_diag.size() - 1; ++i){
        jacobi(i, i + 1) = b_diag(i + 1);
        jacobi(i + 1, i) = b_diag(i + 1);
    }

    // eigval, eigvec = np.linalg.eig(jacobi)
    //VectorXd eivals = A.selfadjointView<Eigen::Lower>().eigenvalues();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenresult(jacobi);
    VectorXd abcissas = eigenresult.eigenvalues();
    VectorXd weights = moments(0) * eigenresult.eigenvectors().row(0).array().pow(2);

    return pair{abcissas, weights};
}
