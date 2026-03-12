#include <iostream>
#include <Eigen/Dense>

int main() {
    // Test Eigen: build a simple 6x6 complex Hermitian matrix and diagonalize it
    using Mat6 = Eigen::Matrix<std::complex<double>, 6, 6>;

    Mat6 H = Mat6::Zero();
    H.diagonal() << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

    Eigen::SelfAdjointEigenSolver<Mat6> solver(H);

    std::cout << "Eigenvalues:\n" << solver.eigenvalues().transpose() << "\n";
    std::cout << "\nEigen is working correctly.\n";

    // Test OpenMP
    #ifdef _OPENMP
        std::cout << "OpenMP is enabled. Max threads: " << omp_get_max_threads() << "\n";
    #else
        std::cout << "OpenMP is NOT enabled.\n";
    #endif

    return 0;
}