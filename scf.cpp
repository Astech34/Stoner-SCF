#include "scf.h"
#include <cmath>
#include <omp.h>

Eigensystem compute_eigensystem_grid(double S, int grid_size, const Params& p) {
    const int N = grid_size * grid_size;

    Eigensystem result;
    result.evals.resize(N);
    result.evecs.resize(N);

    // Precompute k-points (linspace from -pi to pi)
    std::vector<double> k_lin(grid_size);
    for (int i = 0; i < grid_size; i++)
        k_lin[i] = -M_PI + i * (2.0 * M_PI / (grid_size - 1));

    // Precompute SOC once — it's k-independent
    const Mat6 Hsoc = SOC(p.lam);
    const Mat6 Hubbard = HubbardU(S, p.U);

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N; idx++) {
        const int ix = idx % grid_size;
        const int iy = idx / grid_size;

        const double kx = k_lin[ix];
        const double ky = k_lin[iy];

        // Build Hamiltonian for this k-point
        const Mat6 H = H0(kx, ky, p) + Hsoc + Hubbard;

        // Diagonalize — SelfAdjointEigenSolver is not shared, safe per thread
        Eigen::SelfAdjointEigenSolver<Mat6> solver(H);

        result.evals[idx] = solver.eigenvalues();
        result.evecs[idx] = solver.eigenvectors();
    }

    return result;
}