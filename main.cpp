#include <iostream>
#include <iomanip>
#include <cmath>
#include "hamiltonian.h"
#include "scf.h"

int main() {
    // --- Parameters ---
    Params p;
    p.t1      = 1.0;
    p.t_delta = 0.1;
    p.t2      = 0.1;
    p.lam     = 0.1;
    p.U       = 2.0;

    const double S        = 0.2;   // initial Stoner magnetisation guess
    const double T        = 0.05;  // temperature (in units of t1)
    const double N_target = 5.0;   // target electron filling
    const int    grid     = 200;   // k-mesh size (grid x grid)

    // --- Print parameters ---
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "=== Stoner-SCF ===\n";
    std::cout << "Parameters:\n";
    std::cout << "  t1      = " << p.t1      << "\n";
    std::cout << "  t_delta = " << p.t_delta << "\n";
    std::cout << "  t2      = " << p.t2      << "\n";
    std::cout << "  lam     = " << p.lam     << "\n";
    std::cout << "  U       = " << p.U       << "\n";
    std::cout << "  S       = " << S         << "\n";
    std::cout << "  T       = " << T         << "\n";
    std::cout << "  N       = " << N_target  << "\n";
    std::cout << "  grid    = " << grid << " x " << grid << "\n\n";

    // --- Diagonalise over k-mesh ---
    std::cout << "Computing eigensystem..." << std::flush;
    const Eigensystem sys = compute_eigensystem_grid(S, grid, p);
    std::cout << " done.\n";

    // --- Find chemical potential ---
    std::cout << "Finding chemical potential..." << std::flush;
    const double mu = find_mu(sys, grid, T, N_target);
    std::cout << " done.\n\n";

    if (std::isnan(mu)) {
        std::cerr << "ERROR: Could not find chemical potential — Brent's method failed.\n";
        return 1;
    }

    std::cout << "Chemical potential mu = " << mu << "\n";

    return 0;
}