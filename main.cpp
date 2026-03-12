#include <iostream>
#include <iomanip>
#include "hamiltonian.h"
#include "scf.h"

int main() {
    // --- Parameters ---
    Params p;
    p.t1      = 1.0;
    p.t_delta = 0.1;
    p.t2      = 0.1;
    p.lam     = 0.1;
    p.U       = 4.0;

    const double S0       = 0.2;
    const double alpha    = 0.2;
    const double T        = 0.05;
    const double N_target = 5.0;
    const int    grid     = 200;

    // --- Print parameters ---
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=== Stoner-SCF ===\n";
    std::cout << "Parameters:\n";
    std::cout << "  t1      = " << p.t1      << "\n";
    std::cout << "  t_delta = " << p.t_delta << "\n";
    std::cout << "  t2      = " << p.t2      << "\n";
    std::cout << "  lam     = " << p.lam     << "\n";
    std::cout << "  U       = " << p.U       << "\n";
    std::cout << "  S0      = " << S0        << "\n";
    std::cout << "  alpha   = " << alpha     << "\n";
    std::cout << "  T       = " << T         << "\n";
    std::cout << "  N       = " << N_target  << "\n";
    std::cout << "  grid    = " << grid << " x " << grid << "\n\n";

    // --- Run SCF loop ---
    const double S_final = runSelfCalc(S0, alpha, grid, T, N_target, p);

    return 0;
}