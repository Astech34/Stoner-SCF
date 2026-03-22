#include <iostream>
#include <iomanip>
#include <filesystem>
#include <utility>
#include <fstream>
#include "hamiltonian.h"
#include "scf.h"
#include "sweep.h"

int main() {
    Params p;
    p.t1      = 1.0;
    p.t_delta = 0.1;
    p.t2      = 0.1;
    p.lam     = 0.0;
    p.U       = 0.5;

    const double S0       = 0.2;
    const double alpha    = 0.2;
    const double T        = 0.05;
    const double N_target = 5.0;
    const int    grid     = 200;

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

    std::cout << "=== Stoner-SCF: U sweep ===\n";
    std::cout << "U in [0.000, 4.000], 50 points\n\n";

    run_U_sweep(S0, alpha, grid, T, N_target, 0.0, 4.0, 50, p);

    return 0;
}