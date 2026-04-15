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
    p.lam     = 0.1;   // typical upper-end 3d SOC (~50-100 meV for Co/Ni with t1~0.5 eV)
    p.U       = 2;   // well into ferromagnetic phase
    p.t_perp  = 0.3;   // interlayer hopping (yz and xz only)

    const double S0       = 0.3;
    const double alpha    = 0.2;
    const double T        = 0.05;
    const double N_target = 10.0;  // 5 electrons per layer x 2 layers
    const int    grid     = 200;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=== Stoner-SCF: MCA energy ===\n";
    std::cout << "Parameters:\n";
    std::cout << "  t1      = " << p.t1      << "\n";
    std::cout << "  t_delta = " << p.t_delta << "\n";
    std::cout << "  t2      = " << p.t2      << "\n";
    std::cout << "  lam     = " << p.lam     << "\n";
    std::cout << "  U       = " << p.U       << "\n";
    std::cout << "  t_perp  = " << p.t_perp  << "\n";
    std::cout << "  S0      = " << S0        << "\n";
    std::cout << "  alpha   = " << alpha     << "\n";
    std::cout << "  T       = " << T         << "\n";
    std::cout << "  N       = " << N_target  << "\n";
    std::cout << "  grid    = " << grid << " x " << grid << "\n\n";

    // run_U_sweep(S0, alpha, grid, T, N_target, 0.0, 4.0, 50, p);
    run_MCA_lam_sweep(S0, alpha, grid, T, N_target, 0.01, 0.5, 20, p);

    return 0;
}