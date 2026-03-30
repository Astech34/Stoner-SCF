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
    p.U       = 3.0;   // well into ferromagnetic phase

    const double S0       = 0.3;
    const double alpha    = 0.2;
    const double T        = 0.05;
    const double N_target = 5.0;
    const int    grid     = 200;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=== Stoner-SCF: MCA energy ===\n";
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

    // --- [0,0,1] direction: theta=0, phi=0 ---
    p.theta = 0.0;
    p.phi   = 0.0;
    std::cout << "=== n̂ = [0,0,1] ===\n";
    auto res_z = runSelfCalc(S0, alpha, grid, T, N_target, p);

    // --- [1,1,0] direction: theta=pi/2, phi=pi/4 ---
    p.theta = M_PI / 2.0;
    p.phi   = M_PI / 4.0;
    std::cout << "\n=== n̂ = [1,1,0] ===\n";
    auto res_110 = runSelfCalc(S0, alpha, grid, T, N_target, p);

    // --- MCA energy ---
    const double E_MCA = res_110.E_total - res_z.E_total;
    std::cout << "\n=== MCA Energy ===\n";
    std::cout << "  E([1,1,0]) = " << res_110.E_total << "  (S = " << res_110.S_new << ")\n";
    std::cout << "  E([0,0,1]) = " << res_z.E_total   << "  (S = " << res_z.S_new   << ")\n";
    std::cout << "  E_MCA = E([1,1,0]) - E([0,0,1]) = " << E_MCA << " t1\n";

    // run_U_sweep(S0, alpha, grid, T, N_target, 0.0, 4.0, 50, p);

    return 0;
}