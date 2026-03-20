#include <iostream>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include "hamiltonian.h"
#include "scf.h"

int main() {
    // --- Parameters ---
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

    // --- U sweep ---
    const int    N_points = 50;
    const double U_min    = 0.0;
    const double U_max    = 4.0;

    std::cout << "=== Stoner-SCF: U sweep ===\n";
    std::cout << "U in [" << U_min << ", " << U_max << "], "
              << N_points << " points\n\n";

    // --- Output directories ---
    std::filesystem::create_directories("out");
    std::filesystem::create_directories("out/bs_plots");

    std::ofstream outfile("out/stoner_U_sweep.csv");
    outfile << std::fixed << std::setprecision(6);
    outfile << "U,S_final\n";

    for (int i = 0; i < N_points; i++) {
        p.U = U_min + i * (U_max - U_min) / (N_points - 1);

        std::cout << "--- U = " << p.U << " (" << i+1 << "/" << N_points << ") ---\n";

        const double S_final = runSelfCalc(S0, alpha, grid, T, N_target, p);

        outfile << p.U << "," << S_final << "\n";
        outfile.flush();

        // --- Save band structure for this U ---
        // Zero-pad index so files sort correctly: bs_000.csv, bs_001.csv, ...
        std::ostringstream bs_filename;
        bs_filename << "out/bs_plots/bs_"
                    << std::setw(3) << std::setfill('0') << i
                    << "_U" << std::fixed << std::setprecision(3) << p.U
                    << ".csv";

        save_band_structure(S_final, 300, p, bs_filename.str());

        std::cout << "\n";
    }

    outfile.close();
    std::cout << "Results saved to out/stoner_U_sweep.csv\n";
    return 0;
}