#include "sweep.h"
#include "hamiltonian.h"
#include "scf.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <filesystem>

void run_U_sweep(double S0, double alpha, int grid, double T, double N_target,
                 double U_min, double U_max, int N_points, Params p) {

    std::filesystem::create_directories("out");
    std::filesystem::create_directories("out/bs_plots");
    std::filesystem::create_directories("out/dos");

    std::ofstream outfile("out/stoner_U_sweep.csv");
    outfile << std::fixed << std::setprecision(6);
    outfile << "U,S_final,E_total\n";

    for (int i = 0; i < N_points; i++) {
        p.U = U_min + i * (U_max - U_min) / (N_points - 1);
        std::cout << "--- U = " << p.U << " (" << i+1 << "/" << N_points << ") ---\n";

        auto result        = runSelfCalc(S0, alpha, grid, T, N_target, p);
        const double S_final   = result.S_new;
        const double mu        = result.mu;
        const double E_current = result.E_total;

        outfile << p.U << "," << S_final << "," << E_current << "\n";
        outfile.flush();

        std::ostringstream bs_filename;
        bs_filename << "out/bs_plots/bs_"
                    << std::setw(3) << std::setfill('0') << i
                    << "_U" << std::fixed << std::setprecision(3) << p.U
                    << ".csv";

        std::ostringstream dos_filename;
        dos_filename << "out/dos/dos_"
                     << std::setw(3) << std::setfill('0') << i
                     << "_U" << std::fixed << std::setprecision(3) << p.U
                     << ".csv";

        save_band_structure(S_final, 300, p, bs_filename.str(), mu);
        save_dos(S_final, grid, T, N_target, p, dos_filename.str());
        std::cout << "\n";
    }

    outfile.close();
    std::cout << "Results saved to out/stoner_U_sweep.csv\n";
}