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

void run_MCA_lam_sweep(double S0, double alpha, int grid, double T, double N_target,
                       double lam_min, double lam_max, int N_points, Params p) {

    std::filesystem::create_directories("out");

    std::ofstream outfile("out/mca_lam_sweep.csv");
    outfile << std::fixed << std::setprecision(6);

    // Write parameters as comments so the plot script can annotate the figure
    outfile << "# t1      = " << p.t1      << "\n";
    outfile << "# t_delta = " << p.t_delta << "\n";
    outfile << "# t2      = " << p.t2      << "\n";
    outfile << "# U       = " << p.U       << "\n";
    outfile << "# t_perp  = " << p.t_perp  << "\n";

    outfile << "lam,S_110,E_110,S_001,E_001,E_MCA\n";

    for (int i = 0; i < N_points; i++) {
        p.lam = lam_min + i * (lam_max - lam_min) / (N_points - 1);
        std::cout << "--- lam = " << p.lam << " (" << i+1 << "/" << N_points << ") ---\n";

        // [0,0,1]: theta=0, phi=0
        p.theta = 0.0;
        p.phi   = 0.0;
        auto res_001 = runSelfCalc(S0, alpha, grid, T, N_target, p);

        // [1,1,0]: theta=pi/2, phi=pi/4
        p.theta = M_PI / 2.0;
        p.phi   = M_PI / 4.0;
        auto res_110 = runSelfCalc(S0, alpha, grid, T, N_target, p);

        const double E_MCA = res_110.E_total - res_001.E_total;
        std::cout << "E_MCA = " << E_MCA << "\n\n";

        outfile << p.lam << ","
                << res_110.S_new << "," << res_110.E_total << ","
                << res_001.S_new << "," << res_001.E_total << ","
                << E_MCA << "\n";
        outfile.flush();
    }

    outfile.close();
    std::cout << "Results saved to out/mca_lam_sweep.csv\n";
}