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
                       double lam_min, double lam_max, int N_points, Params p, KanamoriParams kp) {

    std::filesystem::create_directories("out");

    std::ofstream outfile("out/mca_lam_sweep.csv");
    outfile << std::fixed << std::setprecision(6);

    outfile << "# t1      = " << p.t1        << "\n";
    outfile << "# t_delta = " << p.t_delta   << "\n";
    outfile << "# t2      = " << p.t2        << "\n";
    outfile << "# U       = " << kp.U        << "\n";
    outfile << "# U'      = " << kp.U_prime  << "\n";
    outfile << "# J       = " << kp.J        << "\n";
    outfile << "# t_perp  = " << p.t_perp    << "\n";

    outfile << "lam,moment_110,E_110,moment_001,E_001,E_MCA\n";

    auto kanamori_moment = [](const KanamoriResult& kres) {
        double up = 0.0, dn = 0.0;
        for (int m = 0; m < 3; m++) {
            up += kres.rho(m,   m  ).real() + kres.rho(m+6,   m+6  ).real();
            dn += kres.rho(m+3, m+3).real() + kres.rho(m+9, m+9).real();
        }
        return up - dn;
    };

    for (int i = 0; i < N_points; i++) {
        p.lam = lam_min + i * (lam_max - lam_min) / (N_points - 1);
        std::cout << "--- lam = " << p.lam << " (" << i+1 << "/" << N_points << ") ---\n";

        // [0,0,1]: theta=0, phi=0
        p.theta = 0.0;
        p.phi   = 0.0;
        const CalcResult stoner_001 = runSelfCalc(S0, alpha, grid, T, N_target, p);
        const Eigensystem sys_001   = compute_eigensystem_grid(stoner_001.S_new, grid, p);
        const Mat12 rho0_001        = compute_density_matrix(sys_001, stoner_001.mu, T);
        const KanamoriResult kres_001 = runKanamoriSCF(rho0_001, alpha, grid, T, N_target, p, kp);

        // [1,1,0]: theta=pi/2, phi=pi/4
        p.theta = M_PI / 2.0;
        p.phi   = M_PI / 4.0;
        const CalcResult stoner_110 = runSelfCalc(S0, alpha, grid, T, N_target, p);
        const Eigensystem sys_110   = compute_eigensystem_grid(stoner_110.S_new, grid, p);
        const Mat12 rho0_110        = compute_density_matrix(sys_110, stoner_110.mu, T);
        const KanamoriResult kres_110 = runKanamoriSCF(rho0_110, alpha, grid, T, N_target, p, kp);

        const double E_MCA = kres_110.E_total - kres_001.E_total;
        std::cout << "E_MCA = " << E_MCA << "\n\n";

        outfile << p.lam << ","
                << kanamori_moment(kres_110) << "," << kres_110.E_total << ","
                << kanamori_moment(kres_001) << "," << kres_001.E_total << ","
                << E_MCA << "\n";
        outfile.flush();
    }

    outfile.close();
    std::cout << "Results saved to out/mca_lam_sweep.csv\n";
}

void run_delta_V_sweep(double S0, double alpha, int grid, double T, double N_target,
                       double dV_min, double dV_max, int N_points,
                       Params p, KanamoriParams kp) {

    std::filesystem::create_directories("out");

    // Stoner bootstrap — delta_V does not enter this path, so run once
    std::cout << "=== delta_V sweep: Stoner bootstrap ===\n";
    p.delta_V = 0.0;
    const CalcResult stoner = runSelfCalc(S0, alpha, grid, T, N_target, p);
    const Eigensystem sys0  = compute_eigensystem_grid(stoner.S_new, grid, p);
    const Mat12 rho_stoner  = compute_density_matrix(sys0, stoner.mu, T);
    std::cout << "\n";

    std::ofstream outfile("out/delta_V_sweep.csv");
    outfile << std::fixed << std::setprecision(6);
    outfile << "delta_V,n_layer1,n_layer2,delta_n\n";

    for (int i = 0; i < N_points; i++) {
        p.delta_V = dV_min + i * (dV_max - dV_min) / (N_points - 1);
        std::cout << "--- delta_V = " << p.delta_V
                  << " (" << i+1 << "/" << N_points << ") ---\n";

        const KanamoriResult kres = runKanamoriSCF(rho_stoner, alpha, grid, T, N_target, p, kp);
        const Mat12& rho = kres.rho;

        double n1 = 0.0, n2 = 0.0;
        for (int m = 0;  m < 6;  m++) n1 += rho(m, m).real();
        for (int m = 6;  m < 12; m++) n2 += rho(m, m).real();
        const double dn = n1 - n2;

        std::cout << "  n_layer1 = " << n1
                  << "  n_layer2 = " << n2
                  << "  delta_n = "  << dn << "\n\n";

        outfile << p.delta_V << "," << n1 << "," << n2 << "," << dn << "\n";
        outfile.flush();
    }

    outfile.close();
    std::cout << "Results saved to out/delta_V_sweep.csv\n";
}