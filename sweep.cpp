#include "sweep.h"
#include "hamiltonian.h"
#include "scf.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <random>

// Random Hermitian perturbation with Frobenius norm = epsilon.
static Mat12 random_hermitian_perturbation(double epsilon, unsigned seed) {
    std::mt19937 rng(seed);
    std::normal_distribution<double> dist(0.0, 1.0);

    Mat12 B;
    for (int i = 0; i < 12; i++)
        for (int j = 0; j < 12; j++)
            B(i, j) = cd(dist(rng), dist(rng));

    Mat12 H = (B + B.adjoint()) / 2.0;
    H *= epsilon / H.norm();
    return H;
}

// Internal order per layer block: yz=0, xz=1, xy=2 (up) | yz=3, xz=4, xy=5 (dn)
void apply_symmetry_breaking(Mat12& rho, double delta) {
    //PO Phase
    /*
    // Layer 1 spin-up: raise dxy, lower dyz, seed dxy-dxz coherence
    rho(2, 2) += delta;   rho(0, 0) -= delta;
    rho(2, 1) += delta;   rho(1, 2) += delta;
    // Layer 2 spin-up: same orbital breaking
    rho(8, 8) += delta;   rho(6, 6) -= delta;
    rho(8, 7) += delta;   rho(7, 8) += delta;
    // Layer 1 spin-down: add electrons (layer-AF seed)
    rho(3, 3) += delta;   rho(4, 4) += delta;   rho(5, 5) += delta;
    // Layer 2 spin-down: remove electrons (layer-AF seed)
    rho(9,  9)  -= delta; rho(10, 10) -= delta; rho(11, 11) -= delta;
    */

    //SYM Phase
    // Layer 1 spin-up: raise dxy, lower dyz, seed dxy-dxz coherence
    rho(2, 2) -= delta;   rho(0, 0) -= delta;
    rho(2, 1) += delta;   rho(1, 2) += delta;
    // Layer 2 spin-up: raise dxy, lower dyz, seed dxy-dxz coherence
    rho(8, 8) -= delta;   rho(6, 6) -= delta;
    rho(8, 7) += delta;   rho(7, 8) += delta;


}

MCAResult compute_MCA(double S0, double alpha, int grid, double T, double N_target,
                      double delta, Params p, KanamoriParams kp) {
    // Stoner bootstrap at [001] — exchange field always along z
    p.theta = 0.0;
    p.phi   = 0.0;
    std::cout << "=== Stoner bootstrap ===\n";
    const CalcResult stoner = runSelfCalc(S0, alpha, grid, T, N_target, p);
    const Eigensystem sys0  = compute_eigensystem_grid(stoner.S_new, grid, p);
    Mat12 rho0 = compute_density_matrix(sys0, stoner.mu, T);
    apply_symmetry_breaking(rho0, delta);

    // Kanamori at [001]
    std::cout << "\n=== Kanamori SCF: [001] ===\n";
    const KanamoriResult res_001 = runKanamoriSCF(rho0, alpha, grid, T, N_target, p, kp,
                                                  MixerType::LinearDIIS);
    std::cout << "\n=== [001] Occupations ===\n";
    printKanamoriOccupations(res_001);

    // Kanamori at [110], seeded from converged [001] rho
    p.theta = M_PI / 2.0;
    p.phi   = M_PI / 4.0;
    std::cout << "\n=== Kanamori SCF: [110] (seeded from [001]) ===\n";
    
    Mat12 rho110SP = res_001.rho + random_hermitian_perturbation(delta, 12345);  // small random perturbation to break any residual symmetries
    const KanamoriResult res_110 = runKanamoriSCF(rho110SP, alpha, grid, T, N_target, p, kp,
                                                  MixerType::Broyden);
    std::cout << "\n=== [110] Occupations ===\n";
    printKanamoriOccupations(res_110);

    const double E_MCA = res_110.E_total - res_001.E_total;
    std::cout << "\nE_MCA = E[110] - E[001] = " << E_MCA << " eV\n";

    return {res_001, res_110, E_MCA};
}

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

        //save_band_structure(S_final, 300, p, bs_filename.str(), mu);
        save_dos(S_final, grid, T, N_target, p, dos_filename.str());
        std::cout << "\n";
    }

    outfile.close();
    std::cout << "Results saved to out/stoner_U_sweep.csv\n";
}

void run_MCA_lam_sweep(double S0, double alpha, int grid, double T, double N_target,
                       double lam_min, double lam_max, int N_points, double delta,
                       Params p, KanamoriParams kp) {

    std::filesystem::create_directories("out");

    std::ofstream outfile("out/mca_lam_sweep.csv");
    outfile << std::fixed << std::setprecision(6);

    outfile << "# tpi      = " << p.tpi        << "\n";
    outfile << "# tdelta = " << p.tdelta   << "\n";
    outfile << "# t2xy      = " << p.t2xy        << "\n";
    outfile << "# t2yzxz      = " << p.t2yzxz        << "\n";
    outfile << "# tg = " << p.tg << "\n";
    outfile << "# U       = " << kp.U        << "\n";
    outfile << "# U'      = " << kp.U_prime  << "\n";
    outfile << "# J       = " << kp.J        << "\n";
    outfile << "# t_perp  = " << p.t_perp    << "\n";

    outfile << "lam,moment_110,E_110,moment_001,E_001,E_MCA\n";

    auto kanamori_moment = [](const KanamoriResult& kres) {
        double up = 0.0, dn = 0.0;
        for (int m = 0; m < 3; m++) {
            up += kres.rho(m,   m  ).real() + kres.rho(m+6, m+6).real();
            dn += kres.rho(m+3, m+3).real() + kres.rho(m+9, m+9).real();
        }
        return up - dn;
    };

    for (int i = 0; i < N_points; i++) {
        p.lam = lam_min + i * (lam_max - lam_min) / (N_points - 1);
        std::cout << "--- lam = " << p.lam << " (" << i+1 << "/" << N_points << ") ---\n";

        const MCAResult mca = compute_MCA(S0, alpha, grid, T, N_target, delta, p, kp);

        outfile << p.lam << ","
                << kanamori_moment(mca.res_110) << "," << mca.res_110.E_total << ","
                << kanamori_moment(mca.res_001) << "," << mca.res_001.E_total << ","
                << mca.E_MCA << "\n";
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

KanamoriResult runKanamoriSCF_random(unsigned seed, double S0, double alpha, int grid_size,
                                     double T, double N_target,
                                     const Params& p, const KanamoriParams& kp,
                                     double epsilon) {
    // Stoner bootstrap along z for a physically valid starting density matrix
    std::cout << "=== Stoner bootstrap ===\n";
    const CalcResult stoner = runSelfCalc(S0, alpha, grid_size, T, N_target, p);
    const Eigensystem sys0  = compute_eigensystem_grid(stoner.S_new, grid_size, p);
    Mat12 rho0 = compute_density_matrix(sys0, stoner.mu, T);

    rho0 += random_hermitian_perturbation(epsilon, seed);

    std::cout << "=== Kanamori SCF (random seed=" << seed
              << ", epsilon=" << epsilon << ") ===\n";
    return runKanamoriSCF(rho0, alpha, grid_size, T, N_target, p, kp);
}