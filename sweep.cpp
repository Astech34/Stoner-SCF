#include "sweep.h"
#include "hamiltonian.h"
#include "scf.h"
#include "seeds.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <random>

// Random Hermitian perturbation with Frobenius norm = epsilon.
Mat12 random_hermitian_perturbation(double epsilon, unsigned seed) {
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

    // Assuming rho is a 2D array or matrix of std::complex<double>
    // and t2g ordering is: 0=dxy, 1=dxz, 2=dyz
    // Layer stride = 6, Spin stride = 3

    // =====================================================================
    // LAYER 1 - SPIN UP (Indices 0, 1, 2)
    // Goal: Raise dxy, lower dyz, seed complex dxz + i*dyz coherence
    // =====================================================================
    //rho(0, 0) -= delta;  // Raise dxy
    //rho(2, 2) += delta;  // Lower dyz

    // Hermiticity requires rho(j,i) = conj(rho(i,j))
    // Setting dxz-dyz coherence to a complex phase
    rho(1, 2) += std::complex<double>(0.0,  delta); // rho(dxz, dyz) = +i*delta
    rho(2, 1) += std::complex<double>(0.0, -delta); // rho(dyz, dxz) = -i*delta


    // =====================================================================
    // LAYER 2 - SPIN UP (Indices 6, 7, 8)
    // Goal: Raise dxy, lower dyz, seed complex dxz + i*dyz coherence
    // =====================================================================
    //rho(6, 6) += delta;  // Raise dxy
    //rho(8, 8) -= delta;  // Lower dyz

    // Setting dxz-dyz coherence to a complex phase
    rho(7, 8) += std::complex<double>(0.0,  delta); // rho(dxz, dyz) = +i*delta

    rho(8, 7) += std::complex<double>(0.0, -delta); // rho(dyz, dxz) = -i*delta


}

MCAResult compute_MCA(double S0, double alpha, int grid, double T, double N_target,
                      double delta, Params p, KanamoriParams kp, 
                      std::string seed001, std::string seed110) {
    // Stoner bootstrap at [001] — exchange field always along z
    p.theta = 0.0;
    p.phi   = 0.0;

    std::cout << "Seeding from: " << seed001;
    std::cout << "Seeding from: " << seed110;

    /*
    std::cout << "=== Stoner bootstrap ===\n";
    const CalcResult stoner = runSelfCalc(S0, alpha, grid, T, N_target, p);
    const Eigensystem sys0  = compute_eigensystem_grid(stoner.S_new, grid, p);
    Mat12 rho0 = compute_density_matrix(sys0, stoner.mu, T);
    apply_symmetry_breaking(rho0, delta);
    */
    Mat12 loaded_rho = load_density_matrix(seed001);
    // Kanamori at [001]
    std::cout << "\n=== Kanamori SCF: [001] ===\n";
    const KanamoriResult res_001 = runKanamoriSCF(loaded_rho, alpha, grid, T, N_target, p, kp,
                                                  MixerType::LinearDIIS);
    std::cout << "\n=== [001] Occupations ===\n";
    printKanamoriOccupations(res_001, p);
    
    //Save
    //std::string Seed110 = "/home/cmp/Documents/Github/Stoner-SCF/out/LMCA110.csv";
    std::ostringstream filename;
    filename << "out/dmatrices/L100Lam" 
             << std::fixed << std::setprecision(2) << p.lam // Controls decimal places
             << "_density_matrix.csv";
    save_density_matrix(res_001.rho, "/home/cmp/Documents/Github/Stoner-SCF/out/LMCA001.csv");

    // Kanamori at [110], seeded from converged [001] rho
    p.theta = M_PI / 2.0;
    p.phi   = M_PI / 4.0;
    std::cout << "\n=== Kanamori SCF: [110] ===\n";
    loaded_rho = load_density_matrix(seed110);
    
    //Mat12 rho110SP = res_001.rho + random_hermitian_perturbation(delta, 12345);  // small random perturbation to break any residual symmetries
    //const KanamoriResult res_110 = runKanamoriSCF(rho0, alpha, grid, T, N_target, p, kp,
                                                  //MixerType::Broyden);
    const KanamoriResult res_110 = runKanamoriSCF(loaded_rho, alpha, grid, T, N_target, p, kp,
                                                  MixerType::LinearDIIS);
    std::cout << "\n=== [110] Occupations ===\n";
    // For S calculation we want to use the unrotated spin operators.
    p.theta = 0.0;
    p.phi = 0.0;
    printKanamoriOccupations(res_110, p);

    //Save
    std::ostringstream fname110;
    fname110 << "out/dmatrices/L110Lam" 
             << std::fixed << std::setprecision(2) << p.lam // Controls decimal places
             << "_density_matrix.csv";
    save_density_matrix(res_110.rho, "/home/cmp/Documents/Github/Stoner-SCF/out/LMCA110.csv");

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

    outfile << "lam,E_110, S_110, L_110, E_001, S_001, L_001, MCA\n";
                    
    // Initial Diagonal Seed
    /*
    Mat12 rho0 = Mat12::Zero();

    double n0 = 0.5;
    double ndelta = 0.1;
    for (int i = 0; i < 12; ++i) {
    if ((i % 6) < 3) {
        rho0(i, i) = n0 + ndelta;
    } else {
        rho0(i, i) = n0 - ndelta;
    }
    */
    // Seed 001
    //Mat12 StartSeed = Mat12::Zero();
    //const KanamoriResult res_boot = runKanamoriSCF(rho0, alpha, grid, T, N_target, p, kp, MixerType::LinearDIIS);
    //StartSeed = res_boot.rho;
    //save_density_matrix(res_boot.rho, "out/dmatrices/MCAStart001_density_matrix.csv");

    //p.theta = M_PI / 2.0;
    //p.phi   = M_PI / 4.0;

    //const KanamoriResult res_boot1 = runKanamoriSCF(rho0, alpha, grid, T, N_target, p, kp, MixerType::LinearDIIS);
    //StartSeed = res_boot1.rho;
    //save_density_matrix(res_boot.rho, "out/dmatrices/MCAStart110_density_matrix.csv");
    
    for (int i = 0; i < N_points; i++) {

        p.lam = lam_min + i * (lam_max - lam_min) / (N_points - 1);
        std::cout << "--- lam = " << p.lam << " (" << i+1 << "/" << N_points << ") ---\n";

        /*
        std::string Seed001 = "out/dmatrices/MCAStart001_density_matrix.csv";
        std::string Seed110 = "out/dmatrices/MCAStart110_density_matrix.csv";

        std::ostringstream flname001;
        std::ostringstream flname110;
        // Choosing Seed
        if (i > 0){
            flname001 << "out/dmatrices/L100Lam" 
                    << std::fixed << std::setprecision(2) << p.lam // Controls decimal places
                    << "_density_matrix.csv";
            
            flname110 << "out/dmatrices/L110Lam" 
                    << std::fixed << std::setprecision(2) << p.lam // Controls decimal places
                    << "_density_matrix.csv";
            
            Seed001 = flname001.str();
            Seed110 = flname110.str();
        }
        */

        std::string Seed001 = "/home/cmp/Documents/Github/Stoner-SCF/out/LMCA001.csv";
        std::string Seed110 = "/home/cmp/Documents/Github/Stoner-SCF/out/LMCA110.csv";
        std::cout << "Seeding from: " << Seed001;
        std::cout << "Seeding from: " << Seed110;   


        const MCAResult mca = compute_MCA(S0, alpha, grid, T, N_target, delta, p, kp, 
            Seed001, Seed110);

        //001 Result — quantization axis must be z for [001]
        p.theta = 0.0;
        p.phi   = 0.0;
        const auto lmom001 = compute_L_moments(mca.res_001.rho, p);
        const auto [l001R1, l001R2] = lmom001[2];
        const auto smom001 = compute_S_moments(mca.res_001.rho, p);
        const auto [s001R1, s001R2] = smom001[2];

        //Update param for [110]
        //p.theta = M_PI / 2.0;
        //p.phi   = M_PI / 4.0;

        //p.theta = 0;
        //p.phi = 0;

        //110 Result
        const auto lmom110 = compute_L_moments(mca.res_110.rho,p);
        const auto [lx1, lx2] = lmom110[0];
        const auto [ly1, ly2] = lmom110[1];
        const double l110_1 = (lx1 + ly1) / std::sqrt(2.0);
        const double l110_2 = (lx2 + ly2) / std::sqrt(2.0);
        const auto smom110 = compute_S_moments(mca.res_110.rho, p);
        const auto [s110R1, s110R2] = smom110[2];

        //std::cout << "\n\nHello" << s110R1 + s110R2;

        //outfile << "lam,E_110, S_110, L_110, E_001, S_001, L_001, MCA\n";
        outfile << p.lam << ","<< mca.res_110.E_total << ","
                << s110R1 + s110R2 << "," << l110_1 + l110_2 << "," << mca.res_001.E_total << ","
                << s001R1 + s001R2 << "," << l001R1 + l001R2 << ","
                << mca.E_MCA << "\n";
        outfile.flush();
    //}
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

// Number of electron sweep
void run_n_electron_sweep(double alpha, int grid, double T, double N_target,
                      Params p, KanamoriParams kp, int max_iter_start){
    
    // Initial Run and Comparison between seeds
    std::vector<std::string> seeds = {"xy","yz","high_spin","low_spin"};
    std::vector<KanamoriResult> results;
    std::vector<std::string> converged_seeds;
    Mat12 loaded_rho = Mat12::Zero();
    for (const auto& seed : seeds) {
        loaded_rho = make_seed(seed, 0.01, 42); //seed, pertubation strength, random seed of pertubation
        KanamoriResult runResult;
        try {
            runResult = runKanamoriSCF(loaded_rho, alpha, grid, T, N_target, p, kp, MixerType::LinearDIIS, max_iter_start);
            converged_seeds.push_back(seed);
            results.push_back(runResult);
            std::cout << "\n\nConverged: Kanamori SCF for seed " << seed << "\n\n\n";
        } catch (const std::exception& e) {
            std::cerr << "Not converged: Kanamori SCF for seed " << seed << ": " << e.what() << std::endl;
        }
        
    }
    auto minIt = std::min_element(results.begin(), results.end(),
        [](const KanamoriResult& a, const KanamoriResult& b) {
            if (std::isnan(a.E_total)) return false;
            if (std::isnan(b.E_total)) return true;
            return a.E_total < b.E_total;
        });

    double minEnergy = minIt->E_total;
    size_t minIndex = static_cast<size_t>(minIt - results.begin());

    std::cout << "Lowest energy value: " << minEnergy << "\n";
    std::cout << "Lowest energy is seed: " << converged_seeds[minIndex] << "\n";
    
    const double degeneracyTol = 1e-6;

    std::vector<size_t> degenerateIndices;
    for (size_t i = 0; i < results.size(); ++i) {
        if (i == minIndex) continue;
        if (std::isnan(results[i].E_total)) continue;
        if (std::abs(results[i].E_total - minEnergy) < degeneracyTol) {
            // If the seeds converge to the same state then this is okay
            double norm = (results[i].rho - results[minIndex].rho).squaredNorm();
            if (norm > 1e-6) {
                degenerateIndices.push_back(i);
            }
            else{
                std::cout << "Seed " << converged_seeds[i] << " has same state as minimum energy seed.\n";
            }
        }
    }

    if (!degenerateIndices.empty()) {
        std::cout << "Warning: " << degenerateIndices.size()
                  << " other result(s) are degenerate with the minimum "
                  << "(within tolerance " << degeneracyTol << "):\n";
        for (size_t idx : degenerateIndices) {
            std::cout << "  seed " << converged_seeds[idx]
                      << ", energy = " << results[idx].E_total
                      << ", ΔE = " << (results[idx].E_total - minEnergy)
                      << "\n";
        }

        std::cout << "Degenerate seeds found. Canceling sweep.\n";
        return;
    }
    else{
        // If not degenerate, look at yz/xz degeneracy
        KanamoriResult lowestEnergy = results[minIndex];
        std::cout << lowestEnergy.rho(3,3).real() << " (yz) vs " << lowestEnergy.rho(4,4).real() << " (xz)\n";
        if (std::abs(lowestEnergy.rho(3,3).real() - lowestEnergy.rho(4,4).real()) > 1e-6) 
        {
            std::cout << "Cancelling sweep due to yz/xz degeneracy broken.\n";
            return;
        }
        else {
            std::cout << "No yz/xz degeneracy broken. Continuing sweep.\n";
            //save_projected_dos(lowestEnergy.rho, grid, T, N_target, p, kp, "out/projected_dos.csv");
            //std::cout << "Press Enter to continue... or ctr+c to cancel...";
            //std::string dummy;
            //std::getline(std::cin, dummy);
            // Increasing N count
            Mat12 run_rh = lowestEnergy.rho;
            bool isYZDeg = false;
            std::vector<double> MCA_V;
            std::vector<double> Nelec_V;
            std::vector<double> isYZDeg_V;


            std::vector<KanamoriResult> sweepResults;
            // Increasing
            for (int i = 0; i < 10; ++i) {
                double n_electron_target = N_target + (i * 0.1);
                try {
                KanamoriResult runResult = runKanamoriSCF(run_rh, alpha, grid, T, n_electron_target, p, kp, MixerType::LinearDIIS, 3000);

                // Check for YZ degeneracy
                if (std::abs(runResult.rho(3,3).real() - runResult.rho(4,4).real()) < 1e-6) {
                    isYZDeg_V.push_back(1.0);
                }
                else{
                    isYZDeg_V.push_back(0.0);
                }

                // Set seed adiabatically
                run_rh = runResult.rho;

                // 110
                Params pc = p;
                pc.theta = M_PI / 2.0;
                pc.phi   = M_PI / 4.0;
                KanamoriResult runResult110 = runKanamoriSCF(run_rh, alpha, grid, T, n_electron_target, pc, kp, MixerType::LinearDIIS, 3000);

                double MCA_D = runResult110.E_total - runResult.E_total;
                MCA_V.push_back(MCA_D);
                Nelec_V.push_back(n_electron_target);
                
                std::cout << "\n\nConverged: Kanamori SCF for n_electron = " << n_electron_target << "\n\n\n";
            } catch (const std::exception& e) {
                std::cerr << "Not converged: Kanamori SCF for n_electron = " << n_electron_target << ": " << e.what() << std::endl;
            }
            }

            // Decreasing
            run_rh = lowestEnergy.rho;
            for (int i = 0; i < 10; ++i) {
                double n_electron_target = N_target - (i * 0.1);
                try {
                KanamoriResult runResult = runKanamoriSCF(run_rh, alpha, grid, T, n_electron_target, p, kp, MixerType::LinearDIIS, 3000);

                // Check for YZ degeneracy
                if (std::abs(runResult.rho(3,3).real() - runResult.rho(4,4).real()) < 1e-6) {
                    isYZDeg_V.push_back(1.0);
                }
                else{
                    isYZDeg_V.push_back(0.0);
                }

                // Set seed adiabatically
                run_rh = runResult.rho;

                // 110
                Params pc = p;
                pc.theta = M_PI / 2.0;
                pc.phi   = M_PI / 4.0;
                KanamoriResult runResult110 = runKanamoriSCF(run_rh, alpha, grid, T, n_electron_target, pc, kp, MixerType::LinearDIIS, 3000);

                double MCA_D = runResult110.E_total - runResult.E_total;
                MCA_V.push_back(MCA_D);
                Nelec_V.push_back(n_electron_target);
                
                std::cout << "\n\nConverged: Kanamori SCF for n_electron = " << n_electron_target << "\n\n\n";
            } catch (const std::exception& e) {
                std::cerr << "Not converged: Kanamori SCF for n_electron = " << n_electron_target << ": " << e.what() << std::endl;
            }
            }


            if (!MCA_V.empty()) {
                std::cout << "MCA values,";
                for (double mca : MCA_V) {
                    std::cout << mca << ",";
                }
                std::cout << "\n";
            }
            if (!Nelec_V.empty()) {
                std::cout << "N electron values,";
                for (double nelec : Nelec_V) {
                    std::cout << nelec << ",";
                }
                std::cout << "\n";
            }

            // Write CSV
            std::string filename = "out/n_electron_sweep.csv";
            std::cout << "Saving results to " << filename << "\n";
            std::ofstream out(filename);
            if (!out.is_open()) {
                std::cerr << "Failed to open " << filename << " for writing\n";
                return;
            }

            // Header
            out << "MCA,Nelec,isYZDeg\n";

            // Assumes all vectors are the same size
            size_t n = MCA_V.size();
            out << std::setprecision(15); // adjust precision as needed

            for (size_t i = 0; i < n; ++i) {
                out << MCA_V[i] << ","
                    << Nelec_V[i] << ","
                    << isYZDeg_V[i] << "\n";
            }

            out.close();
        }
       
    }
    
}