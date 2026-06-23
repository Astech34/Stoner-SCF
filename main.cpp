#include <iostream>
#include <iomanip>
#include <filesystem>
#include <utility>
#include <fstream>
#include "hamiltonian.h"
#include "scf.h"
#include "sweep.h"
#include "params_io.h"

int main(int argc, char* argv[]) {
    const std::string param_file = (argc > 1) ? argv[1] : "params.in";
    const AllParams ap = load_params(param_file);
    const Params&         p   = ap.p;
    const KanamoriParams& kp  = ap.kp;
    const SCFParams&      scf = ap.scf;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "=== Stoner-SCF: MCA energy ===\n";
    std::cout << "Loaded parameters from '" << param_file << "':\n";
    std::cout << "  tpi       = " << p.tpi       << "\n";
    std::cout << "  tdelta    = " << p.tdelta    << "\n";
    std::cout << "  t2xy      = " << p.t2xy      << "\n";
    std::cout << "  t2yzxz    = " << p.t2yzxz    << "\n";
    std::cout << "  tg        = " << p.tg        << "\n";
    std::cout << "  lam       = " << p.lam       << "\n";
    std::cout << "  U         = " << p.U         << "\n";
    std::cout << "  t_perp    = " << p.t_perp    << "\n";
    std::cout << "  t_perp_xy = " << p.t_perp_xy << "\n";
    std::cout << "  delta_cf1 = " << p.delta_cf1 << "\n";
    std::cout << "  delta_cf2 = " << p.delta_cf2 << "\n";
    std::cout << "  delta_V   = " << p.delta_V   << "\n";
    std::cout << "  theta     = " << p.theta     << " rad\n";
    std::cout << "  phi       = " << p.phi       << " rad\n";
    std::cout << "  S0        = " << scf.S0      << "\n";
    std::cout << "  alpha     = " << scf.alpha   << "\n";
    std::cout << "  T         = " << scf.T       << "\n";
    std::cout << "  N_target  = " << scf.N_target << "\n";
    std::cout << "  grid      = " << scf.grid    << " x " << scf.grid << "\n\n";

    std::cout << "KanamoriParams:\n";
    std::cout << "  U       = " << kp.U       << "\n";
    std::cout << "  U'      = " << kp.U_prime << "\n";
    std::cout << "  J       = " << kp.J       << "\n\n";

    
    // const double delta = 0.01;
    std::cout << "=== Stoner bootstrap ===\n";
    const CalcResult stoner = runSelfCalc(scf.S0, scf.alpha, scf.grid, scf.T, scf.N_target, p);
    const Eigensystem sys0  = compute_eigensystem_grid(stoner.S_new, scf.grid, p);
    Mat12 rho0 = compute_density_matrix(sys0, stoner.mu, scf.T);
    apply_symmetry_breaking(rho0, 0.1);

    // Kanamori at [001]
    std::cout << "\n=== Kanamori SCF: [001] ===\n";
    const KanamoriResult res_001 = runKanamoriSCF(rho0, scf.alpha, scf.grid, scf.T, scf.N_target, p, kp, MixerType::LinearDIIS);
    std::cout << "\n=== [001] Occupations ===\n";
    printKanamoriOccupations(res_001);
    
    


    //const double delta = 0.02;
    //const MCAResult mca = compute_MCA(scf.S0, scf.alpha, scf.grid, scf.T, scf.N_target, delta, p, kp);

    //std::cout << "E_MCA = E[110] - E[001] = " << mca.E_MCA << " eV\n";

    //const double delta = 0.01;
    //run_MCA_lam_sweep(scf.S0, scf.alpha, scf.grid, scf.T, scf.N_target, 0.01, 0.05, 9, delta, p, kp);

    // --- delta_V sweep ---
    //std::cout << "\n=== Stage 3: delta_V sweep (0 -> 0.1) ===\n\n";
    //run_delta_V_sweep(scf.S0, scf.alpha, scf.grid, scf.T, scf.N_target, 0.0, 0.1, 20, p, kp);

    //const KanamoriResult res = runKanamoriSCF_random(4, scf.S0, scf.alpha, scf.grid, scf.T, scf.N_target, p, kp, 0.001);
    //std::cout << "\n=== Random-seed result ===\n";
    //printKanamoriOccupations(res);

    save_band_structure(res_001.rho, 200, p, kp, "out/band_structure.csv", res_001.mu);
    save_projected_dos(res_001.rho, scf.grid, scf.T, scf.N_target, p, kp, "out/projected_dos.csv");
    save_density_matrix(res_001.rho, "out/density_matrix.csv");

    return 0;
}