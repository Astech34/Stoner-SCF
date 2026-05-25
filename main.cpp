#include <iostream>
#include <iomanip>
#include <filesystem>
#include <utility>
#include <fstream>
#include "hamiltonian.h"
#include "scf.h"
#include "sweep.h"

// Internal spin-major orbital order per 6x6 layer block: yz=0, xz=1, xy=2 (up), yz=3, xz=4, xy=5 (dn)
static void printKanamoriOccupations(const KanamoriResult& res) {
    const Mat12& rho = res.rho;

    struct OrbEntry { const char* name; int up; int dn; };
    // Display in dxy/dxz/dyz order, mapped to internal indices
    const OrbEntry orbs[] = {{"dxy", 2, 5}, {"dxz", 1, 4}, {"dyz", 0, 3}};

    double total_up = 0.0, total_dn = 0.0;

    for (int layer = 0; layer < 2; layer++) {
        const int base = layer * 6;
        double layer_up = 0.0, layer_dn = 0.0;

        std::cout << "Layer " << (layer + 1) << "\n";
        std::cout << " Orbital   Spin Up   Spin Dn\n";

        for (const auto& o : orbs) {
            const double up = rho(base + o.up, base + o.up).real();
            const double dn = rho(base + o.dn, base + o.dn).real();
            layer_up += up;
            layer_dn += dn;
            total_up += up;
            total_dn += dn;
            std::cout << "   " << std::left  << std::setw(6) << o.name
                      << "   " << std::right << std::fixed << std::setprecision(4)
                      << up << "    " << dn << "\n";
        }
        std::cout << "   Spin Up  " << layer_up << "\n";
        std::cout << "   Spin Dn  " << layer_dn << "\n";
        std::cout << "   Total    " << (layer_up + layer_dn) << "\n";
        std::cout << "   Moment   " << (layer_up - layer_dn) << "\n\n";
    }

    std::cout << "Total e-    " << (total_up + total_dn) << "\n";
    std::cout << "  Moment    " << (total_up - total_dn) << "\n\n";
    std::cout << "Total energy (eV):    " << res.E_total << "\n";
}

int main() {
    Params p;
    p.t1      = 0.3;
    p.t_delta = 0.05;
    p.t2      = 0.03;
    p.lam     = 0.08;   // typical upper-end 3d SOC (~50-100 meV for Co/Ni with t1~0.5 eV)
    p.U       = 4;   // well into ferromagnetic phase
    p.t_perp    = 0.2;   // interlayer hopping for yz and xz
    p.t_perp_xy = 0.05;   // interlayer hopping for xy (set equal to t_perp per advisor)
    p.delta_cf1 = 0.23;   // tetragonal crystal field layer 1: raises xy above yz/xz
    p.delta_cf2 = -0.25;   // tetragonal crystal field layer 2
    //p.delta_V   = 0.073684;   // staggered layer potential (Kanamori only)
    p.delta_V   = 0.0;   // start with zero staggered layer potential, then sweep in stage 3
    p.theta     = M_PI / 2.0;  // polar angle of spin quantization axis (0 = z-axis / out-of-plane)
    p.phi       = M_PI / 4.0;  // azimuthal angle of spin quantization axis

    const double S0       = 0.3;
    const double alpha    = 0.2;
    const double T        = 0.05;
    const double N_target = 10.0;  // 5 electrons per layer x 2 layers
    const int    grid     = 50;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=== Stoner-SCF: MCA energy ===\n";
    std::cout << "Parameters:\n";
    std::cout << "  t1      = " << p.t1      << "\n";
    std::cout << "  t_delta = " << p.t_delta << "\n";
    std::cout << "  t2      = " << p.t2      << "\n";
    std::cout << "  lam     = " << p.lam     << "\n";
    std::cout << "  U       = " << p.U       << "\n";
    std::cout << "  t_perp    = " << p.t_perp    << "\n";
    std::cout << "  t_perp_xy = " << p.t_perp_xy << "\n";
    std::cout << "  delta_cf1 = " << p.delta_cf1 << "\n";
    std::cout << "  delta_cf2 = " << p.delta_cf2 << "\n";
    std::cout << "  delta_V   = " << p.delta_V   << "\n";
    std::cout << "  theta     = " << p.theta     << " rad\n";
    std::cout << "  phi       = " << p.phi       << " rad\n";
    std::cout << "  S0      = " << S0        << "\n";
    std::cout << "  alpha   = " << alpha     << "\n";
    std::cout << "  T       = " << T         << "\n";
    std::cout << "  N       = " << N_target  << "\n";
    std::cout << "  grid    = " << grid << " x " << grid << "\n\n";

    KanamoriParams kp;
    kp.U       = p.U;
    kp.J       = 0.8;
    kp.U_prime = kp.U - 2*kp.J;  // enforce rotational invariance

    std::cout << "KanamoriParams:\n";
    std::cout << "  U       = " << kp.U       << "\n";
    std::cout << "  U'      = " << kp.U_prime << "\n";
    std::cout << "  J       = " << kp.J       << "\n\n";

    const double delta = 0.01;
    const MCAResult mca = compute_MCA(S0, alpha, grid, T, N_target, delta, p, kp);

    std::cout << "\n=== [001] Occupations ===\n";
    printKanamoriOccupations(mca.res_001);
    std::cout << "=== [110] Occupations ===\n";
    printKanamoriOccupations(mca.res_110);
    std::cout << "E_MCA = E[110] - E[001] = " << mca.E_MCA << " eV\n";

    //run_MCA_lam_sweep(S0, alpha, grid, T, N_target, 0.0, 0.2, 10, delta, p, kp);

    // --- delta_V sweep ---
    //std::cout << "\n=== Stage 3: delta_V sweep (0 -> 0.1) ===\n\n";
    //run_delta_V_sweep(S0, alpha, grid, T, N_target, 0.0, 0.1, 20, p, kp);

    return 0;
}