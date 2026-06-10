#include "hamiltonian.h"
#include "scf.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <iomanip>
#include <string>
#include <vector>
#include <array>

Mat6 H0(double kx, double ky, const Params& p, double delta_cf) {

    const double e_xy = -2*p.tpi*(std::cos(kx) + std::cos(ky)) - 4*p.t2xy * std::cos(kx)*std::cos(ky);
    const double e_yz = -2*p.tpi*std::cos(kx) - 2*p.tdelta*std::cos(ky) - 4*p.t2yzxz * std::cos(kx)*std::cos(ky);
    const double e_xz = -2*p.tdelta*std::cos(kx) - 2*p.tpi*std::cos(ky) - 4*p.t2yzxz * std::cos(kx)*std::cos(ky);

    const double g = 4*p.tg * std::sin(kx) * std::sin(ky);
    
    Mat6 H = Mat6::Zero();

    // Spin-major ordering: spin-up (0,1,2) = yz,xz,xy | spin-down (3,4,5) = yz,xz,xy
    H(0,0) = e_yz - 0.5*delta_cf;               H(1,1) = e_xz - 0.5*delta_cf;               H(2,2) = e_xy + delta_cf;
    H(3,3) = e_yz - 0.5*delta_cf;               H(4,4) = e_xz - 0.5*delta_cf;               H(5,5) = e_xy + delta_cf;

    return H;
}
// Checked
Mat6 SOC(double lam, double theta, double phi) {
    // Orbital angular momentum matrices in t2g basis (yz, xz, xy)
    Eigen::Matrix<cd, 3, 3> Lx, Ly, Lz;

    Lx << cd(0,0),  cd(0,0),  cd(0,0),
          cd(0,0),  cd(0,0),  cd(0,-1),
          cd(0,0),  cd(0,1), cd(0,0);

    Ly << cd(0,0),  cd(0,0),  cd(0,1),
          cd(0,0),  cd(0,0),  cd(0,0),
          cd(0,-1),  cd(0,0),  cd(0,0);

    Lz << cd(0,0),  cd(0,-1),  cd(0,0),
          cd(0,1), cd(0,0),  cd(0,0),
          cd(0,0),  cd(0,0),  cd(0,0);

    // Spin-1/2 matrices in z-quantization basis
    Eigen::Matrix<cd, 2, 2> sx, sy, sz;

    sx << cd(0,0),   cd(0.5,0),
          cd(0.5,0), cd(0,0);

    sy << cd(0,0),   cd(0,-0.5),
          cd(0,0.5), cd(0,0);

    sz << cd(0.5,0), cd(0,0),
          cd(0,0),   cd(-0.5,0);

    // Rotate spin quantization axis to n̂(θ,φ): s' = U†·s·U
    if (theta != 0.0 || phi != 0.0) {
        const double c  = std::cos(theta / 2.0);
        const double s  = std::sin(theta / 2.0);
        const cd     ep = std::exp(cd(0,  phi/2));
        const cd     em = std::exp(cd(0, -phi/2));

        Eigen::Matrix<cd, 2, 2> U;
        U << em * cd(c, 0),  -s * em,
              s * ep,   ep * cd(c, 0);

        sx = U.adjoint() * sx * U;
        sy = U.adjoint() * sy * U;
        sz = U.adjoint() * sz * U;
    }

    // Kronecker products: kron(s', L) gives spin-major 6x6
    Mat6 Hsoc = Mat6::Zero();

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int si = 0; si < 2; si++)
                for (int sj = 0; sj < 2; sj++)
                    Hsoc(3*si+i, 3*sj+j) -= lam * (Lx(i,j) * sx(si,sj)
                                                   + Ly(i,j) * sy(si,sj)
                                                   + Lz(i,j) * sz(si,sj));

    return Hsoc;
}
// Checked
Mat6 HubbardU(double S, const Params& p) {
    // Exchange field always along z: H = -U*S * σ_z ⊗ I₃
    // Angles enter only through the SOC rotation, not the exchange field.
    const double scale = -p.U * S;

    Mat6 H = Mat6::Zero();
    H(0,0) = scale;   H(1,1) = scale;   H(2,2) = scale;
    H(3,3) = -scale;  H(4,4) = -scale;  H(5,5) = -scale;

    return H;
}

Mat6 singleLayer(double kx, double ky, double S, const Params& p) {
    return H0(kx, ky, p, p.delta_cf1) + SOC(p.lam, p.theta, p.phi) + HubbardU(S, p);
}

// Checked
Mat6 T_perp_mat(const Params& p) {
    Mat6 T = Mat6::Zero();
    // yz (0) and xz (1): out-of-plane lobes, standard interlayer hopping
    // xy  (2): in-plane lobes, zero for ideal stacking; nonzero if t_perp_xy is set
    // Applied to both spin-up (0,1,2) and spin-down (3,4,5) blocks
    T(0,0) = p.t_perp;     T(1,1) = p.t_perp;     T(2,2) = p.t_perp_xy;
    T(3,3) = p.t_perp;     T(4,4) = p.t_perp;     T(5,5) = p.t_perp_xy;
    return T;
}

Mat12 staggered_potential(const Params& p) {
    Mat12 V = Mat12::Zero();
    for (int m = 0; m < 6; m++)
        V(m, m) = -p.delta_V;   // layer 1: lower energy
    for (int m = 6; m < 12; m++)
        V(m, m) = +p.delta_V;   // layer 2: higher energy
    return V;
}

Mat12 bilayerHamiltonian(double kx, double ky, double S, const Params& p) {
    // Layer-major ordering:
    //   rows/cols 0-5:  layer 1 (spin-up yz,xz,xy | spin-down yz,xz,xy)
    //   rows/cols 6-11: layer 2 (same ordering)
    const Mat6 Hsoc = SOC(p.lam, p.theta, p.phi);
    const Mat6 Hhub = HubbardU(S, p);
    const Mat6 T    = T_perp_mat(p);

    Mat12 H = Mat12::Zero();
    H.block<6,6>(0,0) = H0(kx, ky, p, p.delta_cf1) + Hsoc + Hhub;
    H.block<6,6>(6,6) = H0(kx, ky, p, p.delta_cf2) + Hsoc + Hhub;
    H.block<6,6>(0,6) = T;
    H.block<6,6>(6,0) = T;  // T is real and diagonal, so T† = T
    return H;
}
// -----------------------------------------------------------------------------
// calculate_hubbard_energy — band energy + Hubbard DC correction for a given S
// (convenience wrapper used for convergence checks or E-vs-S plots)
double calculate_hubbard_energy(double S, int grid_size, double T, double N_target,
                                const Params& p) {
    const Eigensystem sys = compute_eigensystem_grid(S, grid_size, p);
    const double mu = find_mu(sys, T, N_target);

    if (std::isnan(mu)) {
        std::cout << "calculate_hubbard_energy: could not find mu, returning NaN.\n";
        return std::numeric_limits<double>::quiet_NaN();
    }

    return calculate_band_energy(sys, mu, T) + 6.0 * p.U * S * S;
}

// -----------------------------------------------------------------------------
// KanamoriMF — Kanamori mean-field Hamiltonian (12x12 bilayer)
// rho(a,b) = <c†_a c_b>, layer-major / spin-major / orbital-minor ordering
// -----------------------------------------------------------------------------
namespace {
// Compute the 6x6 single-layer Kanamori MF block from the 6x6 density matrix.
// Spin-major ordering: up = indices 0-2 (yz, xz, xy), down = indices 3-5.
//
// Three density-density terms (diagonal in orbital-spin):
//   U    intraorbital:          <n_{m↑}> n̂_{m↓}  +  <n_{m↓}> n̂_{m↑}
//   U'   interorbital opp-spin: <n_{m↑}> n̂_{m'↓} +  <n_{m'↓}> n̂_{m↑}   (m≠m')
//   U'-J interorbital same-spin: <n_{mσ}> n̂_{m'σ} + <n_{m'σ}> n̂_{mσ}   (m<m')
//
// Two off-diagonal terms (exchange H_exc and pair-hopping H_ph), each summed
// over all ordered pairs (m, m') with m≠m', coefficient J.
//
// This is well checked!
Mat6 kanamori_layer(const Mat6& rho, const KanamoriParams& kp) {
    const double U  = kp.U;
    const double Up = kp.U_prime;
    const double J  = kp.J;

    Mat6 H = Mat6::Zero();

    double n_up = 0.0, n_dn = 0.0;
    for (int m = 0; m < 3; m++) {
        n_up += rho(m,   m  ).real();
        n_dn += rho(m+3, m+3).real();
    }

    // Diagonal terms
    for (int m = 0; m < 3; m++) {
        const double n_up_m = rho(m,   m  ).real();
        const double n_dn_m = rho(m+3, m+3).real();

        // Up Diagonal
        H(m,   m  ) += U  * n_dn_m
                     + Up * (n_dn - n_dn_m)
                     + (Up - J) * (n_up - n_up_m);
        
        // Down Diagonal
        H(m+3, m+3) += U  * n_up_m
                     + Up * (n_up - n_up_m)
                     + (Up - J) * (n_dn - n_dn_m);
        
        // Hubbard U Off Diagonal Terms
        H(m+3, m) -= U * rho(m, m+3); // -U <d†_{m↑} d_{m↓}> d†_{m↓} d_{m↑}
        H(m, m+3) -= U * rho(m+3, m); // -U <d†_{m↓} d_{m↑}> d†_{m↑} d_{m↓}
    }

    // Off-diagonal terms: exchange (H_exc) and pair-hopping (H_ph)
    for (int m = 0; m < 3; m++) {
        for (int mp = 0; mp < 3; mp++) {
            if (mp == m) continue;

            // U' off diagonal terms
            H(mp+3, m) -= Up * rho(m, mp+3);
            H(m, mp+3) -= Up * rho(mp+3, m);

            // U'-J off diagonal
            // Up Spin
            H(mp, m) -= 0.5 * (Up - J) * rho(m, mp);
            H(m, mp) -= 0.5 * (Up - J) * rho(mp, m);
            // Down Spin
            H(mp+3, m+3) -= 0.5 * (Up - J) * rho(m+3, mp+3);
            H(m+3, mp+3) -= 0.5 * (Up - J) * rho(mp+3, m+3);

            // H_exc 1: <d†_{m↑} d_{m'↑}> d†_{m'↓} d_{m↓}
            H(mp+3, m+3) += J * rho(m, mp);
            // H_ph  1: <d†_{m↑} d_{m'↑}> d†_{m↓} d_{m'↓}
            H(m+3, mp+3) += J * rho(m, mp);

            // H_exc 2 + H_ph 2: (<d†_{m'↓}d_{m↓}> + <d†_{m↓}d_{m'↓}>) d†_{m↑}d_{m'↑}
            H(m, mp) += J * (rho(mp+3, m+3) + rho(m+3, mp+3));

            // H_exc 3: -<d†_{m↑} d_{m↓}> d†_{m'↓} d_{m'↑}
            H(mp+3, mp) -= J * rho(m, m+3);
            // H_ph  3: -<d†_{m↑} d_{m'↓}> d†_{m↓} d_{m'↑}
            H(m+3, mp) -= J * rho(m, mp+3);

            // H_exc 4: -<d†_{m'↓} d_{m'↑}> d†_{m↑} d_{m↓}
            H(m, m+3) -= J * rho(mp+3, mp);
            // H_ph  4: -<d†_{m↓} d_{m'↑}> d†_{m↑} d_{m'↓}
            H(m, mp+3) -= J * rho(m+3, mp);
        }
    }

    return H;
}
} // anonymous namespace

Mat12 KanamoriMF(const Mat12& rho, const KanamoriParams& kp) {
    Mat12 H = Mat12::Zero();
    H.block<6,6>(0, 0) = kanamori_layer(rho.block<6,6>(0, 0), kp);
    H.block<6,6>(6, 6) = kanamori_layer(rho.block<6,6>(6, 6), kp);
    return H;
}

// Plotting Band Structure
// -----------------------------------------------------------------------------
// save_band_structure — diagonalise H along Γ→X→M→Γ and write to CSV
// Columns: path_index, kx, ky, band_0, band_1, ..., band_5
// -----------------------------------------------------------------------------
void save_band_structure(const Mat12& rho, int n_points, const Params& p, const KanamoriParams& kp,
                         const std::string& filename, double mu) {
    // High-symmetry points for square lattice
    // Γ = (0,0)  X = (π,0)  M = (π,π)
    struct KPoint { double kx, ky; const char* label; };
    const std::vector<KPoint> corners = {
        {0.0,   0.0,  "G"},
        {M_PI,  0.0,  "X"},
        {M_PI,  M_PI, "M"},
        {0.0,   0.0,  "G"}
    };

    struct PathPoint { double kx, ky; double path_coord; std::string label; };
    std::vector<PathPoint> path;
    path.reserve(n_points);

    const int n_seg       = static_cast<int>(corners.size()) - 1;
    const int pts_per_seg = n_points / n_seg;
    double path_coord = 0.0;

    for (int seg = 0; seg < n_seg; seg++) {
        const double kx0 = corners[seg].kx,   ky0 = corners[seg].ky;
        const double kx1 = corners[seg+1].kx, ky1 = corners[seg+1].ky;
        const double seg_len = std::sqrt((kx1-kx0)*(kx1-kx0) + (ky1-ky0)*(ky1-ky0));
        const double d_path  = seg_len / pts_per_seg;

        for (int i = 0; i < pts_per_seg; i++) {
            const double t  = static_cast<double>(i) / pts_per_seg;
            const double kx = kx0 + t * (kx1 - kx0);
            const double ky = ky0 + t * (ky1 - ky0);
            std::string lbl = (i == 0) ? corners[seg].label : "";
            path.push_back({kx, ky, path_coord, lbl});
            path_coord += d_path;
        }
    }
    // Add the final Γ point
    path.push_back({corners.back().kx, corners.back().ky, path_coord, corners.back().label});

    std::ofstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("save_band_structure: could not open " + filename);

    // Header — mu column added
    file << "path_coord,kx,ky,label,mu";
    for (int b = 0; b < 12; b++)
        file << ",band_" << b;
    file << "\n";

    file << std::fixed << std::setprecision(8);

    const Mat6 Hsoc  = SOC(p.lam, p.theta, p.phi);
    const Mat6 Tperp = T_perp_mat(p);

    Mat12 H_kfree = KanamoriMF(rho, kp);
    H_kfree.block<6,6>(0, 0) += Hsoc;
    H_kfree.block<6,6>(6, 6) += Hsoc;
    H_kfree.block<6,6>(0, 6) += Tperp;
    H_kfree.block<6,6>(6, 0) += Tperp;
    H_kfree += staggered_potential(p);

    for (const auto& pt : path) {
        Mat12 H = H_kfree;
        H.block<6,6>(0, 0) += H0(pt.kx, pt.ky, p, p.delta_cf1);
        H.block<6,6>(6, 6) += H0(pt.kx, pt.ky, p, p.delta_cf2);
        Eigen::SelfAdjointEigenSolver<Mat12> solver(H);
        const auto& evals = solver.eigenvalues();

        file << pt.path_coord << "," << pt.kx << "," << pt.ky << ","
             << pt.label << "," << mu;
        for (int b = 0; b < 12; b++)
            file << "," << evals[b];
        file << "\n";
    }

    std::cout << "Band structure written to " << filename
              << " (" << path.size() << " k-points)\n";
}
// -----------------------------------------------------------------------------
// save_dos — compute DOS on a fine energy grid and write to CSV
void save_dos(double S, int grid_size, double T, double N_target,
              const Params& p, const std::string& filename,
              int n_energy_points, double sigma)  {

    // --- Compute eigensystem on k-grid ---
    const Eigensystem sys = compute_eigensystem_grid(S, grid_size, p);
    const double mu = find_mu(sys, T, N_target);

    // --- Energy axis: span all eigenvalues with some padding ---
    const int N_k     = static_cast<int>(sys.evals.size());
    const int N_bands = 12;

    double e_min =  std::numeric_limits<double>::infinity();
    double e_max = -std::numeric_limits<double>::infinity();
    for (const auto& ev : sys.evals) {
        e_min = std::min(e_min, ev.minCoeff());
        e_max = std::max(e_max, ev.maxCoeff());
    }
    const double padding = 5.0 * sigma;
    const double e_start = e_min - padding;
    const double e_end   = e_max + padding;
    const double de      = (e_end - e_start) / (n_energy_points - 1);

    // --- Gaussian broadening ---
    // DOS(E) = (1 / N_k) * sum_{k,n} (1 / (sigma * sqrt(2pi))) * exp(-0.5 * ((E - E_nk) / sigma)^2)
    const double norm = 1.0 / (sigma * std::sqrt(2.0 * M_PI));
    std::vector<double> energy(n_energy_points);
    std::vector<double> dos(n_energy_points, 0.0);

    for (int ie = 0; ie < n_energy_points; ie++) {
        energy[ie] = e_start + ie * de;
        for (int idx = 0; idx < N_k; idx++) {
            for (int b = 0; b < N_bands; b++) {
                const double x = (energy[ie] - sys.evals[idx][b]) / sigma;
                dos[ie] += norm * std::exp(-0.5 * x * x);
            }
        }
        dos[ie] /= static_cast<double>(N_k);
    }

    // --- Write to CSV ---
    std::ofstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("save_dos: could not open " + filename);

    file << std::fixed << std::setprecision(8);
    file << "energy,dos,mu\n";
    for (int ie = 0; ie < n_energy_points; ie++)
        file << energy[ie] << "," << dos[ie] << "," << mu << "\n";

    std::cout << "DOS written to " << filename << "\n";
}
// -----------------------------------------------------------------------------
// save_projected_dos — layer/spin/orbital-resolved DOS for the Kanamori system
// -----------------------------------------------------------------------------
// Builds the same Hamiltonian as the converged Kanamori SCF (KanamoriMF(rho) +
// SOC + interlayer hopping + staggered potential + kinetic H0) and resolves the
// DOS onto each of the 12 basis components. Eigenstate |v_nk> contributes its
// projection weight |v_nk(c)|^2 to component c, Gaussian-broadened in energy:
//   PDOS_c(E) = (1/N_k) sum_{k,n} |v_nk(c)|^2 / (sigma*sqrt(2pi)) * exp(-0.5*((E-E_nk)/sigma)^2)
// Since the eigenvectors are normalised, the 12 projections sum to the total DOS.
// Basis ordering (layer-major / spin-major / orbital-minor):
//   0-2  L1 up (yz,xz,xy)   3-5  L1 dn (yz,xz,xy)
//   6-8  L2 up (yz,xz,xy)   9-11 L2 dn (yz,xz,xy)
void save_projected_dos(const Mat12& rho, int grid_size, double T, double N_target,
                        const Params& p, const KanamoriParams& kp,
                        const std::string& filename,
                        int n_energy_points, double sigma) {

    // --- Compute eigensystem with the converged Kanamori MF Hamiltonian ---
    const Eigensystem sys = compute_eigensystem_kanamori(rho, grid_size, p, kp);
    const double mu = find_mu(sys, T, N_target);

    // --- Energy axis: span all eigenvalues with some padding ---
    const int N_k     = static_cast<int>(sys.evals.size());
    const int N_bands = 12;
    const int N_proj  = 12;

    double e_min =  std::numeric_limits<double>::infinity();
    double e_max = -std::numeric_limits<double>::infinity();
    for (const auto& ev : sys.evals) {
        e_min = std::min(e_min, ev.minCoeff());
        e_max = std::max(e_max, ev.maxCoeff());
    }
    const double padding = 5.0 * sigma;
    const double e_start = e_min - padding;
    const double e_end   = e_max + padding;
    const double de      = (e_end - e_start) / (n_energy_points - 1);

    // --- Gaussian broadening, accumulated per basis component ---
    const double norm = 1.0 / (sigma * std::sqrt(2.0 * M_PI));
    std::vector<double> energy(n_energy_points);
    std::vector<std::array<double, 12>> pdos(n_energy_points);
    for (auto& row : pdos) row.fill(0.0);

    #pragma omp parallel for schedule(static)
    for (int ie = 0; ie < n_energy_points; ie++) {
        const double E = e_start + ie * de;
        energy[ie] = E;

        std::array<double, 12> acc{};
        acc.fill(0.0);
        for (int idx = 0; idx < N_k; idx++) {
            for (int b = 0; b < N_bands; b++) {
                const double x = (E - sys.evals[idx][b]) / sigma;
                const double g = norm * std::exp(-0.5 * x * x);
                const auto   v = sys.evecs[idx].col(b);
                for (int c = 0; c < N_proj; c++)
                    acc[c] += g * std::norm(v[c]);   // std::norm(z) = |z|^2
            }
        }
        for (int c = 0; c < N_proj; c++)
            pdos[ie][c] = acc[c] / static_cast<double>(N_k);
    }

    // --- Write to CSV ---
    std::ofstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("save_projected_dos: could not open " + filename);

    file << std::fixed << std::setprecision(8);
    file << "energy,mu,"
            "L1_up_yz,L1_up_xz,L1_up_xy,"
            "L1_dn_yz,L1_dn_xz,L1_dn_xy,"
            "L2_up_yz,L2_up_xz,L2_up_xy,"
            "L2_dn_yz,L2_dn_xz,L2_dn_xy,total\n";

    for (int ie = 0; ie < n_energy_points; ie++) {
        file << energy[ie] << "," << mu;
        double total = 0.0;
        for (int c = 0; c < 12; c++) {
            file << "," << pdos[ie][c];
            total += pdos[ie][c];
        }
        file << "," << total << "\n";
    }

    std::cout << "Projected DOS written to " << filename << "\n";
}