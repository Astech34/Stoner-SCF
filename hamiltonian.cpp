#include "hamiltonian.h"
#include "scf.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <iomanip>
#include <string>
#include <vector>

Mat6 H0(double kx, double ky, const Params& p) {
    const double e_xy = -2*p.t1*(std::cos(kx) + std::cos(ky)) - 4*p.t2*std::cos(kx)*std::cos(ky);
    const double e_xz = -2*p.t1*std::cos(kx) - 2*p.t_delta*std::cos(ky);
    const double e_yz = -2*p.t_delta*std::cos(kx) - 2*p.t1*std::cos(ky);

    Mat6 H = Mat6::Zero();

    // Spin-major ordering: spin-up (0,1,2) = yz,xz,xy | spin-down (3,4,5) = yz,xz,xy
    H(0,0) = e_yz;  H(1,1) = e_xz;  H(2,2) = e_xy;
    H(3,3) = e_yz;  H(4,4) = e_xz;  H(5,5) = e_xy;

    return H;
}
// Checked
Mat6 SOC(double lam) {
    // Orbital angular momentum matrices in t2g basis (yz, xz, xy)
    Eigen::Matrix<cd, 3, 3> Lx, Ly, Lz;

    Lx << cd(0,0),  cd(0,0),  cd(0,0),
          cd(0,0),  cd(0,0),  cd(0,1),
          cd(0,0),  cd(0,-1), cd(0,0);

    Ly << cd(0,0),  cd(0,0),  cd(0,-1),
          cd(0,0),  cd(0,0),  cd(0,0),
          cd(0,1),  cd(0,0),  cd(0,0);

    Lz << cd(0,0),  cd(0,1),  cd(0,0),
          cd(0,-1), cd(0,0),  cd(0,0),
          cd(0,0),  cd(0,0),  cd(0,0);

    // Spin-1/2 matrices
    Eigen::Matrix<cd, 2, 2> sx, sy, sz;

    sx << cd(0,0), cd(0.5,0),
          cd(0.5,0), cd(0,0);

    sy << cd(0,0),    cd(0,-0.5),
          cd(0,0.5),  cd(0,0);

    sz << cd(0.5,0),  cd(0,0),
          cd(0,0),    cd(-0.5,0);

    // Kronecker products: kron(L, s) gives spin-major 6x6
    Mat6 Hsoc = Mat6::Zero();

    // kron(Lx, sx)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int si = 0; si < 2; si++)
                for (int sj = 0; sj < 2; sj++)
                    Hsoc(3*si+i, 3*sj+j) += lam * Lx(i,j) * sx(si,sj);

    // kron(Ly, sy)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int si = 0; si < 2; si++)
                for (int sj = 0; sj < 2; sj++)
                    Hsoc(3*si+i, 3*sj+j) += lam * Ly(i,j) * sy(si,sj);

    // kron(Lz, sz)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int si = 0; si < 2; si++)
                for (int sj = 0; sj < 2; sj++)
                    Hsoc(3*si+i, 3*sj+j) += lam * Lz(i,j) * sz(si,sj);

    return Hsoc;
}
// Checked
Mat6 HubbardU(double S, const Params& p) {
    // Magnetisation direction n̂ = (sin θ cos φ, sin θ sin φ, cos θ)
    const double nx = std::sin(p.theta) * std::cos(p.phi);
    const double ny = std::sin(p.theta) * std::sin(p.phi);
    const double nz = std::cos(p.theta);

    // H = -U*S * (n̂·σ) ⊗ I₃
    // (n̂·σ) = [ nz·I₃,        (nx-i·ny)·I₃ ]
    //          [ (nx+i·ny)·I₃,  -nz·I₃       ]
    const double scale = -p.U * S;

    Mat6 H = Mat6::Zero();

    // Diagonal spin blocks
    H(0,0) = scale * nz;   H(1,1) = scale * nz;   H(2,2) = scale * nz;
    H(3,3) = -scale * nz;  H(4,4) = -scale * nz;  H(5,5) = -scale * nz;

    // Off-diagonal spin blocks (only non-zero when n̂ has x or y component)
    const cd off_up   = scale * cd(nx, -ny);  // spin-up row, spin-down col
    const cd off_down = scale * cd(nx,  ny);  // spin-down row, spin-up col
    H(0,3) = off_up;   H(1,4) = off_up;   H(2,5) = off_up;
    H(3,0) = off_down; H(4,1) = off_down; H(5,2) = off_down;

    return H;
}

Mat6 singleLayer(double kx, double ky, double S, const Params& p) {
    return H0(kx, ky, p) + SOC(p.lam) + HubbardU(S, p);
}

// Checked
Mat6 T_perp_mat(const Params& p) {
    Mat6 T = Mat6::Zero();
    // yz (0) and xz (1) hop between layers; xy (2) does not (lobes lie in-plane)
    // Applied to both spin-up (0,1) and spin-down (3,4) blocks
    T(0,0) = p.t_perp;  T(1,1) = p.t_perp;
    T(3,3) = p.t_perp;  T(4,4) = p.t_perp;
    return T;
}

Mat12 bilayerHamiltonian(double kx, double ky, double S, const Params& p) {
    // Layer-major ordering:
    //   rows/cols 0-5:  layer 1 (spin-up yz,xz,xy | spin-down yz,xz,xy)
    //   rows/cols 6-11: layer 2 (same ordering)
    const Mat6 H_single = H0(kx, ky, p) + SOC(p.lam) + HubbardU(S, p);
    const Mat6 T        = T_perp_mat(p);

    Mat12 H = Mat12::Zero();
    H.block<6,6>(0,0) = H_single;
    H.block<6,6>(6,6) = H_single;
    H.block<6,6>(0,6) = T;
    H.block<6,6>(6,0) = T;  // T is real and diagonal, so T† = T
    return H;
}

// Plotting Band Structure
// -----------------------------------------------------------------------------
// save_band_structure — diagonalise H along Γ→X→M→Γ and write to CSV
// Columns: path_index, kx, ky, band_0, band_1, ..., band_5
// -----------------------------------------------------------------------------
void save_band_structure(double S, int n_points, const Params& p,
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

    const Mat6 Hsoc = SOC(p.lam);
    const Mat6 Hhub = HubbardU(S, p);
    const Mat6 T    = T_perp_mat(p);

    std::ofstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("save_band_structure: could not open " + filename);

    // Header — mu column added
    file << "path_coord,kx,ky,label,mu";
    for (int b = 0; b < 12; b++)
        file << ",band_" << b;
    file << "\n";

    file << std::fixed << std::setprecision(8);

    for (const auto& pt : path) {
        const Mat6  H_single = H0(pt.kx, pt.ky, p) + Hsoc + Hhub;
        Mat12 H = Mat12::Zero();
        H.block<6,6>(0,0) = H_single;
        H.block<6,6>(6,6) = H_single;
        H.block<6,6>(0,6) = T;
        H.block<6,6>(6,0) = T;
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
    const double mu = find_mu(sys, grid_size, T, N_target);

    // --- Energy axis: span all eigenvalues with some padding ---
    const int N_k     = grid_size * grid_size;
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
// calculate total energy at given S (for convergence checks or plotting E vs S)
double calculate_total_energy(double S, int grid_size, double T, double N_target,
                               const Params& p) {
    const int N_k     = grid_size * grid_size;
    const int N_bands = 12;

    const Eigensystem sys = compute_eigensystem_grid(S, grid_size, p);
    const double mu = find_mu(sys, grid_size, T, N_target);

    if (std::isnan(mu)) {
        std::cout << "calculate_total_energy: could not find mu, returning NaN.\n";
        return std::numeric_limits<double>::quiet_NaN();
    }

    double E_total = 0.0;
    for (const auto& ev : sys.evals) {
        for (int b = 0; b < N_bands; b++) {
            const double x = std::clamp((ev[b] - mu) / T, -500.0, 500.0);
            const double f = 1.0 / (std::exp(x) + 1.0);
            E_total += ev[b] * f;
        }
    }

    return E_total / static_cast<double>(N_k);
}