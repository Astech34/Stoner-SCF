#include "hamiltonian.h"
#include <cmath>
#include <fstream>
#include <iostream>
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
                    Hsoc(2*i+si, 2*j+sj) += lam * Lx(i,j) * sx(si,sj);

    // kron(Ly, sy)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int si = 0; si < 2; si++)
                for (int sj = 0; sj < 2; sj++)
                    Hsoc(2*i+si, 2*j+sj) += lam * Ly(i,j) * sy(si,sj);

    // kron(Lz, sz)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int si = 0; si < 2; si++)
                for (int sj = 0; sj < 2; sj++)
                    Hsoc(2*i+si, 2*j+sj) += lam * Lz(i,j) * sz(si,sj);

    return Hsoc;
}

Mat6 HubbardU(double S, double U) {
    const double shift = U * S;

    Mat6 H = Mat6::Zero();

    // spin-up: negative shift | spin-down: positive shift
    H(0,0) = -shift;  H(1,1) = -shift;  H(2,2) = -shift;
    H(3,3) =  shift;  H(4,4) =  shift;  H(5,5) =  shift;

    return H;
}

Mat6 singleLayer(double kx, double ky, double S, const Params& p) {
    // SOC is k-independent — caller may want to cache this for performance
    return H0(kx, ky, p) + SOC(p.lam) + HubbardU(S, p.U);
}

// Plotting Band Structure
// -----------------------------------------------------------------------------
// save_band_structure — diagonalise H along Γ→X→M→Γ and write to CSV
// Columns: path_index, kx, ky, band_0, band_1, ..., band_5
// -----------------------------------------------------------------------------
void save_band_structure(double S, int n_points, const Params& p,
                         const std::string& filename) {
    // High-symmetry points for square lattice
    // Γ = (0,0)  X = (π,0)  M = (π,π)
    struct KPoint { double kx, ky; const char* label; };
    const std::vector<KPoint> corners = {
        {0.0,   0.0,  "G"},
        {M_PI,  0.0,  "X"},
        {M_PI,  M_PI, "M"},
        {0.0,   0.0,  "G"}
    };

    // Build the full path by linearly interpolating between corners
    struct PathPoint { double kx, ky; double path_coord; std::string label; };
    std::vector<PathPoint> path;
    path.reserve(n_points);

    // Total number of segments
    const int n_seg    = static_cast<int>(corners.size()) - 1;
    const int pts_per_seg = n_points / n_seg;

    double path_coord = 0.0;

    for (int seg = 0; seg < n_seg; seg++) {
        const double kx0 = corners[seg].kx,   ky0 = corners[seg].ky;
        const double kx1 = corners[seg+1].kx, ky1 = corners[seg+1].ky;

        // Segment length in k-space (for a physically meaningful x-axis)
        const double seg_len = std::sqrt((kx1-kx0)*(kx1-kx0) + (ky1-ky0)*(ky1-ky0));
        const double d_path  = seg_len / pts_per_seg;

        // Exclude endpoint — it becomes the start of the next segment
        for (int i = 0; i < pts_per_seg; i++) {
            const double t  = static_cast<double>(i) / pts_per_seg;
            const double kx = kx0 + t * (kx1 - kx0);
            const double ky = ky0 + t * (ky1 - ky0);

            // Label only the corner points
            std::string lbl = (i == 0) ? corners[seg].label : "";
            path.push_back({kx, ky, path_coord, lbl});
            path_coord += d_path;
        }
    }

    // Add the final Γ point
    path.push_back({corners.back().kx, corners.back().ky, path_coord, corners.back().label});

    // Precompute SOC once
    const Mat6 Hsoc = SOC(p.lam);

    // Open CSV
    std::ofstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("save_band_structure: could not open " + filename);

    // Header
    file << "path_coord,kx,ky,label";
    for (int b = 0; b < 6; b++)
        file << ",band_" << b;
    file << "\n";
    file << std::fixed << std::setprecision(8);

    // Diagonalise and write
    for (const auto& pt : path) {
        const Mat6 H = H0(pt.kx, pt.ky, p) + Hsoc + HubbardU(S, p.U);

        Eigen::SelfAdjointEigenSolver<Mat6> solver(H);
        const auto& evals = solver.eigenvalues();

        file << pt.path_coord << "," << pt.kx << "," << pt.ky << ","
             << pt.label;
        for (int b = 0; b < 6; b++)
            file << "," << evals[b];
        file << "\n";
    }

    std::cout << "Band structure written to " << filename
              << " (" << path.size() << " k-points)\n";
}