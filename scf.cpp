#include "scf.h"
#include "hamiltonian.h"
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <utility>
#include <stdexcept>
#include <array>
#include <utility>
#include <sstream>
#include <omp.h>
#include <fstream>
#include <string>

// -----------------------------------------------------------------------------
// Brent's method
// Finds x in [a, b] such that f(x) = 0, given f(a)*f(b) < 0.
// Combines bisection, secant, and inverse quadratic interpolation.
// -----------------------------------------------------------------------------
double brent(std::function<double(double)> f, double a, double b,
             double tol, int max_iter) {
    double fa = f(a), fb = f(b);
    assert(fa * fb < 0.0 && "brent: bracket [a,b] must straddle zero");

    // Ensure |f(b)| <= |f(a)| — b is always the best guess
    if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }

    double c  = a, fc = fa;
    double s  = b, fs = fb;
    double d  = 0.0;
    bool   mflag = true;  // bisection was used last iteration

    for (int i = 0; i < max_iter; i++) {
        if (std::abs(b - a) < tol || std::abs(fs) < std::numeric_limits<double>::epsilon())
            return s;

        if (fa != fc && fb != fc) {
            // Inverse quadratic interpolation
            s = a*fb*fc / ((fa-fb)*(fa-fc))
              + b*fa*fc / ((fb-fa)*(fb-fc))
              + c*fa*fb / ((fc-fa)*(fc-fb));
        } else {
            // Secant method
            s = b - fb * (b - a) / (fb - fa);
        }

        // Conditions to fall back to bisection
        const double mid = (a + b) / 2.0;
        const bool cond1 = !((3*a + b)/4.0 < s && s < b) && !((3*a + b)/4.0 > s && s > b);
        const bool cond2 =  mflag && std::abs(s - b) >= std::abs(b - c) / 2.0;
        const bool cond3 = !mflag && std::abs(s - b) >= std::abs(c - d) / 2.0;
        const bool cond4 =  mflag && std::abs(b - c) < tol;
        const bool cond5 = !mflag && std::abs(c - d) < tol;

        if (cond1 || cond2 || cond3 || cond4 || cond5) {
            s     = mid;
            mflag = true;
        } else {
            mflag = false;
        }

        fs = f(s);
        d  = c;
        c  = b; fc = fb;

        if (fa * fs < 0.0) { b = s; fb = fs; }
        else               { a = s; fa = fs; }

        // Keep b as the best guess
        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }

    throw std::runtime_error("brent: failed to converge within max_iter");
}

// -----------------------------------------------------------------------------
// NTotal — Fermi-Dirac sum over all k-points and bands, normalised per k-point
// -----------------------------------------------------------------------------
double NTotal(const Eigensystem& sys, double mu, double T) {
    const int N_k     = static_cast<int>(sys.evals.size());
    const int N_bands = 12;

    double n = 0.0;
    for (const auto& ev : sys.evals)
        for (int b = 0; b < N_bands; b++) {
            const double x = std::clamp((ev[b] - mu) / T, -500.0, 500.0);
            n += 1.0 / (std::exp(x) + 1.0);
        }
    return n / static_cast<double>(N_k);
}

// -----------------------------------------------------------------------------
// find_mu — mirrors scipy brentq on NTotal(mu) - N_target
// -----------------------------------------------------------------------------
// Check normalization over k vs N particles
double find_mu(const Eigensystem& sys, double T, double N_target) {

    // Find the lowest and highest energy along k points
    // to bound the chemical potential.
    double e_min =  std::numeric_limits<double>::infinity();
    double e_max = -std::numeric_limits<double>::infinity();

    for (const auto& ev : sys.evals) {
        e_min = std::min(e_min, ev.minCoeff());
        e_max = std::max(e_max, ev.maxCoeff());
    }

    // Add some padding to either end 
    const double mu_min = e_min - 5.0 * T;
    const double mu_max = e_max + 5.0 * T;

    try {
        return brent([&](double mu) { return NTotal(sys, mu, T) - N_target; }, mu_min, mu_max);
    } catch (const std::runtime_error&) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

// Checked for kx=ky=0 S=0 t=1 tdelta=0 t2=0 lam=0 U=0 t_per=0.3
Eigensystem compute_eigensystem_grid(double S, int grid_size, const Params& p) {
    const int N = grid_size * grid_size;

    Eigensystem result;
    result.evals.resize(N);
    result.evecs.resize(N);

    // Precompute k-points (linspace from -pi to pi)
    // We exclude the endpoint to avoid duplicate k-points at the BZ boundary,
    // but this is a minor detail
    // Check this
    std::vector<double> k_lin(grid_size);
    for (int i = 0; i < grid_size; i++)
        k_lin[i] = -M_PI + i * (2.0 * M_PI / (grid_size));

    // Precompute k-independent pieces once
    const Mat6 Hsoc = SOC(p.lam, p.theta, p.phi);
    const Mat6 Hhub = HubbardU(S, p);
    const Mat6 T    = T_perp_mat(p);

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N; idx++) {
        const int ix = idx % grid_size;
        const int iy = idx / grid_size;

        const double kx = k_lin[ix];
        const double ky = k_lin[iy];

        // Assemble 12x12 bilayer Hamiltonian (layer-major block structure)
        Mat12 H = Mat12::Zero();
        H.block<6,6>(0,0) = H0(kx, ky, p, p.delta_cf1) + Hsoc + Hhub;
        H.block<6,6>(6,6) = H0(kx, ky, p, p.delta_cf2) + Hsoc + Hhub;
        H.block<6,6>(0,6) = T;
        H.block<6,6>(6,0) = T;

        Eigen::SelfAdjointEigenSolver<Mat12> solver(H);

        result.evals[idx] = solver.eigenvalues();
        result.evecs[idx] = solver.eigenvectors();
    }

    return result;
}

// -----------------------------------------------------------------------------
// calculate_band_energy — sum of occupied eigenvalues per unit cell (no DC correction)
// -----------------------------------------------------------------------------
double calculate_band_energy(const Eigensystem& sys, double mu, double T) {
    const int N_k     = static_cast<int>(sys.evals.size());
    const int N_bands = 12;

    double E = 0.0;
    for (const auto& ev : sys.evals) {
        for (int b = 0; b < N_bands; b++) {
            const double x = std::clamp((ev[b] - mu) / T, -500.0, 500.0);
            const double f = 1.0 / (std::exp(x) + 1.0);
            E += ev[b] * f;
        }
    }
    return E / static_cast<double>(N_k);
}

// -----------------------------------------------------------------------------
// compute_density_matrix — rho(a,b) = (1/N_k) sum_{k,n} f(e_nk - mu) * v_nk(a)* * v_nk(b)
// Each k-point contributes f_n * v.conjugate() * v.transpose() (outer product)
// Uses thread-local accumulators reduced via a critical section
// -----------------------------------------------------------------------------
Mat12 compute_density_matrix(const Eigensystem& sys, double mu, double T) {
    const int N_k = static_cast<int>(sys.evals.size());

    Mat12 rho = Mat12::Zero();

    #pragma omp parallel
    {
        Mat12 rho_local = Mat12::Zero();

        #pragma omp for schedule(static)
        for (int idx = 0; idx < N_k; idx++) {
            for (int band = 0; band < 12; band++) {
                const double x = std::clamp((sys.evals[idx][band] - mu) / T, -500.0, 500.0);
                const double f = 1.0 / (std::exp(x) + 1.0);
                const Eigen::Vector<cd, 12> v = sys.evecs[idx].col(band);
                rho_local.noalias() += f * v.conjugate() * v.transpose(); 
            }
        }

        #pragma omp critical
        rho += rho_local;
    }

    return rho / static_cast<double>(N_k);
}

// Layer-major ordering: L1↑(0-2), L1↓(3-5), L2↑(6-8), L2↓(9-11)
// Gather spin-up and spin-down spinors across both layers from one eigenvector.
cd spin_cross(const Eigen::Ref<const Eigen::Vector<cd, 12>>& col,
              Eigen::Vector<cd, 6>& psi_up,
              Eigen::Vector<cd, 6>& psi_dn)
{
    psi_up.segment<3>(0) = col.segment<3>(0);  // L1 spin-up
    psi_up.segment<3>(3) = col.segment<3>(6);  // L2 spin-up
    psi_dn.segment<3>(0) = col.segment<3>(3);  // L1 spin-down
    psi_dn.segment<3>(3) = col.segment<3>(9);  // L2 spin-down
    return psi_dn.dot(psi_up);                 //  ψ↓† · ψ↑
}

// -----------------------------------------------------------------------------
// compute_Lz_moments — <Lz> per layer from the 12x12 density matrix
// Lz acts on t2g orbitals as: Lz(yz,xz)=i, Lz(xz,yz)=-i, xy diagonal = 0
// <Lz>_layer = -2 * (Im[rho(xz↑,yz↑)] + Im[rho(xz↓,yz↓)])
// -----------------------------------------------------------------------------
/*
std::pair<double, double> compute_Lz_moments(const Mat12& rho) {
    auto lz_layer = [&](int base) {
        return cd(0, 1) * (rho(base + 0, base + 1) - rho(base + 1, base + 0)  // xz↑,yz↑
                         + rho(base + 3, base + 4) - rho(base + 4, base + 3)); // xz↓,yz↓
    };
    return {-lz_layer(0).real(), -lz_layer(6).real()};
}

// <Ly> = i * sum_{l,σ} (ρ_{xy,yz} - ρ_{yz,xy})
// Orbital indices per layer: yz=0, xz=1, xy=2 (spin-up); +3 for spin-down
std::pair<double, double> compute_Ly_moments(const Mat12& rho) {
    auto ly_layer = [&](int base) {
        return cd(0, 1) * (rho(base + 2, base + 0) - rho(base + 0, base + 2)  // xy↑,yz↑
                         + rho(base + 5, base + 3) - rho(base + 3, base + 5)); // xy↓,yz↓
    };
    return {ly_layer(0).real(), ly_layer(6).real()};
}

// <Lx> = i * sum_{l,σ} (ρ_{xz,xy} - ρ_{xy,xz})
std::pair<double, double> compute_Lx_moments(const Mat12& rho) {
    auto lx_layer = [&](int base) {
        return cd(0, 1) * (rho(base + 1, base + 2) - rho(base + 2, base + 1)  // xz↑,xy↑
                         + rho(base + 4, base + 5) - rho(base + 5, base + 4)); // xz↓,xy↓
    };
    return {lx_layer(0).real(), lx_layer(6).real()};
}
*/
std::array<std::pair<double,double>, 3> compute_L_moments(const Mat12& rho, const Params& p){
    const auto [Lx, Ly, Lz] = l_matrices();

    const Eigen::Matrix<cd, 2, 2> I2 = Eigen::Matrix<cd, 2, 2>::Identity();
    const Eigen::Matrix<cd, 6, 6> Lx_op = kron(I2, Lx);
    const Eigen::Matrix<cd, 6, 6> Ly_op = kron(I2, Ly);
    const Eigen::Matrix<cd, 6, 6> Lz_op = kron(I2, Lz);

    auto layer_moment = [&](const Eigen::Matrix<cd, 6, 6>& op, int base) {
        return (rho.block<6,6>(base, base) * op).trace().real();
    };

    return {{
        {layer_moment(Lx_op, 0), layer_moment(Lx_op, 6)},
        {layer_moment(Ly_op, 0), layer_moment(Ly_op, 6)},
        {layer_moment(Lz_op, 0), layer_moment(Lz_op, 6)}
    }};
}

// Spin Moments
// Returns {(Sx_L1, Sx_L2), (Sy_L1, Sy_L2), (Sz_L1, Sz_L2)}
// Operator: I_2 ⊗ s_i ⊗ I_3 — identity on layer, spin matrix, identity on orbital.
// Per-layer value: Tr(rho_layer * kron(s_i, I_3)).
std::array<std::pair<double,double>, 3> compute_S_moments(const Mat12& rho, const Params& p){
    const auto [sx, sy, sz] = spin_matrices(p.theta, p.phi);

    const Eigen::Matrix<cd, 3, 3> I3 = Eigen::Matrix<cd, 3, 3>::Identity();
    const Eigen::Matrix<cd, 6, 6> Sx_op = kron(sx, I3);
    const Eigen::Matrix<cd, 6, 6> Sy_op = kron(sy, I3);
    const Eigen::Matrix<cd, 6, 6> Sz_op = kron(sz, I3);

    auto layer_moment = [&](const Eigen::Matrix<cd, 6, 6>& op, int base) {
        return (rho.block<6,6>(base, base) * op).trace().real();
    };

    std::cout << "Compute S";
    std::cout << p.theta << " " << p.phi << "\n";

    return {{
        {layer_moment(Sx_op, 0), layer_moment(Sx_op, 6)},
        {layer_moment(Sy_op, 0), layer_moment(Sy_op, 6)},
        {layer_moment(Sz_op, 0), layer_moment(Sz_op, 6)}
    }};
}
// -----------------------------------------------------------------------------
// printKanamoriOccupations — orbital occupations, spin/Lz moments, total energy
// -----------------------------------------------------------------------------
static void print_density_matrix(const Mat12& rho) {
    constexpr double tol = 5e-5;

    auto fmt = [&](cd v) -> std::string {
        const bool re_zero = std::abs(v.real()) < tol;
        const bool im_zero = std::abs(v.imag()) < tol;
        if (re_zero && im_zero) return "0";
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(4);
        ss << v.real() << "+" << v.imag() << "i";
        /*
        if (re_zero) {
            ss << v.imag() << "i";
        } else if (im_zero) {
            ss << v.real();
        } else {
            ss << v.real();
            if (v.imag() >= 0.0) ss << "+";
            ss << v.imag() << "i";
        }
            */
        return ss.str();
    };

    for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 12; j++) {
            if (j > 0) std::cout << "  ";
            std::cout << std::setw(15) << std::right << fmt(rho(i, j));
        }
        std::cout << "\n";
    }
}

// Takes a filename parameter with a default value of "density_matrix.txt"
void save_density_matrix(const Mat12& rho, const std::string& filename) {
    // Open the file stream
    std::ofstream outfile(filename);

    // Ensure the file was successfully opened before trying to write
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing.\n";
        return;
    }

    constexpr double tol = 5e-5;

    auto fmt = [&](cd v) -> std::string {
        const bool re_zero = std::abs(v.real()) < tol;
        const bool im_zero = std::abs(v.imag()) < tol;
        
        if (re_zero && im_zero) return "0";
        
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(4);
        
        // Cleaned up formatting to avoid double signs like "+-"
        if (re_zero) {
            ss << v.imag() << "i";
        } else if (im_zero) {
            ss << v.real();
        } else {
            ss << v.real();
            if (v.imag() >= 0.0) ss << "+"; // Only add explicit '+' if positive
            ss << v.imag() << "i";
        }
        
        return ss.str();
    };

    // Write to the file exactly as you were writing to cout
    for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 12; j++) {
            if (j > 0) outfile << "  ";
            outfile << std::setw(15) << std::right << fmt(rho(i, j));
        }
        outfile << "\n";
    }

    // Close the file stream
    outfile.close();
    std::cout << "Density matrix saved to " << filename << "\n";
}

Mat12 load_density_matrix(const std::string& filename) {
    // Initialize an empty matrix of zeros
    Mat12 rho = Mat12::Zero(); 
    std::ifstream infile(filename);

    if (!infile.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for reading.\n";
        return rho; // Returns a zero matrix if the file is missing
    }

    // Helper function to decode strings like "3.5+2.1i", "-4.0i", "1.2", or "0"
    auto parse_complex = [](std::string s) -> std::complex<double> {
        if (s == "0") return {0.0, 0.0};

        bool has_i = (s.back() == 'i');
        if (has_i) s.pop_back(); // Chop off the 'i' so we can parse the numbers

        // If there was no 'i', it's a purely real number
        if (!has_i) return {std::stod(s), 0.0}; 

        // Look for a '+' or '-' separating the real and imaginary parts.
        // Because of the fixed formatting we used earlier, we don't have to 
        // worry about 'e-5' scientific notation throwing off the negative sign.
        size_t sign_pos = s.find_last_of("+-");

        // If the only sign is at the very front (e.g., "-2.0"), it's purely imaginary
        if (sign_pos == 0 || sign_pos == std::string::npos) {
            return {0.0, std::stod(s)};
        }

        // Otherwise, split the string at the sign
        std::string re_str = s.substr(0, sign_pos);
        std::string im_str = s.substr(sign_pos); // Includes the sign!

        return {std::stod(re_str), std::stod(im_str)};
    };

    std::string token;
    for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 12; j++) {
            // "infile >> token" grabs the next block of text separated by spaces/newlines
            if (!(infile >> token)) {
                std::cerr << "Error: File ended early or formatting is broken.\n";
                return rho;
            }
            rho(i, j) = parse_complex(token);
        }
    }

    infile.close();
    return rho;
}

void printKanamoriOccupations(const KanamoriResult& res, const Params& p) {
    const Mat12& rho = res.rho;

    // Orbital Moments
    const auto lmom = compute_L_moments(rho, p);
    const auto [lx1, lx2] = lmom[0];
    const auto [ly1, ly2] = lmom[1];
    const auto [lz1, lz2] = lmom[2];
    const double l110_1 = (lx1 + ly1) / std::sqrt(2.0);
    const double l110_2 = (lx2 + ly2) / std::sqrt(2.0);

    // Spin Moments
    Params pcopy = p;
    pcopy.theta = 0;
    pcopy.phi = 0;
    const auto smom = compute_S_moments(rho, pcopy);
    const auto [sx1, sx2] = smom[0];
    const auto [sy1, sy2] = smom[1];
    const auto [sz1, sz2] = smom[2];
    const double s110_1 = (sx1 + sy1) / std::sqrt(2.0);
    const double s110_2 = (sx2 + sy2) / std::sqrt(2.0);

    struct OrbEntry { const char* name; int up; int dn; };
    const OrbEntry orbs[] = {{"dxy", 2, 5}, {"dxz", 1, 4}, {"dyz", 0, 3}};

    double total_up = 0.0, total_dn = 0.0;
    double tlayer1 = 0.0, tlayer2 = 0.0;
    double sm1 = 0.0, sm2 = 0.0;

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
                      << "   " << std::right << std::fixed << std::setprecision(6)
                      << up << "    " << dn << "\n";
        }
        std::cout << "   Spin Up  " << layer_up << "\n";
        std::cout << "   Spin Dn  " << layer_dn << "\n";
        std::cout << "   Total    " << (layer_up + layer_dn) << "\n";
        std::cout << "   Moment   " << (layer_up - layer_dn) << "\n";
        if (layer == 0) {
            std::cout << "   Lx       " << lx1 << "\n";
            std::cout << "   Ly       " << ly1 << "\n";
            std::cout << "   Lz       " << lz1 << "\n";
            std::cout << "   L[110]   " << l110_1 << "\n";
            std::cout << "   Sx       " << sx1 << "\n";
            std::cout << "   Sy       " << sy1 << "\n";
            std::cout << "   Sz       " << sz1 << "\n";
            std::cout << "   S[110]   " << s110_1 << "\n\n";
            tlayer1 = layer_up + layer_dn;
            sm1 = layer_up - layer_dn;
        }
        if (layer == 1) {
            std::cout << "   Lx       " << lx2 << "\n";
            std::cout << "   Ly       " << ly2 << "\n";
            std::cout << "   Lz       " << lz2 << "\n";
            std::cout << "   L[110]   " << l110_2 << "\n";
            std::cout << "   Sx       " << sx2 << "\n";
            std::cout << "   Sy       " << sy2 << "\n";
            std::cout << "   Sz       " << sz2 << "\n";
            std::cout << "   S[110]   " << s110_2 << "\n\n";
            tlayer2 = layer_up + layer_dn;
            sm2 = layer_up - layer_dn;
        }
    }

    std::cout << "Total e-    " << (total_up + total_dn) << "\n";
    std::cout << "Diff e-     " << (tlayer1 - tlayer2) << "\n";
    std::cout << "  Moment    " << (total_up - total_dn) << "\n";
    std::cout << "Lx Total    " << (lx1 + lx2) << "\n";
    std::cout << "Ly Total    " << (ly1 + ly2) << "\n";
    std::cout << "Lz Total    " << (lz1 + lz2) << "\n";
    std::cout << "L[110] Tot  " << (l110_1 + l110_2) << "\n";
    std::cout << "Sx Total    " << (sx1 + sx2) << "\n";
    std::cout << "Sy Total    " << (sy1 + sy2) << "\n";
    std::cout << "Sz Total    " << (sz1 + sz2) << "\n";
    std::cout << "S[110] Tot  " << (s110_1 + s110_2) << "\n\n";

    std::cout << "Total energy (eV):    " << res.E_total << "\n\n";

    std::cout << "Compare with DFT \n";
    std::cout << "\nOccupations \n";
    std::cout << "          n(xz/yz) down     n(xy) down \n";
    std::cout << "Co1       " << rho(4,4).real() << "            " << rho(5,5).real() << "\n";
    std::cout << "Co2       " << rho(10,10).real() << "            " << rho(11,11).real() << "\n";
    std::cout << "\nMoments \n";
    std::cout << "          Spin moment       Lz moment       L[110] moment\n";
    std::cout << "Co1       " << sm1 << "            " << lz1 << "            " << l110_1 << "\n";
    std::cout << "Co2       " << sm2 << "            " << lz2 << "            " << l110_2 << "\n\n";

    std::cout << "--- Initial density matrix (rho0) ---\n";
    print_density_matrix(res.rho0);
    std::cout << "\n--- Final density matrix (rho) ---\n";
    print_density_matrix(res.rho);
    std::cout << "\n";
}

// -----------------------------------------------------------------------------
// calculateS — one SCF step: diagonalise, find mu, compute new magnetisation
// mirrors the Python calculateS function
// -----------------------------------------------------------------------------
CalcResult calculateS(double S, int grid_size, double T, double N_target, const Params& p) {
    const int N_k     = grid_size * grid_size;
    const int N_bands = 12;

    const Eigensystem sys = compute_eigensystem_grid(S, grid_size, p);
    const double mu = find_mu(sys, T, N_target);
    if (std::isnan(mu)) {
        std::cout << "Could not find mu, skipping update.\n";
        return {S, std::numeric_limits<double>::quiet_NaN()};
    }

    // Exchange field is always along z, so S = (1/2)(n_up - n_dn) in the z basis.
    // Angles enter only through SOC; the self-consistency condition is always z-axis.
    double spin_sum = 0.0;
    for (int idx = 0; idx < N_k; idx++) {
        for (int band = 0; band < N_bands; band++) {
            const double x = std::clamp((sys.evals[idx][band] - mu) / T, -500.0, 500.0);
            const double f = 1.0 / (std::exp(x) + 1.0);

            Eigen::Vector<cd, 6> psi_up, psi_dn;
            spin_cross(sys.evecs[idx].col(band), psi_up, psi_dn);

            spin_sum += f * (psi_up.squaredNorm() - psi_dn.squaredNorm());
        }
    }

    // Hubbard DC correction: +n_orb * U * S² (6 orbitals × 2 layers = 6 sites total)
    const double E_total = calculate_band_energy(sys, mu, T) + 6.0 * p.U * S * S;
    return {spin_sum / (2.0 * N_k), mu, E_total};
}

// -----------------------------------------------------------------------------
// runSelfCalc — self-consistent loop with linear mixing
// S_new = alpha * S_calc + (1 - alpha) * S_current
// -----------------------------------------------------------------------------
CalcResult runSelfCalc(double S0, double alpha, int grid_size,
                                       double T, double N_target, const Params& p) {
    constexpr int    max_iter = 5000;
    constexpr double tol      = 1e-5;

    double S_current = S0;
    double mu_current = std::numeric_limits<double>::quiet_NaN();
    double E_current = std::numeric_limits<double>::quiet_NaN();
    bool converged = false;

    for (int i = 0; i < max_iter; i++) {
        auto result = calculateS(S_current, grid_size, T, N_target, p);
        double S_calc = result.S_new;
        mu_current = result.mu;
        E_current  = result.E_total;

        const double diff = std::abs(S_calc - S_current);

        std::cout << "\rIteration " << i
                  << ", S = "    << std::fixed << std::setprecision(6) << S_current
                  << ", Difference: " << diff << "   " << std::flush;

        if (diff < tol) {
            std::cout << "\nSelf-consistency reached!\n";
            converged = true;
            break;
        }

        S_current = alpha * S_calc + (1.0 - alpha) * S_current;
    }

    if (!converged)
        throw std::runtime_error("runSelfCalc: failed to converge within max_iter");

    std::cout << "Final S: " << S_current << ", mu: " << mu_current << ", E_total: " << E_current << "\n";
    return {S_current, mu_current, E_current};
}

// -----------------------------------------------------------------------------
// Kanamori SCF
// -----------------------------------------------------------------------------

// Double-counting correction to the total energy from the Kanamori MF decoupling.
// Computed per single layer from the 6x6 density matrix block.
// Matches the ΔE expression derived in the LaTeX notes (Section: Complete Kanamori Hamiltonian).
// Make test DC should be only real
double kanamori_dc_layer(const Mat6& rho, const KanamoriParams& kp) {
    const double U  = kp.U;
    const double Up = kp.U_prime;
    const double J  = kp.J;

    double dc = 0.0;

    // -U intraorbital:
    for (int m = 0; m < 3; m++){
        dc -= U * (rho(m, m+3)*rho(m+3, m) - 
                    rho(m,m)*rho(m+3,m+3)).real();
    }

    // -U' interorbital opposite-spin: -U' sum_{m≠m'} <n_{m↑}><n_{m'↓}>
    for (int m = 0; m < 3; m++)
        for (int mp = 0; mp < 3; mp++)
            if (mp != m)
                dc -= Up * (rho(m, mp+3) * rho(mp+3, m) - 
                             rho(m,m) * rho(mp+3, mp+3)).real();
    
    std::cout << "DC after -U and -U': " << dc << "\n";

    // -(U'-J) same-spin: -(U'-J) sum_{m<m', σ} <n_{mσ}><n_{m'σ}>
    for (int m = 0; m < 3; m++)
        for (int mp = m+1; mp < 3; mp++) {
            dc -= (Up - J) * (rho(m, mp)*rho(mp,m)
                            - rho(m,m)*rho(mp,mp)).real();  // σ=↑
            dc -= (Up - J) * (rho(m+3, mp+3)*rho(mp+3,m+3)
                            - rho(m+3,m+3)*rho(mp+3,mp+3)).real(); // σ=↓
        }

    // J exchange + pair-hopping DC: sum over all ordered pairs m≠m'
    for (int m = 0; m < 3; m++)
        for (int mp = 0; mp < 3; mp++) {
            if (mp == m) continue;
            // Exchange: <d†_{m↑}d_{m'↑}><d†_{m'↓}d_{m↓}> - <d†_{m↑}d_{m↓}><d†_{m'↓}d_{m'↑}>
            dc += J * (rho(m, mp) * rho(mp+3, m+3)
                     - rho(m, m+3) * rho(mp+3, mp)).real();
            // Pair-hopping: <d†_{m↑}d_{m'↑}><d†_{m↓}d_{m'↓}> - <d†_{m↑}d_{m'↓}><d†_{m↓}d_{m'↑}>
            dc += J * (rho(m, mp) * rho(m+3, mp+3)
                     - rho(m, mp+3) * rho(m+3, mp)).real();
        }

    return dc;
}

// Eigensystem with Kanamori MF: replaces HubbardU(S) with KanamoriMF(rho).
// k-independent pieces (SOC, T_perp, KanamoriMF) are precomputed outside the k-loop.
Eigensystem compute_eigensystem_kanamori(const Mat12& rho, int grid_size,
                                          const Params& p, const KanamoriParams& kp) {
    
    // Initalizing eigensystem result structure and allocating sizes
    const int N = grid_size * grid_size;
    Eigensystem result;
    result.evals.resize(N);
    result.evecs.resize(N);

    // Linspace for k grid from -pi to pi
    // First Brillioun zone of square lattice
    // Will leave out M_PI which is already counted in -M_PI
    std::vector<double> k_lin(grid_size);
    for (int i = 0; i < grid_size; i++)
        k_lin[i] = -M_PI + i * (2.0 * M_PI / grid_size);

    // Precompute all k-independent contributions into a single 12x12 matrix
    const Mat6 Hsoc  = SOC(p.lam, p.theta, p.phi);
    const Mat6 Tperp = T_perp_mat(p);

    Mat12 H_kfree = KanamoriMF(rho, kp);
    H_kfree.block<6,6>(0, 0) += Hsoc;
    H_kfree.block<6,6>(6, 6) += Hsoc;
    H_kfree.block<6,6>(0, 6) += Tperp;
    H_kfree.block<6,6>(6, 0) += Tperp;
    H_kfree += staggered_potential(p);

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N; idx++) {
        const double kx = k_lin[idx % grid_size];
        const double ky = k_lin[idx / grid_size];

        Mat12 H = H_kfree;
        H.block<6,6>(0, 0) += H0(kx, ky, p, p.delta_cf1);
        H.block<6,6>(6, 6) += H0(kx, ky, p, p.delta_cf2);

        Eigen::SelfAdjointEigenSolver<Mat12> solver(H);
        result.evals[idx] = solver.eigenvalues();
        result.evecs[idx] = solver.eigenvectors();
    }

    return result;
}

// Self-consistent loop converging the full density matrix.
// Convergence: Frobenius norm ||rho_new - rho||_F < tol.
// Total energy = band sum + double-counting correction from converged rho.
// Mixing (selected by `mixer`):
//   LinearDIIS : linear (α) for the first diis_start iterations, then Pulay DIIS.
//   Broyden    : modified Broyden second method (Johnson, PRB 38, 12807 (1988)).
KanamoriResult runKanamoriSCF(const Mat12& rho0, double alpha, int grid_size,
                               double T, double N_target,
                               const Params& p, const KanamoriParams& kp,
                               MixerType mixer) {
    constexpr int    max_iter   = 999999;
    double tol        = kp.tol;
    constexpr int    diis_max   = 8;
    constexpr int    diis_start = 200;

    Mat12 rho = rho0;

    // Frobenius inner product Re<A, B>_F (matches the DIIS residual metric).
    auto frob = [](const Mat12& A, const Mat12& B) {
        return (A.array().conjugate() * B.array()).sum().real();
    };

    // --- LinearDIIS state ---
    std::vector<Mat12> diis_rho;  // rho_new history
    std::vector<Mat12> diis_err;  // residual history: e_i = rho_new_i - rho_i

    constexpr int no_improve_max    = 9999999;
    constexpr int linear_reset_steps = 150;
    double best_diis_diff   = std::numeric_limits<double>::max();
    int    no_improve_count = 0;
    int    linear_remaining = 0;  // counts down during post-stall linear phase

    // --- Broyden state ---
    constexpr int    broyden_max = 5;     // max history length
    constexpr double broyden_w0  = 0.01;  // diagonal regularisation weight
    std::vector<Mat12> broyden_dF;  // normalised residual differences
    std::vector<Mat12> broyden_u;   // corresponding update vectors
    Mat12  F_prev      = Mat12::Zero();
    Mat12  rho_in_prev = Mat12::Zero();
    bool   broyden_have_prev = false;

    for (int i = 0; i < max_iter; i++) {
        const Eigensystem sys = compute_eigensystem_kanamori(rho, grid_size, p, kp);
        const double mu = find_mu(sys, T, N_target);
        if (std::isnan(mu)) {
            std::cout << "runKanamoriSCF: could not find mu at iteration " << i << "\n";
            break;
        }

        const Mat12 rho_new = compute_density_matrix(sys, mu, T);
        const double diff = (rho_new - rho).norm();

        const bool using_diis = (mixer == MixerType::LinearDIIS)
                              && (i >= diis_start) && (linear_remaining == 0);
        const char* tag = (mixer == MixerType::Broyden)
                              ? (broyden_have_prev ? " [Broy]" : "  [mix]")
                              : (using_diis        ? " [DIIS]" : "  [mix]");
        std::cout << "\rIteration " << std::setw(5) << i << tag
                  << ", |Δρ|_F = " << std::scientific << std::setprecision(4) << diff
                  << "   " << std::flush;

        if (diff < tol) {
            const double bandsum = calculate_band_energy(sys, mu, T);
            const double dc = kanamori_dc_layer(rho.block<6,6>(0, 0), kp)
                            + kanamori_dc_layer(rho.block<6,6>(6, 6), kp);
            std::cout << "\nKanamori SCF converged! mu = " << mu
                      << "\n Band Energy = " << bandsum
                      << "\n DC Correction = " << dc
                      << ", E_total = " << bandsum - dc << "\n";
            return {rho0, rho, mu, bandsum - dc};
        }

        if (mixer == MixerType::Broyden) {
            // Modified Broyden second method (Johnson 1988).
            // Residual of the fixed-point map ρ ↦ ρ_new: F = ρ_new − ρ.
            const Mat12 F = rho_new - rho;

            if (!broyden_have_prev) {
                // Bootstrap with a single linear-mixing step.
                rho_in_prev = rho;
                F_prev      = F;
                broyden_have_prev = true;
                rho += alpha * F;
            } else {
                const Mat12  dF_raw = F - F_prev;
                const double nrm    = std::sqrt(frob(dF_raw, dF_raw));

                if (nrm > 0.0) {
                    broyden_dF.push_back(dF_raw / nrm);
                    broyden_u.push_back(alpha * (dF_raw / nrm)
                                        + (rho - rho_in_prev) / nrm);

                    if (static_cast<int>(broyden_dF.size()) > broyden_max) {
                        broyden_dF.erase(broyden_dF.begin());
                        broyden_u.erase(broyden_u.begin());
                    }
                }

                rho_in_prev = rho;
                F_prev      = F;

                const int m = static_cast<int>(broyden_dF.size());
                Mat12 rho_next = rho + alpha * F;  // linear step + Broyden correction

                if (m > 0) {
                    // a_ij = <dF_i, dF_j>, regularised; c_k = <dF_k, F>; γ = (w0²I + a)⁻¹ c
                    Eigen::MatrixXd a = Eigen::MatrixXd::Zero(m, m);
                    Eigen::VectorXd c(m);
                    for (int ii = 0; ii < m; ii++) {
                        for (int jj = ii; jj < m; jj++) {
                            const double aij = frob(broyden_dF[ii], broyden_dF[jj]);
                            a(ii, jj) = aij;
                            a(jj, ii) = aij;
                        }
                        a(ii, ii) += broyden_w0 * broyden_w0;
                        c(ii) = frob(broyden_dF[ii], F);
                    }

                    const Eigen::VectorXd gamma = a.colPivHouseholderQr().solve(c);
                    for (int l = 0; l < m; l++)
                        rho_next.noalias() -= gamma(l) * broyden_u[l];
                }

                rho = rho_next;
            }
        } else if (!using_diis) {
            if (linear_remaining > 0) --linear_remaining;
            rho = alpha * rho_new + (1.0 - alpha) * rho;
        } else {
            // Track improvement; on stall, do linear_reset_steps of linear mixing then retry DIIS.
            if (diff < best_diis_diff) {
                best_diis_diff  = diff;
                no_improve_count = 0;
            } else {
                ++no_improve_count;
            }

            if (no_improve_count >= no_improve_max && diff > 1e-7) {
                std::cout << "\n[DIIS stalled " << no_improve_max
                          << " iters — " << linear_reset_steps << " linear mix steps]\n";
                diis_rho.clear();
                diis_err.clear();
                best_diis_diff   = std::numeric_limits<double>::max();
                no_improve_count = 0;
                linear_remaining = linear_reset_steps;
                rho = alpha * rho_new + (1.0 - alpha) * rho;
                continue;
            }

            // Pulay DIIS
            diis_rho.push_back(rho_new);
            diis_err.push_back(rho_new - rho);

            if (static_cast<int>(diis_rho.size()) > diis_max) {
                diis_rho.erase(diis_rho.begin());
                diis_err.erase(diis_err.begin());
            }

            const int m = static_cast<int>(diis_rho.size());

            if (m < 2) {
                rho = rho_new;
            } else {
                // Build (m+1)×(m+1) Pulay system
                //   [B   -1] [c]   [ 0]
                //   [-1   0] [λ] = [-1]
                // B_ij = Re(<e_i, e_j>_F),  Σc_i = 1 enforced by λ
                Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m + 1, m + 1);
                Eigen::VectorXd b = Eigen::VectorXd::Zero(m + 1);
                b(m) = -1.0;

                for (int ii = 0; ii < m; ii++) {
                    for (int jj = ii; jj < m; jj++) {
                        const double Bij = frob(diis_err[ii], diis_err[jj]);
                        A(ii, jj) = Bij;
                        A(jj, ii) = Bij;
                    }
                    A(ii, m) = -1.0;
                    A(m, ii) = -1.0;
                }

                const Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);

                rho = Mat12::Zero();
                for (int ii = 0; ii < m; ii++)
                    rho.noalias() += c(ii) * diis_rho[ii];
            }
        }
    }

    throw std::runtime_error("runKanamoriSCF: failed to converge within max_iter");
}