#include "scf.h"
#include "hamiltonian.h"
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <utility>
#include <stdexcept>
#include <omp.h>

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
    double e_min =  std::numeric_limits<double>::infinity();
    double e_max = -std::numeric_limits<double>::infinity();

    for (const auto& ev : sys.evals) {
        e_min = std::min(e_min, ev.minCoeff());
        e_max = std::max(e_max, ev.maxCoeff());
    }

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
    const Mat6 Hsoc = SOC(p.lam);
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

    // Magnetisation direction n̂ = (sin θ cos φ, sin θ sin φ, cos θ)
    const double nx = std::sin(p.theta) * std::cos(p.phi);
    const double ny = std::sin(p.theta) * std::sin(p.phi);
    const double nz = std::cos(p.theta);

    // Project the density matrix onto n̂: M = Tr(ρ (n̂·σ)) summed over orbitals and k
    // For each eigenstate: ⟨n̂·σ⟩ = nz(|ψ↑|²-|ψ↓|²) + 2·nx·Re[ψ↑†ψ↓] + 2·ny·Im[ψ↑†ψ↓]
    // spin_sum is actually definition for magnetization
    // to return Stoner S we divide by two (in return statment)
    double spin_sum = 0.0;
    for (int idx = 0; idx < N_k; idx++) {
        for (int band = 0; band < N_bands; band++) {
            const double x = std::clamp((sys.evals[idx][band] - mu) / T, -500.0, 500.0);
            const double f = 1.0 / (std::exp(x) + 1.0);

            Eigen::Vector<cd, 6> psi_up, psi_dn;
            const cd cross = spin_cross(sys.evecs[idx].col(band), psi_up, psi_dn);

            const double proj = nz * (psi_up.squaredNorm() - psi_dn.squaredNorm())
                              + 2.0 * nx * cross.real()
                              + 2.0 * ny * cross.imag();

            spin_sum += f * proj;
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

        std::cout << "Iteration " << i
                  << ", S = "    << std::fixed << std::setprecision(6) << S_current
                  << ", Difference: " << diff << std::endl;

        if (diff < tol) {
            std::cout << "Self-consistency reached!\n";
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
namespace {
double kanamori_dc_layer(const Mat6& rho, const KanamoriParams& kp) {
    const double U  = kp.U;
    const double Up = kp.U_prime;
    const double J  = kp.J;

    double dc = 0.0;

    // -U intraorbital:
    for (int m = 0; m < 3; m++)
        dc -= U * (rho(m, m+3).real()*rho(m+3, m).real() - 
                    rho(m,m).real()*rho(m+3,m+3).real());

    // -U' interorbital opposite-spin: -U' sum_{m≠m'} <n_{m↑}><n_{m'↓}>
    for (int m = 0; m < 3; m++)
        for (int mp = 0; mp < 3; mp++)
            if (mp != m)
                dc -= Up * (rho(m, mp+3).real() * rho(mp+3, m).real() - 
                             rho(m,m).real() * rho(mp+3, mp+3).real());

    // -(U'-J) same-spin: -(U'-J) sum_{m<m', σ} <n_{mσ}><n_{m'σ}>
    for (int m = 0; m < 3; m++)
        for (int mp = m+1; mp < 3; mp++) {
            dc -= (Up - J) * (rho(m, mp).real()*rho(mp,m).real()
                            - rho(m,m).real()*rho(mp,mp).real());  // σ=↑
            dc -= (Up - J) * (rho(m+3, mp+3).real()*rho(mp+3,m+3).real()
                            - rho(m+3,m+3).real()*rho(mp+3,mp+3).real()); // σ=↓
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
} // anonymous namespace

// Eigensystem with Kanamori MF: replaces HubbardU(S) with KanamoriMF(rho).
// k-independent pieces (SOC, T_perp, KanamoriMF) are precomputed outside the k-loop.
Eigensystem compute_eigensystem_kanamori(const Mat12& rho, int grid_size,
                                          const Params& p, const KanamoriParams& kp) {
    const int N = grid_size * grid_size;
    Eigensystem result;
    result.evals.resize(N);
    result.evecs.resize(N);

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
KanamoriResult runKanamoriSCF(const Mat12& rho0, double alpha, int grid_size,
                               double T, double N_target,
                               const Params& p, const KanamoriParams& kp) {
    constexpr int    max_iter = 99999;
    constexpr double tol      = 1e-6;

    Mat12 rho = rho0;

    for (int i = 0; i < max_iter; i++) {
        const Eigensystem sys = compute_eigensystem_kanamori(rho, grid_size, p, kp);
        const double mu = find_mu(sys, T, N_target);
        if (std::isnan(mu)) {
            std::cout << "runKanamoriSCF: could not find mu at iteration " << i << "\n";
            break;
        }

        const Mat12 rho_new = compute_density_matrix(sys, mu, T);
        const double diff = (rho_new - rho).norm();

        std::cout << "Iteration " << i
                  << ", |Δρ|_F = " << std::scientific << std::setprecision(4) << diff << std::endl;

        if (diff < tol) {
            const double bandsum = calculate_band_energy(sys, mu, T);
            const double dc = kanamori_dc_layer(rho.block<6,6>(0, 0), kp)
                            + kanamori_dc_layer(rho.block<6,6>(6, 6), kp);
            std::cout << "Kanamori SCF converged! mu = " << mu
                      << ", E_total = " << bandsum - dc << std::endl;
            return {rho, mu, bandsum - dc};
        }

        rho = alpha * rho_new + (1.0 - alpha) * rho;
    }

    throw std::runtime_error("runKanamoriSCF: failed to converge within max_iter");
}