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
// find_mu — mirrors scipy brentq on NTotal(mu) - N_target
// -----------------------------------------------------------------------------
double find_mu(const Eigensystem& sys, int grid_size, double T, double N_target) {
    // Flatten all eigenvalues into one vector for easy min/max and summation
    const int N_k     = grid_size * grid_size;
    const int N_bands = 6;

    double e_min =  std::numeric_limits<double>::infinity();
    double e_max = -std::numeric_limits<double>::infinity();

    for (const auto& ev : sys.evals) {
        e_min = std::min(e_min, ev.minCoeff());
        e_max = std::max(e_max, ev.maxCoeff());
    }

    const double mu_min = e_min - 5.0 * T;
    const double mu_max = e_max + 5.0 * T;

    // NTotal(mu): Fermi-Dirac sum over all k-points and bands
    auto NTotal = [&](double mu) -> double {
        double n = 0.0;
        for (const auto& ev : sys.evals)
            for (int b = 0; b < N_bands; b++) {
                const double x = std::clamp((ev[b] - mu) / T, -500.0, 500.0);
                n += 1.0 / (std::exp(x) + 1.0);
            }
        return n / static_cast<double>(N_k);
    };

    try {
        return brent([&](double mu) { return NTotal(mu) - N_target; }, mu_min, mu_max);
    } catch (const std::runtime_error&) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

Eigensystem compute_eigensystem_grid(double S, int grid_size, const Params& p) {
    const int N = grid_size * grid_size;

    Eigensystem result;
    result.evals.resize(N);
    result.evecs.resize(N);

    // Precompute k-points (linspace from -pi to pi)
    // We exclude the endpoint to avoid duplicate k-points at the BZ boundary,
    // but this is a minor detail
    std::vector<double> k_lin(grid_size);
    for (int i = 0; i < grid_size; i++)
        k_lin[i] = -M_PI + i * (2.0 * M_PI / (grid_size));

    // Precompute SOC once — it's k-independent
    const Mat6 Hsoc = SOC(p.lam);

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < N; idx++) {
        const int ix = idx % grid_size;
        const int iy = idx / grid_size;

        const double kx = k_lin[ix];
        const double ky = k_lin[iy];

        // Build Hamiltonian for this k-point
        const Mat6 H = H0(kx, ky, p) + Hsoc + HubbardU(S, p.U);

        // Diagonalize — SelfAdjointEigenSolver is not shared, safe per thread
        Eigen::SelfAdjointEigenSolver<Mat6> solver(H);

        result.evals[idx] = solver.eigenvalues();
        result.evecs[idx] = solver.eigenvectors();
    }

    return result;
}

// -----------------------------------------------------------------------------
// calculate_total_energy — compute total energy per unit cell for given S
// -----------------------------------------------------------------------------
double calculate_total_energy(const Eigensystem& sys, int grid_size, double mu, double T) {
    const int N_k     = grid_size * grid_size;
    const int N_bands = 6;

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

// -----------------------------------------------------------------------------
// calculateS — one SCF step: diagonalise, find mu, compute new magnetisation
// mirrors the Python calculateS function
// -----------------------------------------------------------------------------
CalcResult calculateS(double S, int grid_size, double T, double N_target, const Params& p) {
    const int N_k     = grid_size * grid_size;
    const int N_bands = 6;

    const Eigensystem sys = compute_eigensystem_grid(S, grid_size, p);
    const double mu = find_mu(sys, grid_size, T, N_target);
    if (std::isnan(mu)) {
        std::cout << "Could not find mu, skipping update.\n";
        return {S, std::numeric_limits<double>::quiet_NaN()};
    }

    Eigen::Vector<double, 6> densities = Eigen::Vector<double, 6>::Zero();
    for (int idx = 0; idx < N_k; idx++) {
        for (int band = 0; band < N_bands; band++) {
            const double x = std::clamp((sys.evals[idx][band] - mu) / T, -500.0, 500.0);
            const double f = 1.0 / (std::exp(x) + 1.0);
            for (int orb = 0; orb < N_bands; orb++)
                densities[orb] += std::norm(sys.evecs[idx](orb, band)) * f;
        }
    }
    densities /= static_cast<double>(N_k);

    const double n_up = densities.head<3>().sum();
    const double n_dn = densities.tail<3>().sum();
    const double E_total = calculate_total_energy(sys, grid_size, mu, T);

    return {(n_up - n_dn) / 2.0, mu, E_total};
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
                  << ", Difference: " << diff << "\n";

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