#pragma once

#include "hamiltonian.h"
#include <vector>
#include <functional>
#include <cassert>
#include <utility>

struct Eigensystem {
    std::vector<Eigen::Vector<double, 12>> evals; // (grid_size^2, 12)
    std::vector<Mat12>                     evecs; // (grid_size^2, 12x12)
};

struct CalcResult { double S_new; double mu; double E_total; };

// Extract spin-up/down spinors from a bilayer eigenvector and return ψ↑† · ψ↓
cd spin_cross(const Eigen::Ref<const Eigen::Vector<cd, 12>>& col,
              Eigen::Vector<cd, 6>& psi_up,
              Eigen::Vector<cd, 6>& psi_dn);

// Brent's method root finder — f must be continuous and f(a)*f(b) < 0
double brent(std::function<double(double)> f, double a, double b,
             double tol = 1e-10, int max_iter = 100);

// Fermi-Dirac sum over all k-points and bands, normalised per k-point
double NTotal(const Eigensystem& sys, double mu, double T);

// Find chemical potential via Brent on the total electron count
double find_mu(const Eigensystem& sys, double T = 0.05, double N_target = 5.0);

// Compute the eigensystem on a k-grid for a given Stoner parameter S
Eigensystem compute_eigensystem_grid(double S, int grid_size = 200, const Params& p = Params{});

// Sum of occupied eigenvalues per unit cell (no double-counting correction)
double calculate_band_energy(const Eigensystem& sys, double mu, double T);

// Compute the full 12x12 density matrix: rho(a,b) = <c†_a c_b>
// Averaged over the k-grid using Fermi-Dirac weights
Mat12 compute_density_matrix(const Eigensystem& sys, double mu, double T);

// Calculate Stoner parameter S for a given guess, and the corresponding chemical potential
CalcResult  calculateS(double S, int grid_size, double T, double N_target, const Params& p);

// Run the self-consistent loop with linear mixing
CalcResult  runSelfCalc(double S0, double alpha, int grid_size, double T, double N_target, const Params& p);

// Orbital Lz moment per layer: <Lz>_layer = Tr(rho_layer * (Lz ⊗ I_spin))
// Returns {Lz_layer1, Lz_layer2}. Uses t2g Lz: Lz(yz,xz)=i, Lz(xz,yz)=-i.
std::pair<double, double> compute_Lz_moments(const Mat12& rho);

// ---- Kanamori SCF ----

struct KanamoriResult { Mat12 rho0; Mat12 rho; double mu; double E_total; };

// Print orbital occupations, spin moments, Lz moments, and total energy.
void printKanamoriOccupations(const KanamoriResult& res);

// Compute eigensystem with the full Kanamori MF Hamiltonian (density matrix input)
Eigensystem compute_eigensystem_kanamori(const Mat12& rho, int grid_size,
                                         const Params& p = Params{},
                                         const KanamoriParams& kp = KanamoriParams{});

// Self-consistent loop converging the full 12x12 density matrix
// rho0: initial guess (e.g. from a converged Stoner SCF)
// Convergence criterion: Frobenius norm ||rho_new - rho||_F < tol
KanamoriResult runKanamoriSCF(const Mat12& rho0, double alpha, int grid_size,
                               double T, double N_target,
                               const Params& p = Params{},
                               const KanamoriParams& kp = KanamoriParams{});