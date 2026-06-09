#pragma once

#include <Eigen/Dense>
#include <complex>
#include <utility>

using Mat6  = Eigen::Matrix<std::complex<double>, 6, 6>;
using Mat12 = Eigen::Matrix<std::complex<double>, 12, 12>;
using cd    = std::complex<double>;

// Model parameters with sensible defaults
struct Params {
    double tpi      = 1.0;
    double tdelta = 0;
    double t2xy      = 0;
    double t2yzxz      = 0;
    double tg = 0;
    double lam     = 0.1;
    double U       = 0.0;
    double theta   = 0.0;  // polar angle of SOC spin quantization axis (0 = z-axis / out-of-plane)
    double phi     = 0.0;  // azimuthal angle of SOC spin quantization axis
    double t_perp    = 0.0;  // interlayer hopping for yz and xz orbitals
    double t_perp_xy = 0.0;  // interlayer hopping for xy orbital (zero by symmetry for ideal stacking)
    double delta_cf1 = 0.0;  // tetragonal crystal field for layer 1 (positive = xy higher)
    double delta_cf2 = 0.0;  // tetragonal crystal field for layer 2
    double delta_V   = 0.0;  // staggered layer potential: layer 1 shifts by -delta_V, layer 2 by +delta_V
};

// Kinetic hopping term (k-dependent, real diagonal). delta_cf shifts the xy orbital energy.
Mat6 H0(double kx, double ky, const Params& p = Params{}, double delta_cf = 0.0);

// Spin-orbit coupling (k-independent, precompute once).
// theta/phi rotate the spin quantization axis via U†(θ,φ)·s·U(θ,φ) on each spin operator.
Mat6 SOC(double lam = 0.1, double theta = 0.0, double phi = 0.0);

// Hubbard mean-field term (k-independent, off-diagonal in spin for general n̂)
Mat6 HubbardU(double S, const Params& p = Params{});

// Full single-layer Hamiltonian
Mat6 singleLayer(double kx, double ky, double S, const Params& p = Params{});

// Interlayer hopping matrix (6x6, spin-major): t_perp on yz and xz, zero on xy
Mat6 T_perp_mat(const Params& p = Params{});

// Staggered layer potential: diagonal 12x12 with -delta_V on layer 1, +delta_V on layer 2
Mat12 staggered_potential(const Params& p = Params{});

// Full bilayer Hamiltonian (12x12, layer-major block structure)
Mat12 bilayerHamiltonian(double kx, double ky, double S, const Params& p = Params{});

// Band energy + Hubbard DC correction (+6·U·S²) for a given Stoner parameter S
double calculate_hubbard_energy(double S, int grid_size, double T, double N_target,
                                const Params& p);

// Kanamori interaction parameters
struct KanamoriParams {
    double U       = 2.0;  // intraorbital Hubbard repulsion
    double U_prime = 1.0;  // interorbital Hubbard repulsion
    double J       = 0.3;  // Hund's coupling (spin-flip exchange + pair hopping)
};

// Kanamori mean-field Hamiltonian (12x12 bilayer)
// rho(a,b) = <c†_a c_b>, layer-major / spin-major / orbital-minor ordering
// Acts within each layer independently (on-site interaction)
Mat12 KanamoriMF(const Mat12& rho, const KanamoriParams& kp = KanamoriParams{});

// Save band structure along high-symmetry path to CSV (for plotting in Python)
void save_band_structure(const Mat12& rho, int n_points, const Params& p, const KanamoriParams& kp,
                         const std::string& filename, double mu);
                         
// Save DOS on a fine energy grid to CSV (for plotting in Python)
void save_dos(double S, int grid_size, double T, double N_target,
              const Params& p, const std::string& filename,
              int n_energy_points = 500, double sigma = 0.05);

// Save layer/spin/orbital-projected DOS for the converged Kanamori system.
// Each basis component contributes |v_nk(c)|^2, Gaussian-broadened on an energy grid.
// Columns: energy, mu, then the 12 projections (L{1,2}_{up,dn}_{yz,xz,xy}) and total.
void save_projected_dos(const Mat12& rho, int grid_size, double T, double N_target,
                        const Params& p, const KanamoriParams& kp,
                        const std::string& filename,
                        int n_energy_points = 500, double sigma = 0.05);