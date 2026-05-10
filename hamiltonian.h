#pragma once

#include <Eigen/Dense>
#include <complex>
#include <utility>

using Mat6  = Eigen::Matrix<std::complex<double>, 6, 6>;
using Mat12 = Eigen::Matrix<std::complex<double>, 12, 12>;
using cd    = std::complex<double>;

// Model parameters with sensible defaults
struct Params {
    double t1      = 1.0;
    double t_delta = 0.1;
    double t2      = 0.1;
    double lam     = 0.1;
    double U       = 2.0;
    double theta   = 0.0;  // polar angle of magnetisation direction (0 = z-axis)
    double phi     = 0.0;  // azimuthal angle of magnetisation direction
    double t_perp    = 0.0;  // interlayer hopping for yz and xz orbitals
    double t_perp_xy = 0.0;  // interlayer hopping for xy orbital (zero by symmetry for ideal stacking)
    double delta_cf  = 0.0;  // tetragonal crystal field: raises xy above yz/xz (positive = xy higher)
    double delta_V   = 0.0;  // staggered layer potential: layer 1 shifts by -delta_V, layer 2 by +delta_V
};

// Kinetic hopping term (k-dependent, real diagonal)
Mat6 H0(double kx, double ky, const Params& p = Params{});

// Spin-orbit coupling (k-independent, precompute once)
Mat6 SOC(double lam = 0.1);

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

// Save band structure along high-symmetry path to CSV (for plotting in Python)
void save_band_structure(double S, int n_points, const Params& p,
                         const std::string& filename, double mu);
                         
// Save DOS on a fine energy grid to CSV (for plotting in Python)
void save_dos(double S, int grid_size, double T, double N_target,
              const Params& p, const std::string& filename,
              int n_energy_points = 500, double sigma = 0.05);

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