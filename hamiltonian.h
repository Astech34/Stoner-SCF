#pragma once

#include <Eigen/Dense>
#include <complex>

using Mat6 = Eigen::Matrix<std::complex<double>, 6, 6>;
using cd   = std::complex<double>;

// Model parameters with sensible defaults
struct Params {
    double t1      = 1.0;
    double t_delta = 0.1;
    double t2      = 0.1;
    double lam     = 0.1;
    double U       = 2.0;
};

// Kinetic hopping term (k-dependent, real diagonal)
Mat6 H0(double kx, double ky, const Params& p = Params{});

// Spin-orbit coupling (k-independent, precompute once)
Mat6 SOC(double lam = 0.1);

// Hubbard mean-field term (k-independent, real diagonal)
Mat6 HubbardU(double S, double U = 2.0);

// Full single-layer Hamiltonian
Mat6 singleLayer(double kx, double ky, double S, const Params& p = Params{});