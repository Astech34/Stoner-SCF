#pragma once

#include "hamiltonian.h"
#include <vector>
#include <functional>
#include <cassert>
#include <utility>

struct Eigensystem {
    std::vector<Eigen::Vector<double, 6>> evals; // (grid_size^2, 6)
    std::vector<Mat6>                    evecs; // (grid_size^2, 6, 6)
};

// Brent's method root finder — f must be continuous and f(a)*f(b) < 0
double brent(std::function<double(double)> f, double a, double b,
             double tol = 1e-10, int max_iter = 100);

// Find chemical potential via Brent on the total electron count
double find_mu(const Eigensystem& sys, int grid_size = 200,
               double T = 0.05, double N_target = 5.0);

Eigensystem compute_eigensystem_grid(double S, int grid_size = 200, const Params& p = Params{});

// Calculate Stoner parameter S for a given guess, and the corresponding chemical potential
std::pair<double, double> calculateS(double S, int grid_size, double T, double N_target, const Params& p);

// Run the self-consistent loop with linear mixing
std::pair<double, double> runSelfCalc(double S0, double alpha, int grid_size, double T, double N_target, const Params& p);