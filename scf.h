#pragma once

#include "hamiltonian.h"
#include <vector>

struct Eigensystem {
    std::vector<Eigen::Vector<double, 6>> evals; // (grid_size^2, 6)
    std::vector<Mat6>                    evecs; // (grid_size^2, 6, 6)
};

Eigensystem compute_eigensystem_grid(double S, int grid_size = 200, const Params& p = Params{});