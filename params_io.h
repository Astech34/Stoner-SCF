#pragma once

#include "hamiltonian.h"
#include <string>
#include <iosfwd>

struct SCFParams {
    double S0       = 0.3;
    double alpha    = 0.2;
    double T        = 0.05;
    double N_target = 10.0;
    int    grid     = 50;
};

struct AllParams {
    Params         p;
    KanamoriParams kp;
    SCFParams      scf;
    std::string    param_file;
    std::string    rho_out_file;
    // Density-matrix seed (see seeds.h). Empty/"none" disables seeding.
    std::string    seed;
    double         seed_strength = 0.0;   // Frobenius norm of random Hermitian perturbation
    unsigned       seed_rng      = 0;     // RNG seed for that perturbation
};

// Load parameters from a key=value file. Lines starting with '#' are ignored.
// KanamoriParams.U defaults to Params.U; U_prime = U - 2*J unless set explicitly.
AllParams load_params(const std::string& filename);

// Write a human-readable report of all loaded parameters to the given stream.
void write_params(std::ostream& os, const AllParams& ap);
