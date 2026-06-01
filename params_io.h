#pragma once

#include "hamiltonian.h"
#include <string>

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
};

// Load parameters from a key=value file. Lines starting with '#' are ignored.
// KanamoriParams.U defaults to Params.U; U_prime = U - 2*J unless set explicitly.
AllParams load_params(const std::string& filename);
