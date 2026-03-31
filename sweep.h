#pragma once
#include "hamiltonian.h"
#include "scf.h"

void run_U_sweep(double S0, double alpha, int grid, double T, double N_target,
                 double U_min, double U_max, int N_points, Params p);

// Sweep SOC strength lambda and compute MCA energy E([110]) - E([001]) at each value.
// Hopping and U are taken from p and held fixed.
void run_MCA_lam_sweep(double S0, double alpha, int grid, double T, double N_target,
                       double lam_min, double lam_max, int N_points, Params p);