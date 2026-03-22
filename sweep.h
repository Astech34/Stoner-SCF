#pragma once
#include "hamiltonian.h"
#include "scf.h"

void run_U_sweep(double S0, double alpha, int grid, double T, double N_target,
                 double U_min, double U_max, int N_points, Params p);