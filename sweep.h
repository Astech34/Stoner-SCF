#pragma once
#include "hamiltonian.h"
#include "scf.h"

void run_U_sweep(double S0, double alpha, int grid, double T, double N_target,
                 double U_min, double U_max, int N_points, Params p);

// Sweep SOC strength lambda and compute MCA energy E([110]) - E([001]) at each value
// using Kanamori SCF. Stoner SCF seeds each direction independently.
void run_MCA_lam_sweep(double S0, double alpha, int grid, double T, double N_target,
                       double lam_min, double lam_max, int N_points, Params p, KanamoriParams kp);

// Sweep staggered layer potential delta_V and record layer occupation difference.
// Stoner bootstrap runs once (delta_V does not enter that path); Kanamori SCF
// is re-run at each delta_V value. Writes delta_V, n1, n2, dn = n1-n2 to CSV.
void run_delta_V_sweep(double S0, double alpha, int grid, double T, double N_target,
                       double dV_min, double dV_max, int N_points,
                       Params p, KanamoriParams kp);