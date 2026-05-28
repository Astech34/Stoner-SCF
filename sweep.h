#pragma once
#include "hamiltonian.h"
#include "scf.h"

void run_U_sweep(double S0, double alpha, int grid, double T, double N_target,
                 double U_min, double U_max, int N_points, Params p);

struct MCAResult {
    KanamoriResult res_001;
    KanamoriResult res_110;
    double E_MCA;
};

// Compute MCA energy E([110]) - E([001]).
// Converges Kanamori at [001] first (with symmetry-breaking seed delta),
// then seeds [110] from the converged [001] rho to ensure both axes are in the same phase.
MCAResult compute_MCA(double S0, double alpha, int grid, double T, double N_target,
                      double delta, Params p, KanamoriParams kp);

// Sweep SOC strength lambda and compute MCA energy E([110]) - E([001]) at each value.
void run_MCA_lam_sweep(double S0, double alpha, int grid, double T, double N_target,
                       double lam_min, double lam_max, int N_points, double delta,
                       Params p, KanamoriParams kp);

// Sweep staggered layer potential delta_V and record layer occupation difference.
// Stoner bootstrap runs once (delta_V does not enter that path); Kanamori SCF
// is re-run at each delta_V value. Writes delta_V, n1, n2, dn = n1-n2 to CSV.
void run_delta_V_sweep(double S0, double alpha, int grid, double T, double N_target,
                       double dV_min, double dV_max, int N_points,
                       Params p, KanamoriParams kp);

// Run the Kanamori SCF seeded from a Stoner bootstrap + random Hermitian perturbation.
// The Stoner SCF runs at theta=phi=0 to produce a physical rho0, then a random
// Hermitian matrix of Frobenius norm epsilon is added to break symmetry.
// Different seeds explore different basins of attraction.
KanamoriResult runKanamoriSCF_random(unsigned seed, double S0, double alpha, int grid_size,
                                     double T, double N_target,
                                     const Params& p, const KanamoriParams& kp,
                                     double epsilon = 0.05);