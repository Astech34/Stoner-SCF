# Parallelization Notes

## Current State

Only one parallel region exists: the k-point loop in `compute_eigensystem_grid()` (`scf.cpp:139`) uses `#pragma omp parallel for schedule(static)`. Everything else — the SCF loop, parameter sweeps, Brent root-finding, NTotal summation, and magnetization computation — is sequential.

---

## Easy OpenMP Wins

Two loops are currently serial but fully independent and could take an `omp parallel for reduction(+:...)`:

- **`NTotal()` loop** (`scf.cpp:81-92`): sums Fermi-Dirac weights over 40,000 k-points × 12 bands. Called 10–20× per SCF step inside Brent's method.
- **Magnetization computation** (`scf.cpp:217-240`): same structure — spin projections over all (k, band) pairs, accumulated into a scalar.

---

## MPI Opportunities

### Parameter sweep distribution (best target)

The U-sweep (`sweep.cpp:21`) and λ-sweep (`sweep.cpp:71`) are embarrassingly parallel — each parameter point runs a completely independent SCF calculation with no mid-flight communication. The pattern is:

1. `MPI_Scatter` to distribute parameter values across ranks
2. Each rank runs its own `runSelfCalc()` independently
3. `MPI_Gather` to collect scalar results (energy, magnetization) on master

Expected speedup: linear in number of MPI ranks, up to the number of sweep points (50 for U-sweep, 20 for λ-sweep).

### MCA magnetization directions

The two SCF runs in `run_MCA_lam_sweep()` (`sweep.cpp:78-83`) — one for [001], one for [110] — are independent by construction. Two ranks can run them simultaneously and reduce with a single subtraction at the end.

### Hybrid MPI + OpenMP

The natural architecture is MPI across nodes (distributing sweep points) with OpenMP within each rank (parallelizing the k-point loop). This is the standard layout used in DFT codes like VASP and Quantum ESPRESSO.

---

## What Cannot Be Parallelized

- **SCF iterations** (`scf.cpp:260-279`): S_{n+1} depends on S_n — inherently sequential. Cannot be parallelized; only accelerated (see below).
- **Brent root-finding** (`scf.cpp:16-76`): each bracket step depends on the previous one, same constraint.

---

## Anderson Acceleration (non-MPI speedup)

Replacing the linear mixing in the SCF loop with Anderson acceleration (DIIS) typically reduces iteration count from ~20–50 down to ~5–10 — a 3–5× speedup without any additional parallelism. This compounds multiplicatively with MPI and OpenMP gains, so it is worth considering alongside the parallelization work.

---

## Computational Context

For reference, the main bottleneck is `compute_eigensystem_grid()`: diagonalizing 40,000 independent 12×12 complex Hermitian matrices per SCF step, storing ~100 MB of eigenvalues and eigenvectors. This is the operation that currently benefits from OpenMP and would benefit most from GPU acceleration (batched cuSolver) if that path were pursued.