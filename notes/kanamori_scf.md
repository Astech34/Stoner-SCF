# Kanamori SCF Loop

## Overview

The Kanamori SCF converges the full 12×12 single-particle density matrix `ρ` rather than just a scalar magnetisation `S` as in the Stoner loop. This captures orbital polarisation, inter-orbital coherences, and Hund's physics that the scalar Stoner order parameter cannot represent.

---

## Order Parameter

The order parameter is the density matrix:

$$\rho_{ab} = \langle c^\dagger_a c_b \rangle = \frac{1}{N_k} \sum_{\mathbf{k},n} f(\varepsilon_{n\mathbf{k}} - \mu)\, v_{n\mathbf{k},a}^* \, v_{n\mathbf{k},b}$$

where `a, b` run over the 12 basis states (layer-major → spin-major → orbital-minor). The diagonal elements are orbital occupations; off-diagonal elements are inter-orbital and inter-spin coherences.

---

## Self-Consistent Loop (`runKanamoriSCF`)

```
rho ← rho0  (initial guess, typically from converged Stoner SCF)

for each iteration:
    1. Build H(k) from rho  →  compute_eigensystem_kanamori
    2. Find μ via Brent's method on N(μ) = N_target
    3. Compute rho_new from eigensystem + Fermi-Dirac weights
    4. Check convergence: ||rho_new - rho||_F < 1e-6
       → if converged: compute total energy and return
    5. Mix:  rho ← α·rho_new + (1-α)·rho
```

Convergence tolerance is `1e-6` on the Frobenius norm; maximum iterations is 5000.

---

## Hamiltonian Construction (`compute_eigensystem_kanamori`)

At each iteration the 12×12 Hamiltonian is assembled as:

$$H(\mathbf{k}) = H_\text{Kanamori}[\rho] + H_\text{SOC} + T_\perp + V_\text{stag} + H_0(\mathbf{k})$$

The k-independent pieces — `KanamoriMF(rho)`, SOC, `T_perp`, and the staggered potential — are precomputed once into a single `H_kfree` matrix before the OpenMP k-loop. Inside the loop only the kinetic term `H0(k)` is added, with `delta_cf1` for layer 1 and `delta_cf2` for layer 2:

```cpp
Mat12 H = H_kfree;
H.block<6,6>(0,0) += H0(kx, ky, p, p.delta_cf1);  // layer 1
H.block<6,6>(6,6) += H0(kx, ky, p, p.delta_cf2);  // layer 2
```

The k-loop is parallelised with `#pragma omp parallel for`.

---

## Linear Mixing (Relaxation)

Simple linear mixing is used to stabilise convergence:

$$\rho \leftarrow \alpha \cdot \rho_\text{new} + (1-\alpha) \cdot \rho$$

- `α = 1` is a direct substitution (fast but can oscillate or diverge)
- `α < 1` damps oscillations at the cost of slower convergence
- Typical value used: `α = 0.2`

The same mixing scheme is used in the Stoner loop for `S`. For the density matrix the same principle applies but the update is now a full 12×12 matrix mix, so convergence is generally slower than the scalar case.

---

## Total Energy and Double-Counting Correction

On convergence the total energy is:

$$E_\text{total} = E_\text{band} + \Delta E_\text{DC}$$

**Band energy** — weighted sum of eigenvalues:

$$E_\text{band} = \frac{1}{N_k} \sum_{\mathbf{k},n} \varepsilon_{n\mathbf{k}} \, f(\varepsilon_{n\mathbf{k}} - \mu)$$

**Double-counting correction** — subtracts the interaction energy that was already counted inside the eigenvalues. Computed per layer from `kanamori_dc_layer`:

$$\Delta E_\text{DC} = -U\sum_m \langle n_{m\uparrow}\rangle\langle n_{m\downarrow}\rangle
- U'\sum_{m\neq m'}\langle n_{m\uparrow}\rangle\langle n_{m'\downarrow}\rangle
- (U'-J)\sum_{m<m',\sigma}\langle n_{m\sigma}\rangle\langle n_{m'\sigma}\rangle
+ J\sum_{m\neq m'}\left[\langle d^\dagger_{m\uparrow}d_{m'\uparrow}\rangle\langle d^\dagger_{m'\downarrow}d_{m\downarrow}\rangle - \ldots\right]$$

The correction is evaluated using the **converged** `rho` (not `rho_new`) and summed over both layers.

---

## Initialisation

The Kanamori loop is seeded from the converged Stoner SCF result. The Stoner eigensystem at its converged `S` is used to build a `rho0` via `compute_density_matrix`. Small symmetry-breaking perturbations are added to `rho0` in `main.cpp` to allow the Kanamori loop to find states with broken orbital or layer symmetry that the symmetric starting point would not reach.
