# Stoner SCF Model Summary

Before Second Layer

## Overview

A self-consistent mean-field model for itinerant ferromagnetism in a t2g system on a 2D square lattice. The three t2g orbitals are yz, xz, and xy, giving a 6-band model (3 orbitals × 2 spins).

All matrices are 6×6 in **spin-major ordering**: rows/cols 0–2 = spin-up (yz, xz, xy), rows/cols 3–5 = spin-down (yz, xz, xy).

---

## Hamiltonian

$$H(\mathbf{k}) = H_0(\mathbf{k}) + H_\text{SOC} + H_\text{Hubbard}$$

### Kinetic term — H₀(k)

Tight-binding hopping, diagonal in spin and orbital:

$$\varepsilon_{yz}(\mathbf{k}) = -2t_\delta \cos k_x - 2t_1 \cos k_y$$
$$\varepsilon_{xz}(\mathbf{k}) = -2t_1 \cos k_x - 2t_\delta \cos k_y$$
$$\varepsilon_{xy}(\mathbf{k}) = -2t_1(\cos k_x + \cos k_y) - 4t_2 \cos k_x \cos k_y$$

Parameters: `t1` (nearest-neighbour), `t_delta` (anisotropic hopping), `t2` (next-nearest-neighbour).

### Spin-orbit coupling — H_SOC

Atomic SOC in the t2g basis, k-independent:

$$H_\text{SOC} = \lambda \, \mathbf{L} \cdot \mathbf{S} = \lambda (L_x \otimes s_x + L_y \otimes s_y + L_z \otimes s_z)$$

For λ = 1, the full matrix is:

$$H_\text{SOC} = \begin{pmatrix} 0 & \frac{i}{2} & 0 & 0 & 0 & -\frac{1}{2} \\ -\frac{i}{2} & 0 & 0 & 0 & 0 & \frac{i}{2} \\ 0 & 0 & 0 & \frac{1}{2} & -\frac{i}{2} & 0 \\ 0 & 0 & \frac{1}{2} & 0 & -\frac{i}{2} & 0 \\ 0 & 0 & \frac{i}{2} & \frac{i}{2} & 0 & 0 \\ -\frac{1}{2} & -\frac{i}{2} & 0 & 0 & 0 & 0 \end{pmatrix}$$

### Hubbard mean-field term — H_Hubbard

Stoner exchange splitting for a magnetisation of magnitude $S$ pointing along a fixed direction $\hat{\mathbf{n}} = (\sin\theta\cos\phi,\, \sin\theta\sin\phi,\, \cos\theta)$:

$$H_\text{Hubbard} = -US\,(\hat{\mathbf{n}}\cdot\boldsymbol{\sigma})\otimes I_3 = -US \begin{pmatrix} n_z I_3 & (n_x - in_y)I_3 \\ (n_x + in_y)I_3 & -n_z I_3 \end{pmatrix}$$

where `S` is the scalar magnetisation magnitude and `U` is the Hubbard parameter. The direction $\hat{\mathbf{n}}$ is a fixed input (set via `theta` and `phi` in `Params`); only $S$ is updated self-consistently. For $\hat{\mathbf{n}} = \hat{z}$ ($\theta=0$) this reduces to the diagonal form $US\cdot\text{diag}(-1,-1,-1,+1,+1,+1)$. The off-diagonal spin blocks are non-zero whenever $\hat{\mathbf{n}}$ has an $x$ or $y$ component.

---

## Self-Consistent Loop

1. **Guess** an initial magnetisation `S₀`
2. **Diagonalise** H(k) over a uniform k-grid (N × N) across the Brillouin zone
3. **Find μ** by solving `N_total(μ) = N_target` using Brent's method on the Fermi-Dirac sum
4. **Update S** by projecting the density matrix onto $\hat{\mathbf{n}}$:

$$S_\text{new} = \frac{1}{2N_k}\sum_{\mathbf{k},n} f(\varepsilon_{n\mathbf{k}}-\mu)\left[n_z\!\left(|\psi^\uparrow_{n\mathbf{k}}|^2 - |\psi^\downarrow_{n\mathbf{k}}|^2\right) + 2n_x\,\text{Re}[\psi^{\uparrow\dagger}_{n\mathbf{k}}\psi^\downarrow_{n\mathbf{k}}] + 2n_y\,\text{Im}[\psi^{\uparrow\dagger}_{n\mathbf{k}}\psi^\downarrow_{n\mathbf{k}}]\right]$$

where $\psi^\uparrow_{n\mathbf{k}}$ and $\psi^\downarrow_{n\mathbf{k}}$ are the spin-up and spin-down 3-component blocks of the eigenvector. The cross term $\psi^{\uparrow\dagger}\psi^\downarrow$ captures off-diagonal spin coherence needed for $\hat{\mathbf{n}}$ with $x$ or $y$ components. Reduces to $(n_\uparrow - n_\downarrow)/2$ for $\hat{\mathbf{n}}=\hat{z}$.

5. **Mix** to stabilise convergence: $S \leftarrow \alpha S_\text{new} + (1-\alpha) S$
6. **Repeat** until $|S_\text{new} - S| < 10^{-5}$

A converged `S ≠ 0` indicates a ferromagnetic ground state; `S = 0` is paramagnetic.

---

## Default Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `t1`      | 1.0   | Nearest-neighbour hopping |
| `t_delta` | 0.1   | Anisotropic hopping |
| `t2`      | 0.1   | Next-nearest-neighbour hopping |
| `lam`     | 0.1   | SOC strength |
| `U`       | 2.0   | Hubbard interaction |
| `theta`   | 0.0   | Polar angle of magnetisation direction (0 = z-axis) |
| `phi`     | 0.0   | Azimuthal angle of magnetisation direction |
| `T`       | 0.05  | Temperature (in units of `t1`) |
| `N_target`| 5.0   | Electrons per unit cell |
