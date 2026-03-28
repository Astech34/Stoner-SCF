# Stoner SCF Model Summary

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

Stoner exchange splitting, diagonal in orbital, opposite sign for each spin:

$$H_\text{Hubbard} = US \cdot \text{diag}(-1,-1,-1,+1,+1,+1)$$

where `S = (n↑ - n↓) / 2` is the magnetisation and `U` is the Hubbard parameter.

---

## Self-Consistent Loop

1. **Guess** an initial magnetisation `S₀`
2. **Diagonalise** H(k) over a uniform k-grid (N × N) across the Brillouin zone
3. **Find μ** by solving `N_total(μ) = N_target` using Brent's method on the Fermi-Dirac sum
4. **Update S** from the spin-resolved orbital occupancies:

$$S_\text{new} = \frac{n_\uparrow - n_\downarrow}{2}, \qquad n_{\uparrow/\downarrow} = \frac{1}{N_k} \sum_{\mathbf{k},n} |{\langle \uparrow/\downarrow | \psi_{n\mathbf{k}} \rangle}|^2 \, f(\varepsilon_{n\mathbf{k}} - \mu)$$

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
| `T`       | 0.05  | Temperature (in units of `t1`) |
| `N_target`| 5.0   | Electrons per unit cell |
