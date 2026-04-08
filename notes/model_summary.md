# Stoner SCF Model Summary

## Overview

A self-consistent mean-field model for itinerant ferromagnetism in a bilayer t2g system on a 2D square lattice (AA stacking). Each layer contains three t2g orbitals (yz, xz, xy) and two spins, giving a **12-band bilayer model** (3 orbitals × 2 spins × 2 layers). The primary application is computing the **magnetocrystalline anisotropy (MCA) energy** as a function of model parameters.

All single-layer matrices are 6×6 in **spin-major ordering**: rows/cols 0–2 = spin-up (yz, xz, xy), rows/cols 3–5 = spin-down (yz, xz, xy). The full bilayer Hamiltonian is 12×12 in **layer-major ordering**: rows/cols 0–5 = layer 1, rows/cols 6–11 = layer 2.

---

## Hamiltonian

### Single-layer Hamiltonian

$$H_\text{single}(\mathbf{k}) = H_0(\mathbf{k}) + H_\text{SOC} + H_\text{Hubbard}$$

#### Kinetic term — H₀(k)

Tight-binding hopping, diagonal in spin and orbital. The orbital shapes dictate the hopping anisotropy:

- **yz**: lobes point along $y$ → hops easily along $y$ (amplitude $t_1$), weakly along $x$ (amplitude $t_\delta$)
- **xz**: lobes point along $x$ → hops easily along $x$ (amplitude $t_1$), weakly along $y$ (amplitude $t_\delta$)
- **xy**: lobes point diagonally → nearest-neighbour hopping $t_1$ plus next-nearest-neighbour $t_2$ along diagonals

$$\varepsilon_{yz}(\mathbf{k}) = -2t_\delta \cos k_x - 2t_1 \cos k_y$$
$$\varepsilon_{xz}(\mathbf{k}) = -2t_1 \cos k_x - 2t_\delta \cos k_y$$
$$\varepsilon_{xy}(\mathbf{k}) = -2t_1(\cos k_x + \cos k_y) - 4t_2 \cos k_x \cos k_y$$

`t_delta` breaks the yz/xz degeneracy along each axis. `t2` controls the curvature and bandwidth of the xy band.

#### Spin-orbit coupling — H_SOC

Atomic SOC in the t2g basis, k-independent:

$$H_\text{SOC} = \lambda \, \mathbf{L} \cdot \mathbf{S} = \lambda (L_x \otimes s_x + L_y \otimes s_y + L_z \otimes s_z)$$

For λ = 1, the full matrix is:

$$H_\text{SOC} = \begin{pmatrix} 0 & \frac{i}{2} & 0 & 0 & 0 & -\frac{1}{2} \\ -\frac{i}{2} & 0 & 0 & 0 & 0 & \frac{i}{2} \\ 0 & 0 & 0 & \frac{1}{2} & -\frac{i}{2} & 0 \\ 0 & 0 & \frac{1}{2} & 0 & -\frac{i}{2} & 0 \\ 0 & 0 & \frac{i}{2} & \frac{i}{2} & 0 & 0 \\ -\frac{1}{2} & -\frac{i}{2} & 0 & 0 & 0 & 0 \end{pmatrix}$$

SOC is the origin of MCA — it couples the spin degree of freedom to the orbital/lattice structure, making the total energy depend on the magnetisation direction.

#### Hubbard mean-field term — H_Hubbard

Stoner exchange splitting for a magnetisation of magnitude $S$ pointing along a fixed direction $\hat{\mathbf{n}} = (\sin\theta\cos\phi,\, \sin\theta\sin\phi,\, \cos\theta)$:

$$H_\text{Hubbard} = -US\,(\hat{\mathbf{n}}\cdot\boldsymbol{\sigma})\otimes I_3 = -US \begin{pmatrix} n_z I_3 & (n_x - in_y)I_3 \\ (n_x + in_y)I_3 & -n_z I_3 \end{pmatrix}$$

$S$ is the scalar magnetisation magnitude and $U$ is the Hubbard parameter. The direction $\hat{\mathbf{n}}$ is a **fixed input** (set via `theta` and `phi` in `Params`); only $S$ is updated self-consistently. For $\hat{\mathbf{n}} = \hat{z}$ ($\theta=0$) this reduces to the diagonal form $US\cdot\text{diag}(-1,-1,-1,+1,+1,+1)$. The off-diagonal spin blocks are non-zero whenever $\hat{\mathbf{n}}$ has an $x$ or $y$ component.

### Bilayer Hamiltonian

For AA stacking (atoms directly above each other), the full 12×12 Hamiltonian is:

$$H(\mathbf{k}) = \begin{pmatrix} H_\text{single}(\mathbf{k}) & T_\perp \\ T_\perp & H_\text{single}(\mathbf{k}) \end{pmatrix}$$

The two diagonal blocks are identical (AA stacking, same crystal environment in each layer). $T_\perp$ is k-independent (direct vertical hopping).

#### Interlayer hopping — T_⊥

$$T_\perp = \begin{pmatrix} t_\perp & 0 & 0 \\ 0 & t_\perp & 0 \\ 0 & 0 & 0 \end{pmatrix} \otimes I_2$$

Only yz and xz orbitals hop between layers — their lobes extend along $z$ toward the neighbouring layer. The xy orbital has lobes lying entirely in-plane so $t_\perp^{xy} = 0$. Spin is conserved during interlayer hopping ($\otimes I_2$). The bilayer Hamiltonian splits each single-layer band into a bonding/antibonding pair separated by $\pm t_\perp$.

---

## Self-Consistent Loop

The SCF loop finds the magnetisation magnitude $S$ that is consistent with the electronic structure it generates.

1. **Guess** an initial magnetisation `S₀`
2. **Diagonalise** $H(\mathbf{k})$ over a uniform $N \times N$ k-grid across the Brillouin zone, storing all 12 eigenvalues and eigenvectors
3. **Find μ** by solving $N_\text{total}(\mu) = N_\text{target}$ using Brent's method on the Fermi-Dirac sum
4. **Update S** by projecting the density matrix onto $\hat{\mathbf{n}}$:

$$S_\text{new} = \frac{1}{2N_k}\sum_{\mathbf{k},n} f(\varepsilon_{n\mathbf{k}}-\mu)\left[n_z\!\left(|\psi^\uparrow_{n\mathbf{k}}|^2 - |\psi^\downarrow_{n\mathbf{k}}|^2\right) + 2n_x\,\text{Re}[\psi^{\uparrow\dagger}_{n\mathbf{k}}\psi^\downarrow_{n\mathbf{k}}] + 2n_y\,\text{Im}[\psi^{\uparrow\dagger}_{n\mathbf{k}}\psi^\downarrow_{n\mathbf{k}}]\right]$$

For the bilayer, $\psi^\uparrow_{n\mathbf{k}}$ is a 6-component vector gathering spin-up components from both layers (indices 0–2 from layer 1, 6–8 from layer 2), and likewise for $\psi^\downarrow_{n\mathbf{k}}$. The cross term $\psi^{\uparrow\dagger}\psi^\downarrow$ captures off-diagonal spin coherence needed when $\hat{\mathbf{n}}$ has $x$ or $y$ components.

5. **Mix** to stabilise convergence: $S \leftarrow \alpha S_\text{new} + (1-\alpha) S$
6. **Repeat** until $|S_\text{new} - S| < 10^{-5}$

A converged $S \neq 0$ indicates a ferromagnetic ground state; $S = 0$ is paramagnetic. For the bilayer, saturation is $S_\text{max} = 1.0$ (vs 0.5 for single layer) since $S$ sums the magnetisation over both layers.

---

## Magnetocrystalline Anisotropy (MCA)

MCA arises from the interplay of SOC and the exchange splitting. SOC couples spin to the orbital/lattice, making the total energy depend on which direction $\hat{\mathbf{n}}$ the magnetisation points. The MCA energy is:

$$E_\text{MCA} = E(\hat{\mathbf{n}}_1) - E(\hat{\mathbf{n}}_2)$$

where each $E(\hat{\mathbf{n}})$ is the converged SCF total energy:

$$E = \frac{1}{N_k} \sum_{\mathbf{k},n} \varepsilon_{n\mathbf{k}} \, f(\varepsilon_{n\mathbf{k}} - \mu)$$

**Key physics:**
- MCA requires both SOC ($\lambda \neq 0$) and magnetic order ($S \neq 0$) simultaneously — SOC breaks rotational symmetry, exchange splitting makes the system care about direction
- MCA scales as $\lambda^2$ (second-order perturbation theory in SOC)
- The square lattice has C4v symmetry, which constrains and partially cancels MCA contributions — lower symmetry gives larger MCA
- Saturation ($S \to S_\text{max}$) suppresses MCA by reducing the number of states near $E_F$ available for SOC mixing

**Procedure:**
1. Run SCF to convergence with $\hat{\mathbf{n}} = \hat{\mathbf{n}}_1$, record $E_1$
2. Run SCF to convergence with $\hat{\mathbf{n}} = \hat{\mathbf{n}}_2$, record $E_2$
3. $E_\text{MCA} = E_1 - E_2$

Negative $E_\text{MCA}$ means $\hat{\mathbf{n}}_1$ is the easy axis (lower energy).

---

## Default Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `t1`      | 1.0   | Nearest-neighbour hopping (sets energy scale) |
| `t_delta` | 0.1   | Anisotropic hopping (breaks yz/xz symmetry along each axis) |
| `t2`      | 0.1   | Next-nearest-neighbour hopping for xy orbital |
| `lam`     | 0.1   | SOC strength (~upper end for 3d metals with t1~0.5 eV) |
| `U`       | 2.0   | Hubbard interaction |
| `t_perp`  | 0.0   | Interlayer hopping for yz and xz (0 = decoupled single layer) |
| `theta`   | 0.0   | Polar angle of magnetisation direction (0 = z-axis) |
| `phi`     | 0.0   | Azimuthal angle of magnetisation direction |
| `T`       | 0.05  | Temperature (in units of `t1`) |
| `N_target`| 10.0  | Electrons per bilayer unit cell (5 per layer) |
