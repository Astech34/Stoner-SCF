#pragma once

#include "hamiltonian.h"
#include "sweep.h"      // random_hermitian_perturbation
#include <string>
#include <array>
#include <algorithm>
#include <cctype>

// =====================================================================
// Density-matrix seeds for bilayer cobalt (12x12, 4 electrons per layer).
//
// Basis ordering:  layer  x  spin  x  orbital(yz, xz, xy)
//   index = layer*6 + spin*3 + orbital
//     layer  : 0, 1
//     spin   : 0 = up, 1 = down
//     orbital: 0 = yz, 1 = xz, 2 = xy
//
// Each seed places 4 electrons per layer (8 total, matching N_target ~ 8.2)
// and is applied identically to both layers. A random Hermitian perturbation
// of given Frobenius norm can be added to break residual symmetry and let the
// SCF explore neighbouring basins of attraction.
// =====================================================================

enum class SeedType {
    HighSpin,   // S = 1 per layer, t2g orbitally degenerate (up=[1,1,1], dn=[1/3,1/3,1/3])
    LowSpin,    // S = 0 per layer, t2g orbitally degenerate (up=dn=[2/3,2/3,2/3])
    OrbitalYZ,  // high spin, minority electron localised in yz (breaks t2g degeneracy)
    OrbitalXZ,  // high spin, minority electron localised in xz
    OrbitalXY   // high spin, minority electron localised in xy
};

// Build a bare seed density matrix (no perturbation).
inline Mat12 make_seed_bare(SeedType type) {
    Mat12 rho = Mat12::Zero();

    // Per-layer occupations: up[yz,xz,xy], dn[yz,xz,xy].
    std::array<double, 3> up{}, dn{};
    switch (type) {
        case SeedType::HighSpin:
            up = {1.0, 1.0, 1.0};
            dn = {1.0/3.0, 1.0/3.0, 1.0/3.0};
            break;
        case SeedType::LowSpin:
            up = {2.0/3.0, 2.0/3.0, 2.0/3.0};
            dn = {2.0/3.0, 2.0/3.0, 2.0/3.0};
            break;
        case SeedType::OrbitalYZ:
            up = {1.0, 1.0, 1.0};
            dn = {1.0, 0.0, 0.0};
            break;
        case SeedType::OrbitalXZ:
            up = {1.0, 1.0, 1.0};
            dn = {0.0, 1.0, 0.0};
            break;
        case SeedType::OrbitalXY:
            up = {1.0, 1.0, 1.0};
            dn = {0.0, 0.0, 1.0};
            break;
    }

    for (int layer = 0; layer < 2; ++layer) {
        const int base = layer * 6;
        for (int orb = 0; orb < 3; ++orb) {
            rho(base + 0 + orb, base + 0 + orb) = up[orb]; // spin up
            rho(base + 3 + orb, base + 3 + orb) = dn[orb]; // spin down
        }
    }
    return rho;
}

// Build a seed and add a random Hermitian perturbation of Frobenius norm
// `strength` (seeded by `rng_seed`). strength <= 0 returns the bare seed.
inline Mat12 make_seed(SeedType type, double strength = 0.0, unsigned rng_seed = 0) {
    Mat12 rho = make_seed_bare(type);
    if (strength > 0.0)
        rho += random_hermitian_perturbation(strength, rng_seed);
    return rho;
}

// Parse a seed name (case-insensitive). Accepted: high_spin/highspin/high,
// low_spin/lowspin/low, yz, xz, xy. Throws on an unknown name.
inline SeedType parse_seed_type(std::string name) {
    std::transform(name.begin(), name.end(), name.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (name == "high_spin" || name == "highspin" || name == "high") return SeedType::HighSpin;
    if (name == "low_spin"  || name == "lowspin"  || name == "low")  return SeedType::LowSpin;
    if (name == "yz")                                                return SeedType::OrbitalYZ;
    if (name == "xz")                                                return SeedType::OrbitalXZ;
    if (name == "xy")                                                return SeedType::OrbitalXY;
    throw std::invalid_argument("parse_seed_type: unknown seed '" + name + "'");
}

// Convenience: build a seed directly from its name.
inline Mat12 make_seed(const std::string& name, double strength = 0.0, unsigned rng_seed = 0) {
    return make_seed(parse_seed_type(name), strength, rng_seed);
}
