# Documentation.md

## scf.cpp

### calculate_total_energy
Takes in a complete eigensystem which was computed over a grid (int grid_size)
```cpp
double calculate_total_energy(const Eigensystem& sys,double mu, double T)
```
Mathematically this is the energy weighted by the fermi dirac distribution. We sum over bands and over a k grid.

$$E = \sum_{n,\vec k} \varepsilon_{nk} f(\varepsilon_{nk} - \mu)$$

Where the energy vals are contained in the ```Eigensystem sys```

sys.evals is {V1, V2, V3.... Vk} where Vk is a 1D vector of length 12 representing the band and Vk is evaluated at a specific (kx, ky)

```cpp
double calculate_total_energy(const Eigensystem& sys, double mu, double T) {
    const int N_k     = static_cast<int>(sys.evals.size());
    const int N_bands = 12;

    double E_total = 0.0;
    for (const auto& ev : sys.evals) {
        for (int b = 0; b < N_bands; b++) {
            const double x = std::clamp((ev[b] - mu) / T, -500.0, 500.0);
            const double f = 1.0 / (std::exp(x) + 1.0);
            E_total += ev[b] * f;
        }
    }

    return E_total / static_cast<double>(N_k);
}
```

Right now this is weighted by number of k points but I'm unsure if this is the correct weighting. This isn't used anywhere else just to compare so the weight does not matter for MCA comparison.

### runSelfCal
This is the self consistency calc. We have a current guess for our order parameter S0. Then with initial conditions we calculate S_new by plugging our S into the model and diagonalizing then we get an S_calc. If this difference is large enough we start again by mixing these two. With $\alpha = 0$ we get the standard average, but with allowing this parameter we can control the mixing which improves convergance.
$$S_{current} = \alpha \cdot S_{calc} + (1 - \alpha) \cdot S_{current}$$
```cpp
CalcResult runSelfCalc(double S0, double alpha, int grid_size,
                                       double T, double N_target, const Params& p) {
    constexpr int    max_iter = 5000;
    constexpr double tol      = 1e-5;

    double S_current = S0;
    double mu_current = std::numeric_limits<double>::quiet_NaN();
    double E_current = std::numeric_limits<double>::quiet_NaN();
    bool converged = false;

    for (int i = 0; i < max_iter; i++) {
        auto result = calculateS(S_current, grid_size, T, N_target, p);
        double S_calc = result.S_new;
        mu_current = result.mu;
        E_current  = result.E_total;

        const double diff = std::abs(S_calc - S_current);

        std::cout << "Iteration " << i
                  << ", S = "    << std::fixed << std::setprecision(6) << S_current
                  << ", Difference: " << diff << "\n";

        if (diff < tol) {
            std::cout << "Self-consistency reached!\n";
            converged = true;
            break;
        }

        S_current = alpha * S_calc + (1.0 - alpha) * S_current;
    }

    if (!converged)
        throw std::runtime_error("runSelfCalc: failed to converge within max_iter");

    std::cout << "Final S: " << S_current << ", mu: " << mu_current << ", E_total: " << E_current << "\n";
    return {S_current, mu_current, E_current};
}
```