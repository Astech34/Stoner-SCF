# Documentation.md

## scf.cpp

### calculate_total_energy
Takes in a complete eigensystem which was computed over a grid (int grid_size)
```cpp
double calculate_total_energy(const Eigensystem& sys,double mu, double T)
```
Mathematically this is the energy weighted by the fermi dirac distribution. We sum over bands and over a k grid.

$$E = \sum_{n,\vec k} E_{nk} f(\varepsilon_{nk} - \mu)$$

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

Right now this is weighted by number of k points but I'm unsure if this is the correct weighting.