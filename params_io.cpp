#include "params_io.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>

AllParams load_params(const std::string& filename) {
    AllParams ap;
    bool U_prime_explicit = false;
    bool kp_U_explicit    = false;

    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("load_params: cannot open '" + filename + "'");

    auto trim = [](std::string& s) {
        const auto start = s.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) { s.clear(); return; }
        s = s.substr(start, s.find_last_not_of(" \t\r\n") - start + 1);
    };

    std::string line;
    while (std::getline(file, line)) {
        auto hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);

        auto eq = line.find('=');
        if (eq == std::string::npos) continue;

        std::string key = line.substr(0, eq);
        std::string val = line.substr(eq + 1);
        trim(key); trim(val);
        if (key.empty() || val.empty()) continue;

        if (key == "param_file") {
            ap.param_file = val;
            continue;
        }
        if (key == "rho_out_file") {
            ap.rho_out_file = val;
            continue;
        }
        if (key == "seed") {
            ap.seed = val;
            continue;
        }

        const double v = std::stod(val);

        // Params
        if      (key == "tpi")       ap.p.tpi       = v;
        else if (key == "tdelta")    ap.p.tdelta    = v;
        else if (key == "t2xy")      ap.p.t2xy      = v;
        else if (key == "t2yzxz")    ap.p.t2yzxz    = v;
        else if (key == "tg")        ap.p.tg        = v;
        else if (key == "lam")       ap.p.lam       = v;
        else if (key == "U")         { ap.p.U = v; if (!kp_U_explicit) ap.kp.U = v; }
        else if (key == "theta")     ap.p.theta     = v;
        else if (key == "phi")       ap.p.phi       = v;
        else if (key == "t_perp")    ap.p.t_perp    = v;
        else if (key == "t_perp_xy") ap.p.t_perp_xy = v;
        else if (key == "delta_cf1") ap.p.delta_cf1 = v;
        else if (key == "delta_cf2") ap.p.delta_cf2 = v;
        else if (key == "delta_V")   ap.p.delta_V   = v;
        // KanamoriParams (K_U overrides the U→kp.U default)
        else if (key == "K_U")       { ap.kp.U = v; kp_U_explicit = true; }
        else if (key == "J")         ap.kp.J        = v;
        else if (key == "U_prime")   { ap.kp.U_prime = v; U_prime_explicit = true; }
        // SCFParams
        else if (key == "S0")        ap.scf.S0       = v;
        else if (key == "alpha")     ap.scf.alpha    = v;
        else if (key == "T")         ap.scf.T        = v;
        else if (key == "N_target")  ap.scf.N_target = v;
        else if (key == "grid")      ap.scf.grid     = static_cast<int>(v);
        // Convergance
        else if (key == "tol")         ap.kp.tol        = v;
        // Seed (see seeds.h)
        else if (key == "seed_strength") ap.seed_strength = v;
        else if (key == "seed_rng")      ap.seed_rng      = static_cast<unsigned>(v);
        else
            std::cerr << "load_params: unknown key '" << key << "', skipping\n";
    }

    // Enforce U' = U - 2J (rotational invariance) unless explicitly set
    if (!U_prime_explicit)
        ap.kp.U_prime = ap.kp.U - 2.0 * ap.kp.J;

    return ap;
}

void write_params(std::ostream& os, const AllParams& ap) {
    const Params&         p   = ap.p;
    const KanamoriParams& kp  = ap.kp;
    const SCFParams&      scf = ap.scf;

    const std::ios::fmtflags saved_flags = os.flags();
    const std::streamsize    saved_prec  = os.precision();
    os << std::fixed << std::setprecision(4);

    os << "Parameters:\n";
    os << "  tpi       = " << p.tpi       << "\n";
    os << "  tdelta    = " << p.tdelta    << "\n";
    os << "  t2xy      = " << p.t2xy      << "\n";
    os << "  t2yzxz    = " << p.t2yzxz    << "\n";
    os << "  tg        = " << p.tg        << "\n";
    os << "  lam       = " << p.lam       << "\n";
    os << "  U         = " << p.U         << "\n";
    os << "  t_perp    = " << p.t_perp    << "\n";
    os << "  t_perp_xy = " << p.t_perp_xy << "\n";
    os << "  delta_cf1 = " << p.delta_cf1 << "\n";
    os << "  delta_cf2 = " << p.delta_cf2 << "\n";
    os << "  delta_V   = " << p.delta_V   << "\n";
    os << "  theta     = " << p.theta     << " rad\n";
    os << "  phi       = " << p.phi       << " rad\n";
    os << "  S0        = " << scf.S0      << "\n";
    os << "  alpha     = " << scf.alpha   << "\n";
    os << "  T         = " << scf.T       << "\n";
    os << "  N_target  = " << scf.N_target << "\n";
    os << "  grid      = " << scf.grid    << " x " << scf.grid << "\n\n";

    os << "KanamoriParams:\n";
    os << "  U       = " << kp.U       << "\n";
    os << "  U'      = " << kp.U_prime << "\n";
    os << "  J       = " << kp.J       << "\n\n";
    os << "ParamFile   = " << ap.param_file << "\n";
    os << "RhoOutFile  = " << ap.rho_out_file << "\n";
    os << "Seed          = " << ap.seed          << "\n";
    os << "SeedStrength  = " << ap.seed_strength << "\n";
    os << "SeedRng       = " << ap.seed_rng      << "\n\n";

    os.flags(saved_flags);
    os.precision(saved_prec);
}
