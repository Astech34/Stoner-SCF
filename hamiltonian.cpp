#include "hamiltonian.h"
#include <cmath>

Mat6 H0(double kx, double ky, const Params& p) {
    const double e_xy = -2*p.t1*(std::cos(kx) + std::cos(ky)) - 4*p.t2*std::cos(kx)*std::cos(ky);
    const double e_xz = -2*p.t1*std::cos(kx) - 2*p.t_delta*std::cos(ky);
    const double e_yz = -2*p.t_delta*std::cos(kx) - 2*p.t1*std::cos(ky);

    Mat6 H = Mat6::Zero();

    // Spin-major ordering: spin-up (0,1,2) = yz,xz,xy | spin-down (3,4,5) = yz,xz,xy
    H(0,0) = e_yz;  H(1,1) = e_xz;  H(2,2) = e_xy;
    H(3,3) = e_yz;  H(4,4) = e_xz;  H(5,5) = e_xy;

    return H;
}

Mat6 SOC(double lam) {
    // Orbital angular momentum matrices in t2g basis (yz, xz, xy)
    Eigen::Matrix<cd, 3, 3> Lx, Ly, Lz;

    Lx << cd(0,0),  cd(0,0),  cd(0,0),
          cd(0,0),  cd(0,0),  cd(0,1),
          cd(0,0),  cd(0,-1), cd(0,0);

    Ly << cd(0,0),  cd(0,0),  cd(0,-1),
          cd(0,0),  cd(0,0),  cd(0,0),
          cd(0,1),  cd(0,0),  cd(0,0);

    Lz << cd(0,0),  cd(0,1),  cd(0,0),
          cd(0,-1), cd(0,0),  cd(0,0),
          cd(0,0),  cd(0,0),  cd(0,0);

    // Spin-1/2 matrices
    Eigen::Matrix<cd, 2, 2> sx, sy, sz;

    sx << cd(0,0), cd(0.5,0),
          cd(0.5,0), cd(0,0);

    sy << cd(0,0),    cd(0,-0.5),
          cd(0,0.5),  cd(0,0);

    sz << cd(0.5,0),  cd(0,0),
          cd(0,0),    cd(-0.5,0);

    // Kronecker products: kron(L, s) gives spin-major 6x6
    Mat6 Hsoc = Mat6::Zero();

    // kron(Lx, sx)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int si = 0; si < 2; si++)
                for (int sj = 0; sj < 2; sj++)
                    Hsoc(2*i+si, 2*j+sj) += lam * Lx(i,j) * sx(si,sj);

    // kron(Ly, sy)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int si = 0; si < 2; si++)
                for (int sj = 0; sj < 2; sj++)
                    Hsoc(2*i+si, 2*j+sj) += lam * Ly(i,j) * sy(si,sj);

    // kron(Lz, sz)
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int si = 0; si < 2; si++)
                for (int sj = 0; sj < 2; sj++)
                    Hsoc(2*i+si, 2*j+sj) += lam * Lz(i,j) * sz(si,sj);

    return Hsoc;
}

Mat6 HubbardU(double S, double U) {
    const double shift = U * S;

    Mat6 H = Mat6::Zero();

    // spin-up: negative shift | spin-down: positive shift
    H(0,0) = -shift;  H(1,1) = -shift;  H(2,2) = -shift;
    H(3,3) =  shift;  H(4,4) =  shift;  H(5,5) =  shift;

    return H;
}

Mat6 singleLayer(double kx, double ky, double S, const Params& p) {
    // SOC is k-independent — caller may want to cache this for performance
    return H0(kx, ky, p) + SOC(p.lam) + HubbardU(S, p.U);
}