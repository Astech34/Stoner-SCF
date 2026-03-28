#include <gtest/gtest.h>
#include "hamiltonian.h"
#include "scf.h"
#include <cmath>

// -----------------------------------------------------------------------------
// Hamiltonian tests
// -----------------------------------------------------------------------------

// H0 at Gamma (kx=ky=0) with lam=0 and S=0 should be real diagonal.
// At Gamma, cos(kx)=cos(ky)=1, so:
//   e_yz = -2*t_delta - 2*t1
//   e_xz = -2*t1     - 2*t_delta
//   e_xy = -4*t1     - 4*t2
// Spin-major ordering: rows/cols (0,1,2) = spin-up  (yz, xz, xy)
//                                (3,4,5) = spin-down (yz, xz, xy)
TEST(Hamiltonian, H0GammaPointDiagonal) {
    Params p;
    p.t1      = 1.0;
    p.t_delta = 0.1;
    p.t2      = 0.1;
    p.lam     = 0.0;
    p.U       = 0.0;

    const Mat6 H = H0(0.0, 0.0, p);

    const double e_yz = -2*p.t_delta - 2*p.t1;  // -2.2
    const double e_xz = -2*p.t1     - 2*p.t_delta;  // -2.2
    const double e_xy = -4*p.t1     - 4*p.t2;   // -4.4

    const double tol = 1e-12;

    // Matrix is 6x6
    EXPECT_EQ(H.rows(), 6);
    EXPECT_EQ(H.cols(), 6);

    // Diagonal entries match analytic values
    EXPECT_NEAR(H(0,0).real(), e_yz, tol);  // spin-up   yz
    EXPECT_NEAR(H(1,1).real(), e_xz, tol);  // spin-up   xz
    EXPECT_NEAR(H(2,2).real(), e_xy, tol);  // spin-up   xy
    EXPECT_NEAR(H(3,3).real(), e_yz, tol);  // spin-down yz
    EXPECT_NEAR(H(4,4).real(), e_xz, tol);  // spin-down xz
    EXPECT_NEAR(H(5,5).real(), e_xy, tol);  // spin-down xy

    // All off-diagonal elements are zero (no SOC, no mixing)
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            if (i != j)
                EXPECT_NEAR(std::abs(H(i,j)), 0.0, tol);
}

// The full Hamiltonian must be Hermitian for any k and parameters.
TEST(Hamiltonian, FullHamiltonianIsHermitian) {
    Params p;
    p.t1      = 1.0;
    p.t_delta = 0.1;
    p.t2      = 0.1;
    p.lam     = 0.3;
    p.U       = 2.0;

    const double S  = 0.3;
    const double kx = 1.1;
    const double ky = 0.7;

    const Mat6 H    = singleLayer(kx, ky, S, p);
    const Mat6 diff = H - H.adjoint();

    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
}

// HubbardU splits spin-up and spin-down bands by ±U*S on the diagonal.
TEST(Hamiltonian, HubbardUShift) {
    const double S = 0.5;
    const double U = 2.0;

    const Mat6 H   = HubbardU(S, U);
    const double shift = U * S;  // 1.0
    const double tol   = 1e-12;

    // spin-up   (rows 0,1,2): negative shift
    EXPECT_NEAR(H(0,0).real(), -shift, tol);
    EXPECT_NEAR(H(1,1).real(), -shift, tol);
    EXPECT_NEAR(H(2,2).real(), -shift, tol);

    // spin-down (rows 3,4,5): positive shift
    EXPECT_NEAR(H(3,3).real(),  shift, tol);
    EXPECT_NEAR(H(4,4).real(),  shift, tol);
    EXPECT_NEAR(H(5,5).real(),  shift, tol);

    // No off-diagonal entries
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            if (i != j)
                EXPECT_NEAR(std::abs(H(i,j)), 0.0, tol);
}

// -----------------------------------------------------------------------------
// SOC matrix tests
// Spin-major ordering: (0,1,2) = up-yz, up-xz, up-xy | (3,4,5) = dn-yz, dn-xz, dn-xy
//
// Analytic non-zero entries for SOC = lam*(Lx⊗sx + Ly⊗sy + Lz⊗sz):
//   Lz⊗sz:  H(0,1)=+i/2  H(1,0)=-i/2  H(3,4)=-i/2  H(4,3)=+i/2   (yz<->xz, same spin)
//   Lx⊗sx:  H(1,5)=+i/2  H(4,2)=+i/2  H(2,4)=-i/2  H(5,1)=-i/2   (xz<->xy, spin-flip)
//   Ly⊗sy:  H(0,5)=-1/2  H(3,2)=+1/2  H(2,3)=+1/2  H(5,0)=-1/2   (yz<->xy, spin-flip)
//   (all values multiplied by lam)
// -----------------------------------------------------------------------------

// lam=0 => zero matrix
TEST(SOC, ZeroLambdaIsZero) {
    const Mat6 H = SOC(0.0);
    EXPECT_NEAR(H.norm(), 0.0, 1e-12);
}

// SOC must be Hermitian for any lam
TEST(SOC, IsHermitian) {
    const double lam = 0.3;
    const Mat6 H    = SOC(lam);
    const Mat6 diff = H - H.adjoint();
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
}

// Diagonal entries are always zero (L·S has no diagonal in t2g basis)
TEST(SOC, DiagonalIsZero) {
    const Mat6 H   = SOC(0.5);
    const double tol = 1e-12;
    for (int i = 0; i < 6; i++)
        EXPECT_NEAR(std::abs(H(i,i)), 0.0, tol);
}

// Lz⊗sz: same-spin yz<->xz coupling
TEST(SOC, LzSzOrbitalCoupling) {
    const double lam = 1.0;
    const Mat6   H   = SOC(lam);
    const double tol = 1e-12;

    // up-yz <-> up-xz : +i*lam/2
    EXPECT_NEAR(H(0,1).real(),  0.0,      tol);
    EXPECT_NEAR(H(0,1).imag(),  lam/2.0,  tol);
    EXPECT_NEAR(H(1,0).real(),  0.0,      tol);
    EXPECT_NEAR(H(1,0).imag(), -lam/2.0,  tol);

    // dn-yz <-> dn-xz : -i*lam/2
    EXPECT_NEAR(H(3,4).real(),  0.0,      tol);
    EXPECT_NEAR(H(3,4).imag(), -lam/2.0,  tol);
    EXPECT_NEAR(H(4,3).real(),  0.0,      tol);
    EXPECT_NEAR(H(4,3).imag(),  lam/2.0,  tol);
}

// Lx⊗sx: spin-flip xz<->xy coupling
TEST(SOC, LxSxSpinFlip) {
    const double lam = 1.0;
    const Mat6   H   = SOC(lam);
    const double tol = 1e-12;

    // up-xz <-> dn-xy : +i*lam/2
    EXPECT_NEAR(H(1,5).real(),  0.0,      tol);
    EXPECT_NEAR(H(1,5).imag(),  lam/2.0,  tol);
    EXPECT_NEAR(H(5,1).real(),  0.0,      tol);
    EXPECT_NEAR(H(5,1).imag(), -lam/2.0,  tol);

    // dn-xz <-> up-xy : +i*lam/2
    EXPECT_NEAR(H(4,2).real(),  0.0,      tol);
    EXPECT_NEAR(H(4,2).imag(),  lam/2.0,  tol);
    EXPECT_NEAR(H(2,4).real(),  0.0,      tol);
    EXPECT_NEAR(H(2,4).imag(), -lam/2.0,  tol);
}

// Ly⊗sy: spin-flip yz<->xy coupling (real)
TEST(SOC, LySySpinFlip) {
    const double lam = 1.0;
    const Mat6   H   = SOC(lam);
    const double tol = 1e-12;

    // up-yz <-> dn-xy : -lam/2 (real)
    EXPECT_NEAR(H(0,5).real(), -lam/2.0,  tol);
    EXPECT_NEAR(H(0,5).imag(),  0.0,      tol);
    EXPECT_NEAR(H(5,0).real(), -lam/2.0,  tol);
    EXPECT_NEAR(H(5,0).imag(),  0.0,      tol);

    // dn-yz <-> up-xy : +lam/2 (real)
    EXPECT_NEAR(H(3,2).real(),  lam/2.0,  tol);
    EXPECT_NEAR(H(3,2).imag(),  0.0,      tol);
    EXPECT_NEAR(H(2,3).real(),  lam/2.0,  tol);
    EXPECT_NEAR(H(2,3).imag(),  0.0,      tol);
}

// Scales linearly with lam
TEST(SOC, ScalesWithLambda) {
    const Mat6 H1 = SOC(1.0);
    const Mat6 H2 = SOC(2.0);
    EXPECT_NEAR((H2 - 2.0*H1).norm(), 0.0, 1e-12);
}

// -----------------------------------------------------------------------------
// Brent root-finder tests
// -----------------------------------------------------------------------------

// f(x) = x^2 - 2  =>  root at sqrt(2) ~ 1.41421356...
TEST(Brent, SquareRootOf2) {
    auto f = [](double x) { return x*x - 2.0; };
    double root = brent(f, 1.0, 2.0);
    EXPECT_NEAR(root, std::sqrt(2.0), 1e-9);
}

// f(x) = cos(x) - x  (Dottie number)  =>  root ~ 0.739085...
TEST(Brent, CosineFixedPoint) {
    auto f = [](double x) { return std::cos(x) - x; };
    double root = brent(f, 0.0, 1.0);
    EXPECT_NEAR(root, 0.7390851332151607, 1e-9);
}

// f(x) = x^3 - x - 2  =>  root at x ~ 1.5213797...
TEST(Brent, CubicRoot) {
    auto f = [](double x) { return x*x*x - x - 2.0; };
    double root = brent(f, 1.0, 2.0);
    EXPECT_NEAR(root, 1.5213797068045676, 1e-9);
}

// f(x) = exp(x) - 3  =>  root at ln(3) ~ 1.09861...
TEST(Brent, ExponentialRoot) {
    auto f = [](double x) { return std::exp(x) - 3.0; };
    double root = brent(f, 0.0, 2.0);
    EXPECT_NEAR(root, std::log(3.0), 1e-9);
}

// Bracket does not straddle zero — brent() should assert/abort
TEST(Brent, BadBracketAsserts) {
    auto f = [](double x) { return x*x + 1.0; };  // always positive
    EXPECT_DEATH(brent(f, 1.0, 2.0), "");
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
