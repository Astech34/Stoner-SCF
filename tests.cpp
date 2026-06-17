#include <gtest/gtest.h>
#include "hamiltonian.h"
#include "scf.h"
#include <cmath>
#include <iomanip>

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
    p.tpi      = 1.0;
    p.tdelta = 0.1;
    p.t2xy      = 0.1;
    p.t2yzxz      = 0.0;
    p.tg = 0.0;
    p.lam     = 0.0;
    p.U       = 0.0;

    const Mat6 H = H0(0.0, 0.0, p);

    const double e_yz = -2*p.tdelta - 2*p.tpi;  // -2.2
    const double e_xz = -2*p.tpi     - 2*p.tdelta;  // -2.2
    const double e_xy = -4*p.tpi     - 4*p.t2xy;   // -4.4

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
    p.tpi      = 1.0;
    p.tdelta = 0.1;
    p.t2xy      = 0.1;
    p.t2yzxz      = 0.0;
    p.tg = 0.0;
    p.lam     = 0.3;
    p.U       = 2.0;

    const double S  = 0.3;
    const double kx = 1.1;
    const double ky = 0.7;

    const Mat6 H    = singleLayer(kx, ky, S, p);
    const Mat6 diff = H - H.adjoint();

    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
}

// HubbardU splits spin-up and spin-down bands by ±U*S on the diagonal (n̂ = ẑ).
TEST(Hamiltonian, HubbardUShift) {
    const double S = 0.5;
    Params p;
    p.U     = 2.0;
    p.theta = 0.0;  // n̂ = ẑ — off-diagonal blocks vanish
    p.phi   = 0.0;

    const Mat6 H   = HubbardU(S, p);
    const double shift = p.U * S;  // 1.0
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

// HubbardU at theta=0.5, phi=0.1 checked against Mathematica.
// Mathematica gives +US*(n̂·σ)⊗I₃; C++ uses scale = -US, so signs are flipped.
// With US=1, nz=cos(0.5)≈0.877583, nx-iny≈0.47703-0.0478627i:
//   spin-up diagonal:          -nz
//   spin-down diagonal:        +nz
//   upper-right off-diagonal:  -(nx-iny)
//   lower-left off-diagonal:   -(nx+iny)
/*
TEST(Hamiltonian, HubbardAngles) {
    const double S = 0.5;
    Params p;
    p.U     = 2.0;

    for (int i = 0; i < 10; i++){
        for (int j = 0; j < 10; j++){
            p.theta = i/5.0 * M_PI;  // 0 to π
            p.phi   = j/5.0 * 2*M_PI;  // 0 to 2π

            const Mat6 H = HubbardU(S, p);

            // Analytically computed from n̂ components (US = 1.0)
            const double nz  =  std::cos(p.theta);
            const double nx  =  std::sin(p.theta) * std::cos(p.phi);
            const double ny  =  std::sin(p.theta) * std::sin(p.phi);
            const double tol = 1e-12;

            // spin-up   (rows 0,1,2): -nz
            EXPECT_NEAR(H(0,0).real(), -nz, tol);
            EXPECT_NEAR(H(1,1).real(), -nz, tol);
            EXPECT_NEAR(H(2,2).real(), -nz, tol);

            // spin-down (rows 3,4,5): +nz
            EXPECT_NEAR(H(3,3).real(),  nz, tol);
            EXPECT_NEAR(H(4,4).real(),  nz, tol);
            EXPECT_NEAR(H(5,5).real(),  nz, tol);

            // upper-right off-diagonal: -(nx - i*ny)
            EXPECT_NEAR(H(0,3).real(), -nx, tol);
            EXPECT_NEAR(H(0,3).imag(),  ny, tol);
            EXPECT_NEAR(H(1,4).real(), -nx, tol);
            EXPECT_NEAR(H(1,4).imag(),  ny, tol);
            EXPECT_NEAR(H(2,5).real(), -nx, tol);
            EXPECT_NEAR(H(2,5).imag(),  ny, tol);

            // lower-left off-diagonal: -(nx + i*ny)
            EXPECT_NEAR(H(3,0).real(), -nx, tol);
            EXPECT_NEAR(H(3,0).imag(), -ny, tol);
            EXPECT_NEAR(H(4,1).real(), -nx, tol);
            EXPECT_NEAR(H(4,1).imag(), -ny, tol);
            EXPECT_NEAR(H(5,2).real(), -nx, tol);
            EXPECT_NEAR(H(5,2).imag(), -ny, tol);
        }
    }
}
*/

TEST(Hamiltonian, BiLayerSimple){
    Params p;
    p.tpi      = 1.0;
    p.tdelta = 0;
    p.t2xy      = 0;
    p.t2yzxz      = 0;
    p.tg = 0;
    p.lam     = 0;
    p.U       = 2.0;
    p.t_perp  = 0.3;
    // kx=ky=S = 0
    Mat12 H = Mat12::Zero();
    const Mat12 HM = bilayerHamiltonian(0, 0, 0, p);

    H(0,0) = -2.0;  H(1, 1) = -2.0; H(2,2) = -4.0;
    H(3,3) = -2.0;  H(4, 4) = -2.0; H(5,5) = -4.0;
    H(6,6) = -2.0;  H(7, 7) = -2.0; H(8,8) = -4.0;
    H(9,9) = -2.0;  H(10,10) = -2.0; H(11,11) = -4.0;

    H(0,6) = 0.3;  H(1,7) = 0.3; H(2,8) = 0;
    H(3,9) = 0.3;  H(4,10) = 0.3; H(5,11) = 0;

    H(6,0) = 0.3;  H(7,1) = 0.3; H(8,2) = 0;
    H(9,3) = 0.3;  H(10,4) = 0.3; H(11,5) = 0;

    EXPECT_NEAR((HM - H).norm(), 0.0, 1e-12);

}

// For printing Hamiltonian for testing
TEST(Hamiltonian, PrintBiLayerH){
    Params p;
    p.tpi      = 1.0;
    p.tdelta = 0.1;
    p.t2xy      = 0.1;
    p.t2yzxz      = 0.1;
    p.tg = 0.1;
    p.lam     = 0.1;
    p.U       = 2.0;
    p.t_perp  = 0.3;
    p.theta   = 0.0;
    p.phi     = 0.0;

    const double S  = 0.3;
    const double kx = 0.0;
    const double ky = 0.0;

    const Mat12 H = bilayerHamiltonian(kx, ky, S, p);

    //std::cout << "\nBilayer H (kx=" << kx << ", ky=" << ky << ", S=" << S << "):\n"
             // << H << "\n";
}

TEST(Hamiltonian, KFSingle){
    const double U = 5;
    const double J = 0.8;
    const double U_prime = U - 2*J;  // 2.4

    KanamoriParams kp;
    kp.U = U;
    kp.J = J;
    kp.U_prime = U_prime;

    Mat6 rho = Mat6::Zero();
    for (int i = 0; i < 6; ++i) {
        for (int j = i; j < 6; ++j) {
            rho(i, j) = cd(i, j);
            rho(j, i) = std::conj(rho(i, j));
            if (i == j){
                rho(i, j) = cd(i, 0);  // make diagonal entries real and nonzero
            }
        }
    }
    std::cout << "\nInput rho:\n" << rho << "\n";

    Mat6 H = kanamori_layer(rho, kp);

    Mat6 H_expected;
    H_expected << cd(53.4, 0), cd(4.8, 2.6), cd(4.8, 5.2), cd(-2.4, 22.2), cd(-0.8, 16), cd(-1.6, 19.4),
        cd(4.8, -2.6), cd(52.4, 0), cd(3.8, 5.2), cd(-3.4, 13.4), cd(-6.6, 26.4), cd(-5, 20.2),
        cd(4.8, -5.2), cd(3.8, -5.2), cd(51.4, 0), cd(-6.8, 14.2), cd(-7.6, 17.6), cd(-10.8, 30.6),
        cd(-2.4, -22.2), cd(-3.4, -13.4), cd(-6.8, -14.2), cd(33.6, 0), cd(-7.8, 10.4), cd(-7.8, 13),
        cd(-0.8, -16), cd(-6.6, -26.4), cd(-7.6, -17.6), cd(-7.8, -10.4), cd(32.6, 0), cd(-8.8, 13),
        cd(-1.6, -19.4), cd(-5, -20.2), cd(-10.8, -30.6), cd(-7.8, -13), cd(-8.8, -13), cd(31.6, 0);

        
    const Mat6 diff = H - H_expected;
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);

    std::cout << "\nKanamori H:\n" << H << "\n";

}

TEST(Hamiltonian, ComputeDensityMatrix){
    const double U = 5;
    const double J = 0.8;
    const double U_prime = U - 2*J;  // 2.4

    KanamoriParams kp;
    kp.U = U;
    kp.J = J;
    kp.U_prime = U_prime;

    Mat6 rho = Mat6::Zero();
    for (int i = 0; i < 6; ++i) {
        for (int j = i; j < 6; ++j) {
            rho(i, j) = cd(i, j);
            rho(j, i) = std::conj(rho(i, j));
            if (i == j){
                rho(i, j) = cd(i, 0);  // make diagonal entries real and nonzero
            }
        }
    }

    Mat6 H = kanamori_layer(rho, kp);

    Mat12 H12 = Mat12::Zero();
    H12.block<6,6>(0,0) = H;
    H12.block<6,6>(6,6) = H;

    Eigensystem result;
    result.evals.resize(1);
    result.evecs.resize(1);

    Eigen::SelfAdjointEigenSolver<Mat12> solver(H12);

    result.evals[0] = solver.eigenvalues();
    result.evecs[0] = solver.eigenvectors();

    Mat12 rho_computed = compute_density_matrix(result, 1.0, 0.005);

    Mat12 rhoexpec;
    rhoexpec << cd(0.07911880142, 0), cd(0.08404216451, 0.02767789033), cd(0.08045967674, 0.06438557567), cd(0.03562743845, 0.1118671034), cd(-0.004194648594, 0.1333657112), cd(-0.06222607281, 0.1376624437), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0),
     cd(0.08404216451, -0.02767789033), cd(0.09895436846, 0), cd(0.1079902898, 0.04024516266), cd(0.07697857341, 0.1063648718), cd(0.04219925126, 0.1431321236), cd(-0.01794020895, 0.167997188), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0),
     cd(0.08045967674, -0.06438557567), cd(0.1079902898, -0.04024516266), cd(0.1342191963, 0), cd(0.1272668424, 0.0847699627), cd(0.1042650782, 0.1390394783), cd(0.04874676953, 0.1906340461), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0),
     cd(0.03562743845, -0.1118671034), cd(0.07697857341, -0.1063648718), cd(0.1272668424, -0.0847699627), cd(0.1742134986, 0), cd(0.1866786523, 0.0659858562), cd(0.1666221303, 0.1499720742), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0),
     cd(-0.004194648594, -0.1333657112), cd(0.04219925126, -0.1431321236), cd(0.1042650782, -0.1390394783), cd(0.1866786523, -0.0659858562), cd(0.225028788, 0), cd(0.2353481837, 0.09759221246), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0),
     cd(-0.06222607281, -0.1376624437), cd(-0.01794020895, -0.167997188), cd(0.04874676953, -0.1906340461), cd(0.1666221303, -0.1499720742), cd(0.2353481837, -0.09759221246), cd(0.2884653473, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0),
     cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0.07911880142, 0), cd(0.08404216451, 0.02767789033), cd(0.08045967674, 0.06438557567), cd(0.03562743845, 0.1118671034), cd(-0.004194648594, 0.1333657112), cd(-0.06222607281, 0.1376624437),
     cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0.08404216451, -0.02767789033), cd(0.09895436846, 0), cd(0.1079902898, 0.04024516266), cd(0.07697857341, 0.1063648718), cd(0.04219925126, 0.1431321236), cd(-0.01794020895, 0.167997188),
     cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0.08045967674, -0.06438557567), cd(0.1079902898, -0.04024516266), cd(0.1342191963, 0), cd(0.1272668424, 0.0847699627), cd(0.1042650782, 0.1390394783), cd(0.04874676953, 0.1906340461),
     cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0.03562743845, -0.1118671034), cd(0.07697857341, -0.1063648718), cd(0.1272668424, -0.0847699627), cd(0.1742134986, 0), cd(0.1866786523, 0.0659858562), cd(0.1666221303, 0.1499720742),
     cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(-0.004194648594, -0.1333657112), cd(0.04219925126, -0.1431321236), cd(0.1042650782, -0.1390394783), cd(0.1866786523, -0.0659858562), cd(0.225028788, 0), cd(0.2353481837, 0.09759221246),
     cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(-0.06222607281, -0.1376624437), cd(-0.01794020895, -0.167997188), cd(0.04874676953, -0.1906340461), cd(0.1666221303, -0.1499720742), cd(0.2353481837, -0.09759221246), cd(0.2884653473, 0);

    
    const Mat12 diff = rho_computed - rhoexpec;
    //std::cout << "\ndiff:\n" << diff << "\n";
    // decreased tolerance due to LAPAC diagonalization
    EXPECT_NEAR(diff.norm(), 0.0, 1e-9);
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

TEST(SOC, kron){
    Eigen::Matrix<cd, 3, 3> A;
    A << 1, 2, 3, 4, 5, 6, 
         7, 8, 9;
    Eigen::Matrix<cd, 2, 2> B;
    B << 1, 2, 3, 4;
    Mat6 C = kron(A, B);
    Mat6 C_expected;
    C_expected << 1, 2, 2, 4, 3, 6,
                  3, 4, 6, 8, 9, 12,
                  4, 8, 5, 10, 6, 12,
                  12, 16, 15, 20, 18, 24,
                  7, 14, 8, 16, 9, 18,
                  21, 28, 24, 32, 27, 36;
    const Mat6 diff = C - C_expected;
    std::cout << "\nC:\n" << C << "\n";
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);

}
TEST(SOC, SOCAngles){
    Mat6 H1 = SOC(0.5, 0.2, 0.3);
    Mat6 mat;
    mat << 
    // Row 1
    cd(0.0, 0.0), cd(0.0, 0.245017), cd(0.0, 0.0146777), cd(0.0, 0.0), cd(0.0, -0.0496673), cd(0.238834, 0.0724074),
    // Row 2
    cd(0.0, -0.245017), cd(0.0, 0.0), cd(0.0, -0.047449), cd(0.0, 0.0496673), cd(0.0, 0.0), cd(0.0738801, -0.234073),
    // Row 3
    cd(0.0, -0.0146777), cd(0.0, 0.047449), cd(0.0, 0.0), cd(-0.238834, -0.0724074), cd(-0.0738801, 0.234073), cd(0.0, 0.0),
    // Row 4
    cd(0.0, 0.0), cd(0.0, -0.0496673), cd(-0.238834, 0.0724074), cd(0.0, 0.0), cd(0.0, -0.245017), cd(0.0, -0.0146777),
    // Row 5
    cd(0.0, 0.0496673), cd(0.0, 0.0), cd(-0.0738801, -0.234073), cd(0.0, 0.245017), cd(0.0, 0.0), cd(0.0, 0.047449),
    // Row 6
    cd(0.238834, -0.0724074), cd(0.0738801, 0.234073), cd(0.0, 0.0), cd(0.0, 0.0146777), cd(0.0, -0.047449), cd(0.0, 0.0);
    Mat6 diff = H1 - mat;
    //std::cout << "\nH1:\n" << diff << "\n";
    EXPECT_NEAR(diff.norm(), 0.0, 1e-3);
}

// -----------------------------------------------------------------------------
// Testing SCF.h functions
// -----------------------------------------------------------------------------

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

// Compute eigensystem
TEST(SCF, BilayerEigenSystem){
    Params p;
    p.tpi      = 1.0;
    p.tdelta = 0;
    p.t2xy      = 0;
    p.lam     = 0;
    p.U       = 0;
    p.t_perp  = 0.3;

    auto result = compute_eigensystem_grid(0.0, 2, p);

    auto evals = result.evals;

    Eigen::Vector<double, 12> expected_evals;
    expected_evals << 1.7, 1.7, 1.7, 1.7, 2.3, 2.3, 2.3, 2.3, 4, 4, 4, 4;
    
    EXPECT_NEAR((evals[0] - expected_evals).norm(), 0.0, 1e-12);

}

// Compute Total Energy single k
TEST(SCF, TotalEnergy){
    Eigensystem sys;
    Eigen::Vector<double, 12> ev;
    ev << 3, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    sys.evals = {ev};

    double mu = 0.5;
    double T = 0.05;
    double E1 = 0.000045397868796;

    double E_total = calculate_band_energy(sys, mu, T);

    //std::cout << std::setprecision(17) << "E_total = " << E_total << "\n";
    //std::cout << std::setprecision(17) << "E expected = " << (E1) << "\n";
    EXPECT_NEAR(E_total, E1, 1e-12);

}

TEST(SCF, TotalEnergyTwoK){
    Eigensystem sys;
    Eigen::Vector<double, 12> ev1;
    Eigen::Vector<double, 12> ev2;
    ev1 << 3, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    ev2 << 6, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    sys.evals = {ev1, ev2};

    double E1 = 0.000045397868796;
    double E2 = 0.000045397868796;

    double mu = 0.5;
    double T = 0.05;

    double E_total = calculate_band_energy(sys, mu, T);
    //std::cout << std::setprecision(17) << "E_total = " << E_total << "\n";
    //std::cout << std::setprecision(17) << "E expected = " << (E1 + E2)/2.0 << "\n";
    EXPECT_NEAR(E_total, (E1 + E2)/(2.0), 1e-12);

}

TEST(SCF, TestNTotal){
    Eigensystem sys;
    Eigen::Vector<double, 12> ev1;
    Eigen::Vector<double, 12> ev2;
    ev1 << 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6;
    ev2 << 3, 1, 3, 4, 5, 6, 3, 1, 3, 4, 5, 6;
    sys.evals = {ev1, ev2};

    double Ncalc = 0.0000907957374984;
    double Ntest = NTotal(sys, 0.5, 0.05);

    std::cout << std::setprecision(17) << "N_diff = " << Ntest-Ncalc << "\n";
    std::cout << std::setprecision(17) << "N_diff * 1e12 = " << (Ntest-Ncalc) * 1e12 << "\n";
    EXPECT_NEAR(Ntest, Ncalc, 1e-12);
}

TEST(SCF, CheckMu){
    Params p;
    p.tpi      = 1.0;
    p.tdelta = 0.1;
    p.t2xy      = 0.1;
    p.lam     = 0.3;
    p.U       = 2.0;

    Eigensystem sys = compute_eigensystem_grid(0.2, 100, p);

    double mu = find_mu(sys, 0.05, 10.0);
    double T = 0.05;

    double Nt = NTotal(sys, mu, T);

    EXPECT_NEAR(Nt, 10.0, 1e-12);

}

TEST(SCF, SpinCross){
    Eigen::Vector<cd, 12> col;
    col << cd(0,1), cd(2,0), cd(0,3),   // L1 spin-up   (indices 0-2)
           cd(4,0), cd(5,0), cd(6,0),   // L1 spin-down (indices 3-5)
           cd(7,0), cd(8,0), cd(9,0),   // L2 spin-up   (indices 6-8)
           cd(10,0),cd(11,0),cd(12,0);  // L2 spin-down (indices 9-11)

    Eigen::Vector<cd, 6> psi_up, psi_dn;
    cd cross = spin_cross(col, psi_up, psi_dn);

    // Check spinor extraction
    EXPECT_EQ(psi_up(0), cd(0,1));   // L1↑
    EXPECT_EQ(psi_up(3), cd(7,0));   // L2↑
    EXPECT_EQ(psi_dn(0), cd(4,0));   // L1↓
    EXPECT_EQ(psi_dn(3), cd(10,0));  // L2↓

    // cross = ψ↑† · ψ↓
    EXPECT_NEAR(cross.real(), 276, 1e-12);
    EXPECT_NEAR(cross.imag(), -22, 1e-12);
}

// calculateS with U=0 and lam=0: no exchange and no SOC means the spin-up and
// spin-down blocks of the Hamiltonian are identical.  Every eigenstate carries
// equal ↑ and ↓ weight, so the spin projection sums to zero exactly regardless
// of the initial guess S or the chosen direction n̂.
TEST(SCF, CalculateS_ZeroExchangeGivesZeroMagnetization) {
    Params p;
    p.tpi      = 1.0;
    p.tdelta = 0.1;
    p.t2xy      = 0.1;
    p.lam     = 0.0;   // no SOC
    p.U       = 0.0;   // no exchange
    p.t_perp  = 0.1;
    p.theta   = 0.3;   // arbitrary, should not matter
    p.phi     = 1.1;

    const CalcResult r = calculateS(0.4, 10, 0.05, 10.0, p);

    EXPECT_NEAR(r.S_new, 0.0, 1e-12);
}

// calculateS with lam=0: the Stoner Hamiltonian has exact SU(2) spin-rotation
// symmetry when SOC is absent — H0 is spin-diagonal with equal ↑/↓ dispersions,
// so the Hubbard term for any n̂ is unitarily equivalent to the n̂=ẑ case.
// S_new must therefore be independent of (theta, phi).
TEST(SCF, CalculateS_RotationalSymmetryNoSOC) {
    Params p;
    p.tpi      = 1.0;
    p.tdelta = 0.1;
    p.t2xy      = 0.1;
    p.lam     = 0.0;   // no SOC → exact spin-rotation symmetry
    p.U       = 2.0;
    p.t_perp  = 0.0;

    const double S_in = 0.3;
    const double T    = 0.05;
    const double N    = 10.0;
    const int    grid = 10;

    p.theta = 0.0;         p.phi = 0.0;
    const double S_z    = calculateS(S_in, grid, T, N, p).S_new;

    p.theta = M_PI / 2.0; p.phi = 0.0;
    const double S_x    = calculateS(S_in, grid, T, N, p).S_new;

    p.theta = M_PI / 4.0; p.phi = M_PI / 3.0;
    const double S_diag = calculateS(S_in, grid, T, N, p).S_new;

    EXPECT_NEAR(S_x,    S_z, 1e-10);
    EXPECT_NEAR(S_diag, S_z, 1e-10);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
