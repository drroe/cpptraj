#ifndef INC_DIAG_3X3_H
#define INC_DIAG_3X3_H
#include <complex>
/*! \file Diag_3x3.h
 *  \brief Analytical diagonalization of 3x3 symmetric matrices.
 *  \author Original code by Andrew Simmonett, NIH.
 *  \author DRR : Remove C++11 requirement, change input/output arguments and eigenvector ordering.
 */

/// Diagonalize input 3x3 symmetric matrix.
/**
  * \param mat Input: flat 3x3 symmetric row-major matrix of dimension 9. Output: eigenvectors in columns (unsorted).
  * \param vec Output eigenvalues (unsorted).
  */
template <typename Real>
void Diag_3x3(Real* mat, Real* val) {
    typedef std::complex<Real> Complex;
    // Real, symmetric matrix version (upper triangle input)
    Real m00 = mat[0];
    Real m01 = mat[1];
    Real m02 = mat[2];
    Real m11 = mat[4];
    Real m12 = mat[5];
    Real m22 = mat[8];

    // Intermediates
    const Real one = static_cast<Real>(1);
    const Real two = static_cast<Real>(2);
    const Real three = static_cast<Real>(3);
    const Real four = static_cast<Real>(4);
    //const Real six = static_cast<Real>(6);
    const Real cbrtTwo = cbrt(two);
    const Real cbrtFour = cbrt(four);
    const Real cbrtEight = cbrtTwo * cbrtFour;
    const Real cbrtSixteen = cbrtFour * cbrtFour;
    const Real sqrtThree = sqrt(three);
    const Real oneSixth = one / 6;
    const Real oneThird = one / three;
    const Complex sqrtThreeMinusI(sqrtThree, -1);
    const Complex sqrtThreePlusI(sqrtThree, 1);

    Real trace = m00 + m11 + m22;
    Real squares = -m00 * m00 - 3 * m01 * m01 - 3 * m02 * m02 + m00 * m11 - m11 * m11 - 3 * m12 * m12 + m00 * m22 +
                   m11 * m22 - m22 * m22;
    Real squares2 = squares * squares;
    Real cubes = 2 * m00 * m00 * m00 + 9 * m00 * m01 * m01 + 9 * m00 * m02 * m02 - 3 * m00 * m00 * m11 +
                 9 * m01 * m01 * m11 - 18 * m02 * m02 * m11 - 3 * m00 * m11 * m11 + 2 * m11 * m11 * m11 +
                 54 * m01 * m02 * m12 - 18 * m00 * m12 * m12 + 9 * m11 * m12 * m12 - 3 * m00 * m00 * m22 -
                 18 * m01 * m01 * m22 + 9 * m02 * m02 * m22 + 12 * m00 * m11 * m22 - 3 * m11 * m11 * m22 +
                 9 * m12 * m12 * m22 - 3 * m00 * m22 * m22 - 3 * m11 * m22 * m22 + 2 * m22 * m22 * m22;
    Complex cbrtTerm = pow(cubes + sqrt(Complex(cubes * cubes + 4 * squares2 * squares, 0)), oneThird);
    Complex cbrtTerm2 = cbrtTerm * cbrtTerm;
    Complex cbrtTerm3 = cbrtTerm2 * cbrtTerm;
    Complex cbrtTerm4 = cbrtTerm3 * cbrtTerm;

    // Eigenvalues
    Real e0 = Complex((trace + cbrtTerm / cbrtTwo - (squares * cbrtTwo) / cbrtTerm) * oneThird).real();
    Real e1 = Complex(((Complex(0, 1) * cbrtTerm * sqrtThreePlusI) / cbrtTwo +
                       (two * (squares + Complex(0, 1) * sqrtThree * squares)) / (cbrtFour * cbrtTerm) + 2 * trace) *
                      oneSixth)
                  .real();
    Real e2 = Complex((-((cbrtTerm + Complex(0, 1) * sqrtThree * cbrtTerm) / cbrtTwo) +
                       (two * (squares - Complex(0, 1) * sqrtThree * squares)) / (cbrtFour * cbrtTerm) + 2 * trace) *
                      oneSixth)
                  .real();

    // Eigenvectors
    Complex v0denom = (three * cbrtTerm * cbrtTwo *
                 (cbrtTerm2 * m02 + cbrtTerm * cbrtTwo * (m00 * m02 - 2 * m02 * m11 + 3 * m01 * m12 + m02 * m22) -
                  cbrtFour * m02 * squares));
    Real v00 =
        Complex((cbrtTerm4 + cbrtTerm3 * cbrtTwo * (2 * m00 - m11 - m22) +
                 cbrtFour * cbrtTerm2 *
                     (m00 * m00 - 2 * m11 * m11 - 9 * m12 * m12 + 5 * m11 * m22 - 2 * m22 * m22 - m00 * (m11 + m22) -
                      2 * squares) +
                 cbrtEight * cbrtTerm * (-2 * m00 + m11 + m22) * squares + cbrtSixteen * squares2) /
                (three * cbrtTerm * cbrtTwo *
                 (cbrtTerm2 * m02 + cbrtTerm * cbrtTwo * (m00 * m02 - 2 * m02 * m11 + 3 * m01 * m12 + m02 * m22) -
                  cbrtFour * m02 * squares)))
            .real();
    Real v01 = Complex((cbrtTerm2 * m01 + cbrtTerm * cbrtTwo * (m00 * m01 + m01 * m11 + 3 * m02 * m12 - 2 * m01 * m22) -
                        cbrtFour * m01 * squares) /
                       (cbrtTerm2 * m02 + cbrtTerm * cbrtTwo * (m00 * m02 - 2 * m02 * m11 + 3 * m01 * m12 + m02 * m22) -
                        cbrtFour * m02 * squares))
                   .real();
    Real N0 = 1 / sqrt(v00 * v00 + v01 * v01 + 1);

    Real v10 =
        Complex(
            (-(cbrtSixteen * cbrtTerm2 *
               (Complex(0, 2) * cbrtFour *
                    (m00 * m00 - 2 * m11 * m11 - 9 * m12 * m12 + 5 * m11 * m22 - 2 * m22 * m22 - m00 * (m11 + m22)) +
                cbrtTerm2 * sqrtThreeMinusI - cbrtTerm * cbrtTwo * (2 * m00 - m11 - m22) * sqrtThreePlusI)) +
             2 * cbrtFour * cbrtTerm * cbrtTwo *
                 (Complex(0, 4) * cbrtTerm + cbrtTwo * (2 * m00 - m11 - m22) * sqrtThreeMinusI) * squares +
             4 * cbrtFour * sqrtThreePlusI * squares2) /
            (3 * cbrtFour * cbrtTerm * cbrtTwo *
             (cbrtFour * cbrtTerm *
                  (Complex(0, -2) * cbrtTwo * (m00 * m02 - 2 * m02 * m11 + 3 * m01 * m12 + m02 * m22) +
                   cbrtTerm * m02 * sqrtThreePlusI) +
              2 * cbrtTwo * m02 * sqrtThreeMinusI * squares)))
            .real();

    Real v11 = Complex((cbrtFour * cbrtTerm *
                            (Complex(0, -2) * cbrtTwo * (m00 * m01 + m01 * m11 + 3 * m02 * m12 - 2 * m01 * m22) +
                             cbrtTerm * m01 * sqrtThreePlusI) +
                        2 * cbrtTwo * m01 * sqrtThreeMinusI * squares) /
                       (cbrtFour * cbrtTerm *
                            (Complex(0, -2) * cbrtTwo * (m00 * m02 - 2 * m02 * m11 + 3 * m01 * m12 + m02 * m22) +
                             cbrtTerm * m02 * sqrtThreePlusI) +
                        2 * cbrtTwo * m02 * sqrtThreeMinusI * squares))
                   .real();
    Real N1 = 1 / sqrt(v10 * v10 + v11 * v11 + 1);

    Real v20 =
        Complex(
            (-(cbrtSixteen * cbrtTerm2 *
               (Complex(0, -2) * cbrtFour *
                    (m00 * m00 - 2 * m11 * m11 - 9 * m12 * m12 + 5 * m11 * m22 - 2 * m22 * m22 - m00 * (m11 + m22)) -
                cbrtTerm * cbrtTwo * (2 * m00 - m11 - m22) * sqrtThreeMinusI + cbrtTerm2 * sqrtThreePlusI)) +
             2 * cbrtFour * cbrtTerm * cbrtTwo *
                 (Complex(0, -4) * cbrtTerm + cbrtTwo * (2 * m00 - m11 - m22) * sqrtThreePlusI) * squares +
             4 * cbrtFour * sqrtThreeMinusI * squares2) /
            (3 * cbrtFour * cbrtTerm * cbrtTwo *
             (cbrtFour * cbrtTerm *
                  (Complex(0, 2) * cbrtTwo * (m00 * m02 - 2 * m02 * m11 + 3 * m01 * m12 + m02 * m22) +
                   cbrtTerm * m02 * sqrtThreeMinusI) +
              2 * cbrtTwo * m02 * sqrtThreePlusI * squares)))
            .real();

    Real v21 = Complex((cbrtFour * cbrtTerm *
                            (Complex(0, 2) * cbrtTwo * (m00 * m01 + m01 * m11 + 3 * m02 * m12 - 2 * m01 * m22) +
                             cbrtTerm * m01 * sqrtThreeMinusI) +
                        2 * cbrtTwo * m01 * sqrtThreePlusI * squares) /
                       (cbrtFour * cbrtTerm *
                            (Complex(0, 2) * cbrtTwo * (m00 * m02 - 2 * m02 * m11 + 3 * m01 * m12 + m02 * m22) +
                             cbrtTerm * m02 * sqrtThreeMinusI) +
                        2 * cbrtTwo * m02 * sqrtThreePlusI * squares))
                   .real();
    Real N2 = 1 / sqrt(v20 * v20 + v21 * v21 + 1);

    val[0] = e0;
    val[1] = e1;
    val[2] = e2;
    // Store output eigenvectors in columns to be consistent with existing
    // diagonalizers used by cpptraj.
    printf("cbrtTerm= %g+%gi  squares= %g  N0= %g  v0denom= %g+%gi v00= %g  v01= %g  N1= %g  v10= %g  v11= %g  N2= %g  v20= %g  v21= %g\n",
            cbrtTerm.real(), cbrtTerm.imag(), squares, N0, v0denom.real(), v0denom.imag(), v00, v01, N1, v10, v11, N2, v20, v21);
    mat[0] = N0 * v00;
    mat[3] = N0 * v01;
    mat[6] = N0;
    mat[1] = N1 * v10;
    mat[4] = N1 * v11;
    mat[7] = N1;
    mat[2] = N2 * v20;
    mat[5] = N2 * v21;
    mat[8] = N2;
}

/// Diagonalize input 3x3 symmetric matrix. Sort by eigenvalues.
template <typename Real>
void Diag_3x3_Sort(Real* M_, Real* Eval) {
  Diag_3x3<Real>( M_, Eval );

  int i1_ = 0;
  int i2_ = 0;
  int i3_ = 0;
  if (Eval[0] > Eval[1] && Eval[0] > Eval[2]) { // 0 is max
    if (Eval[1] > Eval[2]) {
      i1_ = 0; i2_ = 1; i3_ = 2;
    } else {
      i1_ = 0; i2_ = 2; i3_ = 1;
    }
  } else if (Eval[1] > Eval[0] && Eval[1] > Eval[2]) { // 1 is max
    if (Eval[0] > Eval[2]) {
      i1_ = 1; i2_ = 0; i3_ = 2;
    } else {
      i1_ = 1; i2_ = 2; i3_ = 0;
    }
  } else if (Eval[0] > Eval[1]) { // 2 is max
    i1_ = 2; i2_ = 0; i3_ = 1;
  } else {
    i1_ = 2; i2_ = 1; i3_ = 0;
  }

  Real Evec[9];
  for (unsigned int i = 0; i < 9; i++)
    Evec[i] = M_[i];
  // Swap Eigenvectors - place them in rows
  M_[0] = Evec[i1_  ];
  M_[1] = Evec[i1_+3];
  M_[2] = Evec[i1_+6];

  M_[3] = Evec[i2_  ];
  M_[4] = Evec[i2_+3];
  M_[5] = Evec[i2_+6];

  M_[6] = Evec[i3_  ];
  M_[7] = Evec[i3_+3];
  M_[8] = Evec[i3_+6];
  
  // Swap eigenvalues
  Real EvalOut[3];
  EvalOut[0] = Eval[i1_];
  EvalOut[1] = Eval[i2_];
  EvalOut[2] = Eval[i3_];
  for (unsigned int i = 0; i < 3; i++)
    Eval[i] = EvalOut[i];
}
  
/// Determine only eigenvalues from diagonalizing input 3x3 symmetric matrix.
template <typename Real>
void Evals_3x3(Real const* mat, Real* val) {
    typedef std::complex<Real> Complex;
    // Real, symmetric matrix version (upper triangle input)
    Real m00 = mat[0];
    Real m01 = mat[1];
    Real m02 = mat[2];
    Real m11 = mat[4];
    Real m12 = mat[5];
    Real m22 = mat[8];

    // Intermediates
    const Real one = static_cast<Real>(1);
    const Real two = static_cast<Real>(2);
    const Real three = static_cast<Real>(3);
    const Real four = static_cast<Real>(4);
    //const Real six = static_cast<Real>(6);
    const Real cbrtTwo = cbrt(two);
    const Real cbrtFour = cbrt(four);
    const Real cbrtEight = cbrtTwo * cbrtFour;
    const Real cbrtSixteen = cbrtFour * cbrtFour;
    const Real sqrtThree = sqrt(three);
    const Real oneSixth = one / 6;
    const Real oneThird = one / three;
    const Complex sqrtThreeMinusI(sqrtThree, -1);
    const Complex sqrtThreePlusI(sqrtThree, 1);

    Real trace = m00 + m11 + m22;
    Real squares = -m00 * m00 - 3 * m01 * m01 - 3 * m02 * m02 + m00 * m11 - m11 * m11 - 3 * m12 * m12 + m00 * m22 +
                   m11 * m22 - m22 * m22;
    Real squares2 = squares * squares;
    Real cubes = 2 * m00 * m00 * m00 + 9 * m00 * m01 * m01 + 9 * m00 * m02 * m02 - 3 * m00 * m00 * m11 +
                 9 * m01 * m01 * m11 - 18 * m02 * m02 * m11 - 3 * m00 * m11 * m11 + 2 * m11 * m11 * m11 +
                 54 * m01 * m02 * m12 - 18 * m00 * m12 * m12 + 9 * m11 * m12 * m12 - 3 * m00 * m00 * m22 -
                 18 * m01 * m01 * m22 + 9 * m02 * m02 * m22 + 12 * m00 * m11 * m22 - 3 * m11 * m11 * m22 +
                 9 * m12 * m12 * m22 - 3 * m00 * m22 * m22 - 3 * m11 * m22 * m22 + 2 * m22 * m22 * m22;
    Complex cbrtTerm = pow(cubes + sqrt(Complex(cubes * cubes + 4 * squares2 * squares, 0)), oneThird);
    Complex cbrtTerm2 = cbrtTerm * cbrtTerm;
    Complex cbrtTerm3 = cbrtTerm2 * cbrtTerm;
    Complex cbrtTerm4 = cbrtTerm3 * cbrtTerm;

    // Eigenvalues
    Real e0 = Complex((trace + cbrtTerm / cbrtTwo - (squares * cbrtTwo) / cbrtTerm) * oneThird).real();
    Real e1 = Complex(((Complex(0, 1) * cbrtTerm * sqrtThreePlusI) / cbrtTwo +
                       (two * (squares + Complex(0, 1) * sqrtThree * squares)) / (cbrtFour * cbrtTerm) + 2 * trace) *
                      oneSixth)
                  .real();
    Real e2 = Complex((-((cbrtTerm + Complex(0, 1) * sqrtThree * cbrtTerm) / cbrtTwo) +
                       (two * (squares - Complex(0, 1) * sqrtThree * squares)) / (cbrtFour * cbrtTerm) + 2 * trace) *
                      oneSixth)
                  .real();

    val[0] = e0;
    val[1] = e1;
    val[2] = e2;
}

/// Determine only eigenvalues from diagonalizing input 3x3 symmetric matrix. Sort.
template <typename Real>
void Evals_3x3_Sort(Real const* mat, Real* val) {
  Evals_3x3<Real>( mat, val );

  Real EvalOut[3];
  if (val[0] > val[1] && val[0] > val[2]) { // 0 is max
    EvalOut[0] = val[0];
    if (val[1] > val[2]) {
      EvalOut[1] = val[1]; EvalOut[2] = val[2];
    } else {
      EvalOut[1] = val[2]; EvalOut[2] = val[1];
    }
  } else if (val[1] > val[0] && val[1] > val[2]) { // 1 is max
    EvalOut[0] = val[1];
    if (val[0] > val[2]) {
      EvalOut[1] = val[0]; EvalOut[2] = val[2];
    } else {
      EvalOut[1] = val[2]; EvalOut[2] = val[0];
    }
  } else if (val[0] > val[1]) { // 2 is max
    EvalOut[0] = val[2]; EvalOut[1] = val[0]; EvalOut[2] = val[1];
  } else {
    EvalOut[0] = val[2]; EvalOut[1] = val[1]; EvalOut[2] = val[0];
  }
  for (unsigned int i = 0; i < 3; i++)
    val[i] = EvalOut[i];
}
#endif
