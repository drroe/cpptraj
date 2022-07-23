#ifndef INC_JACOBI_H
#define INC_JACOBI_H
#include <cmath>
namespace Cpptraj {
namespace Jacobi {

/// Max number of iterations to execute Jacobi algorithm
const int MAX_ITERATIONS_ = 50;

#define ROTATE(ARR,MAJ1,MIN1,MAJ2,MIN2) { \
  dg = ARR[MAJ1 + MIN1]; \
  dh = ARR[MAJ2 + MIN2]; \
  ARR[MAJ1 + MIN1] = dg - ds*(dh+dg*tau); \
  ARR[MAJ2 + MIN2] = dh + ds*(dg-dh*tau); }

/** Diagonalize the matrix using Jacobi method. Eigenvectors are stored in 
  * columns. 
  * \param vecD Output eigenvalues.
  */
template <typename T>
int Diagonalize( T* vecD, T* M_ ) 
{
  // Store this matrix
  T mat[9];
  mat[0] = M_[0];
  mat[1] = M_[1];
  mat[2] = M_[2];
  mat[3] = M_[3];
  mat[4] = M_[4];
  mat[5] = M_[5];
  mat[6] = M_[6];
  mat[7] = M_[7];
  mat[8] = M_[8];
  // Create identity matrix
  M_[0] = 1;
  M_[1] = 0;
  M_[2] = 0;
  M_[3] = 0;
  M_[4] = 1;
  M_[5] = 0;
  M_[6] = 0;
  M_[7] = 0;
  M_[8] = 1;
  // Set vectors B and D equal to diagonal of mat. vector Z is 0.
  T vecB[3], vecZ[3];
  vecB[0] = vecD[0] = mat[0]; 
  vecB[1] = vecD[1] = mat[4]; 
  vecB[2] = vecD[2] = mat[8];
  vecZ[0] = 0; 
  vecZ[1] = 0; 
  vecZ[2] = 0;
  // MAIN LOOP
  T tresh = 0;
  //int nrot = 0;
  for (int i = 0; i < MAX_ITERATIONS_; ++i) {
    // sm = SUM of UPPER RIGHT TRIANGLE
    T sm = fabs(mat[1]) + fabs(mat[2]) + fabs(mat[5]);
    if (sm == 0) return 0;
    
    if (i < 3)
      tresh = 0.2 * sm / 9;
    else
      tresh = 0;
    // BEGIN INNER LOOP OVER UPPER RIGHT TRIANGLE
    T dt;
    //int p3 = 0;
    int ip, p3;
    for ( ip = p3 = 0; ip < 2; ++ip, p3+=3) {
      for ( int iq = ip + 1; iq < 3; ++iq ) {
        int midx = p3 + iq;
        T dg = 100 * fabs(mat[midx]);
        if ( i > 3 && fabs(vecD[ip]) + dg == fabs(vecD[ip]) && 
                      fabs(vecD[iq]) + dg == fabs(vecD[iq]) )
        {
          mat[midx] = 0;
        } else if ( fabs(mat[midx]) > tresh) {
          T dh = vecD[iq] - vecD[ip];
          if (fabs(dh) + dg == fabs(dh))
            dt = mat[p3 + iq] / dh;
          else {
            T theta = 0.5 * dh / mat[midx];
            dt = 1.0 / (fabs(theta)+(T)sqrt(1.0+theta*theta));
            if (theta < 0.0)
              dt = -dt;
          }
          T dc = 1.0 / (T)sqrt(1+dt*dt);
          T ds = dt * dc;
          T tau = ds / (1.0+dc);
          dh = dt * mat[midx];
          vecZ[ip] -= dh;
          vecZ[iq] += dh;
          vecD[ip] -= dh;
          vecD[iq] += dh;
          mat[midx] = 0;
          int j, j3;
          for (j=j3=0; j<=ip-1; j++,j3+=3)
            ROTATE(mat,j3,ip,j3,iq)
          for (int j=ip+1; j<=iq-1; j++)
            ROTATE(mat,p3,j,j*3,iq)
          for (int j=iq+1; j<3; j++)
            ROTATE(mat,p3,j,iq*3,j)

          for (j3=0; j3<9; j3+=3)
            ROTATE(M_,j3,ip,j3,iq)

          //++nrot;
        }
      }
    } // END INNER LOOP OVER UPPER RIGHT TRIANGLE
    vecB[0] += vecZ[0];
    vecD[0] = vecB[0];
    vecZ[0] = 0;
    vecB[1] += vecZ[1];
    vecD[1] = vecB[1];
    vecZ[1] = 0;
    vecB[2] += vecZ[2];
    vecD[2] = vecB[2];
    vecZ[2] = 0;
  }
  //mprintf("Too many iterations in routine!\n"); // FIXME
  return 1;
}

/** Diagonalize the matrix and sort eigenvalues/eigenvectors in 
  * descending order. Eigenvectors will be stored in rows,
  * (V0x, V0y, V0z, V1x, ... V2z).
  * \param EvalOut Output eigenvalues.
  */
template <typename T>
int Diagonalize_Sort(T* EvalOut, T* M_) 
{
  T Eval[3];
  if ( Diagonalize<T>( Eval, M_ ) ) 
  {
    //mprintf("Convergence failed.\n"); // FIXME
    return 1; 
  }
  //printMatrix_3x3("Eigenvector Matrix", Evec);
  unsigned int i1_, i2_, i3_;
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
  //mprintf("EIGENVALUE ORDER (0=high, 3=med, 6=low): %i %i %i\n",i1_,i2_,i3_);

  // Swap Eigenvectors - place them in rows
  T Evec[9];
  for (unsigned int idx = 0; idx != 9; idx++)
    Evec[idx] = M_[idx];
  //Matrix_3x3 Evec(*this);
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
  EvalOut[0] = Eval[i1_];
  EvalOut[1] = Eval[i2_];
  EvalOut[2] = Eval[i3_];

  return 0;
}

/**  The jacobi diagonalization procedure can sometimes result
  *  in eigenvectors which when applied to transform the coordinates
  *  result in a a chiral inversion about the Y axis.  This code catches
  *  this case, reversing the offending eigenvectors.
  *  
  *  NOTE: the idea of rotating the coordinate basis vectors came from 
  *  some code posted to the computational chemistry mailing list 
  *  (chemistry@osc) in a summary of methods to perform principal axis 
  *  alignment...
  *
  * It is expected that the eigenvector matrix has eigenvectors in rows.
  */
/*
int Matrix_3x3::jacobiCheckChirality()
{
  Matrix_3x3 points(*this);
  Matrix_3x3 result;
  //points.Print("POINTS"); 

  // rotate vector three into XZ plane
  result.RotationAroundZ( points[2], points[5] ); // Ev0z, Ev1z
  result *= points;
  //result.Print("POINTS1");

  // rotate vector three into Z axis
  points.RotationAroundY( result[2], result[8] ); 
  points *= result;
  //points.Print("POINTS2");

  // rotate vector one into XZ
  result.RotationAroundZ( points[0], points[3] );
  result *= points;
  //result.Print("POINTS3");

  // rotate vector one into X 
  points.RotationAroundY( result[2], result[0] );
  points *= result;
  //points.Print("POINTS4");

  // has Y changed sign? If so, flip Y eigenvector (row 1) 
  if ( points[4] < 0 ) {
    M_[3] = -M_[3];
    M_[4] = -M_[4];
    M_[5] = -M_[5];
    return 1;
  }
  return 0;
}

// Matrix_3x3::Diagonalize_Sort_Chirality()
int Matrix_3x3::Diagonalize_Sort_Chirality(Vec3& EvalOut, int debug)
{
  if ( Diagonalize_Sort( EvalOut ) )
    return 1;

  // Invert eigenvector signs based on ordering to avoid reflections
  if (i1_ == 0 && i2_ == 2 && i3_ == 1) {
    M_[3] = -M_[3];
    M_[4] = -M_[4];
    M_[5] = -M_[5];
  } else if (i1_ == 2 && i2_ == 0 && i3_ == 1) {
    M_[0] = -M_[0];
    M_[1] = -M_[1];
    M_[2] = -M_[2];
    M_[3] = -M_[3];
    M_[4] = -M_[4];
    M_[5] = -M_[5];
    M_[6] = -M_[6];
    M_[7] = -M_[7];
    M_[8] = -M_[8];
  }

  // Flip Y-vector if necessary 
  if (jacobiCheckChirality( ) && debug>0)
    mprintf("Warning: PRINCIPAL: CHECK CHIRALITY: Y eigenvector sign swapped!\n");
   
  return 0;
}
*/
} // END namespace Jacobi
} // END namespace Cpptraj
#endif
