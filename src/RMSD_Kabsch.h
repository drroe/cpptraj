#ifndef INC_RMSD_KABSCH_H
#define INC_RMSD_KABSCH_H
#include "Constants.h"
#include "Vec3.h"
#include "Matrix_3x3.h"
#include "Frame.h"
#include "Jacobi.h"
namespace Cpptraj {

template <typename T>
void normalize(T* vIn) {
  T b = 1.0 / sqrt(vIn[0]*vIn[0] + vIn[1]*vIn[1] + vIn[2]*vIn[2]);
  vIn[0] *= b;
  vIn[1] *= b;
  vIn[2] *= b;
}

template <class T>
class RMSD_Kabsch {
  public:
    /// Error statuses
    enum ErrType { OK = 0, ERR_DIV_BY_ZERO };
    /// CONSTRUCTOR
    RMSD_Kabsch() : err_(OK), useMass_(false) {}
    /// Calculate RMSD to already-centered reference
    T RMSD_CenteredRef(Vec3& tgtTrans, Matrix_3x3& Rot, Frame const& Tgt, AtomMask const& maskIn)
    {
      // Rotation will occur around geometric center/center of mass
      T Trans[3];
      Trans[0] = 0;
      Trans[1] = 0;
      Trans[2] = 0;
      T total_mass = 0;
      std::vector<T> selected_mass;
      selected_mass.reserve( maskIn.Nselected() );
      std::vector<T> selected_tgt;
      selected_tgt.reserve( maskIn.Nselected() * 3 );
      if (useMass_) {
        for (AtomMask::const_iterator at = maskIn.begin(); at != maskIn.end(); ++at) {
          total_mass += Tgt.Mass(*at);
          selected_mass.push_back( Tgt.Mass(*at) );
          const double* xyz = Tgt.XYZ( *at );
          Trans[0] += (T)xyz[0] * Tgt.Mass(*at);
          Trans[1] += (T)xyz[1] * Tgt.Mass(*at);
          Trans[2] += (T)xyz[2] * Tgt.Mass(*at);
          selected_tgt.push_back( (T)xyz[0] );
          selected_tgt.push_back( (T)xyz[1] );
          selected_tgt.push_back( (T)xyz[2] );
        }
      } else {
        total_mass = (T)maskIn.Nselected();
        selected_mass.assign( maskIn.Nselected(), 1 );
        for (AtomMask::const_iterator at = maskIn.begin(); at != maskIn.end(); ++at) {
          const double* xyz = Tgt.XYZ( *at );
          Trans[0] += (T)xyz[0];
          Trans[1] += (T)xyz[1];
          Trans[2] += (T)xyz[2];
          selected_tgt.push_back( (T)xyz[0] );
          selected_tgt.push_back( (T)xyz[1] );
          selected_tgt.push_back( (T)xyz[2] );
        }
      }
      if ( total_mass < Constants::SMALL ) {
        err_ = ERR_DIV_BY_ZERO;
        return -1;
      }
      Trans[0] /= total_mass;
      Trans[1] /= total_mass;
      Trans[2] /= total_mass;
      Trans[0] = -Trans[0];
      Trans[1] = -Trans[1];
      Trans[2] = -Trans[2];

      // Use Kabsch algorithm to calculate optimum rotation matrix.
      // U = [(RtR)^.5][R^-1]
      T mwss = 0.0;
      T rot[9];
      for (unsigned int idx = 0; idx != 9; idx++)
        rot[idx] = 0;
      // Calculate covariance matrix of Coords and Reference (R = Xt * Ref)
      typename std::vector<T>::const_iterator mass = selected_mass.begin();
      for (unsigned int i = 0; i < selected_tgt.size(); i += 3, ++mass)
      {
        T xt = selected_tgt[i  ] + Trans[0];
        T yt = selected_tgt[i+1] + Trans[1];
        T zt = selected_tgt[i+2] + Trans[2];
        T xr = selected_ref_[i  ];
        T yr = selected_ref_[i+1];
        T zr = selected_ref_[i+2];
        T atom_mass = *mass;
        mwss += atom_mass * ( (xt*xt)+(yt*yt)+(zt*zt)+(xr*xr)+(yr*yr)+(zr*zr) );
        // Calculate the Kabsch matrix: R = (rij) = Sum(yni*xnj)
        rot[0] += atom_mass*xt*xr;
        rot[1] += atom_mass*xt*yr;
        rot[2] += atom_mass*xt*zr;

        rot[3] += atom_mass*yt*xr;
        rot[4] += atom_mass*yt*yr;
        rot[5] += atom_mass*yt*zr;

        rot[6] += atom_mass*zt*xr;
        rot[7] += atom_mass*zt*yr;
        rot[8] += atom_mass*zt*zr;
      }
      mwss *= 0.5;    // E0 = 0.5*Sum(xn^2+yn^2)
      // Caclulate Kabsch matrix multiplied by its transpose: RtR
      T Evector[9];
      Evector[0] = rot[0]*rot[0] + rot[1]*rot[1] + rot[2]*rot[2];
      Evector[1] = rot[0]*rot[3] + rot[1]*rot[4] + rot[2]*rot[5];
      Evector[2] = rot[0]*rot[6] + rot[1]*rot[7] + rot[2]*rot[8];
      Evector[3] = rot[3]*rot[0] + rot[4]*rot[1] + rot[5]*rot[2];
      Evector[4] = rot[3]*rot[3] + rot[4]*rot[4] + rot[5]*rot[5];
      Evector[5] = rot[3]*rot[6] + rot[4]*rot[7] + rot[5]*rot[8];
      Evector[6] = rot[6]*rot[0] + rot[7]*rot[1] + rot[8]*rot[2];
      Evector[7] = rot[6]*rot[3] + rot[7]*rot[4] + rot[8]*rot[5];
      Evector[8] = rot[6]*rot[6] + rot[7]*rot[7] + rot[8]*rot[8];
      // Diagonalize
      T Eigenvalue[3];
      if ( Cpptraj::Jacobi::Diagonalize_Sort<T>( Eigenvalue, Evector ) ) {
        // FIXME report error
        return 0;
      }
      // a3 = a1 x a2
      Evector[6] = (Evector[1]*Evector[5]) - (Evector[2]*Evector[4]); 
      Evector[7] = (Evector[2]*Evector[3]) - (Evector[0]*Evector[5]); 
      Evector[8] = (Evector[0]*Evector[4]) - (Evector[1]*Evector[3]);
      // Evector dot transpose rot: b = R . ak
      T b[9];
      b[0] = Evector[0]*rot[0] + Evector[1]*rot[3] + Evector[2]*rot[6];
      b[1] = Evector[0]*rot[1] + Evector[1]*rot[4] + Evector[2]*rot[7];
      b[2] = Evector[0]*rot[2] + Evector[1]*rot[5] + Evector[2]*rot[8];
      normalize<T>(b);
      b[3] = Evector[3]*rot[0] + Evector[4]*rot[3] + Evector[5]*rot[6];
      b[4] = Evector[3]*rot[1] + Evector[4]*rot[4] + Evector[5]*rot[7];
      b[5] = Evector[3]*rot[2] + Evector[4]*rot[5] + Evector[5]*rot[8];
      normalize<T>(b+3);
      b[6] = Evector[6]*rot[0] + Evector[7]*rot[3] + Evector[8]*rot[6];
      b[7] = Evector[6]*rot[1] + Evector[7]*rot[4] + Evector[8]*rot[7];
      b[8] = Evector[6]*rot[2] + Evector[7]*rot[5] + Evector[8]*rot[8];
      normalize<T>(b+6);
      // b3 = b1 x b2
      T cp[3];
      cp[0] = (b[1]*b[5]) - (b[2]*b[4]); 
      cp[1] = (b[2]*b[3]) - (b[0]*b[5]); 
      cp[2] = (b[0]*b[4]) - (b[1]*b[3]);
      T sig3;
      if ( (cp[0]*b[6] + cp[1]*b[7] + cp[2]*b[8]) < 0.0 )
        sig3 = -1.0;
      else
        sig3 = 1.0;
      b[6] = cp[0];
      b[7] = cp[1];
      b[8] = cp[2];
      // U has the best rotation
      T U[9]; 
      U[0] = (Evector[0]*b[0]) + (Evector[3]*b[3]) + (Evector[6]*b[6]);  
      U[1] = (Evector[1]*b[0]) + (Evector[4]*b[3]) + (Evector[7]*b[6]);
      U[2] = (Evector[2]*b[0]) + (Evector[5]*b[3]) + (Evector[8]*b[6]);

      U[3] = (Evector[0]*b[1]) + (Evector[3]*b[4]) + (Evector[6]*b[7]);
      U[4] = (Evector[1]*b[1]) + (Evector[4]*b[4]) + (Evector[7]*b[7]);
      U[5] = (Evector[2]*b[1]) + (Evector[5]*b[4]) + (Evector[8]*b[7]);

      U[6] = (Evector[0]*b[2]) + (Evector[3]*b[5]) + (Evector[6]*b[8]);
      U[7] = (Evector[1]*b[2]) + (Evector[4]*b[5]) + (Evector[7]*b[8]);
      U[8] = (Evector[2]*b[2]) + (Evector[5]*b[5]) + (Evector[8]*b[8]);

      // E=E0-sqrt(mu1)-sqrt(mu2)-sig3*sqrt(mu3) 
      T rms_return = mwss - sqrt(fabs(Eigenvalue[0])) 
                          - sqrt(fabs(Eigenvalue[1]))
                          - (sig3*sqrt(fabs(Eigenvalue[2])));
      if (rms_return<0) {
        //mprinterr("RMS returned is <0 before sqrt, setting to 0 (%f)\n", rms_return);
        rms_return = 0.0;
      } else
        rms_return = sqrt((2.0*rms_return)/total_mass);
      //DEBUG
      //printRotTransInfo(U,Trans);
      //fprintf(stdout,"RMS is %lf\n",rms_return);
      return rms_return;
    }
    // -------------------------------------------
    /// Set reference and reference translation
    int SetRmsRef(Frame const& Ref, AtomMask const& maskIn) {
      refTrans_[0] = 0;
      refTrans_[1] = 0;
      refTrans_[2] = 0;
      T total_mass = 0;
      selected_ref_.reserve( maskIn.Nselected() * 3 );
      if (useMass_) {
        for (AtomMask::const_iterator at = maskIn.begin(); at != maskIn.end(); ++at) {
          total_mass += Ref.Mass(*at); // TODO save ref mass?
          const double* xyz = Ref.XYZ( *at );
          refTrans_[0] += (T)xyz[0] * Ref.Mass(*at);
          refTrans_[1] += (T)xyz[1] * Ref.Mass(*at);
          refTrans_[2] += (T)xyz[2] * Ref.Mass(*at);
        }
      } else {
        total_mass = (T)maskIn.Nselected();
        for (AtomMask::const_iterator at = maskIn.begin(); at != maskIn.end(); ++at) {
          const double* xyz = Ref.XYZ( *at );
          refTrans_[0] += (T)xyz[0];
          refTrans_[1] += (T)xyz[1];
          refTrans_[2] += (T)xyz[2];
        }
      }
      if ( total_mass < Constants::SMALL ) {
        err_ = ERR_DIV_BY_ZERO;
        return 1;
      }
      refTrans_[0] /= total_mass;
      refTrans_[1] /= total_mass;
      refTrans_[2] /= total_mass;
      for (AtomMask::const_iterator at = maskIn.begin(); at != maskIn.end(); ++at) {
        const double* xyz = Ref.XYZ( *at );
        selected_ref_.push_back( (T)xyz[0] - refTrans_[0] );
        selected_ref_.push_back( (T)xyz[1] - refTrans_[1] );
        selected_ref_.push_back( (T)xyz[2] - refTrans_[2] );
      }
      return 0;
    }
  private:
    ErrType err_;
    bool useMass_; ///< If true, mass weight the rmsd
    std::vector<T> selected_ref_; ///< Hold selected reference atoms, translated to the origin
    T refTrans_[3]; ///< Hold translation from origin to reference
};
} // END namespace Cpptraj
#endif
