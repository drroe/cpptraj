#include "FrameDist_RMSD.h"
#include "Constants.h"
#include "CpptrajStdio.h"
#include "Diag_3x3.h"
#include <cmath>

using namespace Cpptraj;

/** Set reference, center on origin */
/*void FrameDist_RMSD::SetRef(Frame const& frameIn, AtomMask const& maskIn)
{
  Ref_.SetupFrame(maskIn.Nselected());
  Ref_.SetCoordinates(frameIn, maskIn);
  refTrans_ = Ref_.CenterOnOrigin( useMass_ );
}*/

/// Normalize given vector
static inline void normalize(double* vIn) {
  double b = 1.0 / sqrt(vIn[0]*vIn[0] + vIn[1]*vIn[1] + vIn[2]*vIn[2]);
  vIn[0] *= b;
  vIn[1] *= b;
  vIn[2] *= b;
} 

/** Calculate best-fit RMSD */
double FrameDist_RMSD::RMSD_Fit(Frame& Tgt,
                                Frame const& Ref_,
                                Matrix_3x3& Rot, Vec3& Trans)
//const
{
  t_total_.Start();
  // Sanity check
  if (Tgt.Natom() != Ref_.Natom()) {
    mprinterr("Internal Error: FrameDist_RMSD::RMSD_Fit(): Tgt #atom (%i) != Ref #atom (%i)\n",
              Tgt.Natom(), Ref_.Natom());
    return -1;
  }

  t_trans_.Start();
  double total_mass = 0.0; //, cp[3], sig3, b[9];
  // Rotation will occur around geometric center/center of mass
  Trans.Zero();
  if (useMass_) {
    for (int iat = 0; iat != Tgt.Natom(); iat++)
    {
      const double* XYZ = Tgt.XYZ(iat);
      double mass = Tgt.Mass(iat);
      total_mass += mass;
      Trans[0] += (XYZ[0] * mass);
      Trans[1] += (XYZ[1] * mass);
      Trans[2] += (XYZ[2] * mass);
    }
  } else {
    total_mass = (double)Tgt.Natom();
    for (int iat = 0; iat != Tgt.Natom(); iat++)
    {
      const double* XYZ = Tgt.XYZ(iat);
      Trans[0] += XYZ[0];
      Trans[1] += XYZ[1];
      Trans[2] += XYZ[2];
    }
  }
  if (total_mass < Constants::SMALL) {
    mprinterr("Error: FrameDist_RMSD::RMSD_Fit(): Divide by zero.\n");
    return -1;
  }
  Trans[0] /= total_mass;
  Trans[1] /= total_mass;
  Trans[2] /= total_mass;
  //mprintf("  FRAME COM1: %f %f %f\n",Trans[0],Trans[1],Trans[2]); //DEBUG

  // Shift to common COM
  Trans.Neg();
//  Tgt.Translate(Trans);
  //for (int i = 0; i < natom_; i++) {
  //  mprinterr("\tSHIFTED FRAME %i: %f %f %f\n",i,X_[i*3],X_[i*3+1],X_[i*3+2]); //DEBUG
  //  mprinterr("\tSHIFTED REF %i  : %f %f %f\n",i,Ref.X_[i*3],Ref.X_[i*3+1],Ref.X_[i*3+2]); //DEBUG
  //}
  t_trans_.Stop();
  t_covar_.Start();
  // Use Kabsch algorithm to calculate optimum rotation matrix.
  // U = [(RtR)^.5][R^-1]
  double mwss = 0.0;
  Matrix_3x3 rot(0.0);
  // Calculate covariance matrix of Coords and Reference (R = Xt * Ref)
  if (useMass_) {
    for (int iat = 0; iat != Tgt.Natom(); iat++)
    {
      const double* txyz = Tgt.XYZ(iat);
      const double* rxyz = Ref_.XYZ(iat);
      double xt = txyz[0] + Trans[0];
      double yt = txyz[1] + Trans[1];
      double zt = txyz[2] + Trans[2];
      double xr = rxyz[0];
      double yr = rxyz[1];
      double zr = rxyz[2];

      // Use atom_mass to hold mass for this atom if specified
      double atom_mass = Tgt.Mass(iat);

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
  } else {
    for (int iat = 0; iat != Tgt.Natom(); iat++)
    {
      const double* txyz = Tgt.XYZ(iat);
      const double* rxyz = Ref_.XYZ(iat);
      double xt = txyz[0] + Trans[0];
      double yt = txyz[1] + Trans[1];
      double zt = txyz[2] + Trans[2];
      double xr = rxyz[0];
      double yr = rxyz[1];
      double zr = rxyz[2];

      mwss += ( (xt*xt)+(yt*yt)+(zt*zt)+(xr*xr)+(yr*yr)+(zr*zr) );
      // Calculate the Kabsch matrix: R = (rij) = Sum(yni*xnj)
      rot[0] += xt*xr;
      rot[1] += xt*yr;
      rot[2] += xt*zr;

      rot[3] += yt*xr;
      rot[4] += yt*yr;
      rot[5] += yt*zr;

      rot[6] += zt*xr;
      rot[7] += zt*yr;
      rot[8] += zt*zr;
    }
  }
  mwss *= 0.5;    // E0 = 0.5*Sum(xn^2+yn^2) 
  //DEBUG
  //mprintf("ROT1:\n%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n",
  //          rot[0],rot[1],rot[2],rot[3],rot[4],rot[5],rot[6],rot[7],rot[8]);
  //mprinterr("MWSS: %lf\n",mwss);
  t_covar_.Stop();
  // calculate Kabsch matrix multiplied by its transpose: RtR
  t_matmult_.Start();
  Matrix_3x3 Evector = rot.TransposeMult(rot); 
  t_matmult_.Stop();
  t_diag_.Start();
  // Diagonalize
  Vec3 Eigenvalue;
  //if (Evector.Diagonalize_Sort( Eigenvalue )) return 0;
  Diag_3x3_Sort<double>( Evector.Dptr(), Eigenvalue.Dptr() );
  t_diag_.Stop();
  t_findrot_.Start();
  //Evector.Print("Evector1");
  //Eigenvalue.Print("Evalue1");
  // a3 = a1 x a2
  Evector[6] = (Evector[1]*Evector[5]) - (Evector[2]*Evector[4]); 
  Evector[7] = (Evector[2]*Evector[3]) - (Evector[0]*Evector[5]); 
  Evector[8] = (Evector[0]*Evector[4]) - (Evector[1]*Evector[3]);
  // Evector dot transpose rot: b = R . ak
  double bmat[9];
  bmat[0] = Evector[0]*rot[0] + Evector[1]*rot[3] + Evector[2]*rot[6];
  bmat[1] = Evector[0]*rot[1] + Evector[1]*rot[4] + Evector[2]*rot[7];
  bmat[2] = Evector[0]*rot[2] + Evector[1]*rot[5] + Evector[2]*rot[8];
  normalize(bmat);
  bmat[3] = Evector[3]*rot[0] + Evector[4]*rot[3] + Evector[5]*rot[6];
  bmat[4] = Evector[3]*rot[1] + Evector[4]*rot[4] + Evector[5]*rot[7];
  bmat[5] = Evector[3]*rot[2] + Evector[4]*rot[5] + Evector[5]*rot[8];
  normalize(bmat+3);
  bmat[6] = Evector[6]*rot[0] + Evector[7]*rot[3] + Evector[8]*rot[6];
  bmat[7] = Evector[6]*rot[1] + Evector[7]*rot[4] + Evector[8]*rot[7];
  bmat[8] = Evector[6]*rot[2] + Evector[7]*rot[5] + Evector[8]*rot[8];
  normalize(bmat+6);
  // b3 = b1 x b2
  double cp[3];
  cp[0] = (bmat[1]*bmat[5]) - (bmat[2]*bmat[4]); 
  cp[1] = (bmat[2]*bmat[3]) - (bmat[0]*bmat[5]); 
  cp[2] = (bmat[0]*bmat[4]) - (bmat[1]*bmat[3]);
  double sig3;
  if ( (cp[0]*bmat[6] + cp[1]*bmat[7] + cp[2]*bmat[8]) < 0.0 )
    sig3 = -1.0;
  else
    sig3 = 1.0;
  bmat[6] = cp[0];
  bmat[7] = cp[1];
  bmat[8] = cp[2];
  // Rot has the best rotation 
  Rot[0] = (Evector[0]*bmat[0]) + (Evector[3]*bmat[3]) + (Evector[6]*bmat[6]);  
  Rot[1] = (Evector[1]*bmat[0]) + (Evector[4]*bmat[3]) + (Evector[7]*bmat[6]);
  Rot[2] = (Evector[2]*bmat[0]) + (Evector[5]*bmat[3]) + (Evector[8]*bmat[6]);

  Rot[3] = (Evector[0]*bmat[1]) + (Evector[3]*bmat[4]) + (Evector[6]*bmat[7]);
  Rot[4] = (Evector[1]*bmat[1]) + (Evector[4]*bmat[4]) + (Evector[7]*bmat[7]);
  Rot[5] = (Evector[2]*bmat[1]) + (Evector[5]*bmat[4]) + (Evector[8]*bmat[7]);

  Rot[6] = (Evector[0]*bmat[2]) + (Evector[3]*bmat[5]) + (Evector[6]*bmat[8]);
  Rot[7] = (Evector[1]*bmat[2]) + (Evector[4]*bmat[5]) + (Evector[7]*bmat[8]);
  Rot[8] = (Evector[2]*bmat[2]) + (Evector[5]*bmat[5]) + (Evector[8]*bmat[8]);

  // E=E0-sqrt(mu1)-sqrt(mu2)-sig3*sqrt(mu3) 
  double rms_return = mwss - sqrt(fabs(Eigenvalue[0])) 
                           - sqrt(fabs(Eigenvalue[1]))
                           - (sig3*sqrt(fabs(Eigenvalue[2])));
  if (rms_return<0) {
    //mprinterr("RMS returned is <0 before sqrt, setting to 0 (%f)\n", rms_return);
    rms_return = 0.0;
  } else
    rms_return = sqrt((2.0*rms_return)/total_mass);
  t_findrot_.Stop();
  //DEBUG
  //printRotTransInfo(U,Trans);
  //fprintf(stdout,"RMS is %lf\n",rms_return);
/*  t_move_.Start();
  // Tgt is already the origin. Rotate and translate it.
  for (int i = 0; i < Tgt.size(); i+=3) {
    double* XYZ = Tgt.xAddress() + i;
    double x = XYZ[0] + Trans[0];
    double y = XYZ[1] + Trans[1];
    double z = XYZ[2] + Trans[2];
    XYZ[0] = x*Rot[0] + y*Rot[1] + z*Rot[2] + refTrans_[0];
    XYZ[1] = x*Rot[3] + y*Rot[4] + z*Rot[5] + refTrans_[1];
    XYZ[2] = x*Rot[6] + y*Rot[7] + z*Rot[8] + refTrans_[2];
  }
  t_move_.Stop();*/
  t_total_.Stop();
  return rms_return;
}

void FrameDist_RMSD::PrintTiming() const {
  t_total_.WriteTiming(1, "New RMS total");
  t_trans_.WriteTiming(2, "Trans  ", t_total_.Total());
  t_covar_.WriteTiming(2, "Covar  ", t_total_.Total());
  t_matmult_.WriteTiming(2, "MatMult", t_total_.Total());
  t_diag_.WriteTiming(2, "Diag   ", t_total_.Total());
  t_findrot_.WriteTiming(2, "FindRot", t_total_.Total());
  //t_move_.WriteTiming(2, "Move   ", t_total_.Total());
}
