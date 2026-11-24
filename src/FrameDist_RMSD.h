#ifndef INC_FRAMEDIST_RMSD_H
#define INC_FRAMEDIST_RMSD_H
#include "Frame.h"
namespace Cpptraj {
class FrameDist_RMSD {
  public:
    FrameDist_RMSD() : useMass_(false) {}

    void SetUseMass(bool b) { useMass_ = b; }

//    void SetRef(Frame const&, AtomMask const&);

//    double RMSD_Fit(Frame&, AtomMask const&, Matrix_3x3&, Vec3&) const;
    double RMSD_Fit(Frame&, AtomMask const&, Frame const&, Matrix_3x3&, Vec3&) const;
  private:
//    Frame Ref_;
    Vec3 refTrans_; ///< Translation from origin to original reference center
    bool useMass_;
};
} 
#endif
