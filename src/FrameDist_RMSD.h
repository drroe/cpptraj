#ifndef INC_FRAMEDIST_RMSD_H
#define INC_FRAMEDIST_RMSD_H
#include "Frame.h"
#include "Timer.h"
namespace Cpptraj {
class FrameDist_RMSD {
  public:
    FrameDist_RMSD() : useMass_(false) {}

    void SetUseMass(bool b) { useMass_ = b; }

//    void SetRef(Frame const&, AtomMask const&);

//    double RMSD_Fit(Frame&, AtomMask const&, Matrix_3x3&, Vec3&) const;
//    double RMSD_Fit(Frame&, AtomMask const&, Frame const&, Matrix_3x3&, Vec3&);// const;
    double RMSD_Fit(Frame&, Frame const&, Matrix_3x3&, Vec3&);// const;

    void PrintTiming() const;
  private:
//    Frame Ref_;
    Vec3 refTrans_; ///< Translation from origin to original reference center
    bool useMass_;

    Timer t_total_;
    Timer t_trans_;
    Timer t_covar_;
    Timer t_matmult_;
    Timer t_diag_;
    Timer t_findrot_;
    //Timer t_move_;
};
} 
#endif
