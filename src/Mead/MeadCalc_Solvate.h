#ifndef INC_MEAD_MEADCALC_SOLVATE_H
#define INC_MEAD_MEADCALC_SOLVATE_H
#ifdef HAS_MEAD
#include "MeadCalc.h"
class DataSet_3D;
namespace Cpptraj {
namespace Mead {
class MeadOpts;
class MeadGrid;
class MeadCalc_Solvate : public MeadCalc {
  public:
    MeadCalc_Solvate();
    int Solvate(double&, MeadOpts const&, MeadGrid const&, DataSet_3D*) const;
};
}
}
#endif
#endif
