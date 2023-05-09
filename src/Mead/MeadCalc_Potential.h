#ifndef INC_MEAD_MEADCALC_POTENTIAL_H
#define INC_MEAD_MEADCALC_POTENTIAL_H
#ifdef HAS_MEAD
#include "MeadCalc.h"
#include <vector>
class DataSet_Vector_Scalar;
class Vec3;
namespace Cpptraj {
namespace Mead {
class MeadOpts;
class MeadGrid;

class MeadCalc_Potential : public MeadCalc {
  public:
    MeadCalc_Potential();

    int Potential(DataSet_Vector_Scalar&, MeadOpts const&,
                  MeadGrid const&, std::vector<Vec3> const&) const;
};
}
}
#endif
#endif
