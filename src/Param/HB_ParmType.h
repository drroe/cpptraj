#ifndef INC_PARAM_HBPARMTYPE_H
#define INC_PARAM_HBPARMTYPE_H
#include <vector>
namespace Cpptraj {
namespace Param {
/// Hold LJ 10-12 hbond params
class HB_ParmType {
  public:
    HB_ParmType() : asol_(0), bsol_(0), hbcut_(0) {}
    HB_ParmType(double a, double b, double c) :
                           asol_(a), bsol_(b), hbcut_(c) {}
    inline double Asol()  const { return asol_;  }
    inline double Bsol()  const { return bsol_;  }
    inline double HBcut() const { return hbcut_; }
    void SetAsol(double a)  { asol_ = a;  }
    void SetBsol(double b)  { bsol_ = b;  }
    void SetHBcut(double h) { hbcut_ = h; }
  private:
    double asol_;
    double bsol_;
    double hbcut_;
};
typedef std::vector<HB_ParmType> HB_ParmArray;
}
}
#endif
