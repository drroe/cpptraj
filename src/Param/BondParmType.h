#ifndef INC_PARAM_BONDPARMTYPE_H
#define INC_PARAM_BONDPARMTYPE_H
#include <vector>
#include "../FPcompare.h"
namespace Cpptraj {
namespace Param {
using namespace FPcompare;
/// Hold bond parameters
class BondParmType {
  public:
    BondParmType() : rk_(0), req_(0) {}
    BondParmType(double rk, double req) : rk_(rk), req_(req) {}
    inline double Rk()  const { return rk_;  }
    inline double Req() const { return req_; }
    inline void SetRk(double rk)   { rk_ = rk;   }
    inline void SetReq(double req) { req_ = req; }
    bool operator==(const BondParmType& rhs) const {
      return ( FEQ<double>(rk_,  rhs.rk_ ) &&
               FEQ<double>(req_, rhs.req_) );
    }
    bool operator!=(const BondParmType& rhs) const {
      return ( FNE<double>(rk_,  rhs.rk_ ) ||
               FNE<double>(req_, rhs.req_) );
    }
    bool operator<(const BondParmType& rhs) const {
      if (*this != rhs) {
        if (FEQ<double>(rk_, rhs.rk_)) {
          return (req_ < rhs.req_);
        } else return (rk_ < rhs.rk_);
      } else
        return false;
    }
  private:
    double rk_;
    double req_;
};
typedef std::vector<BondParmType> BondParmArray;
}
}
#endif
