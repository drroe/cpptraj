#ifndef INC_PARAM_HOOKESLAWTYPE_H
#define INC_PARAM_HOOKESLAWTYPE_H
#include <vector>
#include "../FPcompare.h"
namespace Cpptraj {
namespace Param {
using namespace FPcompare;
/// Hold parameters for a Hooke's Law term
class HookesLawType {
  public:
    /// CONSTRUCTOR
    HookesLawType() : rk_(0), req_(0) {}
    /// CONSTRUCTOR - force constant, equilibrium value.
    HookesLawType(double rk, double req) : rk_(rk), req_(req) {}
    /// \return Force constant.
    double Rk()  const { return rk_;  }
    /// \return Equilibrium value (where force is zero).
    double Req() const { return req_; }
    /// Set the force constant.
    void SetRk(double rk)   { rk_ = rk;   }
    /// Set the equilibrium value.
    void SetReq(double req) { req_ = req; }
    /// \return True if given term is same as this one.
    bool operator==(const HookesLawType& rhs) const {
      return ( FEQ<double>(rk_,  rhs.rk_ ) &&
               FEQ<double>(req_, rhs.req_) );
    }
    /// \return True if given term is different than this one.
    bool operator!=(const HookesLawType& rhs) const {
      return ( FNE<double>(rk_,  rhs.rk_ ) ||
               FNE<double>(req_, rhs.req_) );
    }
    /// \return True if this term is less than given term.
    bool operator<(const HookesLawType& rhs) const {
      if (*this != rhs) {
        if (FEQ<double>(rk_, rhs.rk_)) {
          return (req_ < rhs.req_);
        } else return (rk_ < rhs.rk_);
      } else
        return false;
    }
  private:
    double rk_;  ///< Force constant.
    double req_; ///< Equilibrium value.
};
/// Array of Hooke's Law parameters
typedef std::vector<HookesLawType> HookesLawArray;
}
}
#endif
