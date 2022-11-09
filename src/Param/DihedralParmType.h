#ifndef INC_PARAM_DIHEDRALPARMTYPE_H
#define INC_PARAM_DIHEDRALPARMTYPE_H
#include <vector>
#include "../FPcompare.h"
namespace Cpptraj {
namespace Param {
using namespace FPcompare;
/// Hold dihedral parameters (Fourier series)
class DihedralParmType {
  public:
    /// CONSTRUCTOR
    DihedralParmType() : pk_(0), pn_(0), phase_(0), scee_(0), scnb_(0) {}
    /// CONSTRUCTOR - force constant, periodicity, phase, elec scale, nonbond scale
    DihedralParmType(double k, double n, double p, double e, double b) :
                         pk_(k), pn_(n), phase_(p), scee_(e), scnb_(b) {}
    /// CONSTRUCTOR - force constant, phase
    DihedralParmType(double k, double p) :
                         pk_(k), pn_(0), phase_(p), scee_(0), scnb_(0) {}
    /// CONSTRUCTOR - force constant, periodicity, phase
    DihedralParmType(double k, double n, double p) :
                         pk_(k), pn_(n), phase_(p), scee_(0), scnb_(0) {}
    /// \return Force constant
    double Pk()    const { return pk_;    }
    double& Pk()         { return pk_;    }
    double Pn()    const { return pn_;    }
    double Phase() const { return phase_; }
    double SCEE()  const { return scee_;  }
    double SCNB()  const { return scnb_;  }

    void SetPk(double k)        { pk_ = k;       }
    void SetPn(double n)        { pn_ = n;       }
    void SetPhase(double p)     { phase_ = p;    }
    void SetSCEE(double s)      { scee_ = s;     }
    void SetSCNB(double s)      { scnb_ = s;     }
    bool operator==(DihedralParmType const& rhs) const {
      return ( FEQ(pk_, rhs.pk_) &&
               FEQ(pn_, rhs.pn_) &&
               FEQ(phase_, rhs.phase_) &&
               FEQ(scee_, rhs.scee_) &&
               FEQ(scnb_, rhs.scnb_) );
    }
    bool operator<(DihedralParmType const& rhs) const {
      if (pk_ == rhs.pk_) {
        if (pn_ == rhs.pn_) {
          if (phase_ == rhs.phase_) {
            if (scee_ == rhs.scee_) {
              return (scnb_ < rhs.scnb_);
            } else return (scee_ < rhs.scee_);
          } else return (phase_ < rhs.phase_);
        } else return (pn_ < rhs.pn_);
      } else return (pk_ < rhs.pk_);
    }
  private:
    double pk_;
    double pn_;
    double phase_;
    double scee_;
    double scnb_;
};
/// Array of dihedral parameters
typedef std::vector<DihedralParmType> DihedralParmArray;
}
}
#endif
