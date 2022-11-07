#ifndef INC_TOP_BONDTYPE_H
#define INC_TOP_BONDTYPE_H
#include <vector>
#include "NbodyTerm.h"
namespace Cpptraj {
namespace Top {
/// Hold bonded atom indices and parameter index
class BondType : public NbodyTerm {
  public:
    /// CONSTRUCTOR
    BondType() { at_[0]= -1; at_[1] = -1; }
    /// CONSTRUCTOR - 2 atom indices
    BondType(int a1, int a2) { at_[0] = a1; at_[1] = a2; }
    /// CONSTRUCTOR - 2 atom indices, parameter index
    BondType(int a1, int a2, int idx) : NbodyTerm(idx) { at_[0] = a1; at_[1] = a2; }
    /// COPY CONSTRUCTOR
    BondType(BondType const& rhs) : NbodyTerm(rhs) { at_[0] = rhs.at_[0]; at_[1] = rhs.at_[1]; }
    /// ASSIGNMENT
    BondType& operator=(BondType const& rhs) {
      if (this == &rhs) return *this;
      SetIdx(rhs.Idx());
      at_[0] = rhs.at_[0];
      at_[1] = rhs.at_[1];
      return *this;
    }

    int A1()  const { return at_[0];  }
    int A2()  const { return at_[1];  }

    int At(int idx) const { return at_[idx]; }
    int Nat() const { return 2; }
  private:
    int at_[2];
};
/// Array of bonds
typedef std::vector<BondType> BondArray;
}
}
#endif
