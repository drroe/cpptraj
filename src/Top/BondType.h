#ifndef INC_TOP_BONDTYPE_H
#define INC_TOP_BONDTYPE_H
#include <vector>
#include "NbodyTerm.h"
namespace Cpptraj {
namespace Top {
/// Hold bonded atom indices and parameter index
class BondType : public NbodyTerm {
  public:
    BondType() { at_[0]= 0; at_[1] = 0; }
    BondType(int a1, int a2, int idx) : NbodyTerm(idx) { at_[0] = a1; at_[1] = a2; }
    BondType(BondType const& rhs) : NbodyTerm(rhs) { at_[0] = rhs.at_[0]; at_[1] = rhs.at_[1]; }
    inline int A1()  const { return at_[0];  }
    inline int A2()  const { return at_[1];  }
    int At(int idx) const { return at_[idx]; }
    int Nat() const { return 2; }
/*    bool operator<(const BondType& rhs) const {
      if (a1_ == rhs.a1_) {
        return (a2_ < rhs.a2_);
      } else return (a1_ < rhs.a1_);
    }*/
  private:
    int at_[2];
};
typedef std::vector<BondType> BondArray;
}
}
#endif
