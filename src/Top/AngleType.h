#ifndef INC_TOP_ANGLETYPE_H
#define INC_TOP_ANGLETYPE_H
#include <vector>
#include "NbodyTerm.h"
namespace Cpptraj {
namespace Top {
/// Hold angle atom indices and parameter index
class AngleType : public NbodyTerm {
  public:
    /// CONSTRUCTOR
    AngleType() { at_[0]= -1; at_[1] = -1; }
    /// CONSTRUCTOR - 3 atom indices
    AngleType(int a1, int a2, int a3) { at_[0] = a1; at_[1] = a2; at_[2] = a3; }
    /// CONSTRUCTOR - 3 atom indices, parameter index
    AngleType(int a1, int a2, int a3, int idx) : NbodyTerm(idx) { at_[0] = a1; at_[1] = a2; at_[2] = a3; }
    /// COPY CONSTRUCTOR
    AngleType(AngleType const& rhs) : NbodyTerm(rhs) { at_[0] = rhs.at_[0]; at_[1] = rhs.at_[1]; at_[2] = rhs.at_[2]; }
    /// ASSIGNMENT
    AngleType& operator=(AngleType const& rhs) {
      if (this == &rhs) return *this;
      SetIdx(rhs.Idx());
      at_[0] = rhs.at_[0];
      at_[1] = rhs.at_[1];
      at_[2] = rhs.at_[2];
      return *this;
    }

    int A1()  const { return at_[0];  }
    int A2()  const { return at_[1];  }
    int A3()  const { return at_[2];  }

    int At(int idx) const { return at_[idx]; }
    int Nat() const { return 3; }
  private:
    int at_[3];
};
/// Array of angles
typedef std::vector<AngleType> AngleArray;
}
}
#endif
