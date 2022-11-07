#ifndef INC_TOP_CMAPTYPE_H
#define INC_TOP_CMAPTYPE_H
#include <vector>
#include "NbodyTerm.h"
namespace Cpptraj {
namespace Top {
/// Hold CMAP atom indices (1-2-3-4, 2-3-4-5)
class CmapType : public NbodyTerm {
  public:
    CmapType() { for (int i = 0; i < 5; i++) at_[i] = 0; }
    CmapType(int a1, int a2, int a3, int a4, int a5, int idx) : NbodyTerm(idx)
      { for (int i = 0; i < 5; i++) at_[i] = 0; }
    int A1()     const { return at_[0];   }
    int A2()     const { return at_[1];   }
    int A3()     const { return at_[2];   }
    int A4()     const { return at_[3];   }
    int A5()     const { return at_[4];   }

    int At(int idx) const { return at_[idx]; }
    int Nat() const { return 5; }
  private:
    int at_[5];
};
/// Hold array of cmap terms
typedef std::vector<CmapType> CmapArray;
}
}
#endif
