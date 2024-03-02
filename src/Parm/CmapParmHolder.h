#ifndef INC_PARM_CMAPPARMHOLDER_H
#define INC_PARM_CMAPPARMHOLDER_H
class CmapGridType;
#include <vector>
namespace Cpptraj {
namespace Parm {
/// Hold CMAP terms
class CmapParmHolder {
  public:
    CmapParmHolder();
  private:
    typedef std::vector<CmapGridType> Carray;
    Carray CMAP_;
};
}
}
#endif
