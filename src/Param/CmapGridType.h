#ifndef INC_PARAM_CMAPGRIDTYPE_H
#define INC_PARAM_CMAPGRIDTYPE_H
#include <vector>
namespace Cpptraj {
namespace Param {
/// Hold CMAP grid parameters
class CmapGridType {
  public:
    CmapGridType() : resolution_(0) {}
    CmapGridType(unsigned int r) : resolution_(r), grid_((size_t)r*(size_t)r, 0.0) {}
    int Resolution()                  const { return resolution_;       }
    std::vector<double> const& Grid() const { return grid_;             }
    int Size()                        const { return (int)grid_.size(); }
    void SetGridPt(int idx, double d)              { grid_[idx] = d;           }
  private:
    int resolution_;           ///< Number of steps along each phi/psi CMAP axis
    std::vector<double> grid_; ///< CMAP grid (size is resolution_*resolution_), column major order
};
}
}
#endif
