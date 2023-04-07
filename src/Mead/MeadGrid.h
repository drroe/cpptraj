#ifndef INC_MEAD_MEADGRID_H
#define INC_MEAD_MEADGRID_H
#include <vector>
#include "../Vec3.h"
// MEAD fwd declares
class FinDiffMethod;
namespace Cpptraj {
namespace Mead {
/// Wrapper around the MEAD FinDiffMethod (grid) class
class MeadGrid {
  public:
    /// Grid centering modes. 
    enum Center_Mode { C_ON_ORIGIN = 0,    ///< Center on coordinate origin
                       C_ON_CENT_OF_INTR,  ///< Center on center of interest (site)
                       C_ON_GEOM_CENT,     ///< Center on molecule center
                       C_SPECIFIED         ///< Center on specified coordinate
                     };
    /// \return Character string corresponding to Center_Mode
    static const char* Center_ModeStr(Center_Mode);

    MeadGrid();
    /// Add grid level to the FDM object with explicit centering
    int AddGrid(int, float, Vec3 const&);
    /// Add grid level to the FDM object with centering type
    int AddGrid(int, float, Center_Mode);
    /// Print to stdout
    void Print() const;
  private:
    static const char* Center_ModeStr_[]; ///< Hold strings corresponding to Center_Mode

    FinDiffMethod* fdm_;

    class GridOpt {
      public:
      GridOpt(int n, float s, Center_Mode c) : npoints_(n), spacing_(s), mode_(c) {}
      GridOpt(int n, float s, Vec3 const& v) : npoints_(n), spacing_(s), mode_(C_SPECIFIED), coord_(v) {}
      int npoints_;      ///< Number of points along an edge of grid
      float spacing_;    ///< Spacing between grid points
      Center_Mode mode_; ///< Grid centering mode
      Vec3 coord_;       ///< Coordinates if centering mode is C_SPECIFIED
    };
    std::vector<GridOpt> levelOpts_; ///< Grid options for each level of grid
};
}
}
#endif
