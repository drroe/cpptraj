#ifndef INC_CPPTRAJ_GPU_H
#define INC_CPPTRAJ_GPU_H
#include "Box.h"
namespace CpptrajGpu {
/// Set the major CUDA compute version number
void SetComputeVersion(int);
/// \return Max block dimensions if using 2D blocks
unsigned int MaxBlockDim_2D();
/// Default floating point type
typedef float FpType;
/// Used to hold box parameters on host
template <class T> class HostBox {
  public:
    /// CONSTRUCTOR: Take box
    HostBox(Box const& box) {
      gpu_box_[0] = (T)box.Param(Box::X);
      gpu_box_[1] = (T)box.Param(Box::Y);
      gpu_box_[2] = (T)box.Param(Box::Z);
      for (int ibox = 0; ibox != 9; ibox++) {
        gpu_ucell_[ibox] = (T)box.UnitCell()[ibox];
        gpu_frac_[ibox]  = (T)box.FracCell()[ibox];
      }
    }
    /// \return Pointer to box lengths
    T const* BoxLengths() const { return gpu_box_; }
    /// \return Pointer to unit cell vecs
    T const* Ucell() const { return gpu_ucell_; }
    /// \return Pointer to fractional cell vecs
    T const* Frac() const { return gpu_frac_; }
  private:
    T gpu_box_[3];
    T gpu_ucell_[9];
    T gpu_frac_[9];
};
} // END namespace CpptrajGpu
#endif
