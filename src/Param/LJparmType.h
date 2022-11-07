#ifndef INC_PARAM_LJPARMTYPE_H
#define INC_PARAM_LJPARMTYPE_H
#include "../FPcompare.h"
#include <vector>
#include <cmath> // sqrt
namespace Cpptraj {
namespace Param {
/// Hold Lennard-Jones radius and well-depth
class LJparmType {
  public:
    LJparmType() : radius_(0.0), depth_(0.0) {}
    LJparmType(double r, double d) : radius_(r), depth_(d) {}
    double Radius() const { return radius_; }
    double Depth()  const { return depth_;  }
    void SetRadius(double r) { radius_ = r; }
    void SetDepth(double d)  { depth_ = d;  }
    /// \return True if radius and well depth match
    bool operator==(LJparmType const& rhs) const {
      return ( FEQ(radius_, rhs.radius_) &&
               FEQ(depth_,  rhs.depth_) );
    }
    bool operator!=(LJparmType const& rhs) const {
      return ( FNE(radius_, rhs.radius_) ||
               FNE(depth_,  rhs.depth_) );
    }
    /// \return true if radius and well depth are less in that order
    bool operator<(LJparmType const& rhs) const {
      if (*this != rhs) {
        if (FEQ(radius_, rhs.radius_))
          return (depth_ < rhs.depth_);
        else
          return (radius_ < rhs.radius_);
      } else
        return false;
    }
    /// \return LJ A/B params using Lorentz-Berthelot rules.
    NonbondType Combine_LB(LJparmType const& rhs) const {
      double dR = radius_ + rhs.radius_;
      double dE = sqrt( depth_ * rhs.depth_ );
      double dR2 = dR * dR;
      double dR6 = dR2 * dR2 * dR2;
      double dER6 = dE * dR6;
      return NonbondType( dER6*dR6, 2.0*dER6 );
    }
  private:
    double radius_;
    double depth_;
};
typedef std::vector<LJparmType> LJparmArray;
}
}
#endif
