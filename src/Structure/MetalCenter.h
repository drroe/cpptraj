#ifndef INC_STRUCTURE_METALCENTER_H
#define INC_STRUCTURE_METALCENTER_H
#include "../AtomMask.h"
class ArgList;
class Frame;
class Topology;
namespace Cpptraj {
namespace Structure {
/// Used for finding and preparing metal centers
class MetalCenter {
  public:
    /// CONSTRUCTOR
    MetalCenter();
    /// Init with arguments
    int InitMetalCenters(ArgList&, int);
    /// Find metal centers
    int FindMetalCenters(Topology const&, Frame const&);
    /// Print Info to stdout
    void PrintMetalCenterInfo() const;
  private:
    AtomMask metalMask_; ///< Mask containing potential metal centers
    AtomMask coordAtomMask_; ///< Mask containing potential coordinating atoms
    double dcut2_;           ///< Distance cutoff in Ang^2
    int debug_;
};
}
}
#endif
