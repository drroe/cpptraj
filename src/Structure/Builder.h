#ifndef INC_STRUCTURE_BUILDER_H
#define INC_STRUCTURE_BUILDER_H
#include <vector>
class Topology;
class Frame;
class ParameterSet;
namespace Cpptraj {
namespace Structure {
class Zmatrix;
/// Used to attach different topology/frame combos using internal coordinates
class Builder {
    typedef std::vector<bool> Barray;
  public:
    /// CONSTRUCTOR
    Builder();
    /// Combine second fragment into first fragment and bond
    int Combine(Topology&, Frame&, Topology const&, Frame const&, int, int) const;
    /// Model the coordinates around a bond given only some coordinates are known
    int ModelCoordsAroundBond(Frame&, Topology const&, int, int, Zmatrix const*, Zmatrix const*, Barray&) const;
    /// Set debug
    void SetDebug(int d) { debug_ = d; }
    /// Set optional parameter set
    void SetParameters(ParameterSet const*);
  private:

    int debug_;
    ParameterSet const* params_;
};
}
}
#endif
