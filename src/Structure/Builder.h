#ifndef INC_STRUCTURE_BUILDER_H
#define INC_STRUCTURE_BUILDER_H
#include <vector>
class Topology;
class Frame;
class ParameterSet;
namespace Cpptraj {
namespace Structure {
class BuildAtom;
class Zmatrix;
/// Used to attach different topology/frame combos using internal coordinates
class Builder {
    typedef std::vector<bool> Barray;
  public:
    /// CONSTRUCTOR
    Builder();
    /// Set debug
    void SetDebug(int d) { debug_ = d; }
    /// Set optional parameter set
    void SetParameters(ParameterSet const*);

    /// Combine second fragment into first fragment and bond
    int Combine(Topology&, Frame&, Topology const&, Frame const&, int, int) const;
    /// Model the coordinates around a bond given only some coordinates are known
    int ModelCoordsAroundBond(Frame&, Topology const&, int, int, Zmatrix const*, Zmatrix const*, Barray&) const;
  private:
    /// Assign a reasonable value for bond distance given 2 atoms whose position may or may not be known
    int AssignLength(double&, int, int, Topology const&, Frame const&, std::vector<bool> const&) const;
    /// Given atoms J and K, attempt to assign a reasonable value for theta for atom I
    int AssignTheta(double&, int, int, int, Topology const&, Frame const&, std::vector<bool> const&) const;
    /// Insert an internal coord into a zmatrix
    int insertIc(Zmatrix&, int, int, int, int, double,
                 Topology const&, Frame const&, std::vector<bool> const&) const;
    /// Assign internal coordinates for atoms I for torsions around J-K-L.
    int AssignICsAroundBond(Zmatrix&, int, int, int,
                           Topology const&, Frame const&, std::vector<bool> const&,
                           BuildAtom const&) const;

    /// Model coordinates around a bond
    int SetupICsAroundBond(Zmatrix&, int, int, Frame const&, Topology const&,
                           std::vector<bool> const&, std::vector<bool> const&,
                           BuildAtom const&, BuildAtom const&) const;

    int debug_;
    ParameterSet const* params_;
};
}
}
#endif
