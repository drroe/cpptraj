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
class InternalCoords;
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
    /// Set optional Zmatrix with current ICs
    void SetZmatrix(Zmatrix const*);

    /// Combine second fragment into first fragment and bond
    int Combine(Topology&, Frame&, Topology const&, Frame const&, int, int) const;
    /// Model the coordinates around a bond given only some coordinates are known
    int ModelCoordsAroundBond(Frame const&, Topology const&, int, int, Zmatrix&, Barray const&) const;
    /// Update the internal coordinates in given Zmatrix with values from Frame/Parameters
    int UpdateICsFromFrame(Zmatrix&, Frame const&, int, Topology const&, Barray const&) const;
  private:
    typedef std::vector<int> Iarray;
    /// Get length parameter for atoms
    int getLengthParam(double&, int, int, Topology const&) const;
    /// Assign a reasonable value for bond distance given 2 atoms whose position may or may not be known
    int AssignLength(double&, int, int, Topology const&, Frame const&, Barray const&) const;
    /// Get angle parameter for atoms.
    int getAngleParam(double&, int, int, int, Topology const&) const;
    /// Given atoms J and K, attempt to assign a reasonable value for theta for atom I
    int AssignTheta(double&, int, int, int, Topology const&, Frame const&, Barray const&) const;
    /// Calculate an internal coordinate for known atoms
    static inline InternalCoords calcKnownAtomIc(int, int, int, int, Frame const&);
    /// Insert an internal coord into a zmatrix
    int insertIc(Zmatrix&, int, int, int, int, double,
                 Topology const&, Frame const&, Barray const&) const;
    /// Assign internal coordinates for atoms I for torsions around J-K-L.
    int AssignICsAroundBond(Zmatrix&, int, int, int,
                           Topology const&, Frame const&, Barray const&,
                           BuildAtom const&) const;

    /// Model coordinates around a bond
    int SetupICsAroundBond(Zmatrix&, int, int, Frame const&, Topology const&,
                           Barray const&, Barray const&,
                           BuildAtom const&, BuildAtom const&) const;
    /// Generate internal coordinates in the same manner as LEaP
    int GenerateInternals(Zmatrix&, Frame const&, Topology const&, Barray const&);
    int debug_;
    ParameterSet const* params_;
    Zmatrix const* currentZmatrix_; ///< Any existing internal coordinates
};
}
}
#endif
