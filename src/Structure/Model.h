#ifndef INC_STRUCTURE_MODEL_H
#define INC_STRUCTURE_MODEL_H
#include <vector>
//#incl ude "StructureEnum.h"
class Topology;
class Frame;
//template <typename T> class ParmHolder;
//class AtomType;
namespace Cpptraj {
namespace Structure {
class BuildAtom;
class InternalCoords;
/// Routines to create a model for missing bond/angle/torsion parameters
class Model {
  public:
    Model() : debug_(0) {}
    /// Set debug level
    void SetDebug(int d) { debug_ = d; }
    /// Assign a reasonable value for bond distance given 2 atoms whose position may or may not be known
    int AssignLength(double&, int, int, Topology const&, Frame const&, std::vector<bool> const&) const;
    /// Given atoms J and K, attempt to assign a reasonable value for theta for atom I
    int AssignTheta(double&, int, int, int, Topology const&, Frame const&, std::vector<bool> const&) const;
    /// Given atoms J K and L, attempt to assign a reasonable value for phi for atom I
    //int AssignPhi(std::vector<InternalCoords>&, int, int, int, Topology const&, Frame const&, std::vector<bool> const&, BuildAtom const&) const;
    int AssignPhi(double&, int, int, int, int, Topology const&, Frame const&, std::vector<bool> const&, BuildAtom const&) const;
    /// Assign phi around bond
    //int AssignPhiAroundBond(int, int, Topology const&, Frame const&, std::vector<bool> const&, ParmHolder<AtomType> const&) const;
  private:
    int debug_;
};
} // END namespace Structure
} // END namespace Cpptraj
#endif
