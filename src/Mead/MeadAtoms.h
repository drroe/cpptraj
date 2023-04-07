#ifndef INC_MEAD_MEADATOMS_H
#define INC_MEAD_MEADATOMS_H
// Fwd declares
class Topology;
class Frame;
// MEAD fwd declares
class AtomSet;
namespace Cpptraj {
namespace Mead {
/// Wrapper around Mead AtomSet
class MeadAtoms {
  public:
    /// Different radii sets TODO combine with Traj_PDBfile
    enum Radii_Mode { GB = 0, PARSE, VDW };
    /// CONSTRUCTOR
    MeadAtoms();
    /// DESTRUCTOR
    ~MeadAtoms();
    /// Setup AtomSet from top/frame
    int SetupAtoms(Topology const&, Frame const&, Radii_Mode);
    /// \return True if atom set is allocated
    bool IsSetup() const { return atomset_ != 0; }
  private:
    AtomSet* atomset_;
    Radii_Mode rmode_;
};
}
}
#endif
