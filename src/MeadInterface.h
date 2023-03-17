#ifndef INC_MEADINTERFACE_H
#define INC_MEADINTERFACE_H
// Fwd declares
class Vec3;
class Topology;
class Frame;
// MEAD fwd declares
class FinDiffMethod;
class AtomSet;
namespace Cpptraj {
/// Class to interface with libmead.a
class MeadInterface {
  public:
    /// Different radii sets TODO combine with Traj_PDBfile
    enum Radii_Mode { GB = 0, PARSE, VDW };

    /// CONSTRUCTOR
    MeadInterface();
    /// DESTRUCTOR
    ~MeadInterface();
    /// Add a grid to the FDM object
    int AddGrid(int, float, Vec3 const&);
    /// Setup AtomSet from top/frame
    int SetupAtoms(Topology const&, Frame const&, Radii_Mode);
    /// Print info to stdout
    void Print() const;

    /// \return True if finite difference method is allocated.
    bool HasFDM() const { return fdm_ != 0; }
    /// \return True if atom set is allocated
    bool HasAtoms() const { return atomset_ != 0; } 
  private:
    FinDiffMethod* fdm_;
    AtomSet* atomset_;
};
}
#endif
