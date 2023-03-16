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
    /// CONSTRUCTOR
    MeadInterface();
    /// DESTRUCTOR
    ~MeadInterface();
    /// Add a grid to the FDM object
    int AddGrid(int, float, Vec3 const&);
    /// Setup AtomSet from top/frame
    int SetupAtoms(Topology const&, Frame const&);
    /// Print info to stdout
    void Print() const;
  private:
    FinDiffMethod* fdm_;
    AtomSet* atomset_;
};
}
#endif
