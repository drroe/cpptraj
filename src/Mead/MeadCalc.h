#ifndef INC_MEAD_MEADCALC_H
#define INC_MEAD_MEADCALC_H
#include "../Timer.h"
// Fwd declares
class Topology;
class Frame;
// MEAD fwd declares
class AtomSet;
namespace MEAD {
class Atom;
}
namespace Cpptraj {
namespace Mead {
/// Base Class to interface with libmead.a
class MeadCalc {
  public:
    /// Different radii sets TODO combine with Traj_PDBfile
    enum Radii_Mode { GB = 0, PARSE, VDW };

    /// CONSTRUCTOR
    MeadCalc();
    /// DESTRUCTOR - virtual since inherited
    virtual ~MeadCalc();

    /// Setup AtomSet from top/frame
    int SetupAtoms(Topology const&, Frame const&, Radii_Mode);
    /// Access the internal timer
    Timer& TotalTime() { return t_total_; }
    /// Print info to stdout
    void Print() const;
    /// Set MEAD verbosity level
    void MeadVerbosity(int) const;

    /// \return True if atom set is allocated
    bool HasAtoms() const { return atomset_ != 0; }
  protected:
    AtomSet const& InternalAtomset() const { return *atomset_; }
  private:
    void set_at_from_top(MEAD::Atom&, Topology const&, Frame const&, int) const;

    AtomSet* atomset_;
    Radii_Mode rmode_; ///< Which radii set to use
    Timer t_total_;      ///< Total time
};
} /* END namespace Mead */
} /* END namespace Cpptraj */
#endif
