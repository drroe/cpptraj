#ifndef INC_MEAD_MEADINTERFACE_H
#define INC_MEAD_MEADINTERFACE_H
#include <vector>
#include <string>
// Fwd declares
class Vec3;
class Topology;
class Frame;
class DataSet_Vector_Scalar;
class DataSet_3D;
// MEAD fwd declares
class AtomSet;
namespace MEAD {
class Atom;
}
class OutPotat;
class AtomChargeSet;

namespace Cpptraj {
namespace Structure {
class TitrationData;
}
namespace Mead {
class MeadOpts;
class MeadGrid;
class MultiFlexResults;
/// Class to interface with libmead.a
class MeadInterface {
  public:
    /// Different radii sets TODO combine with Traj_PDBfile
    enum Radii_Mode { GB = 0, PARSE, VDW };

    /// CONSTRUCTOR
    MeadInterface();
    /// DESTRUCTOR
    ~MeadInterface();

    /// Setup AtomSet from top/frame
    int SetupAtoms(Topology const&, Frame const&, Radii_Mode);
    /// Print info to stdout
    void Print() const;
    /// Set MEAD verbosity level
    void MeadVerbosity(int) const;

    /// \return True if atom set is allocated
    bool HasAtoms() const { return atomset_ != 0; }

    /// Run potential calc
    int Potential(DataSet_Vector_Scalar&, MeadOpts const&, MeadGrid const&, std::vector<Vec3> const&) const;
    /// Run solvate calc
    int Solvate(double&, MeadOpts const&, MeadGrid const&, DataSet_3D*) const;
    /// Run multiflex calc
    int MultiFlex(MultiFlexResults const&, MeadOpts const&, MeadGrid const&, MeadGrid const&,
                  Topology const&, Frame const&, Structure::TitrationData const&) const;
  private:
    class TitrationCalc;

    inline void set_at_from_top(MEAD::Atom&, Topology const&, Frame const&, int) const;

    void printAtomPotentials(Topology const&, Frame const&, OutPotat*, AtomChargeSet*) const;

    int createModelCompounds(AtomChargeSet&, AtomChargeSet&, AtomChargeSet const&, int, Topology const&, Frame const&) const;

    int setup_titration_calcs(std::vector<TitrationCalc>&, AtomChargeSet&,
                              Topology const&, Frame const&, Structure::TitrationData const&) const;

    AtomSet* atomset_;
    Radii_Mode rmode_; ///< Which radii set to use
};
} /* END namespace Mead */
} /* END namespace Cpptraj */
#endif
