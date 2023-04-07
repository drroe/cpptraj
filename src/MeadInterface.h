#ifndef INC_MEADINTERFACE_H
#define INC_MEADINTERFACE_H
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
class MEADexcept;
namespace MEAD {
class Atom;
}
class OutPotat;
class AtomChargeSet;
namespace Cpptraj {
class MultiFlexResults;
namespace Structure {
class TitrationData;
}
namespace Mead {
class MeadOpts;
class MeadGrid;
}
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
    int Potential(DataSet_Vector_Scalar&, Mead::MeadOpts const&, Mead::MeadGrid const&, std::vector<Vec3> const&) const;
    /// Run solvate calc
    int Solvate(double&, Mead::MeadOpts const&, Mead::MeadGrid const&, DataSet_3D*) const;
    /// Run multiflex calc
    int MultiFlex(MultiFlexResults const&, Mead::MeadOpts const&, Mead::MeadGrid const&, Mead::MeadGrid const&,
                  Topology const&, Frame const&, Structure::TitrationData const&, Radii_Mode) const;
  private:
    class TitrationCalc;

    static int ERR(const char*, MEADexcept&);

    static inline void set_at_from_top(MEAD::Atom&, Topology const&, Frame const&, int, Radii_Mode);

    static inline void printAtomPotentials(Topology const&, Frame const&, OutPotat*, AtomChargeSet*);

    static int createModelCompounds(AtomChargeSet&, AtomChargeSet&, AtomChargeSet const&, int, Topology const&, Frame const&, Radii_Mode);

    static int setup_titration_calcs(std::vector<TitrationCalc>&, AtomChargeSet&,
                                     Topology const&, Frame const&,
                                     Structure::TitrationData const&, Radii_Mode);

    AtomSet* atomset_;
};
} /* END namespace Cpptraj */
#endif
