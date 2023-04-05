#ifndef INC_MEADINTERFACE_H
#define INC_MEADINTERFACE_H
#include <vector>
// Fwd declares
class Vec3;
class Topology;
class Frame;
class DataSet_Vector_Scalar;
class DataSet_3D;
// MEAD fwd declares
class FinDiffMethod;
class AtomSet;
class MEADexcept;
namespace MEAD {
class Atom;
}
class OutPotat;
class AtomChargeSet;
namespace Cpptraj {
namespace Structure {
class TitrationData;
}
/// Class to interface with libmead.a
class MeadInterface {
  public:
    /// Different radii sets TODO combine with Traj_PDBfile
    enum Radii_Mode { GB = 0, PARSE, VDW };
    /// Grid centering modes. 
    enum GridCenter_Mode { C_ON_ORIGIN = 0, C_ON_CENT_OF_INTR, C_ON_GEOM_CENT };

    /// \return Character string corresponding to GridCenter_Mode
    static const char* GridCenter_ModeStr(GridCenter_Mode);

    /// CONSTRUCTOR
    MeadInterface();
    /// DESTRUCTOR
    ~MeadInterface();
    /// Add a grid to the FDM object with explicit centering
    int AddGrid(int, float, Vec3 const&);
    /// Add a grid to the FDM object with centering type
    int AddGrid(int, float, GridCenter_Mode);
    /// Setup AtomSet from top/frame
    int SetupAtoms(Topology const&, Frame const&, Radii_Mode);
    /// Print info to stdout
    void Print() const;
    /// Set MEAD verbosity level
    void MeadVerbosity(int) const;

    /// \return True if finite difference method is allocated.
    bool HasFDM() const { return fdm_ != 0; }
    /// \return True if atom set is allocated
    bool HasAtoms() const { return atomset_ != 0; }

    /// Run potential calc
    int Potential(DataSet_Vector_Scalar&, double, double, std::vector<Vec3> const&) const;
    /// Run solvate calc
    int Solvate(double&, double, double, double, double, double, double, double, DataSet_3D*) const;
    /// Run multiflex calc
    int MultiFlex(double, double, double, double, double, Topology const&, Frame const&, Structure::TitrationData const&, Radii_Mode) const;
  private:
    class TitrationCalc;

    static const char* GridCenter_ModeStr_[];

    static int ERR(const char*, MEADexcept&);

    static inline void set_at_from_top(MEAD::Atom&, Topology const&, Frame const&, int, Radii_Mode);

    static inline void printAtomPotentials(Topology const&, Frame const&, OutPotat*, AtomChargeSet*);

    static int createModelCompounds(AtomChargeSet&, AtomChargeSet&, int, Topology const&, Frame const&, Radii_Mode);

    static int setup_titration_calcs(std::vector<TitrationCalc>&, AtomChargeSet&,
                                     Topology const&, Frame const&,
                                     Structure::TitrationData const&, Radii_Mode);

    FinDiffMethod* fdm_;
    FinDiffMethod* mgm_; ///< Hold grid for model (MULTIFLEX only) TODO wrap in class?
    AtomSet* atomset_;
};
}
#endif
