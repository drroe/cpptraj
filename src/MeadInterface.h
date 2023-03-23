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
  private:
    static int ERR(const char*, MEADexcept&);

    FinDiffMethod* fdm_;
    AtomSet* atomset_;
};
}
#endif
