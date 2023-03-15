#ifndef INC_MEADINTERFACE_H
#define INC_MEADINTERFACE_H
// Fwd declares
class Vec3;
// MEAD fwd declares
class FinDiffMethod;
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
    /// Print info to stdout
    void Print() const;
  private:
    FinDiffMethod* fdm_;
};
}
#endif
