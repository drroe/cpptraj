#ifndef INC_MEADINTERFACE_H
#define INC_MEADINTERFACE_H
class FinDiffMethod;
namespace Cpptraj {
/// Class to interface with libmead.a
class MeadInterface {
  public:
    /// CONSTRUCTOR
    MeadInterface();
    /// DESTRUCTOR
    ~MeadInterface();
  private:
    FinDiffMethod* fdm_;
};
}
#endif
