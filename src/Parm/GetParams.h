#ifndef INC_PARM_GETPARAMS_H
#define INC_PARM_GETPARAMS_H
#include <vector>
class Atom;
class AtomType;
class HB_ParmType;
class NonbondParmType;
class NonbondType;
template<typename Type> class ParmHolder;
namespace Cpptraj {
namespace Parm {
/// Used to get parameters from a Topology
class GetParams {
  public:
    /// CONSTRUCTOR
    GetParams();
    /// Get nonbonded parameters
    void GetLJAtomTypes(ParmHolder<AtomType>&,
                        ParmHolder<NonbondType>&,
                        ParmHolder<NonbondType>&,
                        ParmHolder<double>&,
                        ParmHolder<HB_ParmType>&,
                        std::vector<Atom> const&,
                        NonbondParmType const&) const;
  private:
    int debug_;
};
}
}
#endif
