#ifndef INC_HB_HBCALC_H
#define INC_HB_HBCALC_H
#include "../AtomMask.h"
#include "../PairList.h"
class Atom;
class Topology;
namespace Cpptraj {
namespace HB {
/// Main driver for hydrogen bond calculation
class HbCalc {
  public:
    /// CONSTRUCTOR
    HbCalc();

    int SetupPairlistAtomMask(Topology const&);
  private:
    PairList pairList_; ///< Pair list for atoms involved in hydrogen bond calc

    AtomMask generalMask_; ///< Mask of atoms to potentially calculate hydrogen bonds for

    AtomMask plMask_; ///< Mask selecting atoms to go into the pairlist

    /// Different atom types
    enum Type { HYDROGEN = 0, DONOR, ACCEPTOR, BOTH };

    typedef std::vector<Type> Tarray;

    Tarray plTypes_; ///< Type of each atom in plMask_

    static inline bool IsFON( Atom const& );
};
}
}
#endif
