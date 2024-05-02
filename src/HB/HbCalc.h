#ifndef INC_HB_HBCALC_H
#define INC_HB_HBCALC_H
#include "Site.h"
#include "../AtomMask.h"
#include "../PairList.h"
class ArgList;
class Atom;
class Box;
class Topology;
class Frame;
namespace Cpptraj {
namespace HB {
/// Main driver for hydrogen bond calculation
class HbCalc {
  public:
    /// CONSTRUCTOR
    HbCalc();

    int InitHbCalc(ArgList&, int);

    int SetupHbCalc(Topology const&, Box const&);

    void PrintHbCalcOpts() const;

    int RunCalc_PL(Frame const&);
  private:
    /// Different atom types
    enum Type { HYDROGEN = 0, DONOR, ACCEPTOR, BOTH, VDONOR, VACCEPTOR, VBOTH };
    /// Strings for different atom types
    static const char* TypeStr_[];

    typedef std::vector<Type> Tarray;
    typedef std::vector<Site> Sarray;
    typedef std::vector<int> Iarray;

    /// \return True if Atom element is F, O, or N
    static inline bool IsFON( Atom const& );

    int setupPairlistAtomMask(Topology const&);

    PairList pairList_;    ///< Pair list for atoms involved in hydrogen bond calc
    AtomMask generalMask_; ///< Mask of atoms to potentially calculate hydrogen bonds for
    AtomMask plMask_;      ///< Mask selecting atoms to go into the pairlist
    Tarray plTypes_;       ///< Type of each atom in plMask_
    Sarray Both_;          ///< HB donor/acceptor sites followed by donor-only
    Iarray Acceptor_;      ///< HB acceptor sites
    double dcut2_;         ///< Heavy atom distance cutoff (Ang) squared
};
}
}
#endif
