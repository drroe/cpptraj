#ifndef INC_HB_HBCALC_H
#define INC_HB_HBCALC_H
#include "HbData.h"
#include "../AtomMask.h"
#include "../PairList.h"
class ArgList;
class Atom;
class Box;
class DataFileList;
class Frame;
class Topology;
namespace Cpptraj {
namespace HB {
/// Main driver for hydrogen bond calculation
class HbCalc {
  public:
    /// CONSTRUCTOR
    HbCalc();

    int InitHbCalc(ArgList&, DataSetList*, DataFileList&, int);

    int SetupHbCalc(Topology const&, Box const&);

    void PrintHbCalcOpts() const;

    int RunCalc_PL(Frame const&);

    void FinishHbCalc();
  private:
    /// Different atom types
    enum Type { DONOR=0, ACCEPTOR, BOTH, VDONOR, VACCEPTOR, VBOTH, UNKNOWN };
    /// Strings for different atom types
    static const char* TypeStr_[];

    typedef std::vector<Type> Tarray;
    typedef std::vector<std::string> Sarray;
    typedef std::vector<int> Iarray;
    typedef std::vector<Iarray> Xarray;

    /// \return True if atom is F, O, or N
    static inline bool IsFON( Atom const& );
    /// \return True if interaction is valid between given types
    static inline bool validInteraction(Type, Type);
    /// Set up pair list for generalMask_ and given Topology
    int setupPairlistAtomMask(Topology const&);
    /// Calculate hydrogen bonds between donor site and acceptor
    void CalcSiteHbonds(int, double, int, Iarray const&, int, Frame const&, int&, int);
    /// Calculate a potentially imaged angle
    double Angle(const double*, const double*, const double*, Box const&) const;
    /// Calculate hydrogen bonds between two atoms
    void CalcHbonds(int, double, int, int, Frame const&, int&, int);

    PairList pairList_;    ///< Pair list for atoms involved in hydrogen bond calc
    AtomMask generalMask_; ///< Mask of atoms to potentially calculate hydrogen bonds for
    AtomMask plMask_;      ///< Mask selecting atoms to go into the pairlist
    Tarray plTypes_;       ///< Type of each atom in plMask_
    Sarray plNames_;       ///< Name of each atom in plMask_
    Xarray plHatoms_;      ///< Indices of any hydrogens bonded to each atom in plMask_
    double dcut2_;         ///< Heavy atom distance cutoff (Ang) squared
    double acut_;          ///< Angle cutoff in radians
    HbData hbdata_;        ///< Hold hydrogen bond calculation data.
};
}
}
#endif
