#ifndef INC_HB_HBCALC_H
#define INC_HB_HBCALC_H
#include "HbData.h"
#include "../AtomMask.h"
#include "../PairList.h"
#ifdef TIMER
# include "../Timer.h"
#endif
class ArgList;
class Atom;
class Box;
class DataFileList;
class Frame;
class Topology;
namespace Cpptraj {
namespace HB {
class Hbond;
/// Main driver for hydrogen bond calculation
class HbCalc {
  public:
    /// CONSTRUCTOR
    HbCalc();
    /// Read hydrogen bond calc options
    int InitHbCalc(ArgList&, DataSetList*, DataFileList&, int);
    /// Set up hydrogen bond calc for given Topology; pair list is set up from box
    int SetupHbCalc(Topology const&, Box const&);
    /// Print hydrogen bond calc options to stdout
    void PrintHbCalcOpts() const;
    /// Run hydrogen bond calc on a frame
    int RunCalc_PL(Frame const&, int, int);
    /// Finalize the hydrogen bond calculation
    void FinishHbCalc();
    /// Set hydrogen bond calculation debug level
    void SetDebug(int);
  private:
    /// Different atom types
    enum Type { DONOR=0, ACCEPTOR, BOTH, VDONOR, VACCEPTOR, VBOTH, UNKNOWN };
    /// Strings for different atom types
    static const char* TypeStr_[];

    typedef std::vector<Type> Tarray;
//    typedef std::vector<std::string> Sarray;
    typedef std::vector<int> Iarray;
    typedef std::vector<Iarray> Xarray;

    /// \return True if atom is F, O, or N
    static inline bool IsFON( Atom const& );
    /// \return True if interaction is valid between given types
    static inline bool validInteraction(Type, Type);
    /// Set up pair list for generalMask_ and given Topology
    int setupPairlistAtomMask(Topology const&);
    /// Calculate hydrogen bonds between solute donor site and solute acceptor
    void calc_UU_Hbonds(int, double, int, Iarray const&, int, Frame const&, int&, int);
    /// Calculate hydrogen bonds between solute and solvent
    void calc_UV_Hbonds(int, double, int, Iarray const&, int, Frame const&, int&, bool, int);

    /// Calculate a potentially imaged angle
    double Angle(const double*, const double*, const double*, Box const&) const;
    /// Calculate hydrogen bonds between two atoms
    void CalcHbonds(int, double, int, int, Frame const&, int&, int);

    PairList pairList_;    ///< Pair list for atoms involved in hydrogen bond calc
    AtomMask generalMask_; ///< Mask of atoms to potentially calculate hydrogen bonds for
    AtomMask plMask_;      ///< Mask selecting atoms to go into the pairlist
    Tarray plTypes_;       ///< Type of each atom in plMask_
//    Sarray plNames_;       ///< Name of each atom in plMask_
    Xarray plHatoms_;      ///< Indices of any hydrogens bonded to each atom in plMask_
    double dcut2_;         ///< Heavy atom distance cutoff (Ang) squared
    double acut_;          ///< Angle cutoff in radians
    HbData hbdata_;        ///< Hold hydrogen bond calculation data.
#   ifdef TIMER
    Timer t_action_;       ///< Total time taken by RunCalc_PL
    Timer t_hbcalc_;       ///< Time taken by pairlist loop in RunCalc_PL
#   endif
#   ifdef _OPENMP
    typedef std::vector<Hbond> Harray;
    std::vector<Harray> thread_HBs_; ///< Hold hbonds found by each thread each frame.
#   endif
};
}
}
#endif
