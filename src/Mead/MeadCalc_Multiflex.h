#ifndef INC_MEAD_MEADCALC_MULTIFLEX_H
#define INC_MEAD_MEADCALC_MULTIFLEX_H
#include "MeadCalc.h"
#include "../Timer.h"
// Fwd declares
// MEAD fwd declares
class OutPotat;
class AtomChargeSet;

namespace Cpptraj {
namespace Structure {
class SiteData;
class TitratableSite;
}
namespace Mead {
class MeadOpts;
class MeadGrid;
class MultiFlexResults;
class MeadCalc_Multiflex : public MeadCalc {
  public:
    /// Different radii sets TODO combine with Traj_PDBfile
    enum Radii_Mode { GB = 0, PARSE, VDW };

    /// CONSTRUCTOR
    MeadCalc_Multiflex();

    /// Run multiflex calc
    int MultiFlex(MultiFlexResults&, MeadOpts const&, MeadGrid const&, MeadGrid const&,
                  Topology const&, Frame const&, Structure::SiteData const&, int);

    /// Access timer
    Timer const& TotalTime() const { return t_total_; }
  private:
    class TitrationCalc;

    //void printAtomPotentials(Topology const&, Frame const&, OutPotat*, AtomChargeSet*) const;

    int createModelCompounds(AtomChargeSet&, AtomChargeSet&, AtomChargeSet const&, int, Topology const&, Frame const&) const;

    int setup_titration_site_calc(std::vector<TitrationCalc>&, AtomChargeSet&,
                                  Topology const&, Frame const&,
                                   Structure::TitratableSite const&, int) const;

    int setup_titration_calcs_by_site(std::vector<TitrationCalc>&, AtomChargeSet&,
                              Topology const&, Frame const&, Structure::SiteData const&) const;

  Timer t_total_; ///< Total time
};
} /* END namespace Mead */
} /* END namespace Cpptraj */
#endif
