#include "Protonator.h"
#include "../ArgList.h"
#include "../CpptrajState.h"
#include "../CpptrajStdio.h"
#include "../Mead/MultiFlexResults.h"

using namespace Cpptraj::Structure;
/** CONSTRUCTOR */
Protonator::Protonator() :
  site_intrinsic_pKas_(0),
  site_site_matrix_(0),
  n_mc_steps_(0),
  start_pH_(0),
  stop_pH_(0),
  pH_increment_(0),
  min_wint_(0),
  fract_toler_(0),
  mcmode_(MC_FULL),
  logfile_(0)
{}

/** Set up options */
int Protonator::SetupProtonator(CpptrajState& State, ArgList& argIn,
                                Cpptraj::Mead::MultiFlexResults const& results)
{
  site_intrinsic_pKas_ = results.PkIntSet();
  site_site_matrix_ = results.SiteSiteMatrixSet();

  return 0;
}

/** Print options */
void Protonator::PrintOptions() const {
  if (site_intrinsic_pKas_ != 0) mprintf("\tSite intrinsic pKa set: %s\n", site_intrinsic_pKas_->legend());
}
