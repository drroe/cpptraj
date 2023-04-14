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
  if (site_intrinsic_pKas_ == 0) {
    mprinterr("Internal Error: Site intrinsic pKa data set is missing.\n");
    return 1;
  }
  site_site_matrix_ = results.SiteSiteMatrixSet();
  if (site_site_matrix_ == 0) {
    mprinterr("Internal Error: Site-site interaction matrix data set is missing.\n");
    return 1;
  }
  n_mc_steps_ = argIn.getKeyInt("nmcsteps", 10000);
  start_pH_ = argIn.getKeyDouble("startph", 5.0);
  stop_pH_ = argIn.getKeyDouble("stopph", 10.0);
  pH_increment_ = argIn.getKeyDouble("phincr", 0.5);
  min_wint_ = argIn.getKeyDouble("minwint", 2.0);
  fract_toler_ = argIn.getKeyDouble("fracttol", 0.000001);
  std::string mcstr = argIn.GetStringKey("mcmode");
  if (mcstr.empty())
    mcmode_ = MC_REDUCED;
  else {
    if (mcstr == "full" || mcstr == "0") //  0 for backwards compat.
      mcmode_ = MC_FULL;
    else if (mcstr == "reduced" || mcstr == "1") // 1 for backwards compat.
      mcmode_ = MC_REDUCED;
    else { // CLUSTER not yet implemented
      mprinterr("Error: Bad value for 'mcmode': %s\n", mcstr.c_str());
      return 1;
    }
  }
  logfile_ = State.DFL().AddCpptrajFile( argIn.GetStringKey("mclog"),
                                         "MC log file",
                                         DataFileList::TEXT,
                                         true );

  return 0;
}

/** Print options */
void Protonator::PrintOptions() const {
  if (site_intrinsic_pKas_ != 0) mprintf("\tSite intrinsic pKa set: %s\n", site_intrinsic_pKas_->legend());
  if (site_site_matrix_ != 0) mprintf("\tSite-site interaction matrix set: %s\n", site_site_matrix_->legend());
  mprintf("\tCalculating titration curves from pH %g to %g in increments of %g\n",
          start_pH_, stop_pH_, pH_increment_);
  mprintf("\tCutoff for site-site interactions in pH units is %g\n", min_wint_);
  mprintf("\tTolerance for reduced sites is %g\n", fract_toler_);
  static const char* MCMODESTR[] = { "full", "reduced", "cluster" };
  mprintf("\tMC mode: %s\n", MCMODESTR[mcmode_]);
  if (logfile_ != 0) mprintf("\tLog output to '%s'\n", logfile_->Filename().full());
}
