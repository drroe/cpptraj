#include "Protonator.h"
#include "../ArgList.h"
#include "../CpptrajState.h"
#include "../CpptrajStdio.h"
#include "../DataSet_1D.h"
#include "../DataSet_2D.h"
#include "../Mead/MultiFlexResults.h"

using namespace Cpptraj::Structure;
/** CONSTRUCTOR */
Protonator::Protonator() :
  site_intrinsic_pKas_(0),
  site_site_matrix_(0),
  site_qunprot_(0),
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
  site_qunprot_ = results.QunprotSet();
  if (site_qunprot_ == 0) {
    mprinterr("Internal Error: Site unprotonated state charge data set is missing.\n");
    return 1;
  }
  // TODO check that matrix rows/cols and # sites match?
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

/** Calculate titration curves using MC */
int Protonator::CalcTitrationCurves() const {
  // Calculate the ground energy of the system (no protons)
  unsigned int maxsite = site_intrinsic_pKas_->Size();
  DataSet_1D const& qunprot = static_cast<DataSet_1D const&>( *site_qunprot_ );
  DataSet_2D const& wint = static_cast<DataSet_2D const&>( *site_site_matrix_ );
  // TODO can we assume symmetric?
  double grounde = 0;
  for (unsigned int i = 0; i < maxsite; i++)
    for (unsigned int j = 0; j < maxsite; j++)
      grounde += (qunprot.Dval(i) * qunprot.Dval(j) * wint.GetElement(i,j));
  grounde *= 0.5;
  mprintf("DEBUG: grounde = %g\n", grounde);

  return 0;
}
