#include "Protonator.h"
#include "../ArgList.h"
#include "../CpptrajState.h"
#include "../CpptrajStdio.h"
#include "../DataSet_1D.h"
#include "../DataSet_2D.h"
#include "../Random.h"
#include "../Mead/MultiFlexResults.h"
#include <cmath> // log

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
  beta_(0),
  iseed_(0),
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
  // Beta is 1/kT times Coulombs constant in units of kcal*Ang/mol*e- squared
  double temperature = argIn.getKeyDouble("temperature", -1);
  if (temperature > 0)
    beta_ = (1 / (Constants::GASK_KCAL * temperature)) * Constants::ELECTOAMBER * Constants::ELECTOAMBER;
  else
    // For backwards compat. use value in mcti.f
    // I'm guessing this is the original math:
    // (1/(0.001987*300)) * 18.2223 * 18.2223
    beta_ = 557.04;

  iseed_ = argIn.getKeyInt("iseed", 0);

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
  mprintf("\tValue for converting from charge to kcal/mol: beta= %g\n", beta_);
  mprintf("\tRNG seed: %i\n", iseed_);
}

/** Assign random protonation state to each site.
  * \return total number of protonated sites.
  */
int Protonator::assign_random_state(Iarray& SiteIsProtonated, Random_Number& rng)
const
{
  // FIXME calling set seed here to match mcti init subroutine. Should not always be done.
  // When fixed, rng can be const&
  int ntotal = 0;
  rng.rn_set( iseed_ );
  for (Iarray::iterator it = SiteIsProtonated.begin(); it != SiteIsProtonated.end(); ++it)
  {
    if (rng.rn_gen() > 0.5) {
      *it = 1;
      ntotal++;
    } else {
      *it = 0;
    }
  }
  return ntotal;
}

/** \return Change in energy upon changes in protonation at site flip. */
double Protonator::mc_deltae(Iarray const& SiteIsProtonated, int iflip,
                             unsigned int maxsite,
                             DataSet_2D const& wint, DataSet_1D const& qunprot,
                             Darray const& self)
{
  int dp = 1 - (2*SiteIsProtonated[iflip]);
  double deltae = 0;

  for (unsigned int j=0; j < maxsite; j++) {
    deltae = deltae + ((dp * (qunprot.Dval(j) + SiteIsProtonated[j]))) * wint.GetElement(iflip, j);
  }
  deltae = deltae + (dp * self[iflip]);

  return deltae;
}

/** Perform monte carlo sampling for a given pH value.
  * Determine average values of the protonation for each site and
  * the correlation functions used to compute the energy.
  */
int Protonator::perform_MC_at_pH(double pH, Iarray const& SiteIsProtonated,
                                 DataSet_1D const& pkint,
                                 DataSet_2D const& wint, DataSet_1D const& qunprot,
                                 Random_Number const& rng)
const
{
  unsigned int maxsite = pkint.Size();
  // Compute the self energy to protonate each site at this pH
  Darray self;
  self.reserve( maxsite );
  double fac = -(log(10.0) / beta_);
  mprintf("DEBUG: MC self E:");
  for (unsigned int i = 0; i < maxsite; i++) {
    self[i] = fac * (pkint.Dval(i) - pH);
    mprintf(" %E", self[i]);
  }
  mprintf("\n");
  // Thermalization: do n_mc_steps_ of Monte Carlo to approach equilibrium
  //                 before data is taken.
  // NOTE: This currently does not allow a 2 site transition.
  int maxsteps = n_mc_steps_ * maxsite;
  for (int i = 0; i < maxsteps; i++) {
    // Choose site
    int iflip = (int)((double)(maxsite-1) * rng.rn_gen());
    // Choose whether to flip
  }

  return 0;
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

  // Calculate pairs.
  // Pairs of strongly interacting sites are allowed to simultaneously
  // change their protonation states.
  typedef std::pair<int,int> StatePair;
  typedef std::vector<StatePair> PairArray;
  PairArray pairs;
  double min_g = min_wint_ * log(10.0) / beta_;
  mprintf("DEBUG: min_g = %g\n", min_g);
  for (unsigned int i = 0; i < maxsite; i++) {
    for (unsigned int j = i+1; j < maxsite; j++) {
      if (wint.GetElement(i,j) > min_g)
        pairs.push_back( StatePair(i,j) );
    }
  }
  mprintf("DEBUG: %zu pairs.\n", pairs.size());
  // Count the number of pH values
  int nph = (int)((stop_pH_ - start_pH_) / pH_increment_) + 1;
  std::vector<double> pH_values;
  pH_values.reserve( nph );
  double pH = start_pH_;
  while (pH <= stop_pH_) {
    pH_values.push_back( pH );
    pH += pH_increment_;
  }
  logfile_->Printf("maxsite= %u  #pH vals= %zu  nph= %i\n", maxsite, pH_values.size(), nph);

  Iarray SiteIsProtonated(maxsite, 0);

  // Loop over pH values
  for (Darray::const_iterator ph = pH_values.begin(); ph != pH_values.end(); ++ph)
  {
    mprintf("DEBUG: pH %g\n", *ph);
    Random_Number rng; // TODO use 1 overall rng
    // Assign initial protonation
    int nprotonated = assign_random_state(SiteIsProtonated, rng);
    mprintf("Initial states (%i total):", nprotonated);
    for (Iarray::const_iterator it = SiteIsProtonated.begin(); it != SiteIsProtonated.end(); ++it)
      mprintf(" %i", *it);
    mprintf("\n");
    // Do the MC trials at this pH
    int err =  perform_MC_at_pH(*ph, SiteIsProtonated,
                                static_cast<DataSet_1D const&>( *site_intrinsic_pKas_ ),
                                wint, qunprot, rng);
    if (err != 0) {
      mprinterr("Error: Could not perform MC at pH %g\n", *ph);
      return 1;
    }
  } // END loop over pH values
  

  return 0;
}
