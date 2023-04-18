#include "Protonator.h"
#include "../ArgList.h"
#include "../CpptrajState.h"
#include "../CpptrajStdio.h"
#include "../DataSet_1D.h"
#include "../DataSet_2D.h"
#include "../DataSet_string.h"
#include "../Random.h"
#include "../Mead/MultiFlexResults.h"
#include <cmath> // log

using namespace Cpptraj::Structure;
/** CONSTRUCTOR */
Protonator::Protonator() :
  site_intrinsic_pKas_(0),
  site_site_matrix_(0),
  site_qunprot_(0),
  site_names_(0),
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
  site_names_ = results.SiteNamesSet();
  if (site_names_ == 0) {
    mprinterr("Internal Error: Site name data set is missing.\n");
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
  if (site_qunprot_ != 0) mprintf("\tSite unprotonated charge set: %s\n", site_qunprot_->legend());
  if (site_names_ != 0) mprintf("\tSite names set: %s\n", site_names_->legend());
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

// -----------------------------------------------------------------------------
/** Class used to track/manipulate protonation state at each site. */
class Protonator::StateArray {
  public:
    StateArray() : nprotonated_(0) {}
    /// CONSTRUCTOR - Number of sites
    StateArray(unsigned int nsites) : prot_(nsites, 0), nprotonated_(0) {}
    /// \return Total # of sites that are protonated
    int Nprotonated() const { return nprotonated_; }
    /// \return Protonation state at site (1 = protonated, 0 = unprotonated)
    int operator[](int idx) const { return prot_[idx]; }
    /// \return Total energy based on current protonation states
    double Esum(Darray const&, DataSet_1D const&, DataSet_2D const&) const;
    /// Print state to log file
    void PrintState(CpptrajFile*) const;
    /// Assign state to string
    void StateStr(std::string&) const;

    /// Assign random protonation state to each site.
    void AssignRandom(Random_Number const& rng);
    /// Flip protonation state at site
    void FlipProt(int);
    /// Assign incoming state to this state, return true if different
    bool AssignFromState(StateArray const&);

  private:
    typedef std::vector<int> Iarray;
    Iarray prot_;     ///< Hold protonation state for each site
    int nprotonated_; ///< Total # of sites that are protonated
};

/** Assign random protonation state to each site.
  */
void Protonator::StateArray::AssignRandom(Random_Number const& rng)
{
  nprotonated_ = 0;
  for (Iarray::iterator it = prot_.begin(); it != prot_.end(); ++it)
  {
    if (rng.rn_gen() > 0.5) {
      *it = 1;
      nprotonated_++;
    } else {
      *it = 0;
    }
  }
}

/** Flip protonation state at site, update count. */
void Protonator::StateArray::FlipProt(int site) {
  if (prot_[site] == 0) {
    prot_[site] = 1;
    nprotonated_++;
  } else {
    prot_[site] = 0;
    nprotonated_--;
  }
}

/// Calculate total energy based on current protonation states
double Protonator::StateArray::Esum(Darray const& self, DataSet_1D const& qunprot,
                                    DataSet_2D const& wint)
const
{
  unsigned int maxsite = self.size();
  // TODO check sizes?
  double esum = 0;
  for (unsigned int j = 0; j < maxsite; j++) {
    esum = esum + (prot_[j] * self[j]);
    for (unsigned int i = 0; i < maxsite; i++) {
      esum = esum + 0.5 * ( ((qunprot.Dval(i) + prot_[i]) * (qunprot.Dval(j) + prot_[j])) -
                            (qunprot.Dval(i) * qunprot.Dval(j)) ) * wint.GetElement(i,j);
    }
  }
  return esum;
}

/** Assign protonation state to this state array.
  * \return true if this state was different than incoming state.
  */
bool Protonator::StateArray::AssignFromState(StateArray const& rhs) {
  bool isDifferent = false;
  if (prot_.empty()) {
    isDifferent = true;
    prot_ = rhs.prot_;
    nprotonated_ = rhs.nprotonated_;
  } else {
    if (prot_.size() != rhs.prot_.size()) { // sanity check
      mprinterr("Internal Error:  Protonator::StateArray::AssignFromState: Sizes do not match.\n");
      return true;
    }
    for (unsigned int idx = 0; idx != rhs.prot_.size(); idx++) {
      if (prot_[idx] != rhs.prot_[idx])
        isDifferent = true;
      prot_[idx] = rhs.prot_[idx];
    }
  }
  return isDifferent;
}

/** Print state to log file */
void Protonator::StateArray::PrintState(CpptrajFile* logfile) const {
  static const char* fmt = "%12i";
  logfile->Printf(fmt, nprotonated_);
  for (Iarray::const_iterator it = prot_.begin(); it != prot_.end(); ++it)
    logfile->Printf(fmt, *it);
  logfile->Printf("\n");
}

/** Assign state string. */
void Protonator::StateArray::StateStr(std::string& stateChar) const {
  stateChar.resize( prot_.size() );
  for (unsigned int idx = 0; idx != prot_.size(); idx++) {
    if (prot_[idx] == 0)
      stateChar[idx] = '0';
    else
      stateChar[idx] = '1';
  }
}

// -----------------------------------------------------------------------------

/** \return Change in energy upon changes in protonation at site flip. */
double Protonator::mc_deltae(StateArray const& SiteIsProtonated, int iflip,
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

/** Perform a monte carlo step.
  * Consider pairs of strongly interacting residues as separate "sites"
  * that can undergo change (i.e. a two-site transition).
  */
void Protonator::mc_step(double& energy,
                        unsigned int maxsite, Darray const& self,
                        PairArray const& pairs,
                        StateArray& SiteIsProtonated,
                        DataSet_2D const& wint, DataSet_1D const& qunprot,
                        Random_Number const& rng)
const
{
  for (unsigned int i = 0; i < maxsite + pairs.size(); i++)
  {
    // Choose site
    // FIXME use rng.rn_num_interval(0, maxsite+npairs)
    int iflip = (int)((double)(maxsite+pairs.size()) * rng.rn_gen()); // FIXME THIS IS FOR TEST ONLY
    // 2 site or 1 site
    if ((unsigned int)iflip >= maxsite) {
      // 2 site transition
      iflip -= (int)maxsite;
    } else {
      // 1 site transition
      double de = mc_deltae(SiteIsProtonated, iflip, maxsite, wint, qunprot, self);
      // Standard Metropolis criterion.
      // if de < 0 change
      // if de >= 0 change with probability exp(-beta*de)
      if ( de < 0 ) {
        SiteIsProtonated.FlipProt(iflip);
        energy = energy + de;
      } else if ( exp(-(beta_*de)) > rng.rn_gen() ) {
        SiteIsProtonated.FlipProt(iflip);
        energy = energy + de;
      }
    }
  } // END MC loop
}
  

/** Perform monte carlo sampling for a given pH value.
  * Determine average values of the protonation for each site and
  * the correlation functions used to compute the energy.
  */
int Protonator::perform_MC_at_pH(double pH, StateArray& SiteIsProtonated,
                                 DataSet_1D const& pkint,
                                 DataSet_2D const& wint, DataSet_1D const& qunprot,
                                 Random_Number const& rng,
                                 PairArray const& pairs)
const
{
  unsigned int maxsite = pkint.Size();
  // Compute the self energy to protonate each site at this pH
  Darray self;
  self.reserve( maxsite );
  double fac = -(log(10.0) / beta_);
  //mprintf("DEBUG: MC self E:");
  for (unsigned int i = 0; i < maxsite; i++) {
    self.push_back( fac * (pkint.Dval(i) - pH) );
    logfile_->Printf("%26.16E", self[i]);
  }
  logfile_->Printf("\n");
  // Thermalization: do n_mc_steps_ of Monte Carlo to approach equilibrium
  //                 before data is taken.
  // NOTE: This currently does not allow a 2 site transition.
  int maxsteps = n_mc_steps_ * maxsite;
  for (int i = 0; i < maxsteps; i++) {
    // Choose site FIXME use rng.rn_num_interval(0, maxsite)
    int iflip = (int)((double)maxsite * rng.rn_gen()); // FIXME THIS IS FOR TEST ONLY
    // Choose whether to flip
    double de = mc_deltae(SiteIsProtonated, iflip, maxsite, wint, qunprot, self);
    if ( de < 0 )
      SiteIsProtonated.FlipProt(iflip);
    else if ( exp(-(beta_*de)) > rng.rn_gen() )
      SiteIsProtonated.FlipProt(iflip);
  }
  //mprintf("DEBUG: After thermal:");
  //for (unsigned int i = 0; i < maxsite; i++)
  //  mprintf(" %i", SiteIsProtonated[i]);
  //mprintf("\n");
  SiteIsProtonated.PrintState(logfile_);
  double energy = SiteIsProtonated.Esum(self, qunprot, wint);
  logfile_->Printf("Energy after thermal= %16.8f\n", energy);
  // Statistics loop: take 1 MC step and calulate averages.
  StateArray tempProt;
  std::string stateChar;
  stateChar.resize( maxsite );
  for (int mct = 0; mct < n_mc_steps_; mct++) {
    mc_step(energy,
            maxsite, self,
            pairs,
            SiteIsProtonated,
            wint, qunprot,
            rng);
    // Is the overall state different?
    bool isDifferentState = tempProt.AssignFromState( SiteIsProtonated );
    // DEBUG: Print if different
    if (isDifferentState) {
      SiteIsProtonated.StateStr( stateChar );
      logfile_->Printf("%12i State : %s Energy : %g\n", mct, stateChar.c_str(), energy);
    }

  }

  return 0;
}

/** Calculate titration curves using MC */
int Protonator::CalcTitrationCurves() const {
  // Calculate the ground energy of the system (no protons)
  DataSet_1D const& pkint = static_cast<DataSet_1D const&>( *site_intrinsic_pKas_ );
  DataSet_1D const& qunprot = static_cast<DataSet_1D const&>( *site_qunprot_ );
  DataSet_2D const& wint = static_cast<DataSet_2D const&>( *site_site_matrix_ );
  DataSet_string const& siteNames = static_cast<DataSet_string const&>( *site_names_ );
  unsigned int maxsite = pkint.Size();
  // TODO can we assume symmetric?
  double grounde = 0;
  for (unsigned int i = 0; i < maxsite; i++)
    for (unsigned int j = 0; j < maxsite; j++)
      grounde += (qunprot.Dval(i) * qunprot.Dval(j) * wint.GetElement(i,j));
  grounde *= 0.5;
  logfile_->Printf("grounde = %16.8f\n", grounde);

  // Calculate pairs.
  // Pairs of strongly interacting sites are allowed to simultaneously
  // change their protonation states.
  PairArray pairs;
  double min_g = min_wint_ * log(10.0) / beta_;
  logfile_->Printf("min_g = %16.8f\n", min_g);
  for (unsigned int i = 0; i < maxsite; i++) {
    for (unsigned int j = i+1; j < maxsite; j++) {
      if (wint.GetElement(i,j) > min_g)
        pairs.push_back( StatePair(i,j) );
    }
  }
  logfile_->Printf("npairs = %6zu\n", pairs.size());
  // Count the number of pH values
  int nph = (int)((stop_pH_ - start_pH_) / pH_increment_) + 1;
  std::vector<double> pH_values;
  pH_values.reserve( nph );
  double pH = start_pH_;
  while (pH <= stop_pH_) {
    pH_values.push_back( pH );
    pH += pH_increment_;
  }
  mprintf("DEBUG: maxsite= %u  #pH vals= %zu  nph= %i\n", maxsite, pH_values.size(), nph);
  logfile_->Printf("%6u%6i\n", maxsite, nph);
  for (unsigned int i = 0; i < maxsite; i++) {
    char siteType;
    if (qunprot.Dval(i) < 0)
      siteType = 'A';
    else
      siteType = 'C';
    logfile_->Printf("%10.5f %c %s\n", pkint.Dval(i), siteType, siteNames[i].c_str());
  }

  StateArray SiteIsProtonated(maxsite);

  // Loop over pH values
  for (Darray::const_iterator ph = pH_values.begin(); ph != pH_values.end(); ++ph)
  {
    Random_Number rng; // TODO use 1 overall rng
    // Assign initial protonation
    // FIXME calling set seed here to match mcti init subroutine. Should not always be done.
    rng.rn_set( iseed_ );
    SiteIsProtonated.AssignRandom( rng );
    logfile_->Printf("pH= %6.2f\n", *ph);
    SiteIsProtonated.PrintState(logfile_);
    //mprintf("Initial states (%i total):", SiteIsProtonated.Nprotonated());
    //for (Iarray::const_iterator it = SiteIsProtonated.begin(); it != SiteIsProtonated.end(); ++it)
    //  mprintf(" %i", *it);
    //mprintf("\n");
    // Do the MC trials at this pH
    int err =  perform_MC_at_pH(*ph, SiteIsProtonated,
                                pkint, wint, qunprot, rng, pairs);
    if (err != 0) {
      mprinterr("Error: Could not perform MC at pH %g\n", *ph);
      return 1;
    }
  } // END loop over pH values
  

  return 0;
}
