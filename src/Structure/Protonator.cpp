#include "Protonator.h"
#include "../ArgList.h"
#include "../CpptrajState.h"
#include "../CpptrajStdio.h"
#include "../DataSet_1D.h"
#include "../DataSet_2D.h"
#include "../DataSet_string.h"
#include "../Random.h"
#include "../StringRoutines.h"
#include "../Mead/MultiFlexResults.h"
#include <cmath> // log, fabs, sqrt

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

/** For testing, read data from .pkint and .g files. */
int Protonator::read_files(CpptrajState& State, std::string const& prefix)
{
  CpptrajFile infile;
  // Pkint - site_intrinsic_pKas_, site_qunprot_, site_names_
  std::string pkintfile = prefix + ".pkint";
  if (infile.OpenRead( pkintfile )) {
    mprinterr("Error: Opening file %s\n", pkintfile.c_str());
    return 1;
  }
  site_intrinsic_pKas_ = State.DSL().AddSet( DataSet::DOUBLE, MetaData(pkintfile, "pkint") );
  site_qunprot_ = State.DSL().AddSet(DataSet::DOUBLE, MetaData(pkintfile, "qunprot") );
  site_names_ = State.DSL().AddSet(DataSet::STRING, MetaData(pkintfile, "names" ));
  const char* line = infile.NextLine();
  int isite = 0;
  while (line != 0) {
    ArgList argline(line, " ");
    double pkint = convertToDouble(argline[0]);
    double q;
    if (argline[1] == "A")
      q = -1;
    else
      q = 0;
    site_intrinsic_pKas_->Add( isite, &pkint );
    site_qunprot_->Add( isite, &q );
    site_names_->Add( isite, argline[2].c_str() );
    isite++;
    line = infile.NextLine();
  }
  infile.CloseFile();
  // G - site_site_matrix_
  std::string gfile = prefix + ".g";
  if (infile.OpenRead( gfile )) {
    mprinterr("Error: Opening file %s\n", gfile.c_str());
    return 1;
  }
  site_site_matrix_ = State.DSL().AddSet(DataSet::MATRIX_DBL, gfile);
  DataSet_2D& mat = static_cast<DataSet_2D&>( *site_site_matrix_ );
  mat.AllocateTriangle( isite );
  line = infile.NextLine();
  while (line != 0) {
    ArgList argline(line, " ");
    int ii = convertToInteger(argline[0]);
    int jj = convertToInteger(argline[1]);
    if (ii < jj)
      mat.SetElement(ii-1, jj-1, convertToDouble(argline[2]));
    line = infile.NextLine();
  }
  // DEBUG print matrix
  for (int ii = 0; ii < isite; ii++)
    for (int jj = 0; jj < isite; jj++)
      mprintf("%i %i %g\n", ii+1, jj+1, mat.GetElement(ii, jj));

  return 0;
}

/** Set up options */
int Protonator::SetupProtonator(CpptrajState& State, ArgList& argIn,
                                Cpptraj::Mead::MultiFlexResults const& results)
{
  std::string mcprefix = argIn.GetStringKey("mcprefix");
  if (!mcprefix.empty()) {
    mprintf("\tLoading previous MEAD results for MC using prefix '%s'\n", mcprefix.c_str());
    if (read_files(State, mcprefix)) return 1;
  } else {
    site_intrinsic_pKas_ = results.PkIntSet();
    site_site_matrix_ = results.SiteSiteMatrixSet();
    site_qunprot_ = results.QunprotSet();
    site_names_ = results.SiteNamesSet();
  }
  if (site_intrinsic_pKas_ == 0) {
    mprinterr("Internal Error: Site intrinsic pKa data set is missing.\n");
    return 1;
  }
  if (site_site_matrix_ == 0) {
    mprinterr("Internal Error: Site-site interaction matrix data set is missing.\n");
    return 1;
  }
  if (site_qunprot_ == 0) {
    mprinterr("Internal Error: Site unprotonated state charge data set is missing.\n");
    return 1;
  }
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
    /// \return Total # of sites
    unsigned int Nsites() const { return prot_.size(); }
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

/** Attempt a simultaneous change in the protonation state of 2
  * sites, and accept the change based on the standard
  * Metropolis criteria.
  */
void Protonator::mc_pair_flip(double& energy,
                              unsigned int maxsite, Darray const& self,
                              StatePair const& pairIn,
                              StateArray& SiteIsProtonated,
                              DataSet_2D const& wint, DataSet_1D const& qunprot,
                              Random_Number const& rng)
const
{
  unsigned int pair[2];
  pair[0] = pairIn.first;
  pair[1] = pairIn.second;
  // New potential protonation states of sites
  int newprot[2];
  // Change in protonation state of sites in pair
  int dp[2];
  for (int i = 0; i != 2; i++) {
    // FIXME use rng.rn_num_interval(0, 1)
    newprot[i] = (int)(2*rng.rn_gen());
    if (newprot[i] > 1) { // FIXME will not be necessary with above interval
      mprinterr("Error: bad pair flip.\n");
      return;
    }
    // -1 deprotonated, 1 protonated, 0 no change
    dp[i] = newprot[i] - SiteIsProtonated[ pair[i] ];
  }

  if (dp[0] == 0 && dp[1] == 0) return;

  double de = dp[0]*self[ pair[0] ] + dp[1]*self[ pair[1] ];

  // Change in self energy
  for (unsigned int i = 0; i < maxsite; i++) {
    if (i != pair[0] && i != pair[1]) {
      de = de + dp[0] * (qunprot.Dval(i) + SiteIsProtonated[i]) * wint.GetElement(i, pair[0]) +
                dp[1] * (qunprot.Dval(i) + SiteIsProtonated[i]) * wint.GetElement(i, pair[1]);
    }
  }

  // Change in interaction energy
  de = de + wint.GetElement(pair[0], pair[1]) * (
                 (qunprot.Dval(pair[0]) + newprot[0])*
                 (qunprot.Dval(pair[1]) + newprot[1]) -
                 (qunprot.Dval(pair[0]) + SiteIsProtonated[pair[0]])*
                 (qunprot.Dval(pair[1]) + SiteIsProtonated[pair[1]]) );

  // Check for flip
  if (de < 0) {
    if (dp[0] != 0) SiteIsProtonated.FlipProt(pair[0]);
    if (dp[1] != 0) SiteIsProtonated.FlipProt(pair[1]);
    energy = energy + de;
  } else if ( exp(-(beta_*de)) > rng.rn_gen() ) {
    if (dp[0] != 0) SiteIsProtonated.FlipProt(pair[0]);
    if (dp[1] != 0) SiteIsProtonated.FlipProt(pair[1]);
    energy = energy + de;
  }
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
      mc_pair_flip(energy, maxsite, self, pairs[iflip], SiteIsProtonated,
                   wint, qunprot, rng);
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
// -----------------------------------------------
/** Class for holding correlation function results from MC. */
class Protonator::MC_Corr {
  public:
    typedef std::vector<Darray> D2array;
    typedef std::vector<int> Iarray;
    typedef std::vector<Iarray> I2array;

    MC_Corr() : eave_(0), seave_(0), taumax_(100) {}

    void Init(unsigned int);

    void Update(double energy, StateArray const&);

    void Finish(int, unsigned int);

    double get_err(Darray const&, int) const;

    double eave_; ///< Average energy
    double seave_;
    int taumax_; ///< Time for correlation function
    Darray cenergy_; ///< Correlation function for energy [taumax+1]
    Darray se_; // [taumax+1]
    Darray eold_; // [taumax+1]
    Darray aveprot_; ///< Average total protonation [maxsite+1]
    Darray sqave_; // mcti save [maxsite+1]
    Darray prot_error_; // [maxsite+1]
    D2array corr_; // mcti c [maxsite+1][taumax+1]
    D2array scorr_; // mcti s [maxsite+1][taumax+1]
    I2array iold_;  // [maxsite+1][taumax+1]
};

void Protonator::MC_Corr::Init(unsigned int maxsite) {
  eave_ = 0;
  seave_ = 0;
  int taumax1 = taumax_+1;
  cenergy_.assign(taumax1, 0);
  se_.assign(     taumax1, 0);
  eold_.assign(   taumax1, 0);
  unsigned int maxsite1 = maxsite  + 1;
  aveprot_.assign(    maxsite1, 0);
  sqave_.assign(      maxsite1, 0);
  prot_error_.assign( maxsite1, 0);
  corr_.assign(       maxsite1, Darray(taumax1, 0));
  scorr_.assign(      maxsite1, Darray(taumax1, 0));
  iold_.assign(       maxsite1, Iarray(taumax1, 0));
}

static inline void update_prot(double& aveprot, double& sqave,
                               std::vector<double>& scorr,
                               std::vector<int>& iold, int prot, int taumax)
{
  aveprot = aveprot + prot;
  if (aveprot < 0) {
    mprintf("Warning: Average protonation is negative: avg= %g  prot= %i\n", aveprot, prot);
  }
  sqave = sqave + prot;
  scorr[0] = scorr[0] + (prot * prot);
  for (int tau = 1; tau <= taumax; tau++)
    scorr[tau] = scorr[tau] + (prot * iold[tau]);
  for (int tau = taumax; tau > 1; tau--)
    iold[tau] = iold[tau-1];
  iold[1] = prot;
}

void Protonator::MC_Corr::Update(double energy, StateArray const& SiteIsProtonated)
{
   // Update the correlation functions
    eave_ = eave_ + energy;
    seave_ = seave_ + energy;
    se_[0] = se_[0] + (energy * energy);
    for (int tau = 1; tau <= taumax_; tau++) {
      se_[tau] = se_[tau] + (energy * eold_[tau]);
    }
    for (int tau = taumax_; tau > 1; tau--) {
      eold_[tau] = eold_[tau - 1];
    }
    eold_[1] = energy;
  // Update total protonation
  update_prot(aveprot_[0], sqave_[0], scorr_[0], iold_[0], SiteIsProtonated.Nprotonated(), taumax_);
  // Update protonation for each site

  for (unsigned int site = 0; site != SiteIsProtonated.Nsites(); site++) {
    unsigned int idx = site + 1;
    update_prot(aveprot_[idx], sqave_[idx], scorr_[idx], iold_[idx], SiteIsProtonated[site], taumax_);
  }
}

/** Compute the standard deviation from given correlation fxn. If the
  * correlation does not fall to 1/10 the value of c[0] then assign 999
  * to the error. */
double Protonator::MC_Corr::get_err(Darray const& cc, int mcstepsIn)
const
{
  int ict = -1;
  bool flag = false;
  double error = 0;
  double ctol = 0.1*cc[0];

  // If cc[0] is already small, the correlation is 0
  if (fabs(cc[0]) < 0.00001) {
    ict = 0;
    flag = true;
  }

  for (int i = 1; i <= taumax_; i++) {
    if ( cc[i] < ctol && !flag ) {
      ict = i;
      flag = true;
    }
  }
  if (!flag)
    error = 999;
  else
    // mcsteps / ict = num. of statistically independent samples
    error = sqrt( (cc[0]*ict) / (double)mcstepsIn );
  return error;
}

void Protonator::MC_Corr::Finish(int mcstepsIn, unsigned int maxsiteIn) {
  double mcsteps = (double)mcstepsIn;
  eave_ = eave_ / mcsteps;
  seave_ = (seave_/mcsteps) * (seave_/mcsteps);
  for (int tau = 0; tau <= taumax_; tau++)
    cenergy_[tau] = ((1.0 / (mcsteps - (double)tau)) * se_[tau]) - seave_;

  unsigned int maxsite1 = maxsiteIn+1;
  for (unsigned int idx = 0; idx < maxsite1; idx++) {
    aveprot_[idx] = aveprot_[idx] / mcsteps;
    sqave_[idx]   = (sqave_[idx]/mcsteps) * (sqave_[idx]/mcsteps);
    for (int tau = 0; tau <= taumax_; tau++) {
      corr_[idx][tau] = ((1.0 / (mcsteps - (double)tau)) * scorr_[idx][tau]) - sqave_[idx];
    }
    prot_error_[idx] = get_err(corr_[idx], mcstepsIn);
  }
}


// -----------------------------------------------

/** Perform monte carlo sampling for a given pH value.
  * Determine average values of the protonation for each site and
  * the correlation functions used to compute the energy.
  */
int Protonator::perform_MC_at_pH(double pH, StateArray& SiteIsProtonated,
                                 MC_Corr& corr,
                                 DataSet_1D const& pkint,
                                 DataSet_2D const& wint, DataSet_1D const& qunprot,
                                 Random_Number const& rng,
                                 PairArray const& pairs)
const
{
  unsigned int maxsite = pkint.Size();
  // Initialize correlation variables
  corr.Init(maxsite);
  // Compute the self energy to protonate each site at this pH
  Darray self;
  self.reserve( maxsite );
  double fac = -(log(10.0) / beta_);
  //mprintf("DEBUG: MC self E:");
  for (unsigned int i = 0; i < maxsite; i++) {
    self.push_back( fac * (pkint.Dval(i) - pH) );
    //logfile_->Printf("%26.16E", self[i]);
    logfile_->Printf("%16.8f\n", self[i]);
  }
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
      logfile_->Printf("%12i State : %s\n", mct+1, stateChar.c_str());
      logfile_->Printf("Energy : %16.8f\n", energy);
    }
    corr.Update(energy, SiteIsProtonated);
  }
  corr.Finish(n_mc_steps_, maxsite);

  // DEBUG
  for (int tau = 0; tau <= corr.taumax_; tau++) {
    logfile_->Printf("Cenergy %6i%12.5E\n", tau, corr.cenergy_[tau]);
    logfile_->Printf("Senergy %6i%12.5E\n", tau, corr.se_[tau]);
  }
  for (unsigned int j = 0; j <= maxsite; j++) {
    for (int tau = 0; tau <= corr.taumax_; tau++) {
      logfile_->Printf("Corr %6i%6i%12.5f\n", j, tau, corr.corr_[j][tau]);
      logfile_->Printf("Scorr %6i%6i%12.5f\n", j, tau, corr.scorr_[j][tau]);
      logfile_->Printf("aveprot %6i%12.5f\n", j, corr.aveprot_[j]);
      logfile_->Printf("save %6i%12.5f\n", j, corr.sqave_[j]);
    }
  }
  for (unsigned int j = 0; j <= maxsite; j++)
    logfile_->Printf("Proterror %6i%12.5f\n", j, corr.prot_error_[j]);
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
    logfile_->Printf("ph= %6.2f\n", *ph);
    SiteIsProtonated.PrintState(logfile_);
    //mprintf("Initial states (%i total):", SiteIsProtonated.Nprotonated());
    //for (Iarray::const_iterator it = SiteIsProtonated.begin(); it != SiteIsProtonated.end(); ++it)
    //  mprintf(" %i", *it);
    //mprintf("\n");
    // Do the MC trials at this pH
    MC_Corr corr;
    int err =  perform_MC_at_pH(*ph, SiteIsProtonated, corr,
                                pkint, wint, qunprot, rng, pairs);
    if (err != 0) {
      mprinterr("Error: Could not perform MC at pH %g\n", *ph);
      return 1;
    }
  } // END loop over pH values
  

  return 0;
}
