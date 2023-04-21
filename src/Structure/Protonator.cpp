#include "Protonator.h"
#include "../ArgList.h"
#include "../CpptrajState.h"
#include "../CpptrajStdio.h"
#include "../DataSet_1D.h"
#include "../DataSet_2D.h"
#include "../DataSet_string.h"
#include "../DataSet_double.h"
#include "../DataSet_MatrixDbl.h"
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
  n_reduced_mc_steps_(0),
  start_pH_(0),
  stop_pH_(0),
  pH_increment_(0),
  min_wint_(0),
  fract_toler_(0),
  beta_(0),
  iseed_(0),
  mcmode_(MC_FULL),
  logfile_(0),
  pkoutfile_(0)
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
  n_reduced_mc_steps_ = argIn.getKeyInt("redsteps", 10000);
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
  pkoutfile_ = State.DFL().AddCpptrajFile( argIn.GetStringKey("pkout"),
                                           "PKout file",
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
  if (pkoutfile_ != 0) mprintf("\tPKout file: '%s'\n", pkoutfile_->Filename().full());
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
// -----------------------------------------------------------------------------
/** Class for holding correlation function results from MC. */
class Protonator::MC_Corr {
  public:
    typedef std::vector<Darray> D2array;
    typedef std::vector<Iarray> I2array;

    MC_Corr() : eave_(0), seave_(0), taumax_(100) {}

    void Init(unsigned int);

    void Update(double energy, StateArray const&);

    void Finish(int, unsigned int);

    void PrintDebug(CpptrajFile*) const;

    void PrintSiteStddev(CpptrajFile*) const;
    /// \return Average protonation for site (mcti aveprot)
    double SiteAvgProtonation(int idx) const { return aveprot_[idx+1]; }
    /// \return SD of correlation function (mcti prot_error) for site
    double SiteStddev(int idx) const { return prot_error_[idx+1]; }

    /// Set avg. protonation (mcti aveprot) and SD of corr. fxn (mcti prot_error) from reduced
    void SetFromReduced(MC_Corr const&, Iarray const&);
  private:
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

/** Init MC corr fxns. */
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

/** Used to update protonation corr fxns */
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

/** Update corr fxns during MC iterations. */
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

/** Finish corr fxns after MC iterations. */
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

/** Print debug info to stdout. */
void Protonator::MC_Corr::PrintDebug(CpptrajFile* logfile) const {
  for (int tau = 0; tau <= taumax_; tau++) {
    //logfile->Printf("Cenergy %6i%12.5E\n", tau, cenergy_[tau]);
    logfile->Printf("Senergy %6i%12.5f\n", tau, se_[tau]);
  }
  for (unsigned int j = 0; j < corr_.size(); j++) {
    for (int tau = 0; tau <= taumax_; tau++) {
      //logfile->Printf("Corr %6i%6i%12.5f\n", j, tau, corr_[j][tau]);
      logfile->Printf("Scorr %6i%6i%12.5f\n", j, tau, scorr_[j][tau]);
      logfile->Printf("aveprot %6i%12.5f\n", j, aveprot_[j]);
      logfile->Printf("save %6i%12.5f\n", j, sqave_[j]);
    }
  }
}

/** Print std dev to stdout. */
void Protonator::MC_Corr::PrintSiteStddev(CpptrajFile* logfile) const {
  for (unsigned int j = 0; j < prot_error_.size(); j++)
    logfile->Printf("Proterror %6i%12.5f\n", j, prot_error_[j]);
}

/** Set avg prot and error from reduced set */
void Protonator::MC_Corr::SetFromReduced(MC_Corr const& r_corr, Iarray const& ridx_to_site) {
  for (unsigned int ir = 0; ir != ridx_to_site.size(); ir++) {
    int isite = ridx_to_site[ir];
    double error = r_corr.SiteStddev(ir);
    if (error < 0.00000001) {
      aveprot_[isite+1] = r_corr.SiteAvgProtonation(ir);
      prot_error_[isite+1] = error;
    } else {
      double w1 = 1.0 / (SiteStddev(isite)*SiteStddev(isite));
      double w2 = 1.0 / (error * error);
      aveprot_[isite+1] = (w1*SiteAvgProtonation(isite) + w2*r_corr.SiteAvgProtonation(ir))/(w1+w2);
      prot_error_[isite+1] = sqrt(1.0/(w1+w2));
    }
  }
  aveprot_[0] = 0;
  prot_error_[0] = 0;
  for (unsigned int idx = 1; idx < aveprot_.size(); idx++) {
    aveprot_[0] = aveprot_[0] + aveprot_[idx];
    prot_error_[0] = prot_error_[0] + prot_error_[idx]*prot_error_[idx];
  }
  prot_error_[0] = sqrt( prot_error_[0] );
}

// -----------------------------------------------------------------------------

/** Perform monte carlo sampling for a given pH value.
  * Determine average values of the protonation for each site and
  * the correlation functions used to compute the energy.
  */
int Protonator::perform_MC_at_pH(double pH, StateArray& SiteIsProtonated,
                                 MC_Corr& corr,
                                 int mcsteps,
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
  // Thermalization: do mcsteps of Monte Carlo to approach equilibrium
  //                 before data is taken.
  // NOTE: This currently does not allow a 2 site transition.
  int maxsteps = mcsteps * maxsite;
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
  if (maxsite > 0) {
    for (int mct = 0; mct < mcsteps; mct++) {
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
  }
  corr.Finish(mcsteps, maxsite);

  // DEBUG
  corr.PrintDebug( logfile_ );

  //for (unsigned int j = 0; j <= maxsite; j++)
  //  logfile_->Printf("Proterror %6i%12.5f\n", j, corr.prot_error_[j]);
  return 0;
}

/** Reduce number of sites in MC calculation to only those which have
  * some factional protonation.
  */
int Protonator::reduce_sites(DataSet_double&   r_pkint, DataSet_1D&       r_qunprot, DataSet_2D&       r_Wint,
                             Iarray& ridx_to_site,
                             DataSet_1D const& pkint,   DataSet_1D const& qunprot,   DataSet_2D const& wint,
                             MC_Corr const& corr)
const
{
  int r_maxsite = 0;
  int nskipped = 0;
  for (unsigned int isite = 0; isite != pkint.Size(); isite++)
  {
    double aveprot = corr.SiteAvgProtonation(isite);
    logfile_->Printf("REDUCE %6u%16.8f%16.8f\n", isite+1, aveprot, fract_toler_);
    if (aveprot > 1 - fract_toler_ ||
        aveprot < fract_toler_)
    {
      nskipped++;
    } else {
      ridx_to_site.push_back( isite );
      double dval = qunprot.Dval(isite);
      r_qunprot.Add( r_maxsite, &dval );
      dval = pkint.Dval(isite);
      r_pkint.Add( r_maxsite, &dval );
      r_maxsite++;
    }
  }
  logfile_->Printf("N REDUCED SITES: %6i\n", r_maxsite);
  // Fill interaction matrix
  if (r_maxsite > 0)
    r_Wint.AllocateTriangle( r_maxsite );
  for (int j = 0; j < r_maxsite; j++) {
    int jsite = ridx_to_site[j];
    for (int i = j + 1; i < r_maxsite; i++) {
      int isite = ridx_to_site[i];
      r_Wint.SetElement(i, j, wint.GetElement(isite, jsite));
    }
  }
  // Adjust r_pkint for fixed charges
  for (int ir = 0; ir < r_maxsite; ir++) {
    int isite = ridx_to_site[ir];
    for (unsigned int jsite = 0; jsite < pkint.Size(); jsite++) {
      bool site_fixed = true;
      for (int jr = 0; jr < r_maxsite; jr++) {
        if ((unsigned int)ridx_to_site[jr] == jsite) site_fixed = false;
      }
      if (site_fixed) {
        r_pkint[ir] = r_pkint[ir]
                      - (qunprot.Dval(jsite) + corr.SiteAvgProtonation(jsite))
                      * ( (beta_*wint.GetElement(isite,jsite))/log(10.0) );
      }
    }
  }
  for (int ir = 0; ir < r_maxsite; ir++) {
    int isite = ridx_to_site[ir];
    logfile_->Printf("RSITE %6i%6i%12.5f%12.5f\n", ir+1, isite+1, r_pkint[ir], r_qunprot.Dval(ir));
  }
  for (int ir = 0; ir < r_maxsite; ir++) {
    for (int j = 0; j < r_maxsite; j++) {
      logfile_->Printf("RWINT %6i%6i%12.5f\n", ir+1, j+1, r_Wint.GetElement(ir, j));
    }
  }

  return 0;
}

/** Get pairs of strongly interacting sites. */
Protonator::PairArray Protonator::get_pairs(double min_g, unsigned int maxsite, DataSet_2D const& wint) const {
  PairArray pairs;
  for (unsigned int i = 0; i < maxsite; i++) {
    for (unsigned int j = i+1; j < maxsite; j++) {
      if (wint.GetElement(i,j) > min_g)
        pairs.push_back( StatePair(i,j) );
    }
  }
  return pairs;
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
  double min_g = min_wint_ * log(10.0) / beta_;
  logfile_->Printf("min_g = %16.8f\n", min_g);
  PairArray pairs = get_pairs(min_g, maxsite, wint);
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
  // Allocate arrays for storing titration data
  Darray pkhalf(maxsite, -999); // TODO is -999 a "safe" value?
  std::vector<Darray> pklist(maxsite, Darray(nph, 0));

  StateArray SiteIsProtonated(maxsite);

  Random_Number rng;
  rng.rn_set( iseed_ );
  // Loop over pH values
  for (Darray::const_iterator ph = pH_values.begin(); ph != pH_values.end(); ++ph)
  {
    // Assign initial protonation
    SiteIsProtonated.AssignRandom( rng );
    logfile_->Printf("ph= %6.2f\n", *ph);
    SiteIsProtonated.PrintState(logfile_);
    //mprintf("Initial states (%i total):", SiteIsProtonated.Nprotonated());
    //for (Iarray::const_iterator it = SiteIsProtonated.begin(); it != SiteIsProtonated.end(); ++it)
    //  mprintf(" %i", *it);
    //mprintf("\n");
    // Do the MC trials at this pH
    MC_Corr corr;
    int err =  perform_MC_at_pH(*ph, SiteIsProtonated, corr, n_mc_steps_,
                                pkint, wint, qunprot, rng, pairs);
    if (err != 0) {
      mprinterr("Error: Could not perform MC at pH %g\n", *ph);
      return 1;
    }
    corr.PrintSiteStddev( logfile_ );

    if (mcmode_ != MC_FULL) {
      // Reduced site techniques
      // Allocate arrays for reduced sites
      DataSet_double r_pkint;
      DataSet_double r_qunprot;
      DataSet_MatrixDbl r_Wint;
      Iarray ridx_to_site;
      if (reduce_sites(r_pkint, r_qunprot, r_Wint, ridx_to_site,
                       pkint, qunprot, wint, corr))
      {
        mprinterr("Error: Site reduction failed.\n");
        return 1;
      }

      if (mcmode_ == MC_REDUCED) {
        // Reduced monte carlo
        PairArray r_pairs = get_pairs(min_g, maxsite, wint);
        StateArray r_prot(r_pkint.Size());
        // Assign initial protonation for reduced set of sites
        r_prot.AssignRandom( rng );
        // Do MC trials at this pH for reduced set of sites
        MC_Corr r_corr;
        err = perform_MC_at_pH(*ph, r_prot, r_corr, n_reduced_mc_steps_,
                               r_pkint, r_Wint, r_qunprot, rng, r_pairs);
        if (err != 0) {
          mprinterr("Error: Could not perform reduced MC at pH %g\n", *ph);
          return 1;
        }
        // Set avg prot and error from reduced
        corr.SetFromReduced(r_corr, ridx_to_site);
        for (unsigned int ir = 0; ir != r_pkint.Size(); ir++) {
          int isite = ridx_to_site[ir];
          logfile_->Printf("RRESULT %6i%12.5f%12.5f\n",isite+1, corr.SiteAvgProtonation(isite), corr.SiteStddev(isite));
        }
        // Calc pKhalf, save titration data in array
        long int phidx = ph - pH_values.begin();
        for (unsigned int isite = 0; isite < maxsite; isite++) {
          pklist[isite][phidx] = corr.SiteAvgProtonation(isite);
          if (phidx > 0) {
            double diff1 = pklist[isite][phidx-1] - 0.5;
            double diff2 = pklist[isite][phidx  ] - 0.5;
            if (diff1 > 0 && diff2 <= 0) {
              pkhalf[isite] = *ph + ((*ph - pH_values[phidx-1]) * (diff2 / (diff1 - diff2)));
            }
          }
        }
        // END mcmode_ == REDUCED
      } else {
        mprinterr("Internal Error; Invalid MC mode.\n");
        return 1;
      } 
    }
  } // END loop over pH values

  // Write pkout file
  if (pkoutfile_ != 0) {
    // Pk and pkhalf for each site
    pkoutfile_->Printf("%7s      %5s          %8s   %5s\n",  "Residue", "pKint", "pK_(1/2)", "     ");
    for (unsigned int isite = 0; isite < maxsite; isite++) {
      if (pkhalf[isite] > -99)
        pkoutfile_->Printf("%-13s%7.3f       %7.3f\n", siteNames[isite].c_str(), pkint.Dval(isite), pkhalf[isite]);
      else {
        if (pklist[isite][0] < 0.5)
          pkoutfile_->Printf("%-13s%7.3f       %c%8.3f\n", siteNames[isite].c_str(), pkint.Dval(isite), '<', pH_values.front());
        else
          pkoutfile_->Printf("%-13s%7.3f       %c%8.3f\n", siteNames[isite].c_str(), pkint.Dval(isite), '>', pH_values.back());
      }
    }
    // Curve for each site
    for (unsigned int isite = 0; isite < maxsite; isite++) {
      pkoutfile_->Printf(" Site %s\n", siteNames[isite].c_str());
      for (int phidx = 0; phidx < nph; phidx++) {
        double dval = pklist[isite][phidx];
        if (dval > 0.05 && dval < 0.95)
          pkoutfile_->Printf("%10.5f %9.5f       %s\n", pH_values[phidx], dval, siteNames[isite].c_str());
      }
    }

  } // END pkout write 

  return 0;
}
