#ifndef STRUCTURE_PROTONATOR_MC_CORR_H
#define STRUCTURE_PROTONATOR_MC_CORR_H
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

#endif
