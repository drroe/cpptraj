#ifndef INC_STRUCTURE_PROTONATOR_STATEARRAY_H
#define INC_STRUCTURE_PROTONATOR_STATEARRAY_H
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

#endif
