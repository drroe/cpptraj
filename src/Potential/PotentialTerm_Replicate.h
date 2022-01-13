#ifndef INC_POTENTIALTERM_REPLICATE_H
#define INC_POTENTIALTERM_REPLICATE_H
#include <vector>
#include "PotentialTerm.h"
#include "MdOpts.h"
/// Hold individual terms for replicates
class PotentialTerm_Replicate : public PotentialTerm {
  public:
    PotentialTerm_Replicate();
    ~PotentialTerm_Replicate();

    int InitTerm(MdOpts const&);
    int SetupTerm(Topology const&, Box const&, CharMask const&, EnergyArray&);
    void CalcForce(Frame&, CharMask const&) const;
  private:
    /// Clear terms
    void clearTerms();
    /// Clear masks
    void clearMasks();
    /// Clear energy arrays
    void clearEarrays();
    /// Add a term for replicate with given options
    int addRepTerm(MdOpts const&);

    typedef std::vector<PotentialTerm*> Parray;
    typedef std::vector<CharMask*> Carray;
    typedef std::vector<EnergyArray*> Earray;
    Parray REPTERM_; ///< Individual energy Terms, one for each replicate
    Carray REPMASK_; ///< Selected atoms for each replicate
    Earray REPENE_;  ///< Energy arrays, one for each replicate
    MdOpts opts_;    ///< Options for each term
    double* ene_;    ///< Pointer into overall energy array
};
#endif
