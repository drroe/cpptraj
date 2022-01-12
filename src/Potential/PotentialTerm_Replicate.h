#ifndef INC_POTENTIALTERM_REPLICATE_H
#define INC_POTENTIALTERM_REPLICATE_H
#include <vector>
#include "PotentialTerm.h"
#include "MdOpts.h"
class Topology;
/// Hold individual terms for replicates
class PotentialTerm_Replicate : public PotentialTerm {
  public:
    PotentialTerm_Replicate();
    ~PotentialTerm_Replicate();

    int InitTerm(MdOpts const&);
    int SetupTerm(Topology const&, Box const&, CharMask const&, EnergyArray&);
  private:
    /// Add a term for replicate with given options
    int addRepTerm(MdOpts const&);

    typedef std::vector<PotentialTerm*> Parray;
    typedef std::vector<Topology*> Tarray;
    typedef std::vector<EnergyArray*> Earray;
    Parray REPTERM_; ///< Terms, one for each replicate
    Tarray REPTOPS_; ///< Topologies, one for each replicate
    Earray REPENE_;  ///< Energy arrays, one for each replicate
    MdOpts opts_;    ///< Options for each term
};
#endif
