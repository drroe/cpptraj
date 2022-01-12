#ifndef INC_POTENTIALTERM_REPLICATE_H
#define INC_POTENTIALTERM_REPLICATE_H
#include <vector>
#include "PotentialTerm.h"
/// Hold individual terms for replicates
class PotentialTerm_Replicate : public PotentialTerm {
  public:
    PotentialTerm_Replicate();
    ~PotentialTerm_Replicate();
  private:
    typedef std::vector<PotentialTerm*> Parray;
    Parray REPTERM_; ///< Terms, one for each replicate
};
#endif
