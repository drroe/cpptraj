#include "PotentialTerm_Replicate.h"

/** CONSTRUCTOR */
PotentialTerm_Replicate::PotentialTerm_Replicate() :
  PotentialTerm(REPLICATE)
{}

/** DESTRUCTOR */
PotentialTerm_Replicate::~PotentialTerm_Replicate() {
  for (Parray::iterator it = REPTERM_.begin(); it != REPTERM_.end(); ++it)
    delete *it;
}
