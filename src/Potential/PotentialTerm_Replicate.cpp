#include "PotentialTerm_Replicate.h"
#include "PotentialTerm_OpenMM.h"
#include "EnergyArray.h"
#include "../Topology.h"
#include "../CpptrajStdio.h"
#include "../StringRoutines.h"
#include "../CharMask.h"

/** CONSTRUCTOR */
PotentialTerm_Replicate::PotentialTerm_Replicate() :
  PotentialTerm(REPLICATE)
{}

/** DESTRUCTOR */
PotentialTerm_Replicate::~PotentialTerm_Replicate() {
  for (Parray::iterator it = REPTERM_.begin(); it != REPTERM_.end(); ++it)
    delete *it;
  for (Tarray::iterator it = REPTOPS_.begin(); it != REPTOPS_.end(); ++it)
    delete *it;
  for (Earray::iterator it = REPENE_.begin(); it != REPENE_.end(); ++it)
    delete *it;
}

/** Add term with given options */
int PotentialTerm_Replicate::addRepTerm(MdOpts const& opts) {
  PotentialTerm_OpenMM* pterm = new PotentialTerm_OpenMM();

  if (pterm->InitTerm( opts )) {
    delete pterm;
    return 1;
  }

  REPTERM_.push_back( (PotentialTerm*)pterm );
  return 0;
}

/** Init openmm options. */
int PotentialTerm_Replicate::InitTerm(MdOpts const& opts) {
  opts_ = opts;
  if (addRepTerm( opts_ )) return 1;
  return 0;
}

/** Set up openmm terms. */
int PotentialTerm_Replicate::SetupTerm(Topology const& topIn, Box const& boxIn,
                                    CharMask const& maskIn, EnergyArray& earrayIn)
{
  // FIXME this is the big kludge. Split off into n different topologies. Yuck.
  REPTOPS_.clear();
  REPENE_.clear();
  for (unsigned int rep = 0; rep != topIn.Replicates().size(); rep++) {
    // Atoms selected in replica
    CharMask cmask;
    for (int idx = topIn.Replicates()[rep].RepUnit().Front();
             idx != topIn.Replicates()[rep].RepUnit().Back(); ++idx)
      cmask.AddAtom(maskIn.AtomInCharMask( idx ));
    // Atoms in replica
    AtomMask tempMask( topIn.Replicates()[rep].RepUnit().Front(),
                       topIn.Replicates()[rep].RepUnit().Back() );
    mprintf("DEBUG: rep %i, %i selected.\n", rep, tempMask.Nselected());
    Topology* repTop = topIn.modifyStateByMask( tempMask );
    repTop->Brief( std::string("Replicate " + integerToString(rep)).c_str() );
    REPTOPS_.push_back( repTop );
    REPENE_.push_back( new EnergyArray() );
    if (rep > 0) {
      if (addRepTerm( opts_ )) {
        mprinterr("Error: Initializing replicate %i\n", rep);
        return 1;
      }
    }
    if (REPTERM_.back()->SetupTerm( *(REPTOPS_.back()), boxIn, cmask, *(REPENE_.back()))) {
      mprinterr("Error: Setting up replicate %i\n", rep);
      return 1;
    }
  }

  return 1; // FIXME
}
