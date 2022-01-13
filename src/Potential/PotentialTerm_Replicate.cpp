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

/** Clear terms. */
void PotentialTerm_Replicate::clearTerms() {
  for (Parray::iterator it = REPTERM_.begin(); it != REPTERM_.end(); ++it)
    delete *it;
}

/** Clear masks. */
void PotentialTerm_Replicate::clearMasks() {
  for (Carray::iterator it = REPMASK_.begin(); it != REPMASK_.end(); ++it)
    delete *it;
}

/** Clear energy arrays. */
void PotentialTerm_Replicate::clearEarrays() {
  for (Earray::iterator it = REPENE_.begin(); it != REPENE_.end(); ++it)
    delete *it;
}

/** DESTRUCTOR */
PotentialTerm_Replicate::~PotentialTerm_Replicate() {
  clearTerms();
  clearMasks();
  clearEarrays();
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
  clearTerms();
  opts_ = opts;
  if (addRepTerm( opts_ )) return 1;
  return 0;
}

/** Set up openmm terms. */
int PotentialTerm_Replicate::SetupTerm(Topology const& topIn, Box const& boxIn,
                                    CharMask const& maskIn, EnergyArray& earrayIn)
{
  clearMasks();
  clearEarrays();
  // Create a CharMask with selected atoms from each replicate
  for (unsigned int rep = 0; rep != topIn.Replicates().size(); rep++) {
    REPMASK_.push_back( new CharMask() );
    CharMask& repCmask = *(REPMASK_.back());
    for (int topidx = topIn.Replicates()[rep].RepUnit().Front();
             topidx != topIn.Replicates()[rep].RepUnit().Back(); ++topidx)
    {
      repCmask.AddAtom(maskIn.AtomInCharMask( topidx ));
    }
    mprintf("DEBUG: Replicate %i, %i selected.\n", rep, repCmask.Nselected());
    if (rep > 0) {
      // Already added the term for replicate 0 in InitTerm()
      if (addRepTerm( opts_ )) {
        mprinterr("Error: Initializing replicate %i\n", rep);
        return 1;
      }
    }
    // Create energy array for replicate
    REPENE_.push_back( new EnergyArray() );
    // Set up the replicate term
    if (REPTERM_.back()->SetupTerm( topIn, boxIn, repCmask, *(REPENE_.back()))) {
      mprinterr("Error: Setting up replicate %i\n", rep);
      return 1;
    }
  } // END loop over replicates

  // Set up overall energy term
  ene_ = earrayIn.AddType( EnergyArray::E_OPENMM );

  return 0;
}

/** Calculate force from openmm */
void PotentialTerm_Replicate::CalcForce(Frame& frameIn, CharMask const& maskIn) const
{
  for (unsigned int rep = 0; rep != REPTERM_.size(); ++rep) {
    CalcForce(frameIn, *(REPMASK_[rep]));
    mprintf("DEBUG: Replicate %i energy= %f \n", REPENE_[rep]->Ene(EnergyArray::E_OPENMM));
    *ene_ += REPENE_[rep]->Ene(EnergyArray::E_OPENMM);
  }
}
    
