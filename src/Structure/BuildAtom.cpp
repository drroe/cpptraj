#include "BuildAtom.h"
#include "Chirality.h"
#include "../Constants.h"
#include "../CpptrajStdio.h"
#include "../Frame.h"
#include "../Topology.h"
#include "../TorsionRoutines.h"
#include <algorithm> // std::sort

/** Determine chirality and priority around an atom. */
int Cpptraj::Structure::BuildAtom::DetermineChirality(int atnum, Topology const& topIn,
                                                      Frame const& frameIn, int debugIn)
{
  priority_.clear();
  ctype_ = Cpptraj::Structure::SetPriority(priority_, atnum, topIn, frameIn, debugIn);
  if (ctype_ == CHIRALITY_ERR)
    return 1;
  return 0;
}

/** Determine only priority around an atom. */
int Cpptraj::Structure::BuildAtom::SetPriority(int atnum, Topology const& topIn,
                                               Frame const& frameIn, int debugIn)
{
  priority_.clear();
  if (Cpptraj::Structure::SetPriority(priority_, atnum, topIn, frameIn, debugIn) == CHIRALITY_ERR)
    return 1;
  return 0;
}
