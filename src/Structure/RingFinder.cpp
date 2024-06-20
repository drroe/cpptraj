#include "RingFinder.h"
#include "../AtomMask.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
RingFinder::RingFinder()
{}

/** Init ring finder. */
int RingFinder::InitRingFinder(ArgList& argIn) {
  return 0;
}

/** Set up ring finder for topology. */
int RingFinder::SetupRingFinder(Topology const& topIn) {
  rings_.clear();

  return 0;
}
