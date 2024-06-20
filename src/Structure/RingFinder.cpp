#include "RingFinder.h"
#include "../AtomMask.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
RingFinder::RingFinder()
{}

/** Init ring finder. */
int RingFinder::InitRingFinder(ArgList& argIn) {
  return 0;
}

static inline void visitAtom( int at, int previousAt, int res, Topology const& topIn, std::vector<bool>& Visited, int dist )
{
  if (Visited[at]) {
    mprintf("Already visited %s %i\n", topIn.AtomMaskName(at).c_str(), dist);
  } else {
    Visited[at] = true;
    Atom const& currentAt = topIn[at];
    for (Atom::bond_iterator bat = currentAt.bondbegin(); bat != currentAt.bondend(); ++bat)
    {
      if (*bat != previousAt && topIn[*bat].Nbonds() > 1 && topIn[*bat].ResNum() == res) {
        visitAtom( *bat, at, res, topIn, Visited, dist+1 );
      }
    }
  }
}

/** Set up ring finder for topology. */
int RingFinder::SetupRingFinder(Topology const& topIn) {
  rings_.clear();

  std::vector<bool> Visited(topIn.Natom(), false); // TODO only have visited for residue atoms
  // Do not look for rings that span residues.
  for (int res = 0; res < topIn.Nres(); res++) {
    Residue const& currentRes = topIn.Res(res);
    int idx = 0;
    for (int at = currentRes.FirstAtom(); at != currentRes.LastAtom(); at++, idx++)
    {
      Atom const& currentAt = topIn[at];
      if (currentAt.Nbonds() > 1) {
        mprintf("Starting at %s\n", topIn.AtomMaskName(at).c_str());
        Visited.assign( topIn.Natom(), false );
        Visited[at] = true;
        for (Atom::bond_iterator bat = currentAt.bondbegin(); bat != currentAt.bondend(); ++bat)
        {
          if (topIn[*bat].Nbonds() > 1 && topIn[*bat].ResNum() == res) {
            visitAtom( *bat, at, res, topIn, Visited, 1 );
          }
        }
      }
    }
  }

  return 0;
}
