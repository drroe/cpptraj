#include "GenerateAngles.h"
#include "../CpptrajStdio.h"
#include "../ParameterTypes.h"
#include "../Topology.h"
#include <algorithm> // std::sort

static inline void enumerateAngles(int at1, int at2, Topology const& topIn) {
  Atom const& A2 = topIn[at2];
  if (A2.Nbonds() > 1) {
    for (Atom::bond_iterator bat = A2.bondbegin(); bat != A2.bondend(); ++bat)
    {
      if (*bat != at1 && at1 < *bat) {
        mprintf("%s - %s - %s\n",
                topIn.AtomMaskName(at1).c_str(),
                topIn.AtomMaskName(at2).c_str(),
                topIn.AtomMaskName(*bat).c_str());
      }
    }
  }
}

int Cpptraj::Structure::GenerateAngles(Topology& topIn) {
  if (topIn.Nbonds() < 1) {
    mprintf("Warning: No bonds in '%s', no angles to generate.\n", topIn.c_str());
    return 0;
  }
  // Create a combined bonds array
  BondArray allBonds;
  allBonds.reserve(topIn.Nbonds());

  for (BondArray::const_iterator it = topIn.Bonds().begin(); it != topIn.Bonds().end(); ++it)
    allBonds.push_back( *it );
  for (BondArray::const_iterator it = topIn.BondsH().begin(); it != topIn.BondsH().end(); ++it)
    allBonds.push_back( *it );
  std::sort( allBonds.begin(), allBonds.end() );

  mprintf("DEBUG: Sorted bonds:\n");
  for (BondArray::const_iterator it = allBonds.begin(); it != allBonds.end(); ++it)
    mprintf("\t%s - %s\n",
            topIn.AtomMaskName(it->A1()).c_str(),
            topIn.AtomMaskName(it->A2()).c_str());

  // Angles
  AngleArray angles;
  AngleArray anglesH;
  for (BondArray::const_iterator it = allBonds.begin(); it != allBonds.end(); ++it)
  {
    // Forward direction. A1-A2-X
    enumerateAngles( it->A1(), it->A2(), topIn );
    // Reverse direction. A2-A1-X
    enumerateAngles( it->A2(), it->A1(), topIn );
  }
  return 0;
}
