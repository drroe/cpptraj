#include "GenerateAngles.h"
#include "../CpptrajStdio.h"
#include "../ParameterTypes.h"
#include "../Topology.h"
#include <algorithm> // std::sort

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

  return 0;
}
