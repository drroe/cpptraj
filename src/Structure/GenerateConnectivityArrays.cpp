#include "GenerateConnectivityArrays.h"
#include "../ParameterTypes.h"
#include "../Atom.h"
#include "../Residue.h"
#include "../CpptrajStdio.h"

/** From atom connectivity, generate a bond array in the same order as LEaP. */ // TODO use in GenerateBAT
BondArray Cpptraj::Structure::GenerateBondArray(std::vector<Residue> const& residues,
                                                std::vector<Atom> const& atoms)
{
  BondArray out;
  // BONDS
  int bidx = 0;
  for (std::vector<Residue>::const_iterator res = residues.begin(); res != residues.end(); ++res)
  {
    for (int iat = res->LastAtom()-1; iat >= res->FirstAtom(); iat--)
    {
      Atom const& At = atoms[iat];
      for (Atom::bond_iterator bat = At.bondbegin(); bat != At.bondend(); ++bat)
      {
        if (iat < *bat) {
          mprintf("DEBUG: BOND  i= %i  %i - %i (%i %i)\n",  bidx++, iat+1, *bat+1, iat*3, *bat*3);
          out.push_back( BondType(iat, *bat, -1) );
        }
        //else
        //  mprintf("DEBUG: X    i= %i  %i - %i (%i %i)\n",   bidx++, iat+1, *bat+1, iat*3, *bat*3);
      }
    }
  }
  return out;
}


