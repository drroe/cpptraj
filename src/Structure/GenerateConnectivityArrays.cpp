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

/** From atom connectiviy, generate an angle array in the same order as LEaP. */ // TODO use in GenerateBAT
AngleArray Cpptraj::Structure::GenerateAngleArray(std::vector<Residue> const& residues,
                                                  std::vector<Atom> const& atoms)
{
  AngleArray out;
  // ANGLES TODO combine above
  int aidx = 0;
  for (std::vector<Residue>::const_iterator res = residues.begin(); res != residues.end(); ++res)
  {
    for (int iat1 = res->LastAtom()-1; iat1 >= res->FirstAtom(); iat1--)
    {
      Atom const& At1 = atoms[iat1];
      for (int bidx1 = 0; bidx1 < At1.Nbonds(); bidx1++) {
        int iat2 = At1.Bond(bidx1);
        Atom const& At2 = atoms[iat2];
        for (int bidx2 = 0; bidx2 < At2.Nbonds(); bidx2++) {
          int iat3 = At2.Bond(bidx2);
          if (iat1 < iat3) {
            mprintf("DEBUG: ANGLE  i= %i  %i - %i - %i (%i %i %i)\n", aidx++, iat1+1, iat2+1, iat3+1, iat1*3, iat2*3, iat3*3);
            out.push_back( AngleType(iat1, iat2, iat3, -1) );
          }
        }
      }
    }
  }
  return out;
}
