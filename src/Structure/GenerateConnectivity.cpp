#include "GenerateConnectivity.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"
#include <algorithm> // std::sort

static inline void enumerateAngles(int& aidx, int at1, int at2, Topology& topIn) {
  Atom const& A2 = topIn[at2];
  if (A2.Nbonds() > 1) {
    for (Atom::bond_iterator bat = A2.bondbegin(); bat != A2.bondend(); ++bat)
    {
      if (*bat != at1 && at1 < *bat) {
        topIn.AddAngle(at1, at2, *bat);
        mprintf("DEBUG: ANGLE  i= %i  %i - %i - %i (%i %i %i)\n",  aidx++, at1+1, at2+1, *bat+1, at1*3, at2*3, *bat*3);
        //mprintf("%s - %s - %s\n",
        //        topIn.AtomMaskName(at1).c_str(),
        //        topIn.AtomMaskName(at2).c_str(),
        //        topIn.AtomMaskName(*bat).c_str());
      }
    }
  }
}

static inline void enumerateDihedrals(int at1, int at2, Topology& topIn) {
  Atom const& A1 = topIn[at1];
  Atom const& A2 = topIn[at2];
  if (A1.Nbonds() > 1 && A2.Nbonds() > 1) {
    for (Atom::bond_iterator bat1 = A1.bondbegin(); bat1 != A1.bondend(); ++bat1)
    {
      if (*bat1 != at2) {
        for (Atom::bond_iterator bat2 = A2.bondbegin(); bat2 != A2.bondend(); ++bat2)
        {
          if (*bat2 != at1) {
            // LEaP convention appears to be first atom less than last atom
            if (*bat1 < *bat2)
              topIn.AddDihedral(*bat1, at1, at2, *bat2);
            else
              topIn.AddDihedral(*bat2, at2, at1, *bat1);
            //mprintf("%s - %s - %s - %s\n",
            //        topIn.AtomMaskName(*bat1).c_str(),
            //        topIn.AtomMaskName(at1).c_str(),
            //        topIn.AtomMaskName(at2).c_str(),
            //        topIn.AtomMaskName(*bat2).c_str());
          }
        }
      }
    }
  }
}

/* Set bonds, angles, and dihedral arrays for a given topology based on
 * a sorted and combined version of the current bond arrays.
 */
int Cpptraj::Structure::GenerateAngleTorsionArrays(Topology& topIn) {
  if (topIn.Nbonds() < 1) {
    mprintf("Warning: No bonds in '%s', no angles to generate.\n", topIn.c_str());
    return 0;
  }
  bool has_angles = topIn.Nangles() > 0;
  if (has_angles) {
    mprintf("Warning: Topology '%s' already has angle information.\n", topIn.c_str());
    //return 0;
  }
  bool has_dihedrals = topIn.Ndihedrals() > 0;
  if (has_dihedrals) {
    mprintf("Warning: Topology '%s' already has dihedral information.\n", topIn.c_str());
    //return 0;
  }
  if (has_angles && has_dihedrals) return 0;
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

  // Angles and Torsiions
  int aidx = 0;
  for (BondArray::const_iterator it = allBonds.begin(); it != allBonds.end(); ++it)
  {
    if (!has_angles) {
      // Forward direction. A1-A2-X
      enumerateAngles( aidx, it->A1(), it->A2(), topIn );
      // Reverse direction. A2-A1-X
      enumerateAngles( aidx, it->A2(), it->A1(), topIn );
    }
    if (!has_dihedrals) {
      // Dihedrals
      enumerateDihedrals( it->A1(), it->A2(), topIn );
    }
  }
  return 0;
}

/* Set bonds, angles, and dihedral arrays for a given topology based on 
 * current atom connectivity.
 * This is done in the same manner as LEaP, which goes in increasing
 * residue order but decreasing atom index.
 */
int Cpptraj::Structure::GenerateBondAngleTorsionArrays(Topology& topIn) {
  if (topIn.Nbonds() < 1) {
    mprintf("Warning: No bonds in '%s', no angles to generate.\n", topIn.c_str());
    return 0;
  }
  bool has_angles = topIn.Nangles() > 0;
  if (has_angles) {
    mprintf("Warning: Topology '%s' already has angle information.\n", topIn.c_str());
    //return 0;
  }
  bool has_dihedrals = topIn.Ndihedrals() > 0;
  if (has_dihedrals) {
    mprintf("Warning: Topology '%s' already has dihedral information.\n", topIn.c_str());
    //return 0;
  }
  if (has_angles && has_dihedrals) return 0;
  // Clear existing bond information TODO clear angles and dihedrals?
  topIn.ClearBondArrays();

  // BONDS
  int bidx = 0;
  for (int ires = 0; ires < topIn.Nres(); ires++)
  {
    Residue const& res = topIn.Res(ires);
    for (int iat = res.LastAtom()-1; iat >= res.FirstAtom(); iat--)
    {
      Atom const& At = topIn[iat];
      for (Atom::bond_iterator bat = At.bondbegin(); bat != At.bondend(); ++bat)
      {
        if (iat < *bat) {
          mprintf("DEBUG: BOND  i= %i  %i - %i (%i %i)\n",  bidx++, iat+1, *bat+1, iat*3, *bat*3);
          topIn.AddToBondArrays( BondType(iat, *bat, -1) );
        }
        //else
        //  mprintf("DEBUG: X    i= %i  %i - %i (%i %i)\n",   bidx++, iat+1, *bat+1, iat*3, *bat*3);
      }
    }
  }

  // ANGLES TODO combine above
  int aidx = 0;
  for (int ires = 0; ires < topIn.Nres(); ires++)
  {
    Residue const& res = topIn.Res(ires);
    for (int iat1 = res.LastAtom()-1; iat1 >= res.FirstAtom(); iat1--)
    {
      Atom const& At1 = topIn[iat1];
      for (int bidx1 = 0; bidx1 < At1.Nbonds(); bidx1++) {
        int iat2 = At1.Bond(bidx1);
        Atom const& At2 = topIn[iat2];
        for (int bidx2 = 0; bidx2 < At2.Nbonds(); bidx2++) {
          int iat3 = At2.Bond(bidx2);
          if (iat1 < iat3) {
            mprintf("DEBUG: ANGLE  i= %i  %i - %i - %i (%i %i %i)\n", aidx++, iat1+1, iat2+1, iat3+1, iat1*3, iat2*3, iat3*3);
            topIn.AddAngle(iat1, iat2, iat3);
          }
        }
      }
    }
  }

  // TORSIONS TODO combine above
  int didx = 0;
  for (int ires = 0; ires < topIn.Nres(); ires++)
  {
    Residue const& res = topIn.Res(ires);
    for (int iat1 = res.LastAtom()-1; iat1 >= res.FirstAtom(); iat1--)
    {
      Atom const& At1 = topIn[iat1];
      for (int bidx1 = 0; bidx1 < At1.Nbonds(); bidx1++) {
        int iat2 = At1.Bond(bidx1);
        Atom const& At2 = topIn[iat2];
        for (int bidx2 = 0; bidx2 < At2.Nbonds(); bidx2++) {
          int iat3 = At2.Bond(bidx2);
          if (iat3 != iat1) {
            Atom const& At3 = topIn[iat3];
            for (int bidx3 = 0; bidx3 < At3.Nbonds(); bidx3++) {
              int iat4 = At3.Bond(bidx3);
              if (iat4 != iat2 && iat1 < iat4) {
                mprintf("DEBUG: DIHEDRAL  i= %i  %i - %i - %i - %i (%i %i %i %i)\n", didx++, iat1+1, iat2+1, iat3+1, iat4+1, iat1*3, iat2*3, iat3*3, iat4*3);
                topIn.AddDihedral( iat1, iat2, iat3, iat4 );
              }
            }
          }
        }
      }
    }
  }

  return 0;
}
