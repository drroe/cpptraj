#include "GenerateAngles.h"
#include "../CpptrajStdio.h"
#include "../GuessAtomHybridization.h"
#include "../ParameterTypes.h"
#include "../Topology.h"
#include <algorithm> // std::sort

static inline void enumerateAngles(int at1, int at2, Topology& topIn) {
  Atom const& A2 = topIn[at2];
  if (A2.Nbonds() > 1) {
    for (Atom::bond_iterator bat = A2.bondbegin(); bat != A2.bondend(); ++bat)
    {
      if (*bat != at1 && at1 < *bat) {
        topIn.AddAngle(at1, at2, *bat);
        mprintf("%s - %s - %s\n",
                topIn.AtomMaskName(at1).c_str(),
                topIn.AtomMaskName(at2).c_str(),
                topIn.AtomMaskName(*bat).c_str());
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
            mprintf("%s - %s - %s - %s\n",
                    topIn.AtomMaskName(*bat1).c_str(),
                    topIn.AtomMaskName(at1).c_str(),
                    topIn.AtomMaskName(at2).c_str(),
                    topIn.AtomMaskName(*bat2).c_str());
          }
        }
      }
    }
  }
}

        

int Cpptraj::Structure::GenerateAngles(Topology& topIn) {
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

  // Angles
  AngleArray angles;
  AngleArray anglesH;
  for (BondArray::const_iterator it = allBonds.begin(); it != allBonds.end(); ++it)
  {
    if (!has_angles) {
      // Forward direction. A1-A2-X
      enumerateAngles( it->A1(), it->A2(), topIn );
      // Reverse direction. A2-A1-X
      enumerateAngles( it->A2(), it->A1(), topIn );
    }
    if (!has_dihedrals) {
      // Dihedrals
      enumerateDihedrals( it->A1(), it->A2(), topIn );
    }
  }
  return 0;
}

/** Try to determine impropers for topology. */
int Cpptraj::Structure::GenerateImpropers(Topology& topIn) {
  for (int iat = 0; iat != topIn.Natom(); iat++) {
    Atom const& AJ = topIn[iat];
    if (AJ.Nbonds() == 3) { // TODO only 3 atoms OK?
      // FIXME pass in atom types
      AtomType::HybridizationType hybrid = Cpptraj::GuessAtomHybridization( AJ, topIn );
      if (hybrid == AtomType::SP2) {
        mprintf("DEBUG: Potential improper center: %s\n", topIn.AtomMaskName(iat).c_str());
      }
    }
  }
  return 0;
}
