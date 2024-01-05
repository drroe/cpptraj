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

/** From atom connectivity, generate a dihedral array in the same order as LEaP. */ // TODO use in GenerateBAT
DihedralArray Cpptraj::Structure::GenerateDihedralArray(std::vector<Residue> const& residues,
                                                        std::vector<Atom> const& atoms)
{
  DihedralArray out;
  // TORSIONS TODO combine above
  int didx = 0;
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
          if (iat3 != iat1) {
            Atom const& At3 = atoms[iat3];
            for (int bidx3 = 0; bidx3 < At3.Nbonds(); bidx3++) {
              int iat4 = At3.Bond(bidx3);
              if (iat4 != iat2 && iat1 < iat4) {
                mprintf("DEBUG: DIHEDRAL  i= %i  %i - %i - %i - %i (%i %i %i %i)\n", didx++, iat1+1, iat2+1, iat3+1, iat4+1, iat1*3, iat2*3, iat3*3, iat4*3);
                out.push_back( DihedralType( iat1, iat2, iat3, iat4, -1 ) );
              }
            }
          }
        }
      }
    }
  }
  return out;
}

/** Try to order an improper the same way that LEaP does.
  * LEaP has wild card names first, followed by atom types
  * in alphabetical order.
  */
static void order_improper_atoms(int* indices, std::vector<Atom> const& atoms)
{
  if (atoms[indices[0]].Type() > atoms[indices[1]].Type()) std::swap( indices[0], indices[1] );
  if (atoms[indices[1]].Type() > atoms[indices[2]].Type()) std::swap( indices[1], indices[2] );
  if (atoms[indices[0]].Type() > atoms[indices[1]].Type()) std::swap( indices[0], indices[1] );
  if (atoms[indices[1]].Type() > atoms[indices[2]].Type()) std::swap( indices[1], indices[2] );
}

// DEBUG
static inline void printName(Atom const& AJ) {
  mprintf(" :%i@%s", AJ.ResNum()+1, AJ.Name().Truncated().c_str());
}

/** From atom connectivity, generate an improper array in the same order as LEaP.
  * No attempt is made to determine if this is an sp2 center; that is done
  * during parameterization.
  */
DihedralArray Cpptraj::Structure::GenerateImproperArray(std::vector<Residue> const& residues,
                                                        std::vector<Atom> const& atoms)
{
  DihedralArray out;
  int iidx = 0;
  for (std::vector<Residue>::const_iterator res = residues.begin(); res != residues.end(); ++res)
  {
    for (int iat3 = res->LastAtom()-1; iat3 >= res->FirstAtom(); iat3--)
    {
      Atom const& AJ = atoms[iat3];
      if (AJ.Nbonds() >= 3) {
        for (int bidx0 = 0; bidx0 < AJ.Nbonds(); bidx0++) {
          for (int bidx1 = bidx0 + 1; bidx1 < AJ.Nbonds(); bidx1++) {
            for (int bidx2 = bidx1 + 1; bidx2 < AJ.Nbonds(); bidx2++) {
              int iat1 = AJ.BondIdxArray()[bidx0];
              int iat2 = AJ.BondIdxArray()[bidx1];
              int iat4 = AJ.BondIdxArray()[bidx2];
              mprintf("DEBUG: IMPROPER  i= %i  %i - %i - %i - %i (%i %i %i %i)\n", iidx++, iat1+1, iat2+1, iat3+1, iat4+1, iat1*3, iat2*3, iat3*3, iat4*3);
              int indices[3];
              indices[0] = iat1;
              indices[1] = iat2;
              indices[2] = iat4;
              order_improper_atoms(indices, atoms);
              out.push_back( DihedralType(indices[0], indices[1], iat3, indices[2], DihedralType::BOTH) );
              // DEBUG
              mprintf("DEBUG:\tOriginal order :");
              printName(atoms[iat1]);
              printName(atoms[iat2]);
              printName(atoms[iat3]);
              printName(atoms[iat4]);
              mprintf("\nDEBUG:\tLeap order     :");
              printName(atoms[indices[0]]);
              printName(atoms[indices[1]]);
              printName(atoms[iat3]);
              printName(atoms[indices[2]]);
              mprintf("\n");
            }
          }
        } // END outer loop over bond indices
      }
    } // END loop over residue atoms
  } // END loop over residues
  return out;
}
