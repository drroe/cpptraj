#include "GenerateImpropers.h"
#include "../CpptrajStdio.h"
#include "../Atom.h"
#include "../AtomType.h"
#include "../GuessAtomHybridization.h"
#include "../ParameterTypes.h"
#include "../ParameterHolders.h"
#include <algorithm> // std::swap

/** Try to order an improper the same way that LEaP does.
  * LEaP has wild card names first, followed by atom types
  * in alphabetical order. The third atom is always the central
  * atom.
  */
static int order_improper_atoms(int* leapdih, int centralAt, std::vector<Atom> const& topIn)
{
  Atom const& Ak = topIn[centralAt];
  std::vector<int> const& bonds = Ak.BondIdxArray();
  int indices[3];
  indices[0] = bonds[0];
  indices[1] = bonds[1];
  indices[2] = bonds[2];
  if (topIn[indices[0]].Type() > topIn[indices[1]].Type()) std::swap( indices[0], indices[1] );
  if (topIn[indices[1]].Type() > topIn[indices[2]].Type()) std::swap( indices[1], indices[2] );
  if (topIn[indices[0]].Type() > topIn[indices[1]].Type()) std::swap( indices[0], indices[1] );
  if (topIn[indices[1]].Type() > topIn[indices[2]].Type()) std::swap( indices[1], indices[2] );
  leapdih[0] = indices[0];
  leapdih[1] = indices[1];
  leapdih[2] = centralAt;
  leapdih[3] = indices[2];
  return 0;
}

// DEBUG
static inline void printName(Atom const& AJ) {
  mprintf(" :%i@%s", AJ.ResNum()+1, AJ.Name().Truncated().c_str());
}

/** Try to determine impropers for topology. */ // TODO option for charmm improper
DihedralArray Cpptraj::Structure::GenerateImproperArray(std::vector<Atom> const& atoms,
                                                        ParmHolder<AtomType> const& AT)
{
  DihedralArray out;
  long int atBegin = (long int)(atoms.size()) - 1;
  for (long int iat = atBegin; iat >= 0; iat--) {
    Atom const& AJ = atoms[iat];
    if (AJ.Nbonds() == 3) { // TODO only 3 atoms OK?
      AtomType::HybridizationType hybrid = AtomType::UNKNOWN_HYBRIDIZATION; 
      bool found;
      AtomType atype = AT.FindParam(TypeNameHolder(AJ.Type()), found);
      if (found)
        hybrid = atype.Hybridization();
      if (hybrid == AtomType::UNKNOWN_HYBRIDIZATION)
        hybrid = Cpptraj::GuessAtomHybridization( AJ, atoms );
      if (hybrid == AtomType::SP2) {
        mprintf("DEBUG: Potential improper center:");
        printName(AJ);
        mprintf("\n");
        int leapdih[4];
        order_improper_atoms(leapdih, iat, atoms);
        mprintf("DEBUG: Bond order");
        printName(AJ);
        mprintf(": %i %i %i %i", 
                AJ.BondIdxArray()[0]+1,
                AJ.BondIdxArray()[1]+1,
                iat+1,
                AJ.BondIdxArray()[2]+1);
        printName(atoms[AJ.BondIdxArray()[0]]);
        printName(atoms[AJ.BondIdxArray()[1]]);
        printName(AJ);
        printName(atoms[AJ.BondIdxArray()[2]]);
        mprintf("\n");
        mprintf("DEBUG: Leap order: %i %i %i %i\n", leapdih[0]+1, leapdih[1]+1, leapdih[2]+1, leapdih[3]+1);
        out.push_back( DihedralType(leapdih[0], leapdih[1], leapdih[2], leapdih[3], DihedralType::BOTH) );
      } else if (hybrid == AtomType::UNKNOWN_HYBRIDIZATION) {
        mprintf("Warning: When searching for impropers could not determine hybridization of");
        printName(AJ);
        mprintf("\n");
      }
    }
  }
  return out;
}
