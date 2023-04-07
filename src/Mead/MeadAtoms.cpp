#include "MeadAtoms.h"
#include "../Topology.h"
#include "../Frame.h"
#include "../../mead/AtomSet.h"

using namespace Cpptraj::Mead;
/** CONSTRUCTOR */
MeadAtoms::MeadAtoms() :
  atomset_(0),
  rmode_(GB)
{}

/** DESTRUCTOR */
MeadAtoms::~MeadAtoms() {
  if (atomset_ != 0) delete atomset_;
}

/** Setup an AtomSet from Frame and Topology. */
/*int MeadAtoms::SetupAtoms(Topology const& topIn, Frame const& frameIn, Radii_Mode radiiMode)
{
  // Sanity checking
  if (topIn.Natom() != frameIn.Natom()) {
    mprinterr("Internal Error: MeadInterface::SetupAtoms(): Top '%s' has %i atoms, frame has %i atoms.\n",
              topIn.c_str(), topIn.Natom(), frameIn.Natom());
    return 1;
  }

  if (atomset_ != 0) delete atomset_;
  atomset_ = new AtomSet();

  bool has_radii = false;

  for (int aidx = 0; aidx < topIn.Natom(); aidx++)
  {
    MEAD::Atom at;

    set_at_from_top(at, topIn, frameIn, aidx, radiiMode);
    if (at.rad > 0)
      has_radii = true;
    try {
      atomset_->insert( at );
    }
    catch (MEADexcept& e) {
      return ERR("SetupAtoms()", e);
    }
  }
  if (!has_radii) {
    mprinterr("Error: No radii set for topology '%s'\n", topIn.c_str());
    return 1;
  }
  
  return 0;
}*/

