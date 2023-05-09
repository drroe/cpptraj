#ifdef HAS_MEAD
#include "MeadCalc.h"
#include "MeadError.h"
#include "../Frame.h"
#include "../Topology.h"
#include "../CpptrajStdio.h"
// MEAD includes
#include "../../mead/MEADexcept.h"
#include "../../mead/AtomSet.h"

using namespace Cpptraj::Mead;

/** CONSTRUCTOR */
MeadCalc::MeadCalc() :
  atomset_(0),
  rmode_(GB)
{ }

/** DESTRUCTOR */
MeadCalc::~MeadCalc() {
  if (atomset_ != 0) delete atomset_;
}

/** Set MEAD Atom from Topology Atom. */
void MeadCalc::set_at_from_top(MEAD::Atom& at, Topology const& topIn, Frame const& frameIn, int aidx)
const
{
    Atom const& thisAtom = topIn[aidx];
    at.atname.assign( thisAtom.Name().Truncated() );
    int rnum = thisAtom.ResNum();
    Residue const& thisRes = topIn.Res(rnum);
    at.resname.assign( thisRes.Name().Truncated() );
    at.resnum = rnum;
    // NOTE: Not using chain ID here. Since we are using the residue index
    //       (and not original res#) this is enough to distinguish all residues.
    //if (thisRes.HasChainID())
    //  at.chainid.assign( 1, thisRes.ChainId() );
    const double* xyz = frameIn.XYZ(aidx);
    at.coord.x = xyz[0];
    at.coord.y = xyz[1];
    at.coord.z = xyz[2];
    at.charge = thisAtom.Charge();
    switch (rmode_) {
      case MeadCalc::GB : at.rad = thisAtom.GBRadius(); break;
      case MeadCalc::PARSE : at.rad = thisAtom.ParseRadius(); break;
      case MeadCalc::VDW   : at.rad = topIn.GetVDWradius(aidx); break;
    }
}

/** Setup an AtomSet from Frame and Topology. */
int MeadCalc::SetupAtoms(Topology const& topIn, Frame const& frameIn, Radii_Mode radiiMode)
{
  rmode_ = radiiMode;
  // Sanity checking
  if (topIn.Natom() != frameIn.Natom()) {
    mprinterr("Internal Error: MeadCalc::SetupAtoms(): Top '%s' has %i atoms, frame has %i atoms.\n",
              topIn.c_str(), topIn.Natom(), frameIn.Natom());
    return 1;
  }

  if (atomset_ != 0) delete atomset_;
  atomset_ = new AtomSet();

  bool has_radii = false;

  for (int aidx = 0; aidx < topIn.Natom(); aidx++)
  {
    MEAD::Atom at;
    set_at_from_top(at, topIn, frameIn, aidx);
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
}
#endif
