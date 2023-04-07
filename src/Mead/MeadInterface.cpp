#include "MeadInterface.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"
#include "../Frame.h"
#include "../../mead/MEADexcept.h"

/** Print MEAD error message. */
int Cpptraj::Mead::ERR(const char* fxn, MEADexcept& e) {
  mprinterr("Error: MEAD error in '%s': '%s' '%s' '%s'\n", fxn,
              e.get_error1().c_str(),
              e.get_error2().c_str(),
              e.get_error3().c_str());
  return 1;
}

/** Set MEAD Atom from Topology Atom. */
/*void Cpptraj::Mead::set_at_from_top(MEAD::Atom& at, Topology const& topIn, Frame const& frameIn, int aidx, Radii_Mode radiiMode)
{
    Atom const& thisAtom = topIn[aidx];
    at.atname.assign( thisAtom.Name().Truncated() );
    int rnum = thisAtom.ResNum();
    Residue const& thisRes = topIn.Res(rnum);
    at.resname.assign( thisRes.Name().Truncated() );
    at.resnum = thisRes.OriginalResNum();
    if (thisRes.HasChainID())
      at.chainid.assign( 1, thisRes.ChainId() );
    const double* xyz = frameIn.XYZ(aidx);
    at.coord.x = xyz[0];
    at.coord.y = xyz[1];
    at.coord.z = xyz[2];
    at.charge = thisAtom.Charge();
    switch (radiiMode) {
      case MeadInterface::GB : at.rad = thisAtom.GBRadius(); break;
      case MeadInterface::PARSE : at.rad = thisAtom.ParseRadius(); break;
      case MeadInterface::VDW   : at.rad = topIn.GetVDWradius(aidx); break;
    }
}*/

