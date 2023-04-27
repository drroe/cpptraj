#include "MeadCalc.h"
#include "MeadError.h"
#include "../Frame.h"
#include "../Topology.h"
#include "../CpptrajStdio.h"
// MEAD includes
#include "../../mead/MEADexcept.h"
#include "../../mead/AtomSet.h"
/*
#inc lude "../../mead/ChargeDist.h"
#inc lude "../../mead/AtomChargeSet.h"
#inc lude "../../mead/DielectricEnvironment.h"
#inc lude "../../mead/DielByAtoms.h"
#inc lude "../../mead/ElectrolyteEnvironment.h"
#inc lude "../../mead/ElectrolyteByAtoms.h"
#inc lude "../../mead/FinDiffElstatPot.h"
#inc lude "../../mead/Potat.h"
*/
using namespace Cpptraj::Mead;

/** CONSTRUCTOR */
MeadCalc::MeadCalc() :
  atomset_(0),
  rmode_(GB),
  t_total_("Mead Total")
{
  t_total_.AddSubTimer(Timer("SetupAtoms"));  // 0
/*  t_total_.AddSubTimer(Timer("MultiFlex "));  // 1
  t_total_[1].AddSubTimer(Timer("Setup sites")); // 1,0
  t_total_[1].AddSubTimer(Timer("MAC1       ")); // 1,1
  t_total_[1].AddSubTimer(Timer("MAC2       ")); // 1,2
  t_total_[1].AddSubTimer(Timer("MOD1       ")); // 1,3
  t_total_[1].AddSubTimer(Timer("MOD2       ")); // 1,4*/
}

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
    at.resnum = thisRes.OriginalResNum();
    if (thisRes.HasChainID())
      at.chainid.assign( 1, thisRes.ChainId() );
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
  t_total_[0].Start();
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
  t_total_[0].Stop();
  
  return 0;
}

/** Print debug info. */
void MeadCalc::Print() const {
  //std::cout << *fdm_;
  //std::cout << *mgm_;
}

/** Set verbosity of underlying MEAD library. */
void MeadCalc::MeadVerbosity(int i) const {
  blab1pt = &cnull;
  blab2pt = &cnull;
  blab3pt = &cnull;
  if (i >= 1)
    blab1pt = &std::cout;
  if (i >= 2)
    blab2pt = &std::cout;
  if (i >= 3)
    blab3pt = &std::cout;
  if (i != 0)
    mprintf("Info: MEAD verbosity set to %i\n", i);
}

