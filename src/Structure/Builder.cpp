#include "Builder.h"
#include "BuildAtom.h"
#include "Zmatrix.h"
#include "../CpptrajStdio.h"
#include "../Frame.h"
#include "../Topology.h"
#include <algorithm> // std::copy

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Builder::Builder() :
  debug_(0)
{}

/** Combine two units. Fragment 1 will be merged into Fragment 0 and bonded. */
int Builder::Combine(Topology&       frag0Top, Frame&       frag0frm,
                     Topology const& frag1Top, Frame const& frag1frm,
                     int bondAt0, int bondAt1)
const
{
  mprintf("DEBUG: Calling COMBINE between %s and %s\n",
          frag0Top.AtomMaskName(bondAt0).c_str(),
          frag1Top.AtomMaskName(bondAt1).c_str());
  int natom0 = frag0Top.Natom();
  int newNatom = natom0 + frag1Top.Natom();

  // Determine which "direction" we will be combining the fragments.
  // Make atA belong to the smaller fragment. atB fragment will be "known".
  // Ensure atB index is what it will be after fragments are combined.
  Barray posKnown( newNatom, false );
  int atA, atB;
  Zmatrix zmatrixA;
  if (frag0Top.HeavyAtomCount() < frag1Top.HeavyAtomCount()) {
    // Fragment 1 is larger
    atA = bondAt0;
    atB = bondAt1 + natom0;
    for (int at = frag0Top.Natom(); at != newNatom; at++)
      posKnown[at] = true;
    zmatrixA.SetFromFrameAndConnect( frag0frm, frag0Top );
  } else {
    // Fragment 0 is larger or equal
    atA = bondAt1 + natom0;
    atB = bondAt0;
    for (int at = 0; at != natom0; at++)
     posKnown[at] = true;
    zmatrixA.SetFromFrameAndConnect( frag1frm, frag1Top );
    zmatrixA.OffsetIcIndices( natom0 );
  }

  // Combine fragment1 into fragment 0 topology
  Topology& combinedTop = frag0Top;
  combinedTop.AppendTop( frag1Top );
  // Combined fragment1 into fragment 0 coords.
  // Need to save the original coords in frame0 since SetupFrameV does not preserve.
  double* tmpcrd0 = new double[natom0*3];
  std::copy( frag0frm.xAddress(), frag0frm.xAddress()+frag0frm.size(), tmpcrd0 );
  frag0frm.SetupFrameV( combinedTop.Atoms(), CoordinateInfo(frag0frm.BoxCrd(), false, false, false));
  std::copy( tmpcrd0, tmpcrd0+natom0*3, frag0frm.xAddress() );
  std::copy( frag1frm.xAddress(), frag1frm.xAddress()+frag1frm.size(), frag0frm.xAddress()+natom0*3 );
  Frame& CombinedFrame = frag0frm;
  delete[] tmpcrd0;

  int chiralityDebug;
  if (debug_ < 1)
    chiralityDebug = 0;
  else
    chiralityDebug = debug_ - 1;
  // Get the chirality around each atom before the bond is added.
  BuildAtom AtomA;
  if (combinedTop[atA].Nbonds() > 2) {
    if (AtomA.DetermineChirality(atA, combinedTop, CombinedFrame, chiralityDebug)) return 1;
  }
  if (debug_ > 0)
    mprintf("DEBUG:\tAtom %4s chirality before bonding is %6s\n", combinedTop.AtomMaskName(atA).c_str(), chiralStr(AtomA.Chirality()));
  BuildAtom AtomB;
  if (combinedTop[atB].Nbonds() > 2) {
    if (AtomB.DetermineChirality(atB, combinedTop, CombinedFrame, chiralityDebug)) return 1;
  }
  if (debug_ > 0)
    mprintf("DEBUG:\tAtom %4s chirality before bonding is %6s\n", combinedTop.AtomMaskName(atB).c_str(), chiralStr(AtomB.Chirality()));

  // Create the bond
  if (debug_ > 0)
    mprintf("DEBUG: Bonding atom %s to %s\n", combinedTop.AtomMaskName(atA).c_str(), combinedTop.AtomMaskName(atB).c_str());
  combinedTop.AddBond( atA, atB ); // TODO pseudo-parameter?
  // // Regenerate the molecule info FIXME should Topology just do this?
  if (combinedTop.DetermineMolecules()) return 1;

  // Determine new priorities around atoms that were just bonded
  if (combinedTop[atA].Nbonds() > 2) {
    if (AtomA.SetPriority(atA, combinedTop, CombinedFrame, chiralityDebug)) return 1;
  }
  if (combinedTop[atB].Nbonds() > 2) {
    if (AtomB.SetPriority(atB, combinedTop, CombinedFrame, chiralityDebug)) return 1;
  }

  // Generate Zmatrix only for ICs involving bonded atoms
  Zmatrix bondZmatrix;

  // Note which atoms already have an IC in zmatrixA
  std::vector<bool> hasIC(combinedTop.Natom(), false);
  for (Zmatrix::const_iterator it = zmatrixA.begin(); it != zmatrixA.end(); ++it)
    hasIC[it->AtI()] = true;

  //bondZmatrix.SetDebug( debug_ );
  bondZmatrix.SetDebug( 1 ); // FIXME
  if (bondZmatrix.SetupICsAroundBond(atA, atB, CombinedFrame, combinedTop, posKnown, hasIC, AtomA, AtomB)) {
    mprinterr("Error: Zmatrix setup for ICs around %s and %s failed.\n",
              combinedTop.AtomMaskName(atA).c_str(),
              combinedTop.AtomMaskName(atB).c_str());
    return 1;
  }
  mprintf("DEBUG: CONVERT WITH bondZmatrix.\n");
  //if (debug_ > 0)
    bondZmatrix.print(&combinedTop);
  if (bondZmatrix.SetToFrame( CombinedFrame, posKnown )) {
    mprinterr("Error: Conversion from bondZmatrix to Cartesian coords failed.\n");
    return 1;
  }
  mprintf("DEBUG: CONVERT WITH zmatrixA.\n");
  zmatrixA.print( &combinedTop );
  if (zmatrixA.SetToFrame( CombinedFrame, posKnown )) {
    mprinterr("Error: Conversion from zmatrixA to Cartesian coords failed.\n");
    return 1;
  }

  return 0;
}

/// \return The total number of atoms whose position is known for the specified residue
static inline int known_count(int ires, Topology const& topIn, std::vector<bool> const& hasPosition)
{
  int count = 0;
  for (int at = topIn.Res(ires).FirstAtom(); at != topIn.Res(ires).LastAtom(); at++)
    if (hasPosition[at])
      count++;
  return count;
}

/** Model the coordinates around a bond */
int Builder::ModelCoordsAroundBond(Frame& frameIn, Topology const& topIn, int bondAt0, int bondAt1,
                                   Zmatrix const* zmatrix0, Zmatrix const* zmatrix1, Barray& hasPosition)
const
{
  mprintf("DEBUG: Model coords around bond %s - %s\n", topIn.AtomMaskName(bondAt0).c_str(), topIn.AtomMaskName(bondAt1).c_str());
  int res0 = topIn[bondAt0].ResNum();
  int res1 = topIn[bondAt1].ResNum();
  if (res0 == res1) {
    mprinterr("Internal Error: ModelCoordsAroundBond(): Atoms are in the same residue.\n");
    return 1;
  }
  // Determine which "direction" we will be combining the fragments.
  // Make atA belong to the less-known fragment. atB fragment will be "known".
  int known0 = known_count(res0, topIn, hasPosition);
  int known1 = known_count(res1, topIn, hasPosition);
  int atA, atB;
  Zmatrix const* zA;
  if (known0 < known1) {
    // Fragment 1 is better-known 
    atA = bondAt0;
    atB = bondAt1;
    zA  = zmatrix0;
    //zB  = zmatrix1;
  } else {
    // Fragment 0 is better or equally known 
    atA = bondAt1;
    atB = bondAt0;
    zA  = zmatrix1;
    //zB  = zmatrix0;
  }

  mprintf("DEBUG: More well-known atom: %s\n", topIn.AtomMaskName(atB).c_str());
  mprintf("DEBUG: Less well-known atom: %s  has zmatrix= %i\n", topIn.AtomMaskName(atA).c_str(), (int)(zA != 0));

  int chiralityDebug;
  if (debug_ < 1)
    chiralityDebug = 0;
  else
    chiralityDebug = debug_ - 1;

  // Get the chirality around each atom before the bond is added.
  // Determine priorities
  BuildAtom AtomA;
  if (topIn[atA].Nbonds() > 2) {
    if (AtomA.DetermineChirality(atA, topIn, frameIn, chiralityDebug)) return 1;
  }
  if (debug_ > 0)
    mprintf("DEBUG:\tAtom %4s chirality %6s\n", topIn.AtomMaskName(atA).c_str(), chiralStr(AtomA.Chirality()));
  BuildAtom AtomB;
  if (topIn[atB].Nbonds() > 2) {
    if (AtomB.DetermineChirality(atB, topIn, frameIn, chiralityDebug)) return 1;
  }
  if (debug_ > 0)
    mprintf("DEBUG:\tAtom %4s chirality %6s\n", topIn.AtomMaskName(atB).c_str(), chiralStr(AtomB.Chirality()));

  // Generate Zmatrix only for ICs involving bonded atoms
  Zmatrix bondZmatrix;

  // Note which atoms already have an IC in zmatrix A
  std::vector<bool> hasIC(topIn.Natom(), false);
  if (zA != 0) {
    for (Zmatrix::const_iterator it = zA->begin(); it != zA->end(); ++it)
      hasIC[it->AtI()] = true;
  }

  bondZmatrix.SetDebug( debug_ );
  if (bondZmatrix.SetupICsAroundBond(atA, atB, frameIn, topIn, hasPosition, hasIC, AtomA, AtomB)) {
    mprinterr("Error: Zmatrix setup for ICs around %s and %s failed.\n",
              topIn.AtomMaskName(atA).c_str(),
              topIn.AtomMaskName(atB).c_str());
    return 1;
  }
  if (debug_ > 0)
    bondZmatrix.print(&topIn);
  if (bondZmatrix.SetToFrame( frameIn, hasPosition )) {
    mprinterr("Error: Conversion from bondZmatrix to Cartesian coords failed.\n");
    return 1;
  }

  return 0;
}

