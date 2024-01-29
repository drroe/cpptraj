#include "Builder.h"
#include "BuildAtom.h"
#include "Model.h"
#include "Zmatrix.h"
#include "../Constants.h"
#include "../CpptrajStdio.h"
#include "../DistRoutines.h"
#include "../Frame.h"
#include "../GuessAtomHybridization.h"
#include "../Topology.h"
#include "../TorsionRoutines.h"
#include <algorithm> // std::copy

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Builder::Builder() :
  debug_(0),
  params_(0)
{}

/** Set optional parameter set. */
void Cpptraj::Structure::Builder::SetParameters(ParameterSet const* paramsIn) {
  if (paramsIn == 0) {
    mprinterr("Internal Error: Builder::SetParmaters called with null set.\n");
    return;
  }
  params_ = paramsIn;
}

// -----------------------------------------------------------------------------
/** Assign reasonable value for bond distance. */
int Cpptraj::Structure::Builder::AssignLength(double& dist, int ai, int aj, Topology const& topIn, Frame const& frameIn, std::vector<bool> const& atomPositionKnown)
const
{
  if (atomPositionKnown[ai] && atomPositionKnown[aj])
    dist = sqrt( DIST2_NoImage( frameIn.XYZ(ai), frameIn.XYZ(aj) ) );
  else {
    // One or both positions unknown. Use estimated bond length or parameters.
    bool foundParam = false;
    if (params_ != 0 && topIn[ai].HasType() && topIn[aj].HasType()) {
      TypeNameHolder btypes(2);
      btypes.AddName( topIn[ai].Type() );
      btypes.AddName( topIn[aj].Type() );
      ParmHolder<BondParmType>::const_iterator it = params_->BP().GetParam( btypes );
      if (it != params_->BP().end()) {
        dist = it->second.Req();
        foundParam = true;
        mprintf("DEBUG: Found bond parameter for %s (%s) - %s (%s): req=%g rk=%g\n",
                topIn.AtomMaskName(ai).c_str(), *(topIn[ai].Type()),
                topIn.AtomMaskName(aj).c_str(), *(topIn[aj].Type()),
                it->second.Req(), it->second.Rk());
      }
    }
    if (!foundParam)
      dist = Atom::GetBondLength( topIn[ai].Element(), topIn[aj].Element() );
  }
  return 0;
}

/** Attempt to assign a reasonable value for theta internal coordinate for
  * atom i given that atoms j and k have known positions.
  */
int Cpptraj::Structure::Builder::AssignTheta(double& theta, int ai, int aj, int ak, Topology const& topIn, Frame const& frameIn, std::vector<bool> const& atomPositionKnown)
const
{
  // Figure out hybridization and chirality of atom j.
  if (debug_ > 0)
    mprintf("DEBUG: AssignTheta for atom j : %s\n", topIn.AtomMaskName(aj).c_str());
  if (atomPositionKnown[ai] && atomPositionKnown[aj] && atomPositionKnown[ak])
  {
    theta = CalcAngle(frameIn.XYZ(ai), frameIn.XYZ(aj), frameIn.XYZ(ak));
    return 0;
  }

  Atom const& AJ = topIn[aj];
  if (debug_ > 0) {
    mprintf("DEBUG:\t\tI %s Nbonds: %i\n", topIn[ai].ElementName(), topIn[ai].Nbonds());
    mprintf("DEBUG:\t\tJ %s Nbonds: %i\n", AJ.ElementName(), AJ.Nbonds());
    mprintf("DEBUG:\t\tK %s Nbonds: %i\n", topIn[ak].ElementName(), topIn[ak].Nbonds());
  }
  // Sanity check
  if (AJ.Nbonds() < 2) {
    mprinterr("Internal Error: AssignTheta() called for atom J %s with fewer than 2 bonds.\n", topIn.AtomMaskName(aj).c_str());
    return 1;
  }
  AtomType::HybridizationType hybrid = AtomType::UNKNOWN_HYBRIDIZATION;
  if (params_ != 0) {
    ParmHolder<AtomType>::const_iterator it = params_->AT().GetParam( TypeNameHolder(AJ.Type()) );
    if (it != params_->AT().end())
      hybrid = it->second.Hybridization();
  }
  if (hybrid == AtomType::UNKNOWN_HYBRIDIZATION)
    hybrid = GuessAtomHybridization(AJ, topIn.Atoms());
/*
  HybridizationType hybrid = UNKNOWN_HYBRIDIZATION;
  // Handle specific elements
  switch (AJ.Element()) {
    case Atom::CARBON :
      switch (AJ.Nbonds()) {
        case 2 : hybrid = SP; break;
        case 3 : hybrid = SP2; break;
        case 4 : hybrid = SP3; break;
      }
      break;
    case Atom::NITROGEN :
      switch (AJ.Nbonds()) {
        case 2 : hybrid = SP2; break;
        case 3 :
          // Check for potential SP2. If only 1 of the bonded atoms is
          // hydrogen, assume SP2. TODO actually check for aromaticity.
          int n_hydrogens = 0;
          for (Atom::bond_iterator bat = AJ.bondbegin(); bat != AJ.bondend(); ++bat)
            if (topIn[*bat].Element() == Atom::HYDROGEN)
              n_hydrogens++;
          if (n_hydrogens == 1)
            hybrid = SP2;
          else
            hybrid = SP3;
          break;
      }
      break;
    case Atom::OXYGEN :
      switch (AJ.Nbonds()) {
        case 2 : hybrid = SP3; break;
      }
      break;
    case Atom::SULFUR :
      switch (AJ.Nbonds()) {
        case 2 : hybrid = SP3; break;
      }
      break;
    default: hybrid = UNKNOWN_HYBRIDIZATION; break;
  }*/
  // Fill in what values we can for known atoms
/*  std::vector<double> knownTheta( AJ.Nbonds() );
  int knownIdx = -1;
  for (int idx = 0; idx < AJ.Nbonds(); idx++) {
    int atnum = AJ.Bond(idx);
    if (atnum != ak && atomPositionKnown[atnum]) {
      knownTheta[idx] = CalcAngle(frameIn.XYZ(atnum),
                                  frameIn.XYZ(aj),
                                  frameIn.XYZ(ak));
      mprintf("DEBUG:\t\tKnown theta for %s = %g\n", topIn.AtomMaskName(atnum).c_str(), knownTheta[idx]*Constants::RADDEG);
      if (knownIdx == -1) knownIdx = idx; // FIXME handle more than 1 known
    }
  }
  if (knownIdx == -1) {*/
    //mprintf("DEBUG:\t\tNo known theta.\n");
  if (hybrid == AtomType::UNKNOWN_HYBRIDIZATION) {
    // Assign a theta based on number of bonds 
    switch (AJ.Nbonds()) {
      case 4 : hybrid = AtomType::SP3; break;
      case 3 : hybrid = AtomType::SP2; break;
      case 2 : hybrid = AtomType::SP; break;
      default : mprinterr("Internal Error: AssignTheta(): Unhandled # bonds for %s (%i)\n", topIn.AtomMaskName(aj).c_str(), AJ.Nbonds()); return 1;
    }
  }
  // Assign a theta based on hybridization
  switch (hybrid) {
    case AtomType::SP3 : theta = 109.5 * Constants::DEGRAD; break;
    case AtomType::SP2 : theta = 120.0 * Constants::DEGRAD; break;
    case AtomType::SP  : theta = 180.0 * Constants::DEGRAD; break;
    default : mprinterr("Internal Error: AssignTheta(): Unhandled hybridization for %s (%i)\n", topIn.AtomMaskName(aj).c_str(), AJ.Nbonds()); return 1;
  }/*
  } else {
    theta = knownTheta[knownIdx]; // TODO just use above guess via hybrid?
  }*/

  return 0;
}

// -----------------------------------------------------------------------------
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
  if (SetupICsAroundBond(bondZmatrix, atA, atB, CombinedFrame, combinedTop, posKnown, hasIC, AtomA, AtomB)) {
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
  if (SetupICsAroundBond(bondZmatrix, atA, atB, frameIn, topIn, hasPosition, hasIC, AtomA, AtomB)) {
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

/** Given two bonded atoms A and B, where B has a depth of at least 2
  * (i.e., it is possible to have B be atom J where we can form J-K-L),
  * set up a complete set of internal coordinates involving A and B in the
  * direction of atom A. This means all internal coordinates with A and B
  * as I and J (should be only 1), as J and K, and as K and L.
  */
int Builder::SetupICsAroundBond(Zmatrix& zmatrix,
                                int atA, int atB, Frame const& frameIn, Topology const& topIn,
                                std::vector<bool> const& atomPositionKnown,
                                std::vector<bool> const& hasICin,
                                BuildAtom const& AtomA, BuildAtom const& AtomB)
const
{
  if (debug_ > 0)
    mprintf("DEBUG: SetupICsAroundBond: atA= %s (%i)  atB= %s (%i) total # atoms %i\n",
            topIn.AtomMaskName(atA).c_str(), atA+1,
            topIn.AtomMaskName(atB).c_str(), atB+1,
            topIn.Natom());
  zmatrix.clear();

  //Barray hasIC( topIn.Natom(), false );
  Barray hasIC = hasICin;
  unsigned int nHasIC = 0;
  for (std::vector<bool>::const_iterator it = hasIC.begin(); it != hasIC.end(); ++it) {
    if (*it) {
      nHasIC++;
      mprintf("DEBUG:\tAtom %s already has an IC.\n", topIn.AtomMaskName(it-hasIC.begin()).c_str());
    }
  }
  // Mark known atoms as already having IC
  for (std::vector<bool>::const_iterator it = atomPositionKnown.begin();
                                         it != atomPositionKnown.end(); ++it)
  {
    //mprintf("DEBUG: MARKING KNOWN ATOMS. %li\n", it - atomPositionKnown.begin());
    if (*it) { //MARK( it - atomPositionKnown.begin(), hasIC, nHasIC );
      hasIC[it - atomPositionKnown.begin()] = true;
      nHasIC++;
    }
  }

  // First, make sure atom B as a bond depth of at least 2.
  // Choose K and L atoms given atA is I and atB is J.
  Atom const& AJ = topIn[atB];
  typedef std::pair<int,int> Apair;
  std::vector<Apair> KLpairs;
  for (Atom::bond_iterator kat = AJ.bondbegin(); kat != AJ.bondend(); ++kat)
  {
    if (*kat != atA) {
      //mprintf("DEBUG: kat= %s\n", topIn.AtomMaskName(*kat).c_str());
      Atom const& AK = topIn[*kat];
      for (Atom::bond_iterator lat = AK.bondbegin(); lat != AK.bondend(); ++lat)
      {
        if (*lat != atB && *lat != atA) {
          //mprintf("DEBUG: lat= %s\n", topIn.AtomMaskName(*lat).c_str());
          KLpairs.push_back( Apair(*kat, *lat) );
        }
      }
    }
  }
  if (debug_ > 0) {
    for (std::vector<Apair>::const_iterator it = KLpairs.begin();
                                            it != KLpairs.end(); ++it)
      mprintf("DEBUG:\t\tKL pair %s - %s\n", topIn.AtomMaskName(it->first).c_str(),
              topIn.AtomMaskName(it->second).c_str());
  }
  if (KLpairs.empty()) {
    mprinterr("Error: SetFromFrameAroundBond(): Could not find an atom pair bonded to atom %s\n",
              topIn.AtomMaskName(atB).c_str());
    return 1;
  }
  // TODO be smarter about how K and L are selected?
  double maxMass = topIn[KLpairs[0].first].Mass() + topIn[KLpairs[0].second].Mass();
  unsigned int maxIdx = 0;
  for (unsigned int idx = 1; idx < KLpairs.size(); idx++) {
    double sumMass = topIn[KLpairs[idx].first].Mass() + topIn[KLpairs[idx].second].Mass();
    if (sumMass > maxMass) {
      maxMass = sumMass;
      maxIdx = idx;
    }
  }
  int atk0 = KLpairs[maxIdx].first;
  int atl0 = KLpairs[maxIdx].second;
  int modelDebug = 0;
  if (debug_ > 0) {
    mprintf("DEBUG: Chosen KL pair: %s - %s\n",topIn.AtomMaskName(atk0).c_str(),
              topIn.AtomMaskName(atl0).c_str());
    modelDebug = debug_ - 1;
  }
  Cpptraj::Structure::Model model;
  //model.SetDebug( modelDebug ); //FIXME
  model.SetDebug( 1 );
  // FIXME
  mprintf("DEBUG: ------------------------------------------------\n");
  //std::vector<InternalCoords> tmpic;
  // I J: Set up ICs for X atB K L
  if (model.AssignICsAroundBond(zmatrix, atB, atk0, atl0, topIn, frameIn, atomPositionKnown, AtomB)) {
    mprinterr("Error: AssignICsAroundBond (I J) failed.\n");
    return 1;
  }
  // J K: Set up ICs for X atA atB K
  if (model.AssignICsAroundBond(zmatrix, atA, atB, atk0, topIn, frameIn, atomPositionKnown, AtomA)) {
    mprinterr("Error: AssignICsAroundBond (J K) failed.\n");
    return 1;
  }
  // K L: Set up ICs for X iat atA atB
  Atom const& AJ1 = topIn[atA];
  for (Atom::bond_iterator iat = AJ1.bondbegin(); iat != AJ1.bondend(); ++iat)
  {
    if (*iat != atB) {
      // Only do this if one or more of the atoms bonded to iat does not
      // already have an IC.
      unsigned int needsKLic = 0;
      for (Atom::bond_iterator bat = topIn[*iat].bondbegin(); bat != topIn[*iat].bondend(); ++bat) {
        if ( *bat != atA && !hasIC[*bat] ) {
          mprintf("DEBUG:\tAtom %s around %s has no IC.\n", topIn.AtomMaskName(*bat).c_str(), topIn.AtomMaskName(*iat).c_str());
          needsKLic++;
        }
      }
      if (needsKLic > 0) {
        BuildAtom AtomC;
        if (topIn[*iat].Nbonds() > 2) {
          if (AtomC.DetermineChirality(*iat, topIn, frameIn, modelDebug)) return 1;
        }
        mprintf("DEBUG: K L IC needed for %s.\n", topIn.AtomMaskName(*iat).c_str());
        if (model.AssignICsAroundBond(zmatrix, *iat, atA, atB, topIn, frameIn, atomPositionKnown, AtomC)) {
          mprinterr("Error: AssignICsAroundBond (K L) failed.\n");
          return 1;
        }
      } else {
        mprintf("DEBUG: K L IC not needed for %s.\n", topIn.AtomMaskName(*iat).c_str());
      }
    }
  }
  // Mark/Print ICs
  for (Zmatrix::const_iterator it = zmatrix.begin(); it != zmatrix.end(); ++it) {
    it->printIC( topIn ); // DEBUG
    //MARK( it->AtI(), hasIC, nHasIC );
  }
  mprintf("DEBUG: END AssignICsAroundBond ------------------------\n");
  // FIXME
  // ---- I J: Set dist, theta, phi for atA atB K L internal coord ---
/*  if (debug_ > 0)
    mprintf("DEBUG: IC (i j) %i - %i - %i - %i\n", atA+1, atB+1, atk0+1, atl0+1);
  double newDist = 0;
  if (model.AssignLength(newDist, atA, atB, topIn, frameIn, atomPositionKnown)) {
    mprinterr("Error: length (i j) assignment failed.\n");
    return 1;
  }
  if (debug_ > 0) mprintf("DEBUG:\t\tnewDist= %g\n", newDist);
  double newTheta = 0;
  if (model.AssignTheta(newTheta, atA, atB, atk0, topIn, frameIn, atomPositionKnown)) {
    mprinterr("Error: theta (i j) assignment failed.\n");
    return 1;
  }
  if (debug_ > 0) mprintf("DEBUG:\t\tnewTheta = %g\n", newTheta*Constants::RADDEG);
  double newPhi = 0;
  if (model.AssignPhi(newPhi, atA, atB, atk0, atl0, topIn, frameIn,
                                           atomPositionKnown, AtomB))
  {
    mprinterr("Error: phi (i j) assignment failed.\n");
    return 1;
  }
  if (debug_ > 0) mprintf("DEBUG:\t\tnewPhi = %g\n", newPhi*Constants::RADDEG);
  IC_.push_back(InternalCoords( atA, atB, atk0, atl0, newDist, newTheta*Constants::RADDEG, newPhi*Constants::RADDEG ));
  if (debug_ > 0) {
    mprintf("DEBUG: MODEL I J IC: ");
    IC_.back().printIC(topIn);
  }
  MARK( atA, hasIC, nHasIC );*/
/*
  // ----- J K: Set up ICs for X atA atB K ---------------------------
  Atom const& AJ1 = topIn[atA];
  int ati = -1;
  //Atom const& AK1 = topIn[atB];
  //Atom const& AL1 = topIn[atk0];
  for (Atom::bond_iterator iat = AJ1.bondbegin(); iat != AJ1.bondend(); ++iat)
  {
    if (*iat != atB) {
      if (ati == -1) ati = *iat;
/      // Set bond dist
      if (model.AssignLength(newDist, *iat, atA, topIn, frameIn, atomPositionKnown)) {
        mprinterr("Error: length (j k) assignment failed.\n");
        return 1;
      }
      // Set theta for I atA atB
      newTheta = 0;
      if (model.AssignTheta(newTheta, *iat, atA, atB, topIn, frameIn, atomPositionKnown)) {
        mprinterr("Error: theta (j k) assignment failed.\n");
        return 1;
      }
      if (debug_ > 0)
        mprintf("DEBUG:\t\tnewTheta = %g\n", newTheta*Constants::RADDEG);
      // Set phi for I atA atB K
      newPhi = 0;
      if (model.AssignPhi(newPhi, *iat, atA, atB, atk0, topIn, frameIn,
                                               atomPositionKnown, AtomA))
      {
        mprinterr("Error: phi (j k) assignment failed.\n");
        return 1;
      }
      if (debug_ > 0)
        mprintf("DEBUG:\t\tnewPhi = %g\n", newPhi*Constants::RADDEG);
      IC_.push_back(InternalCoords( *iat, atA, atB, atk0, newDist, newTheta*Constants::RADDEG, newPhi*Constants::RADDEG ));
      if (debug_ > 0) {
        mprintf("DEBUG: MODEL J K IC: ");
        IC_.back().printIC(topIn);
      }
      MARK( *iat, hasIC, nHasIC );
      // ----- K L: Set up ICs for X iat atA atB ---------------------
      Atom const& AJ2 = topIn[*iat];
      for (Atom::bond_iterator i2at = AJ2.bondbegin(); i2at != AJ2.bondend(); ++i2at)
      {
        if (*i2at != atA && *i2at != atB && !hasIC[*i2at]) {
          mprintf("DEBUG: MODEL K L IC: %s %s %s %s %i %i %i %i\n",
                  topIn.AtomMaskName(*i2at).c_str(),
                  topIn.AtomMaskName(*iat).c_str(),
                  topIn.AtomMaskName(atA).c_str(),
                  topIn.AtomMaskName(atB).c_str(),
                  *i2at+1, *iat+1, atA+1, atB+1);
          if (model.AssignLength(newDist, *i2at, *iat, topIn, frameIn, atomPositionKnown)) {
            mprinterr("Error: length (k l) assignment failed.\n");
            return 1;
          }
          mprintf("DEBUG: K L distance= %g\n", newDist);
          //newTheta = CalcAngle( frameIn.XYZ(*i2at), frameIn.XYZ(*iat), frameIn.XYZ(atA) );
          if (model.AssignTheta(newTheta, *iat, atA, atB, topIn, frameIn, atomPositionKnown)) {
            mprinterr("Error: theta (k l) assignment failed.\n");
            return 1;
          }
          mprintf("DEBUG: K L angle= %g\n", newTheta*Constants::RADDEG);
          // Set phi for X iat atA atB
          BuildAtom AtomC;
          if (topIn[*iat].Nbonds() > 2) {
            if (AtomC.DetermineChirality(*iat, topIn, frameIn, modelDebug)) return 1;
          }
          newPhi = 0;
          //model.SetDebug( 1 ); // FIXME
          if (model.AssignPhi(newPhi, *i2at, *iat, atA, atB, topIn, frameIn,
                                                   atomPositionKnown, AtomC))
          {
            mprinterr("Error: phi (k l) assignment failed.\n");
            return 1;
          }
          mprintf("DEBUG: K L Phi = %g\n", newPhi*Constants::RADDEG);
          IC_.push_back(InternalCoords( *i2at, *iat, atA, atB, newDist, newTheta*Constants::RADDEG, newPhi*Constants::RADDEG ));
          mprintf("DEBUG: MODEL K L IC: ");
          IC_.back().printIC(topIn);
          MARK( *i2at, hasIC, nHasIC );
          // Trace from atA *iat *i2at outwards
          //ToTrace.push_back(atA);
          //ToTrace.push_back(*iat);
          //ToTrace.push_back(*i2at);
          //if (traceMol(atA, *iat, *i2at, frameIn, topIn, topIn.Natom(), nHasIC, hasIC)) return 1;

        }
      } 
    }
  }
*/

/*
  // Handle remaining atoms.
  if (AJ1.Nbonds() > 1) {
    if (AJ1.Nbonds() == 2) {
      if (debug_ > 0) mprintf("DEBUG: 2 bonds to %s.\n", topIn.AtomMaskName(atA).c_str());
      if (traceMol(atB, atA, ati, frameIn, topIn, topIn.Natom(), nHasIC, hasIC))
        return 1;
    } else {
      // 3 or more bonds
      std::vector<int> const& priority = AtomA.Priority();
      int at0 = -1;
      int at1 = -1;
      std::vector<int> remainingAtoms;
      if (debug_ > 0) mprintf("DEBUG: %i bonds to %s\n", AJ1.Nbonds(), topIn.AtomMaskName(atA).c_str());
      for (std::vector<int>::const_iterator it = priority.begin(); it != priority.end(); ++it) {
        if (*it != atB) {
          if (debug_ > 0) mprintf("DEBUG:\t\t%s\n", topIn.AtomMaskName(*it).c_str());
          if (at0 == -1)
            at0 = *it;
          else if (at1 == -1)
            at1 = *it;
          else
            remainingAtoms.push_back( *it );
        }
      }
      // at0 atA at1
      if (traceMol(at1, atA, at0, frameIn, topIn, topIn.Natom(), nHasIC, hasIC))
        return 1;
      // at1 atA, at0
      if (traceMol(at0, atA, at1, frameIn, topIn, topIn.Natom(), nHasIC, hasIC))
        return 1;
      // Remaining atoms.
      for (std::vector<int>::const_iterator it = remainingAtoms.begin(); it != remainingAtoms.end(); ++it) {
        if (traceMol(at0, atA, *it, frameIn, topIn, topIn.Natom(), nHasIC, hasIC))
          return 1;
      }
    }
  }
*/
  return 0;
}

