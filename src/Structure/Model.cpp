#include "Model.h"
#include "BuildAtom.h"
//#incl ude "InternalCoords.h"
//#incl ude "Chirality.h"
#include "../CpptrajStdio.h"
#include "../GuessAtomHybridization.h"
#include "../Topology.h"
#include "../TorsionRoutines.h"
#include "../DistRoutines.h"
#include "../Constants.h"
#include <cmath>

/** Assign reasonable value for bond distance. */
int Cpptraj::Structure::Model::AssignLength(double& dist, int ai, int aj, Topology const& topIn, Frame const& frameIn, std::vector<bool> const& atomPositionKnown)
const
{
  if (atomPositionKnown[ai] && atomPositionKnown[aj])
    dist = sqrt( DIST2_NoImage( frameIn.XYZ(ai), frameIn.XYZ(aj) ) );
  else
    // One or both positions unknown. Use bond length. TODO use parameters
    dist = Atom::GetBondLength( topIn[ai].Element(), topIn[aj].Element() );
  return 0;
}

/** Attempt to assign a reasonable value for theta internal coordinate for
  * atom i given that atoms j and k have known positions.
  */
int Cpptraj::Structure::Model::AssignTheta(double& theta, int ai, int aj, int ak, Topology const& topIn, Frame const& frameIn, std::vector<bool> const& atomPositionKnown)
const
{
  // Figure out hybridization and chirality of atom j.
  if (debug_ > 0)
    mprintf("DEBUG: AssignTheta for atom j : %s\n", topIn.AtomMaskName(aj).c_str());

  //enum HybridizationType { SP = 0, SP2, SP3, UNKNOWN_HYBRIDIZATION };

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
  // TODO pass in atom type info?
  AtomType::HybridizationType hybrid = GuessAtomHybridization(AJ, topIn.Atoms());
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

/// Recursive function to return depth from an atom along bonds
static int atom_depth(int& depth,
                      int at, Topology const& topIn, std::vector<bool>& visited, int maxdepth)
{
  if (depth == maxdepth) return 0;
  depth++;
  visited[at] = true;
  int depthFromHere = 1;
  for (Atom::bond_iterator bat = topIn[at].bondbegin(); bat != topIn[at].bondend(); ++bat)
  {
    if (!visited[*bat])
      depthFromHere += atom_depth( depth, *bat, topIn, visited, maxdepth );
  }
  return depthFromHere;
}

static inline double wrap360(double phi) {
  if (phi > Constants::PI)
    return phi - Constants::TWOPI;
  else if (phi < -Constants::PI)
    return phi + Constants::TWOPI;
  else
    return phi;
}

/** Attempt to assign reasonable values for phi internal coordinates for atoms
  * bonded to atom j given that atoms j k and l have known positions.
  *   j - k
  *  /     \
  * i       l
  */
//int Cpptraj::Structure::Model::AssignPhi(std::vector<InternalCoords>& IC, int aj, int ak, int al,
int Cpptraj::Structure::Model::AssignPhi(double& phi, int ai, int aj, int ak, int al,
                                         Topology const& topIn, Frame const& frameIn,
                                         std::vector<bool> const& atomPositionKnown,
                                         BuildAtom const& AtomJ)
const
{
  Atom const& AJ = topIn[aj];
  // If atom J has only 1 bond this is not needed.
  if (AJ.Nbonds() < 2) return 0;

  if (debug_ > 0) mprintf("DEBUG:\t\tNbonds: %i\n", AJ.Nbonds());
  // If atom J only has 2 bonds, ai-aj-ak-al is the only possibility.
  if (AJ.Nbonds() < 3) {
    if (debug_ > 0)
      mprintf("DEBUG:\t\tFewer than 3 bonds. Setting phi to -180.\n");
    double currentPhi = -180 * Constants::DEGRAD;
    for (int idx = 0; idx < AJ.Nbonds(); idx++) {
      if (AJ.Bond(idx) != ak) {
        //IC.push_back( InternalCoords(AJ.Bond(idx), aj, ak, al, 0, 0, currentPhi) );
        phi = currentPhi;
        break;
       }
    }
    return 0;
  }
  // Figure out hybridization and chirality of atom j.
  if (debug_ > 0)
    mprintf("DEBUG: AssignPhi for atom j : %s\n", topIn.AtomMaskName(aj).c_str());
//  std::vector<int> priority;
//  int chiralDebug = debug_;
//  if (chiralDebug > 0)
//    chiralDebug--;
//  ChiralType chirality = SetPriority(priority, aj, topIn, frameIn, chiralDebug);


  // TODO check that atom i actually ends up on the list?
  std::vector<int> const& priority = AtomJ.Priority();
  ChiralType chirality = AtomJ.Chirality();
  if (debug_ > 0) {
    mprintf("DEBUG: Original chirality around J %s is %s\n", topIn.AtomMaskName(aj).c_str(), chiralStr(chirality));
    mprintf("DEBUG:\t\tPriority around J %s(%i) is", 
            topIn.AtomMaskName(aj).c_str(), (int)atomPositionKnown[aj]);
    for (int idx = 0; idx < AJ.Nbonds(); idx++)
      mprintf(" %s(%i)", topIn.AtomMaskName(priority[idx]).c_str(), (int)atomPositionKnown[priority[idx]]);
    mprintf("\n");
  }

  // Fill in what values we can for known atoms
  std::vector<double> knownPhi( AJ.Nbonds() );
  std::vector<bool> isKnown( AJ.Nbonds(), false );
  int knownIdx = -1;
  double knownInterval = 0;
  bool hasKnownInterval = false;
  for (int idx = 0; idx < AJ.Nbonds(); idx++) {
    int atnum = priority[idx];
    if (atnum != ak && atomPositionKnown[atnum] &&
                       atomPositionKnown[aj] &&
                       atomPositionKnown[ak] &&
                       atomPositionKnown[al])
    {
      knownPhi[idx] = Torsion(frameIn.XYZ(atnum),
                              frameIn.XYZ(aj),
                              frameIn.XYZ(ak),
                              frameIn.XYZ(al));
      isKnown[idx] = true;
      if (debug_ > 0)
        mprintf("DEBUG:\t\tKnown phi for %s (pos=%i) = %g\n", topIn.AtomMaskName(atnum).c_str(), idx, knownPhi[idx]*Constants::RADDEG);
      if (knownIdx == -1) knownIdx = idx; // FIXME handle more than 1 known
      if (idx > 0 && isKnown[idx-1] && isKnown[idx]) {
        knownInterval = wrap360(knownPhi[idx] - knownPhi[idx-1]);
        hasKnownInterval = true;
      }
    }
  }
  if (hasKnownInterval)
    mprintf("DEBUG:\t\tKnown interval = %g\n", knownInterval * Constants::RADDEG);

  // If we have to assign an initial phi, make trans the longer branch
  if (knownIdx == -1) {
    std::vector<bool> visited = atomPositionKnown;
    // TODO: Ensure bonded atoms are not yet visited?
    visited[aj] = true;
    visited[ak] = true;
    //std::vector<int> depth( AJ.Nbonds() );
    int max_depth = 0;
    int max_idx = -1;
    for (int idx = 0; idx < AJ.Nbonds(); idx++) {
      int atnum = priority[idx];
      if (atnum != ak) {
        int currentDepth = 0;
        //depth[idx] = atom_depth(currentDepth, atnum, topIn, visited, 10);
        int depth = atom_depth(currentDepth, atnum, topIn, visited, 10);
        if (debug_ > 0)
          mprintf("DEBUG:\t\tAJ %s depth from %s is %i\n",
                  topIn.AtomMaskName(aj).c_str(), topIn.AtomMaskName(atnum).c_str(), depth);
        //if (knownIdx == -1 && depth[idx] < 3) {
        //  knownIdx = idx;
        //  knownPhi[idx] = 0;
        //}
        if (max_idx == -1 || depth > max_depth) {
          max_depth = depth;
          max_idx = idx;
        }
      }
    }
    mprintf("DEBUG:\t\tLongest depth is for atom %s (%i)\n", topIn.AtomMaskName(priority[max_idx]).c_str(), max_depth);
    knownIdx = max_idx;
    knownPhi[max_idx] = -180 * Constants::DEGRAD;
    isKnown[max_idx] = true;
  }

  // Sanity check
  if (knownIdx < 0) {
    mprinterr("Internal Error: AssignPhi(): knownIdx is < 0\n");
    return 1;
  }

  // The interval will be 360 / (number of bonds - 1)
  double interval = Constants::TWOPI / (AJ.Nbonds() - 1);
  if (hasKnownInterval) {
    if (chirality == IS_UNKNOWN_CHIRALITY) {
      mprintf("DEBUG: Setting chirality from known interval.\n");
      if (knownInterval < 0)
        chirality = IS_S;
      else
        chirality = IS_R;
    } else if (chirality == IS_S) {
      if (knownInterval > 0)
        mprinterr("Error: Detected chirality S does not match known interval %g\n", knownInterval*Constants::RADDEG);
    } else if (chirality == IS_R) {
      if (knownInterval < 0)
        mprinterr("Error: Detected chriality R does not match known interval %g\n", knownInterval*Constants::RADDEG);
    }
  }
  if (chirality == IS_S || chirality == IS_UNKNOWN_CHIRALITY)
    interval = -interval;
 

  if (debug_ > 0) {
    mprintf("DEBUG:\t\tStart phi is %g degrees\n", knownPhi[knownIdx]*Constants::RADDEG);
    mprintf("DEBUG:\t\tInterval is %g, chirality around J is %s\n", interval*Constants::RADDEG, chiralStr(chirality));
  }

  // Forwards from the known index
  double currentPhi = knownPhi[knownIdx];
  for (int idx = knownIdx; idx < AJ.Nbonds(); idx++) {
    int atnum = priority[idx];
    if (atnum != ak) {
      if (isKnown[idx])
        currentPhi = knownPhi[idx];
      else
        currentPhi = wrap360(currentPhi + interval);
      if (atnum == ai) phi = currentPhi;
      //IC.push_back( InternalCoords(atnum, aj, ak, al, 0, 0, currentPhi) );
      if (debug_ > 0)
        mprintf("DEBUG:\t\t\t%s (at# %i) phi= %g\n", topIn.AtomMaskName(atnum).c_str(), atnum+1, currentPhi*Constants::RADDEG);
    }
  }
  // Backwards from the known index
  currentPhi = knownPhi[knownIdx];
  for (int idx = knownIdx - 1; idx > -1; idx--) {
    int atnum = priority[idx];
    if (atnum != ak) {
      if (isKnown[idx])
        currentPhi = knownPhi[idx];
      else
        currentPhi = wrap360(currentPhi - interval);
      if (atnum == ai) phi = currentPhi;
      //IC.push_back( InternalCoords(atnum, aj, ak, al, 0, 0, currentPhi) );
      if (debug_ > 0)
        mprintf("DEBUG:\t\t\t%s (at# %i) phi= %g\n", topIn.AtomMaskName(atnum).c_str(), atnum+1, currentPhi*Constants::RADDEG);
    }
  }
/*
  double startPhi;
  if (knownIdx == -1) {
    startPhi = -180*Constants::DEGRAD;
    if (debug_ > 0) mprintf("DEBUG:\t\tNo known phi. Setting to %g.\n", startPhi*Constants::RADDEG);
    knownIdx = 0;
  } else
    startPhi = knownPhi[knownIdx];

  if (AtomJ.Chirality() == IS_R) {
    startPhi = -startPhi;
    interval = -interval;
  }
    
  // Forward direction
  double currentPhi = startPhi;
  for (int idx = knownIdx; idx < AJ.Nbonds(); idx++) {
    int atnum = priority[idx];
    if (atnum != ak) {
      if (atnum == ai) phi = currentPhi;
      if (debug_ > 0)
        mprintf("DEBUG:\t\t\t%s phi= %g\n", topIn.AtomMaskName(atnum).c_str(), currentPhi*Constants::RADDEG);
      currentPhi += interval;
      currentPhi = wrap360(currentPhi);
    }
  }
  // Reverse direction
  currentPhi = startPhi - interval;
  for (int idx = knownIdx - 1; idx > -1; idx--) {
    int atnum = priority[idx];
    if (atnum != ak) {
      if (atnum == ai) phi = currentPhi;
      if (debug_ > 0)
        mprintf("DEBUG:\t\t\t%s phi= %g\n", topIn.AtomMaskName(atnum).c_str(), currentPhi*Constants::RADDEG);
      currentPhi -= interval;
      currentPhi = wrap360(currentPhi);
    }
  }
*/
  return 0;
}
/*
static inline AtomType::HybridizationType getAtomHybridization(ParmHolder<AtomType> const& AT, Topology const& topIn, int ix)
{
  AtomType::HybridizationType hx;
  ParmHolder<AtomType>::const_iterator itx = AT.GetParam( TypeNameHolder( topIn[ix].Type() ) );
  if (itx == AT.end())
    hx = Cpptraj::GuessAtomHybridization(topIn[ix], topIn.Atoms());
  else
    hx = itx->second.Hybridization();
  return hx;
}
*/
/** Assign phi values around a bond. */
/*int Cpptraj::Structure::Model::AssignPhiAroundBond(int ix, int iy, Topology const& topIn, Frame const& frameIn,
                                                   std::vector<bool> const& atomPositionKnown,
                                                   ParmHolder<AtomType> const& AT)
const
{
  // Order atoms by hybridization; AX > AY
  int ax, ay;
  AtomType::HybridizationType hx = getAtomHybridization(AT, topIn, ix);
  AtomType::HybridizationType hy = getAtomHybridization(AT, topIn, iy);
  if (hx < hy) {
    ax = iy;
    ay = ix;
  } else {
    ax = ix;
    ay = iy;
  }
  mprintf("DEBUG: Assign torsions around %s (%i) - %s (%i)\n",
          topIn.AtomMaskName(ix).c_str(), (int)hx,
          topIn.AtomMaskName(iy).c_str(), (int)hy);
  mprintf("DEBUG: Ordered by hybridization: %s %s\n", 
          topIn.AtomMaskName(ax).c_str(),
          topIn.AtomMaskName(ay).c_str());
 
  return 0;
}*/