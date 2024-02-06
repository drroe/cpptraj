#include "Builder.h"
#include "BuildAtom.h"
#include "GenerateConnectivityArrays.h"
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
  params_(0),
  currentZmatrix_(0),
  currentTop_(0),
  currentFrm_(0),
  hasPosition_(0),
  newZmatrix_(0)
{}

/** Set optional parameter set. */
void Cpptraj::Structure::Builder::SetParameters(ParameterSet const* paramsIn) {
  if (paramsIn == 0) {
    mprinterr("Internal Error: Builder::SetParmaters called with null set.\n");
    return;
  }
  params_ = paramsIn;
}

/** Set optional Zmatrix. */
void Cpptraj::Structure::Builder::SetZmatrix(Zmatrix const* zmatrixIn) {
  if (zmatrixIn == 0) {
    mprinterr("Internal Error: Builder::SetZmatrix called with null set.\n");
    return;
  }
  currentZmatrix_ = zmatrixIn;
}

// -----------------------------------------------------------------------------
/** Get length from parameter set if present.
  * \return 1 if a length parameter was found.
  */
int Cpptraj::Structure::Builder::getLengthParam(double& dist, int ai, int aj, Topology const& topIn)
const
{
  if (params_ != 0 && topIn[ai].HasType() && topIn[aj].HasType()) {
    TypeNameHolder btypes(2);
    btypes.AddName( topIn[ai].Type() );
    btypes.AddName( topIn[aj].Type() );
    ParmHolder<BondParmType>::const_iterator it = params_->BP().GetParam( btypes );
    if (it != params_->BP().end()) {
      dist = it->second.Req();
      mprintf("DEBUG: Found bond parameter for %s (%s) - %s (%s): req=%g rk=%g\n",
              topIn.AtomMaskName(ai).c_str(), *(topIn[ai].Type()),
              topIn.AtomMaskName(aj).c_str(), *(topIn[aj].Type()),
              it->second.Req(), it->second.Rk());
      return 1;
    }
  }
  return 0;
}

/** Assign reasonable value for bond distance. */
int Cpptraj::Structure::Builder::AssignLength(double& dist, int ai, int aj, Topology const& topIn, Frame const& frameIn, Barray const& atomPositionKnown)
const
{
  if (atomPositionKnown[ai] && atomPositionKnown[aj]) {
    dist = sqrt( DIST2_NoImage( frameIn.XYZ(ai), frameIn.XYZ(aj) ) );
    return 0;
  }

  // One or both positions unknown. Use estimated bond length or parameters.
  if (getLengthParam(dist, ai, aj, topIn)) return 0;

  // Default to bond length based on elements
  dist = Atom::GetBondLength( topIn[ai].Element(), topIn[aj].Element() );
  return 0;
}

/** Get angle parameter if present.
  * \return 1 if parameter found.
  */
int Cpptraj::Structure::Builder::getAngleParam(double& theta, int ai, int aj, int ak, Topology const& topIn)
const
{
  if (params_ != 0 &&
      topIn[ai].HasType() &&
      topIn[aj].HasType() &&
      topIn[ak].HasType())
  {
    TypeNameHolder atypes(3);
    atypes.AddName( topIn[ai].Type() );
    atypes.AddName( topIn[aj].Type() );
    atypes.AddName( topIn[ak].Type() );
    ParmHolder<AngleParmType>::const_iterator it = params_->AP().GetParam( atypes );
    if (it != params_->AP().end()) {
      theta = it->second.Teq();
      mprintf("DEBUG: Found angle parameter for %s (%s) - %s (%s) - %s (%s): teq=%g tk=%g\n",
                topIn.AtomMaskName(ai).c_str(), *(topIn[ai].Type()),
                topIn.AtomMaskName(aj).c_str(), *(topIn[aj].Type()),
                topIn.AtomMaskName(ak).c_str(), *(topIn[ak].Type()),
                it->second.Teq()*Constants::RADDEG, it->second.Tk());
      return 1;
    }
  }
  return 0;
}

/** Attempt to assign a reasonable value for theta internal coordinate for
  * atom i given that atoms j and k have known positions.
  */
int Cpptraj::Structure::Builder::AssignTheta(double& theta, int ai, int aj, int ak, Topology const& topIn, Frame const& frameIn, Barray const& atomPositionKnown)
const
{
  if (debug_ > 0)
    mprintf("DEBUG: AssignTheta for atom j : %s\n", topIn.AtomMaskName(aj).c_str());
  Atom const& AJ = topIn[aj];

  // Sanity check
  if (AJ.Nbonds() < 2) {
    mprinterr("Internal Error: AssignTheta() called for atom J %s with fewer than 2 bonds.\n", topIn.AtomMaskName(aj).c_str());
    return 1;
  }

  // Check if all positions are already known
  if (atomPositionKnown[ai] && atomPositionKnown[aj] && atomPositionKnown[ak])
  {
    theta = CalcAngle(frameIn.XYZ(ai), frameIn.XYZ(aj), frameIn.XYZ(ak));
    return 0;
  }

  // Figure out angles from known atoms + ICs
  int nAngles = (AJ.Nbonds() * (AJ.Nbonds()-1)) / 2;
  if (nAngles == 3) {
    mprintf("DEBUG: Expect %i angles around AJ.\n", nAngles);
    std::vector<double> thetaVals;
    int tgtIdx = -1;
    int numKnown = 0;
    for (int idx0 = 0; idx0 < AJ.Nbonds(); idx0++) {
      int atomi = AJ.Bond(idx0);
      for (int idx1 = idx0 + 1; idx1 < AJ.Nbonds(); idx1++) {
        int atomk = AJ.Bond(idx1);
        mprintf("DEBUG: AssignTheta(): Angle %zu atoms %i - %i - %i\n", thetaVals.size(), atomi+1, aj+1, atomk+1);
        if (ai == atomi && ak == atomk)
          tgtIdx = (int)thetaVals.size();
        double ajTheta = 0;
        if (atomPositionKnown[atomi] && atomPositionKnown[aj] && atomPositionKnown[atomk]) {
          ajTheta = CalcAngle(frameIn.XYZ(atomi), frameIn.XYZ(aj), frameIn.XYZ(atomk));
          numKnown++;
          mprintf("DEBUG: AssignTheta(): Known angle centered on atom J: %s - %s - %s = %g\n",
                  topIn.LeapName(atomi).c_str(), 
                  topIn.LeapName(aj).c_str(), 
                  topIn.LeapName(atomk).c_str(),
                  ajTheta * Constants::RADDEG);
        } else if (currentZmatrix_ != 0) { // FIXME faster search?
          for (Zmatrix::const_iterator ic = currentZmatrix_->begin(); ic != currentZmatrix_->end(); ++ic) {
            if (ic->AtI() == atomi && ic->AtJ() == aj && ic->AtK() == atomk) {
              // TODO: Check that repeats are equal?
              ajTheta = ic->Theta() * Constants::DEGRAD;
              numKnown++;
              mprintf("DEBUG: AssignTheta(): IC angle centered on atomJ: %s - %s - %s = %g\n",
                      topIn.LeapName(atomi).c_str(), 
                      topIn.LeapName(aj).c_str(), 
                      topIn.LeapName(atomk).c_str(),
                      ic->Theta());
              break;
            }
          }
        }
        thetaVals.push_back( ajTheta );
      } // END inner loop over bonds to AJ
    } // END outer loop over bonds to AJ
    mprintf("DEBUG: AssignTheta(): thetaVals=");
    for (std::vector<double>::const_iterator it = thetaVals.begin(); it != thetaVals.end(); ++it) {
      mprintf(" %g", *it * Constants::RADDEG);
      if (tgtIdx == it - thetaVals.begin()) mprintf("*");
    }
    mprintf("\n");
    if (tgtIdx != -1 && numKnown >= 2) {
      double sumTheta = 0;
      for (std::vector<double>::const_iterator it = thetaVals.begin(); it != thetaVals.end(); ++it) {
        if (it - thetaVals.begin() != tgtIdx) {
          double tval = *it;
          if (tval < 0)
            tval += Constants::TWOPI;
          else if (tval > Constants::TWOPI)
            tval -= Constants::TWOPI;
          sumTheta += tval;
        }
      }
      theta = Constants::TWOPI - sumTheta;
      mprintf("DEBUG: AssignTheta(): Setting from existing atoms/ICs: %g\n", theta * Constants::RADDEG);
      return 0;
    }
  }

  // See if a parameter is defined for these atom types
  if (getAngleParam(theta, ai, aj, ak, topIn)) return 0;

  // Figure out hybridization and chirality of atom j.
  if (debug_ > 0) {
    mprintf("DEBUG:\t\tI %s Nbonds: %i\n", topIn[ai].ElementName(), topIn[ai].Nbonds());
    mprintf("DEBUG:\t\tJ %s Nbonds: %i\n", AJ.ElementName(), AJ.Nbonds());
    mprintf("DEBUG:\t\tK %s Nbonds: %i\n", topIn[ak].ElementName(), topIn[ak].Nbonds());
  }
  AtomType::HybridizationType hybrid = AtomType::UNKNOWN_HYBRIDIZATION;
  // Check params for hybrid
  if (params_ != 0) {
    ParmHolder<AtomType>::const_iterator it = params_->AT().GetParam( TypeNameHolder(AJ.Type()) );
    if (it != params_->AT().end())
      hybrid = it->second.Hybridization();
  }
  // Guess hybrid if needed
  if (hybrid == AtomType::UNKNOWN_HYBRIDIZATION)
    hybrid = GuessAtomHybridization(AJ, topIn.Atoms());
  // Set from number of bonds if still unknown. This is a pretty crude guess.
  if (hybrid == AtomType::UNKNOWN_HYBRIDIZATION) {
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
  }

  return 0;
}

/** Calculate internal coordinate for atoms i j k l with known positions. */
Cpptraj::Structure::InternalCoords Builder::calcKnownAtomIc(int ai, int aj, int ak, int al, Frame const& frameIn)
{
  double newDist = DIST2_NoImage(frameIn.XYZ(ai), frameIn.XYZ(aj));

  double newTheta = CalcAngle( frameIn.XYZ(ai), frameIn.XYZ(aj), frameIn.XYZ(ak) );

  double newPhi = Torsion( frameIn.XYZ(ai),
                           frameIn.XYZ(aj),
                           frameIn.XYZ(ak),
                           frameIn.XYZ(al) );
  return InternalCoords(ai, aj, ak, al, newDist, newTheta*Constants::RADDEG, newPhi*Constants::RADDEG);
}

/** Insert internal coordinates with bond i-j, angle i-j-k, and torsion i-j-k-l. */
int Cpptraj::Structure::Builder::insertIc(Zmatrix& zmatrix,
                                        int ai, int aj, int ak, int al, double newPhi,
                                        Topology const& topIn, Frame const& frameIn,
                                        Barray const& atomPositionKnown)
const
{
  if (atomPositionKnown[ai]) {
    mprintf("DEBUG:\tAtom position already known for %s, skipping IC.\n", topIn.AtomMaskName(ai).c_str());
    return 0;
  }
  double newDist = 0;
  if (AssignLength(newDist, ai, aj, topIn, frameIn, atomPositionKnown)) {
    mprinterr("Error: AssignLength failed for %s - %s \n",
              topIn.AtomMaskName(ai).c_str(), topIn.AtomMaskName(aj).c_str());
    return 1;
  }
  double newTheta = 0;
  if (AssignTheta(newTheta, ai, aj, ak, topIn, frameIn, atomPositionKnown)) {
    mprinterr("Error: AssignTheta failed for %s - %s - %s\n",
              topIn.AtomMaskName(ai).c_str(),
              topIn.AtomMaskName(aj).c_str(),
              topIn.AtomMaskName(ak).c_str());
    return 1;
  }
  zmatrix.AddIC( InternalCoords(ai, aj, ak, al, newDist, newTheta*Constants::RADDEG, newPhi*Constants::RADDEG) );
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

/// Wrap given value between -PI and PI
static inline double wrap360(double phi) {
  if (phi > Constants::PI)
    return phi - Constants::TWOPI;
  else if (phi < -Constants::PI)
    return phi + Constants::TWOPI;
  else
    return phi;
}

/** Assign internal coordinates for atoms I for torsions around J-K-L. */
int Cpptraj::Structure::Builder::AssignICsAroundBond(Zmatrix& zmatrix,
                                                   int aj, int ak, int al,
                                                   Topology const& topIn, Frame const& frameIn,
                                                   Barray const& atomPositionKnown,
                                                   BuildAtom const& AtomJ)
const
{
  mprintf("DEBUG: AssignICsAroundBond: X - %s - %s - %s  %i - %i - %i\n",
          topIn.AtomMaskName(aj).c_str(),
          topIn.AtomMaskName(ak).c_str(),
          topIn.AtomMaskName(al).c_str(),
          aj+1, ak+1, al+1);
  // Ideally, atoms J K and L should be known
  if (!atomPositionKnown[aj] ||
      !atomPositionKnown[ak] ||
      !atomPositionKnown[al])
  {
    mprintf("Warning: AssignICsAroundBond(): Not all atom positions known.\n"
            "Warning: %i %s (%i) - %i %s (%i) - %i %s (%i)\n",
              aj+1, topIn.AtomMaskName(aj).c_str(), (int)atomPositionKnown[aj],
              ak+1, topIn.AtomMaskName(ak).c_str(), (int)atomPositionKnown[ak],
              al+1, topIn.AtomMaskName(al).c_str(), (int)atomPositionKnown[al]);
    //return 1;
  }
  Atom const& AJ = topIn[aj];
  // If atom J has only 1 bond this is not needed.
  if (AJ.Nbonds() < 2) return 0;
  

  if (debug_ > 0) mprintf("DEBUG:\t\tNbonds: %i\n", AJ.Nbonds());
  // If atom J only has 2 bonds, ai-aj-ak-al is the only possibility.
  if (AJ.Nbonds() < 3) {
    if (debug_ > 0)
      mprintf("DEBUG:\t\tFewer than 3 bonds. Setting phi to -180.\n");
    double newPhi = -180 * Constants::DEGRAD;
    for (int idx = 0; idx < AJ.Nbonds(); idx++) {
      if (AJ.Bond(idx) != ak) {
        int ai = AJ.Bond(idx);
        if (insertIc(zmatrix, ai, aj, ak, al, newPhi, topIn, frameIn, atomPositionKnown)) return 1;
        break;
      }
    }
    return 0;
  } // END only 2 bonds

  // 3 or more bonds.
  std::vector<int> const& priority = AtomJ.Priority();
  ChiralType chirality = AtomJ.Chirality();

  if (debug_ > 0) {
    mprintf("DEBUG:\t\tOriginal chirality around J %s is %s\n", topIn.AtomMaskName(aj).c_str(), chiralStr(chirality));
    mprintf("DEBUG:\t\tPriority around J %s(%i) is", 
            topIn.AtomMaskName(aj).c_str(), (int)atomPositionKnown[aj]);
    for (int idx = 0; idx < AJ.Nbonds(); idx++)
      mprintf(" %s(%i)", topIn.AtomMaskName(priority[idx]).c_str(), (int)atomPositionKnown[priority[idx]]);
    for (int idx = 0; idx < AJ.Nbonds(); idx++)
      mprintf(" %i", priority[idx]);
    mprintf("\n");
  }

  // Get index of atom K in the priority list, starting from least priority.
  int kPriorityIdx = -1;
  for (int idx = AJ.Nbonds()-1; idx > -1; idx--) {
    if (priority[idx] == ak) {
      kPriorityIdx = AJ.Nbonds() - 1 - idx;
      break;
    }
  }
  mprintf("DEBUG:\t\tK priority index is %i\n", kPriorityIdx);
  if (kPriorityIdx < 0) {
    mprinterr("Error: Could not find atom K %s in atom J %s bond list.\n",
              topIn.AtomMaskName(ak).c_str(), topIn.AtomMaskName(aj).c_str());
    return 1;
  }

  // Fill in what values we can for known atoms
  std::vector<double> knownPhi( AJ.Nbonds() );
  Barray isKnown( AJ.Nbonds(), false );
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

  // Check known interval if set
  if (hasKnownInterval) {
    mprintf("DEBUG:\t\tKnown interval = %g\n", knownInterval * Constants::RADDEG);
    if (chirality == IS_UNKNOWN_CHIRALITY) {
      mprintf("DEBUG:\t\tSetting chirality from known interval.\n");
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

  // If still no chirality use the detected orientation
  if (chirality == IS_UNKNOWN_CHIRALITY) {
    chirality = AtomJ.Orientation();
    mprintf("Warning: Unknown chirality around %s; using detected orientation of %s\n",
            topIn.AtomMaskName(aj).c_str(), chiralStr(chirality));
  }

  // Determine the interval
  bool intervalIsSet = false;
  double interval = 0;
  if (params_ != 0) {
    ParmHolder<AtomType>::const_iterator it = params_->AT().GetParam( TypeNameHolder(AJ.Type()) );
    if (it != params_->AT().end()) {
      if (it->second.Hybridization() == AtomType::SP2) {
        interval = 180 * Constants::DEGRAD;
        intervalIsSet = true;
      } else if (it->second.Hybridization() == AtomType::SP3) {
        interval = 120 * Constants::DEGRAD;
        intervalIsSet = true;
      }
    }
    if (intervalIsSet) mprintf("DEBUG:\t\tInterval was set from atom J hybridization.\n");
  }
  if (!intervalIsSet) {
    // The interval will be 360 / (number of bonds - 1)
    interval = Constants::TWOPI / (AJ.Nbonds() - 1);
    mprintf("DEBUG:\t\tInterval was set from number of bonds.\n");
  }

  // Adjust interval based on chirality and priority of the k atom
  if (chirality == IS_S || chirality == IS_UNKNOWN_CHIRALITY)
    interval = -interval;
  if ( (kPriorityIdx%2) == 0 ) {
    mprintf("DEBUG:\t\tFlipping interval based on priority index of %i\n", kPriorityIdx);
    interval = -interval;
  }

  // If there is a known interval, compare it to the determined one.
  if (hasKnownInterval) {
    double deltaInterval = fabs(interval - knownInterval);
    mprintf("DEBUG:\t\tDetermined interval %g, known interval %g, delta %g\n",
            interval*Constants::RADDEG, knownInterval*Constants::RADDEG, deltaInterval*Constants::RADDEG);
    interval = fabs(knownInterval);
  }

  // If we have to assign an initial phi, make trans the longer branch
  if (knownIdx == -1) {
    Barray visited = atomPositionKnown;
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
    knownPhi[max_idx] = interval;// -180 * Constants::DEGRAD;
    isKnown[max_idx] = true;
  }

  // Sanity check
  if (knownIdx < 0) {
    mprinterr("Internal Error: AssignPhi(): knownIdx is < 0\n");
    return 1;
  }

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
      //if (atnum == ai) phi = currentPhi;
      //IC.push_back( InternalCoords(atnum, aj, ak, al, 0, 0, currentPhi) );
      if (insertIc(zmatrix, atnum, aj, ak, al, currentPhi, topIn, frameIn, atomPositionKnown)) return 1;
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
      //if (atnum == ai) phi = currentPhi;
      //IC.push_back( InternalCoords(atnum, aj, ak, al, 0, 0, currentPhi) );
      if (insertIc(zmatrix, atnum, aj, ak, al, currentPhi, topIn, frameIn, atomPositionKnown)) return 1;
      if (debug_ > 0)
        mprintf("DEBUG:\t\t\t%s (at# %i) phi= %g\n", topIn.AtomMaskName(atnum).c_str(), atnum+1, currentPhi*Constants::RADDEG);
    }
  }

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
  Barray hasIC(combinedTop.Natom(), false);
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
/*static inline int known_count(int ires, Topology const& topIn, std::vector<bool> const& hasPosition)
{
  int count = 0;
  for (int at = topIn.Res(ires).FirstAtom(); at != topIn.Res(ires).LastAtom(); at++)
    if (hasPosition[at])
      count++;
  return count;
}*/

/** Model the internal coordinates around a bond between two residues defined
  * by bondAt0 (in residue0) and bondAt1 in (residue1).
  *
  */
int Builder::ModelCoordsAroundBond(Frame const& frameIn, Topology const& topIn, int bondAt0, int bondAt1,
                                   Zmatrix& zmatrix0, Barray const& hasPosition)
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
  // Ensure atA belongs to the less-known fragment. atB fragment will be "known".
//  int known0 = known_count(res0, topIn, hasPosition);
//  int known1 = known_count(res1, topIn, hasPosition);
////  int atA, atB;
//  if (known0 > known1) {
//    mprinterr("Internal Error: ModelCoordsAroundBond(): Residue0 '%s' is more well-known than Residue1 '%s'\n",
//              topIn.TruncResNameNum( topIn[bondAt0].ResNum() ).c_str(),
//              topIn.TruncResNameNum( topIn[bondAt1].ResNum() ).c_str());
//    return 1;
//  }
  int atA = bondAt0;
  int atB = bondAt1;
  Zmatrix* zA = &zmatrix0; // FIXME
//  Zmatrix const* zA;
//  if (known0 < known1) {
//    // Fragment 1 is better-known 
//    atA = bondAt0;
//    atB = bondAt1;
//    zA  = zmatrix0;
//    //zB  = zmatrix1;
//  } else {
//    // Fragment 0 is better or equally known 
//    atA = bondAt1;
//    atB = bondAt0;
//    zA  = zmatrix1;
//    //zB  = zmatrix0;
//  }

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
  Barray hasIC(topIn.Natom(), false);
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
  // Add the remaining ICs FIXME check for duplicates
  for (Zmatrix::const_iterator it = zmatrix0.begin(); it != zmatrix0.end(); ++it)
    bondZmatrix.AddIC( *it );
  if (debug_ > 0)
    bondZmatrix.print(&topIn);
  zmatrix0 = bondZmatrix;
//  if (bondZmatrix.SetToFrame( frameIn, hasPosition )) {
//    mprinterr("Error: Conversion from bondZmatrix to Cartesian coords failed.\n");
//    return 1;
//  }

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
                                Barray const& atomPositionKnown,
                                Barray const& hasICin,
                                BuildAtom const& AtomA, BuildAtom const& AtomB)
const
{
  if (debug_ > 0)
    mprintf("DEBUG: SetupICsAroundBond: atA= %s (%i)  atB= %s (%i) total # atoms %i\n",
            topIn.AtomMaskName(atA).c_str(), atA+1,
            topIn.AtomMaskName(atB).c_str(), atB+1,
            topIn.Natom());
//  zmatrix.clear();

  //Barray hasIC( topIn.Natom(), false );
  Barray hasIC = hasICin;
  unsigned int nHasIC = 0;
  for (Barray::const_iterator it = hasIC.begin(); it != hasIC.end(); ++it) {
    if (*it) {
      nHasIC++;
      mprintf("DEBUG:\tAtom %s already has an IC.\n", topIn.AtomMaskName(it-hasIC.begin()).c_str());
    }
  }
  // Mark known atoms as already having IC
  for (Barray::const_iterator it = atomPositionKnown.begin();
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
  //Cpptraj::Structure::Model model;
  //model.SetDebug( modelDebug ); //FIXME
  //model.SetDebug( 1 );
  // FIXME
  mprintf("DEBUG: ------------------------------------------------\n");
  //std::vector<InternalCoords> tmpic;
  // I J: Set up ICs for X atB K L
  if (AssignICsAroundBond(zmatrix, atB, atk0, atl0, topIn, frameIn, atomPositionKnown, AtomB)) {
    mprinterr("Error: AssignICsAroundBond (I J) failed.\n");
    return 1;
  }
  // J K: Set up ICs for X atA atB K
  if (AssignICsAroundBond(zmatrix, atA, atB, atk0, topIn, frameIn, atomPositionKnown, AtomA)) {
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
        if (AssignICsAroundBond(zmatrix, *iat, atA, atB, topIn, frameIn, atomPositionKnown, AtomC)) {
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

/** For existing torsions, see if all coordinates in that torsion
  * exist. If so, update the IC from the existing coordinates.
  */
int Builder::UpdateICsFromFrame(Zmatrix& zmatrix, Frame const& frameIn, int ires, Topology const& topIn, Barray const& hasPosition)
const
{
  // Update bond/angle values for atoms with no position.
//  for (unsigned int idx = 0; idx != zmatrix.N_IC(); idx++)
//  {
/*    InternalCoords& thisIc = zmatrix.ModifyIC(idx);
    if (hasPosition[thisIc.AtI()] && hasPosition[thisIc.AtJ()]) {
      // Bond exists.
      double dist = DIST2_NoImage( frameIn.XYZ(thisIc.AtI()), frameIn.XYZ(thisIc.AtJ()) );
      mprintf("DEBUG: Replacing existing bond length %g for %s - %s with existing length %g\n",
              thisIc.Dist(), topIn.AtomMaskName(thisIc.AtI()).c_str(),
              topIn.AtomMaskName(thisIc.AtJ()).c_str(), dist);
      thisIc.SetDist( dist );
      if (hasPosition[thisIc.AtK()]) {
        // Angle also exists
        double theta = CalcAngle(frameIn.XYZ(thisIc.AtI()),
                                 frameIn.XYZ(thisIc.AtJ()),
                                 frameIn.XYZ(thisIc.AtK()));
        theta *= Constants::RADDEG;
        mprintf("DEBUG: Replacing existing angle %g for %s - %s - %s with existing angle %g\n",
                thisIc.Theta(), topIn.AtomMaskName(thisIc.AtI()).c_str(),
                topIn.AtomMaskName(thisIc.AtJ()).c_str(), topIn.AtomMaskName(thisIc.AtK()).c_str(),
                theta);
        thisIc.SetTheta( theta );
    }
  }*/
/*
    if (!hasPosition[thisIc.AtI()]) {
      double dist = 0;
      if (getLengthParam( dist, thisIc.AtI(), thisIc.AtJ(), topIn )) {
        mprintf("DEBUG: Replacing existing bond length %g for %s - %s with parameter %g\n",
                thisIc.Dist(), topIn.AtomMaskName(thisIc.AtI()).c_str(),
                topIn.AtomMaskName(thisIc.AtJ()).c_str(), dist);
        thisIc.SetDist( dist );
      }
      double theta = 0;
      if (getAngleParam( theta, thisIc.AtI(), thisIc.AtJ(), thisIc.AtK(), topIn)) {
        theta *= Constants::RADDEG;
        mprintf("DEBUG: Replacing existing angle %g for %s - %s - %s with parameter %g\n",
                thisIc.Theta(), topIn.AtomMaskName(thisIc.AtI()).c_str(),
                topIn.AtomMaskName(thisIc.AtJ()).c_str(), topIn.AtomMaskName(thisIc.AtK()).c_str(),
                theta);
        thisIc.SetTheta( theta );
      }
    }
*/
//  } // END loop over zmatrix ICs
  // Update torsions
  Barray isUsed( zmatrix.N_IC(), false );
  //unsigned int Nused = 0;
  // Get list of bonds 
  BondArray myBonds = GenerateBondArray( std::vector<Residue>(1, topIn.Res(ires)), topIn.Atoms() );
  for (BondArray::const_iterator bnd = myBonds.begin(); bnd != myBonds.end(); ++bnd) {
    if (topIn[bnd->A1()].ResNum() == ires && topIn[bnd->A2()].ResNum() == ires) {
      mprintf("DEBUG: Looking at torsions around: %s - %s\n", topIn.AtomMaskName(bnd->A1()).c_str(), topIn.AtomMaskName(bnd->A2()).c_str());    // Find all ICs that share atoms 1 and 2 (J and K)
      Iarray bondICs;
      bool needsUpdate = false;
      double tDiff = 0;
      for (unsigned int idx = 0; idx != zmatrix.N_IC(); idx++)
      {
        if (!isUsed[idx]) {
          InternalCoords& thisIc = zmatrix.ModifyIC(idx);
          if ( (bnd->A1() == thisIc.AtJ() && bnd->A2() == thisIc.AtK()) ||
               (bnd->A2() == thisIc.AtJ() && bnd->A1() == thisIc.AtK()) )
          {
            // This IC has this bond at the center
            bondICs.push_back( idx );
            isUsed[idx] = true;
            //MARK(idx, isUsed, Nused);
            if (hasPosition[thisIc.AtI()] &&
                hasPosition[thisIc.AtJ()] &&
                hasPosition[thisIc.AtK()] &&
                hasPosition[thisIc.AtL()])
            {
              mprintf("DEBUG:\tMeasuring torsion of fixed atoms: %s - %s - %s - %s\n",
                      topIn.LeapName(thisIc.AtI()).c_str(),
                      topIn.LeapName(thisIc.AtJ()).c_str(),
                      topIn.LeapName(thisIc.AtK()).c_str(),
                      topIn.LeapName(thisIc.AtL()).c_str());
              InternalCoords frameIc = calcKnownAtomIc(thisIc.AtI(), thisIc.AtJ(), thisIc.AtK(), thisIc.AtL(), frameIn);
              double dTorsion = frameIc.Phi() * Constants::DEGRAD;
              double dInternalValue = thisIc.Phi() * Constants::DEGRAD;
              tDiff = (dTorsion - dInternalValue) * Constants::RADDEG;
              mprintf("DEBUG:\tdTorsion= %f  dInternalValue= %f\n", dTorsion, dInternalValue);
              thisIc = frameIc;
              needsUpdate = true;
            } // END all IC coords present
          } // END this IC matches current bond
        } // END IC is not used
      } // END loop searching for ICs matching current bond
      // If any difference was found, shift all of the torsions
      if (needsUpdate) {
        mprintf("DEBUG: Twisting torsions centered on %s - %s by %f degrees\n",
                topIn.LeapName(bnd->A1()).c_str(),
                topIn.LeapName(bnd->A2()).c_str(),
                tDiff);
        for (Iarray::const_iterator it = bondICs.begin(); it != bondICs.end(); ++it)
        {
          InternalCoords& thisIc = zmatrix.ModifyIC(*it);
          double dNew = thisIc.Phi() + tDiff;
          mprintf("DEBUG:\tTwisting torsion for atoms: %s-%s-%s-%s\n",
                  topIn.AtomMaskName(thisIc.AtI()).c_str(),
                  topIn.AtomMaskName(thisIc.AtJ()).c_str(),
                  topIn.AtomMaskName(thisIc.AtK()).c_str(),
                  topIn.AtomMaskName(thisIc.AtL()).c_str());
          mprintf("DEBUG:\t------- From %f to %f\n", thisIc.Phi(), dNew);
          thisIc.SetPhi( dNew );
        }
      } // END ICs need update
    } // END both bond atoms belong to this residue
  } // END loop over bonds
  return 0;
}

// ------------------------------------------------------------------------------
/// Used to track atoms for mock externals
class MockAtom {
  public:
    /// CONSTRUCTOR
    MockAtom() : idx_(-1), pos_(0.0), known_(false) {}
    /// CONSTRUCTOR - index, position
    MockAtom(int idx, Vec3 const& pos) : idx_(idx), pos_(pos), known_(true) {}
    /// CONSTRUCTOR - index, unknown position
    MockAtom(int idx) : idx_(idx), pos_(0.0), known_(false) {}
    /// Set position
    void SetPos(Vec3 const& p) { pos_ = p; known_ = true; }
    /// Set unknown
    void SetUnknown() { pos_ = Vec3(0.0); known_ = false; }

    int Idx()         const { return idx_; }
    Vec3 const& Pos() const { return pos_; }
    bool Known()      const { return known_; }
  private:
    int idx_;    ///< Atom index
    Vec3 pos_;   ///< Atom position
    bool known_; ///< True if atom position is known
};

// -----------------------------------------------
/** Store info for modelling torsions around X-Y */
class Cpptraj::Structure::Builder::TorsionModel {
  public:
    typedef std::vector<MockAtom> Marray;
    /// CONSTRUCTOR
    TorsionModel() : dAbsolute_(0), Xorientation_(0), Yorientation_(0), axHasKnownAtoms_(false), ayHasKnownAtoms_(false) {}
    /// CONSTRUCTOR - AX and AY atom indices
    //TorsionModel(int ax, int ay) : ax_(ax), ay_(ay),  dAbsolute_(0), Xorientation_(0), Yorientation_(0) {}
    /// Initialize torsions around bonded atoms
    int InitTorsion(int, int, Frame const&, Topology const&, std::vector<bool> const&);
    /// Set up torsions around bonded atoms
    int SetupTorsion(AtomType::HybridizationType, AtomType::HybridizationType, Topology const& topIn);
    /// Build mock externals from given internals
    int BuildMockExternals(std::vector<InternalCoords> const&, Topology const&);

    /// \return Value of A-X-Y-D torsion in radians
    double Absolute() const { return dAbsolute_; }
    /// \return Value of orientation around X
    double XOrientation() const { return Xorientation_; }
    /// \return Value of orientation around Y
    double YOrientation() const { return Yorientation_; }
    /// \return Sorted atoms bonded to X excluding Y
    Marray const& SortedAx() const { return sorted_ax_; }
    /// \return Sorted atoms bonded to Y excluding X
    Marray const& SortedAy() const { return sorted_ay_; }
    /// \return Atom X
    MockAtom const& AtX() const { return atX_; }
    /// \return Atom Y
    MockAtom const& AtY() const { return atY_; }
    /// \return True if AX has known bonded atoms
    bool AxHasKnownAtoms() const { return axHasKnownAtoms_; }
    /// \return True if AY has known bonded atoms
    bool AyHasKnownAtoms() const { return ayHasKnownAtoms_; }
  private:
    static int LeapAtomWeight(Atom const&);
    //static inline std::vector<int> SiftBondedAtomsLikeLeap(unsigned int&, Atom const&, std::vector<bool> const&);
    static inline Marray SortBondedAtomsLikeLeap(unsigned int&, Topology const& topIn, Marray const&);
    static inline void swap_heaviest(Marray&, Topology const&);

    MockAtom atX_;              ///< Atom X
    MockAtom atY_;              ///< Atom Y
    Marray sorted_ax_;    ///< Hold the leap-sorted bonded atoms for atom X
    Marray sorted_ay_;    ///< Hold the leap-sorted bonded atoms for atom Y
    double dAbsolute_;    ///< Hold the value of the A-X-Y-D torsion in radians
    double Xorientation_; ///< Orientation around the X atom TODO make an enum?
    double Yorientation_; ///< Orientation around the Y atoms
    bool axHasKnownAtoms_; ///< True if any atom around ax (other than ay) is known
    bool ayHasKnownAtoms_; ///< True if any atom around ay (other than ax) is known
};
    
/** \return the LEaP 'weight' of an atom.
  * Originally used to force the 'heaviest' atoms around a torsion trans to
  * each other. The 'weight' of an atom is defined as its element number,
  * unless the atom is CARBON, then it is 1000, making it the 'heaviest' atom.
  */
int Cpptraj::Structure::Builder::TorsionModel::LeapAtomWeight(Atom const& At)
{
  if ( At.Element() == Atom::CARBON )
    return 1000;
  return At.AtomicNumber();
}

/** Place atoms with known position ahead of atoms with no known position.
  * \param firstUnknownIdx Position in output array of the first unknown atom.
  * \param At The atom to sift.
  * \param hasPosition Array indicating whether atoms have position.
  */
/*std::vector<int> 
  Cpptraj::Structure::Builder::TorsionModel::SiftBondedAtomsLikeLeap(unsigned int& firstUnknownIdx,
                                                                     Atom const& At,
                                                                     std::vector<bool> const& hasPosition)
{
  std::vector<int> out;
  out.reserve( At.Nbonds() );
  for (Atom::bond_iterator bat = At.bondbegin(); bat != At.bondend(); ++bat)
    if (hasPosition[*bat])
      out.push_back( *bat );
  firstUnknownIdx = out.size();
  for (Atom::bond_iterator bat = At.bondbegin(); bat != At.bondend(); ++bat)
    if (!hasPosition[*bat])
      out.push_back( *bat );
  return out;
}*/

static inline void swap_mock_atom(MockAtom& lhs, MockAtom& rhs) {
  MockAtom tmp = lhs;
  lhs = rhs;
  rhs = tmp;
}

void Cpptraj::Structure::Builder::TorsionModel::swap_heaviest(Marray& bondedAtoms, Topology const& topIn)
{
  // Find the index of the heaviest atom
  int iHighest = 0;
  int iPos = 0;
  for (Marray::const_iterator it = bondedAtoms.begin(); it != bondedAtoms.end(); ++it)
  {
    int iWeight = LeapAtomWeight( topIn[it->Idx()] );
    if ( iHighest < iWeight ) {
      iHighest = iWeight;
      iPos = (int)(it - bondedAtoms.begin());
    }
  }
  // If highest weight atom not already in front, swap it there.
  if (iPos != 0) swap_mock_atom( bondedAtoms[0], bondedAtoms[iPos] );
}

/** Order atoms bonded to the given atom in a manner similar to LEaP's
  * zModelOrderAtoms. In that routine, first atoms were sorted into
  * known position > unknown position. Then the heaviest atom in each
  * subgroup was swapped with the first element of that list. 
  * The ignore atom is the index of the atom this atom is bonded to that
  * forms the torsion we are interested in.
  */
std::vector<MockAtom>
  Cpptraj::Structure::Builder::TorsionModel::SortBondedAtomsLikeLeap(unsigned int& firstUnknownIdx,
                                                                     Topology const& topIn,
                                                                     std::vector<MockAtom> const& bondedAtoms)
{
  typedef std::vector<MockAtom> Marray;
  Marray known_out;
  known_out.reserve( bondedAtoms.size() );
  // Sift so that atoms with known position are at the front
  for (Marray::const_iterator it = bondedAtoms.begin(); it != bondedAtoms.end(); ++it)
    if (it->Known())
      known_out.push_back( *it );
  firstUnknownIdx = known_out.size();
  swap_heaviest(known_out, topIn);
  // Unknown atoms go at the back
  Marray unknown_out;
  unknown_out.reserve(bondedAtoms.size() - known_out.size() + 1);
  for (Marray::const_iterator it = bondedAtoms.begin(); it != bondedAtoms.end(); ++it)
    if (!it->Known())
      unknown_out.push_back( *it );
  swap_heaviest(unknown_out, topIn);

  for (Marray::const_iterator it = unknown_out.begin(); it != unknown_out.end(); ++it)
    known_out.push_back( *it );

  return known_out;
}

/** LEaP routine for determining atom chirality.
  * This is done by crossing A to B and then dotting the
  * result with C. TODO use Chirality in BuildAtom?
  * The chirality of the vectors is determined by the sign of
  * the result, which is determined by whether or not C has
  * a component in the direction AxB or in the opposite direction.
  */
static inline double VectorAtomChirality(Vec3 const& Center, Vec3 const& A, Vec3 const& B, Vec3 const& C)
{
  Vec3 vA = A - Center;
  Vec3 vB = B - Center;
  Vec3 vC = C - Center;
  Vec3 vCross = vA.Cross( vB );
  double dot = vCross * vC;
  if (dot > 0)
    return 1.0;
  else if (dot < 0)
    return -1.0;
  return 0.0;
}

/** Assuming atoms have been ordered with SortBondedAtomsLikeLeap,
  * calculate the orientation of the iB atom with respect to the
  * triangle (iA, iX, iY). This orientation will be used by the
  * CreateSpXSpX routines to determine which torsion values to use.
  */
static inline double calculateOrientation(MockAtom const& iX, MockAtom const& iA, MockAtom const& iY, MockAtom const& iB)
{
  double dOrientation = 1.0;
  if (iX.Known() &&
      iA.Known() &&
      iY.Known() &&
      iB.Known())
  {
    dOrientation = VectorAtomChirality( iX.Pos(), iA.Pos(), iY.Pos(), iB.Pos() );
  } else {
    mprinterr("Internal Error: Builder::calculateOrientation not yet set up for unknown positions.\n");
  }
  return dOrientation;
}

/** Initialize model torsion for bonded atoms. */
int Cpptraj::Structure::Builder::TorsionModel::InitTorsion(int ax, int ay,
                                                           Frame const& frameIn,
                                                           Topology const& topIn,
                                                           std::vector<bool> const& hasPosition)
{
  if (hasPosition[ax])
    atX_ = MockAtom(ax, frameIn.XYZ(ax));
  else
    atX_ = MockAtom(ax);
  if (hasPosition[ay])
    atY_ = MockAtom(ay, frameIn.XYZ(ay));
  else
    atY_ = MockAtom(ay);
  Atom const& AX = topIn[ax];
  Atom const& AY = topIn[ay];
  // Create array of AX bonded atoms
  axHasKnownAtoms_ = false;
  sorted_ax_.clear();
  sorted_ax_.reserve( AX.Nbonds() - 1 );
  for (Atom::bond_iterator bat = AX.bondbegin(); bat != AX.bondend(); ++bat) {
    if (*bat != ay) {
      if (hasPosition[*bat]) {
        axHasKnownAtoms_ = true;
        sorted_ax_.push_back( MockAtom(*bat, frameIn.XYZ(*bat)) );
      } else
        sorted_ax_.push_back( MockAtom(*bat) );
    }
  }
  // Create array of AY bonded atoms
  ayHasKnownAtoms_ = false;
  sorted_ay_.clear();
  sorted_ay_.reserve( AY.Nbonds() - 1 );
  for (Atom::bond_iterator bat = AY.bondbegin(); bat != AY.bondend(); ++bat) {
    if (*bat != ax) {
      if (hasPosition[*bat]) {
        ayHasKnownAtoms_ = true;
        sorted_ay_.push_back( MockAtom(*bat, frameIn.XYZ(*bat)) );
      } else
        sorted_ay_.push_back( MockAtom(*bat) );
    }
  }
  return 0;
}

/** Set up model torsion for bonded atoms. */
int Cpptraj::Structure::Builder::TorsionModel::SetupTorsion(AtomType::HybridizationType Hx,
                                                            AtomType::HybridizationType Hy,
                                                            Topology const& topIn)
{
  if (Hx != AtomType::UNKNOWN_HYBRIDIZATION && Hy != AtomType::UNKNOWN_HYBRIDIZATION) {
    if (Hy > Hx) {
      mprinterr("Internal Error: TorsionModel::SetupTorsion() called with AX hybrid > AY hybrid.\n");
      return 1;
    }
  }
  // Sort AX bonds
  unsigned int firstUnknownIdxX = 0;
  sorted_ax_ = SortBondedAtomsLikeLeap(firstUnknownIdxX, topIn, sorted_ax_);
  // Sort AY bonds
  unsigned int firstUnknownIdxY = 0;
  sorted_ay_ = SortBondedAtomsLikeLeap(firstUnknownIdxY, topIn, sorted_ay_);
  // Calculate the chirality around atom X
  Xorientation_ = 0;
  if (Hx == AtomType::SP3) {
    Xorientation_ = calculateOrientation( atX_, sorted_ax_[0], atY_, sorted_ax_[1] );
  }
  // Calculate the chirality around atom Y
  Yorientation_ = 0;
  if (Hy == AtomType::SP3) {
    Yorientation_ = calculateOrientation( atY_, sorted_ay_[0], atX_, sorted_ay_[1] );
  }
  // DEBUG
  Atom const& AX = topIn[atX_.Idx()];
  Atom const& AY = topIn[atY_.Idx()];
  mprintf("Orientation around: %s = %f\n", *(AX.Name()), Xorientation_);
  //for (Atom::bond_iterator bat = AX.bondbegin(); bat != AX.bondend(); ++bat) mprintf(" %s", *(topIn[*bat].Name()));
  //mprintf("}\n");
  for (Marray::const_iterator it = sorted_ax_.begin(); it != sorted_ax_.end(); ++it)
      mprintf("Atom %li: %s\n", it - sorted_ax_.begin(), *(topIn[it->Idx()].Name()));
  mprintf("Orientation around: %s = %f\n", *(AY.Name()), Yorientation_);
  //for (Atom::bond_iterator bat = AY.bondbegin(); bat != AY.bondend(); ++bat) mprintf(" %s", *(topIn[*bat].Name()));
  //mprintf("}\n");
  for (Marray::const_iterator it = sorted_ay_.begin(); it != sorted_ay_.end(); ++it)
      mprintf("Atom %li: %s\n", it - sorted_ay_.begin(), *(topIn[it->Idx()].Name()));
  // Calculate the actual torsion angle between A-X-Y-D
  if (sorted_ax_[0].Known() &&
      atX_.Known() &&
      atY_.Known() &&
      sorted_ay_[0].Known())
  {
    dAbsolute_ = Torsion( sorted_ax_[0].Pos().Dptr(),
                          atX_.Pos().Dptr(),
                          atY_.Pos().Dptr(),
                          sorted_ay_[0].Pos().Dptr() );
  } else {
    dAbsolute_ = 180.0 * Constants::DEGRAD;
  }
  mprintf("DABSOLUTE= %g\n", dAbsolute_);
  return 0;
}

/// Used to find mock atom
static inline std::vector<MockAtom>::iterator find_mock_atom(std::vector<MockAtom>& outerAtoms, int idx)
{
  std::vector<MockAtom>::iterator it = outerAtoms.begin();
  for (; it != outerAtoms.end(); ++it)
    if (it->Idx() == idx) return it;
  return outerAtoms.end();
}

/** Build mock external coordinates around the given torsion using 
  * the given internal coordinates.
  * By definition, the two central atoms will be the same for each
  * IC in iaTorsions.
  */
int Cpptraj::Structure::Builder::TorsionModel::BuildMockExternals(std::vector<InternalCoords> const& iaTorsions,
                                                                  Topology const& topIn) // DEBUG topIn for debug only
{
  if (iaTorsions.empty()) {
    mprinterr("Internal Error: Builder::buildMockExternals() called with no internal torsions.\n");
    return 1;
  }
  mprintf("=======  Started mock coords from: %s\n", topIn.LeapName(iaTorsions.front().AtI()).c_str());

  for (std::vector<InternalCoords>::const_iterator ic = iaTorsions.begin();
                                                   ic != iaTorsions.end(); ++ic)
    mprintf("------- Known torsion: %s - %s - %s - %s  %f\n",
            topIn.LeapName(ic->AtI()).c_str(),
            topIn.LeapName(ic->AtJ()).c_str(),
            topIn.LeapName(ic->AtK()).c_str(),
            topIn.LeapName(ic->AtL()).c_str(),
            ic->Phi());

  // Define coordinates for the central atoms.
  //Vec3 posX(0, 0, 0);
  //Vec3 posY(1, 0, 0);
  atX_.SetPos( Vec3(0, 0, 0) );
  atY_.SetPos( Vec3(1, 0, 0) );
  Vec3 const& posX = atX_.Pos();
  Vec3 const& posY = atY_.Pos();

  // Tell the outer atoms they do not have defined positions.
  for (Marray::iterator it = sorted_ax_.begin(); it != sorted_ax_.end(); ++it)
    it->SetUnknown();
  for (Marray::iterator it = sorted_ay_.begin(); it != sorted_ay_.end(); ++it)
    it->SetUnknown();

//  // Hold info on X-Y outer atoms
  typedef std::vector<MockAtom> Marray;
  Marray outerAtoms;

  // Define outer atoms
  for (std::vector<InternalCoords>::const_iterator ic = iaTorsions.begin();
                                                   ic != iaTorsions.end(); ++ic)
  {
    if (ic == iaTorsions.begin()) {
      // Define first outer atom as being in the XY plane
      outerAtoms.push_back( MockAtom(ic->AtI(), Vec3(1, 1, 0)) );
    } else {
      Marray::iterator mi = find_mock_atom( outerAtoms, ic->AtI() );
      if (mi == outerAtoms.end())
        outerAtoms.push_back( MockAtom(ic->AtI()) );
    }
    Marray::iterator ml = find_mock_atom( outerAtoms, ic->AtL() );
    if (ml == outerAtoms.end())
      outerAtoms.push_back( MockAtom(ic->AtL()) );
  }
  mprintf("DEBUG: Outer atoms:\n");
  for (Marray::const_iterator it = outerAtoms.begin(); it != outerAtoms.end(); ++it)
    mprintf("DEBUG:\t\t%i %4s (%i) {%f %f %f}\n", it->Idx()+1, topIn.AtomMaskName(it->Idx()).c_str(),
            (int)it->Known(), it->Pos()[0], it->Pos()[1], it->Pos()[2]);

  // Loop through the known torsions looking for those that
  // have one position defined, then build coords for the
  // other atom and mark the torsion as used.
  std::vector<bool> used( iaTorsions.size(), false );
  unsigned int nused = 0;
  for (unsigned int idx = 0; idx != iaTorsions.size(); idx++) {
    //bool gotOne = false;
    for (unsigned int jdx = 0; jdx != iaTorsions.size(); jdx++) {
      if (!used[jdx]) {
        InternalCoords const& iInt = iaTorsions[jdx];
        Marray::iterator tmpAt1 = find_mock_atom(outerAtoms, iInt.AtI());
        Marray::iterator tmpAt4 = find_mock_atom(outerAtoms, iInt.AtL());
        if (tmpAt1 == outerAtoms.end()) {
          mprinterr("Internal Error: Builder::buildMockExternals(): Outer atom I %i not found.\n", iInt.AtI()+1);
          return 1;
        }
        if (tmpAt4 == outerAtoms.end()) {
          mprinterr("Internal Error: Builder::buildMockExternals(): Outer atom L %i not found.\n", iInt.AtL()+1);
          return 1;
        }
        Vec3 maPC1, maPC2;
        Marray::iterator tgt = outerAtoms.end();
        Marray::iterator knownAt = outerAtoms.end();
        if (tmpAt4->Known()) {
          tgt     = tmpAt1;
          maPC1   = posX;
          maPC2   = posY;
          knownAt = tmpAt4;
        } else if (tmpAt1->Known()) {
          //gotOne = true; // FIXME needed?
          tgt     = tmpAt4;
          maPC1   = posY;
          maPC2   = posX;
          knownAt = tmpAt1;
        }
        if (tgt != outerAtoms.end()) {
          //gotOne = true;
          mprintf("======= Building mock coord for: %s\n", topIn.LeapName(tgt->Idx()).c_str());
          mprintf("======= Using torsion: %s - %s - %s - %s (p1known= %i, p4known= %i)\n",
                   topIn.LeapName(iInt.AtI()).c_str(),
                   topIn.LeapName(iInt.AtJ()).c_str(),
                   topIn.LeapName(iInt.AtK()).c_str(),
                   topIn.LeapName(iInt.AtL()).c_str(),
                   (int)tmpAt1->Known(), (int)tmpAt4->Known());
          // Now build the coordinate for the target atom
          tgt->SetPos( Zmatrix::AtomIposition(maPC1, maPC2, knownAt->Pos(), 1.0, 90.0, iInt.Phi()) );
          mprintf("ZMatrixAll:  %f,%f,%f\n", tgt->Pos()[0], tgt->Pos()[1], tgt->Pos()[2]);
          used[jdx] = true;
          nused++;
          break;
        }
      } // END IC not yet used
    } // END inner loop over ICs
  } // END outer loop over ICs
  if (nused < used.size()) {
    mprinterr("Error: There are %u torsions left over for mock coords.\n", used.size() - nused);
    return 1;
  }
  // Update the outer atoms for this torsion
  mprintf("DEBUG: Final outer atoms:\n");
  for (Marray::const_iterator it = outerAtoms.begin(); it != outerAtoms.end(); ++it) {
    mprintf("DEBUG:\t\t%i %4s (%i) {%f %f %f}\n", it->Idx()+1, topIn.AtomMaskName(it->Idx()).c_str(),
            (int)it->Known(), it->Pos()[0], it->Pos()[1], it->Pos()[2]);
    Marray::iterator itx = find_mock_atom(sorted_ax_, it->Idx());
    if (itx != sorted_ax_.end()) {
      *itx = *it;
    } else {
      itx = find_mock_atom(sorted_ay_, it->Idx());
      if (itx != sorted_ay_.end()) {
        *itx = *it;
      } else {
        mprinterr("Internal Error: TorsionModel::BuildMockExternals(): Could not update mock atom.\n");
        return 1;
      }
    }
  }

  return 0;
}

// -----------------------------------------------------------------------------
/** Model torsion */
void Builder::ModelTorsion(TorsionModel const& MT, unsigned int iBondX, unsigned int iBondY, double dvalIn)
{
  if (iBondX >= MT.SortedAx().size() ||
      iBondY >= MT.SortedAy().size())
    return;
  mprintf("CALLING ModelTorsion for iBondX=%u iBondY=%u dVal=%g\n",iBondX,iBondY,dvalIn*Constants::RADDEG);
  int aa = MT.SortedAx()[iBondX].Idx();
  int ax = MT.AtX().Idx();
  int ay = MT.AtY().Idx();
  int ad = MT.SortedAy()[iBondY].Idx();
  // Get bond lengths FIXME deal with unknown positions
  double l0, l1;
  if (AssignLength(l0, aa, ax, *currentTop_, *currentFrm_, *hasPosition_)) {
    mprinterr("Error: Could not assign length between %s and %s\n",
              currentTop_->AtomMaskName(aa).c_str(),
              currentTop_->AtomMaskName(ax).c_str());
    return;
  }
  if (AssignLength(l1, ad, ay, *currentTop_, *currentFrm_, *hasPosition_)) {
    mprinterr("Error: Could not assign length between %s and %s\n",
              currentTop_->AtomMaskName(ad).c_str(),
              currentTop_->AtomMaskName(ay).c_str());
    return;
  }
  // Get angles
  double t0, t1;
  if (AssignTheta(t0, aa, ax, ay, *currentTop_, *currentFrm_, *hasPosition_)) {
    mprinterr("Error: Could not assign angle between %s and %s and %s\n",
              currentTop_->AtomMaskName(aa).c_str(),
              currentTop_->AtomMaskName(ax).c_str(),
              currentTop_->AtomMaskName(ay).c_str());
    return;
  }
  if (AssignTheta(t1, ad, ay, ax, *currentTop_, *currentFrm_, *hasPosition_)) {
    mprinterr("Error: Could not assign angle between %s and %s and %s\n",
              currentTop_->AtomMaskName(ad).c_str(),
              currentTop_->AtomMaskName(ay).c_str(),
              currentTop_->AtomMaskName(ax).c_str());
    return;
  }
  // If the coordinates for the atoms are defined then
  // measure the torsion angle between them and use that for
  // the internal.
  MockAtom const& AA = MT.SortedAx()[iBondX];
  MockAtom const& AX = MT.AtX();
  MockAtom const& AY = MT.AtY();
  MockAtom const& AD = MT.SortedAy()[iBondY];

  double phiVal = dvalIn;
  if (AA.Known() &&
      AX.Known() &&
      AY.Known() &&
      AD.Known())
  {
    phiVal = Torsion( AA.Pos().Dptr(),
                      AX.Pos().Dptr(),
                      AY.Pos().Dptr(),
                      AD.Pos().Dptr() );
  } else {
    mprinterr("Internal Error: Need to implement torsion lookup.\n");
  }
  mprintf("++++Torsion INTERNAL: %f to %s - %s - %s - %s\n", phiVal*Constants::RADDEG,
          currentTop_->LeapName(aa).c_str(),
          currentTop_->LeapName(ax).c_str(),
          currentTop_->LeapName(ay).c_str(),
          currentTop_->LeapName(ad).c_str());
  newZmatrix_->AddIC( InternalCoords(aa, ax, ay, ad, l0, t0*Constants::RADDEG, phiVal*Constants::RADDEG) );
  newZmatrix_->AddIC( InternalCoords(ad, ay, ax, aa, l1, t1*Constants::RADDEG, phiVal*Constants::RADDEG) );
}

/** Create torsions around SP3-SP3. */
void Builder::createSp3Sp3Torsions(TorsionModel const& MT) {
  // First twist the torsion so that the AD torsion has
  // the same absolute angle that is measured
  // and twist all the others with it.
  static const double PIOVER3 = Constants::PI / 3.0;
  double dADOffset =         MT.Absolute() - Constants::PI ;
  double d180      =         Constants::PI + dADOffset     ;
  double dm60      =         -PIOVER3      + dADOffset     ;
  double d60       =          PIOVER3      + dADOffset     ;

  if ( MT.XOrientation() > 0.0 ) {
    if ( MT.YOrientation() > 0.0 ) {
      ModelTorsion( MT, 0, 0, d180);
      ModelTorsion( MT, 0, 1, dm60);
      ModelTorsion( MT, 0, 2, d60);
      ModelTorsion( MT, 1, 0, dm60);
      ModelTorsion( MT, 1, 1, d60);
      ModelTorsion( MT, 1, 2, d180);
      ModelTorsion( MT, 2, 0, d60);
      ModelTorsion( MT, 2, 1, d180);
      ModelTorsion( MT, 2, 2, dm60);
    } else {
      ModelTorsion( MT, 0, 0,  d180);
      ModelTorsion( MT, 0, 1,  d60);
      ModelTorsion( MT, 0, 2,  dm60);
      ModelTorsion( MT, 1, 0,  dm60);
      ModelTorsion( MT, 1, 1,  d180);
      ModelTorsion( MT, 1, 2,  d60);
      ModelTorsion( MT, 2, 0,  d60);
      ModelTorsion( MT, 2, 1,  dm60);
      ModelTorsion( MT, 2, 2,  d180);
    }
  } else {
    if ( MT.YOrientation() > 0.0 ) {
      ModelTorsion( MT, 0, 0, d180);
      ModelTorsion( MT, 0, 1, dm60);
      ModelTorsion( MT, 0, 2, d60);
      ModelTorsion( MT, 1, 0, d60);
      ModelTorsion( MT, 1, 1, d180);
      ModelTorsion( MT, 1, 2, dm60);
      ModelTorsion( MT, 2, 0, dm60);
      ModelTorsion( MT, 2, 1, d60);
      ModelTorsion( MT, 2, 2, d180);
    } else {
      ModelTorsion( MT, 0, 0, d180);
      ModelTorsion( MT, 0, 1, d60);
      ModelTorsion( MT, 0, 2, dm60);
      ModelTorsion( MT, 1, 0, d60);
      ModelTorsion( MT, 1, 1, dm60);
      ModelTorsion( MT, 1, 2, d180);
      ModelTorsion( MT, 2, 0, dm60);
      ModelTorsion( MT, 2, 1, d180);
      ModelTorsion( MT, 2, 2, d60);
    }
  }

  return;
}

/** Create torsions around SP3-SP2. */
void Builder::createSp3Sp2Torsions(TorsionModel const& MT) {
  // First twist the torsion so that the AD torsion has
  // the same absolute angle that is measured
  // and twist all the others with it.
  static const double PIOVER3   = Constants::PI / 3.0;
  static const double PIOVER3x2 = PIOVER3 * 2.0;
  double dADOffset =         MT.Absolute() - Constants::PI ;
  double   d180    =         Constants::PI + dADOffset     ;
  double   dm60    =         -PIOVER3      + dADOffset     ;
  double   d60     =          PIOVER3      + dADOffset     ;
  double  dm120    =         -PIOVER3x2    + dADOffset     ;
  double  d120     =          PIOVER3x2    + dADOffset     ;
  double  d0       =                         dADOffset     ;

  if ( MT.XOrientation() > 0.0 ) {
    ModelTorsion( MT, 0, 0, d180);
    ModelTorsion( MT, 0, 1, d0);
    ModelTorsion( MT, 1, 0, dm60);
    ModelTorsion( MT, 1, 1, d120);
    ModelTorsion( MT, 2, 0, d60);
    ModelTorsion( MT, 2, 1, dm120);
  } else {
    ModelTorsion( MT, 0, 0, d180);
    ModelTorsion( MT, 0, 1, d0);
    ModelTorsion( MT, 1, 0, d60);
    ModelTorsion( MT, 1, 1, dm120);
    ModelTorsion( MT, 2, 0, dm60);
    ModelTorsion( MT, 2, 1, d120);
  }

  return;
}

void Builder::createSp2Sp2Torsions(TorsionModel const& MT) {
  // First twist the torsion so that the AD torsion has
  // the same absolute angle that is measured
  // and twist all the others with it.
  double dADOffset =         MT.Absolute() - Constants::PI ;
  double d180      =         Constants::PI + dADOffset     ;
  double d0        =                 dADOffset             ;
  mprintf("In ModelCreateSp2Sp2Torsions, dAbsolute is %g, dADOffset= %g, d180= %g, d0= %g\n",
          MT.Absolute()*Constants::RADDEG,
          dADOffset * Constants::RADDEG,
          d180 * Constants::RADDEG,
          d0 * Constants::RADDEG);

  ModelTorsion( MT, 0, 0, d180 );
  ModelTorsion( MT, 0, 1, d0 );
  ModelTorsion( MT, 1, 0, d0 );
  ModelTorsion( MT, 1, 1, d180 );
  return;
}

/** Find any existing internal coords around ax-ay. */
std::vector<InternalCoords> Builder::getExistingInternals(int ax, int ay) const {
  std::vector<InternalCoords> iTorsions;
  if (currentZmatrix_ != 0) {
    for (Zmatrix::const_iterator it = currentZmatrix_->begin(); it != currentZmatrix_->end(); ++it)
    {
      if (it->AtJ() == ax && it->AtK() == ay) {
        iTorsions.push_back( *it );
      }
    }
  }
  return iTorsions;
}


/** Assign torsions around bonded atoms in manner similar to LEaP's ModelAssignTorsionsAround. */
int Builder::assignTorsionsAroundBond(Zmatrix& zmatrix, int a1, int a2, Frame const& frameIn, Topology const& topIn, Barray const& hasPosition)
{
  // Save addresses of zmatrix, frame, topology, and hasPosition.
  // These are required for the createSpXSpX routines. TODO zero them at the end?
  newZmatrix_ = &zmatrix;
  currentFrm_ = &frameIn;
  currentTop_ = &topIn;
  hasPosition_ = &hasPosition;
  // No need to do this if either atom only has 1 bond.
  if (topIn[a1].Nbonds() < 2 || topIn[a2].Nbonds() < 2)
    return 0;
  // Get atom hybridizations
  AtomType::HybridizationType H1 = AtomType::UNKNOWN_HYBRIDIZATION;
  AtomType::HybridizationType H2 = AtomType::UNKNOWN_HYBRIDIZATION;
  if (params_ != 0) {
    ParmHolder<AtomType>::const_iterator it;
    if (topIn[a1].HasType()) {
      it = params_->AT().GetParam( TypeNameHolder(topIn[a1].Type()) );
      if (it != params_->AT().end())
        H1 = it->second.Hybridization();
    }
    if (topIn[a2].HasType()) {
      it = params_->AT().GetParam( TypeNameHolder(topIn[a2].Type()) );
      if (it != params_->AT().end())
        H2 = it->second.Hybridization();
    }
  }
  if (H1 == AtomType::UNKNOWN_HYBRIDIZATION)
    H1 = GuessAtomHybridization(topIn[a1], topIn.Atoms());
  if (H2 == AtomType::UNKNOWN_HYBRIDIZATION)
    H2 = GuessAtomHybridization(topIn[a2], topIn.Atoms());
  if (H1 == AtomType::UNKNOWN_HYBRIDIZATION)
    mprintf("Warning: No hybridization set for atom %s\n", topIn.AtomMaskName(a1).c_str());
  if (H2 == AtomType::UNKNOWN_HYBRIDIZATION)
    mprintf("Warning: No hybridization set for atom %s\n", topIn.AtomMaskName(a2).c_str());
  // Ensure the hybridization of ax is > ay
  int ax, ay;
  AtomType::HybridizationType Hx, Hy; // DEBUG
  if (H1 == AtomType::UNKNOWN_HYBRIDIZATION || H2 == AtomType::UNKNOWN_HYBRIDIZATION) {
    // Do not try to sort.
    ax = a1;
    ay = a2;
    Hx = H1;
    Hy = H2;
  } else if (H1 < H2) {
    ax = a2;
    ay = a1;
    Hx = H2;
    Hy = H1;
  } else {
    ax = a1;
    ay = a2;
    Hx = H1;
    Hy = H2;
  }
  static const char* hstr[] = { "SP", "SP2", "SP3", "Unknown" };
  mprintf("DEBUG: assignTorsionsAroundBond: AX= %s (%s)  AY= %s (%s)\n",
          topIn.AtomMaskName(ax).c_str(), hstr[Hx],
          topIn.AtomMaskName(ay).c_str(), hstr[Hy]);
  TorsionModel mT;
  if (mT.InitTorsion( ax, ay, frameIn, topIn, hasPosition )) {
    mprinterr("Error: Could not init model torsion.\n");
    return 1;
  }
  // Check if there is at least one atom on either side of the ax-ay pair
  // whose position is known.
  //Atom const& AX = topIn[ax];
  //Atom const& AY = topIn[ay];
  mprintf("bKnownX=%i  bKnownY=%i\n", (int)mT.AxHasKnownAtoms(), (int)mT.AyHasKnownAtoms());
  if (!(mT.AxHasKnownAtoms() && mT.AyHasKnownAtoms())) {
    // Find any existing internal coords around ax-ay
    std::vector<InternalCoords> iTorsions = getExistingInternals(ax, ay);
    if (!iTorsions.empty()) {
      mprintf("Using INTERNALs to fit new torsions around: %s - %s\n",
              topIn.LeapName(ax).c_str(), topIn.LeapName(ay).c_str());
      //for (std::vector<InternalCoords>::const_iterator ic = iTorsions.begin(); ic != iTorsions.end(); ++ic)
      //  ic->printIC( topIn );
      if (mT.BuildMockExternals(iTorsions, topIn)) {
        mprinterr("Error: Building mock externals around %s - %s failed.\n",
                  topIn.AtomMaskName(ax).c_str(), topIn.AtomMaskName(ay).c_str());
        return 1;
      }
      //return 0; // FIXME DEBUG 
    } else {
      mprintf("Completely free in assigning new torsions for: %s - %s\n",
              topIn.LeapName(ax).c_str(), topIn.LeapName(ay).c_str());
      return 0; // FIXME DEBUG 
    }
  } else {
    // Use existing atoms to determine torsions
    mprintf("DEBUG: Using externals to fit new torsions around: %s - %s\n",
            topIn.LeapName(ax).c_str(),
            topIn.LeapName(ay).c_str());
  }

  if (mT.SetupTorsion(Hx, Hy, topIn)) {
    mprinterr("Error: Could not set up torsions around %s - %s\n",
              topIn.LeapName(ax).c_str(),
              topIn.LeapName(ay).c_str());
    return 1;
  } 

  // Build the new internals
  if (Hx == AtomType::SP3 && Hy == AtomType::SP3) {
    mprintf("SP3 SP3\n");
    createSp3Sp3Torsions(mT);
  } else if (Hx == AtomType::SP3 && Hy == AtomType::SP2) {
    mprintf("SP3 SP2\n");
    createSp3Sp2Torsions(mT);
  } else if (Hx == AtomType::SP2 && Hy == AtomType::SP2) {
    mprintf("SP2 SP2\n");
    createSp2Sp2Torsions(mT);
  } else {
    mprinterr("Error: Currently only Sp3-Sp3/Sp3-Sp2/Sp2-Sp2 are supported\n"
              "Error: ---Tried to superimpose torsions for: *-%s-%s-*\n"
              "Error: --- With %s - %s\n"
              "Error: --- Sp0 probably means a new atom type is involved\n"
              "Error: --- which needs to be defined prior to this routine.\n",
              topIn.AtomMaskName(ax).c_str(), topIn.AtomMaskName(ay).c_str(),
              hstr[Hx], hstr[Hy]);
   return 1;
  }

  return 0;
}

/** Generate internal coordinates in the same
  * manner as LEaP's BuildInternalsForContainer/ModelAssignTorsionsAround.
  */
int Builder::GenerateInternals(Zmatrix& zmatrix, Frame const& frameIn, Topology const& topIn, Barray const& hasPosition)
{
  mprintf("DEBUG: ----- Entering Builder::GenerateInternals. -----\n");
  zmatrix.clear();
  // First generate the bond array
  BondArray bonds = GenerateBondArray( topIn.Residues(), topIn.Atoms() );
  // Loop over bonds
  for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
  {
    if (assignTorsionsAroundBond( zmatrix, bnd->A1(), bnd->A2(), frameIn, topIn, hasPosition )) {
      mprinterr("Error Assign torsions around bond %s - %s failed.\n",
                topIn.AtomMaskName(bnd->A1()).c_str(),
                topIn.AtomMaskName(bnd->A2()).c_str());
      return 1;
    }
  }
  //zmatrix.print( &topIn );
  mprintf("DEBUG: ----- Leaving Builder::GenerateInternals. ------\n");
  return 0;
}

/** Build internal coordinates around an atom. */
int Builder::generateAtomInternals(int at, Frame const& frameIn, Topology const& topIn, Barray const& hasPosition)
{
  mprintf( "Building internals for: %s\n", topIn.LeapName(at).c_str());
  Atom const& AtA = topIn[at];
  for (Atom::bond_iterator bat = AtA.bondbegin(); bat != AtA.bondend(); ++bat) {
    Atom const& AtB = topIn[*bat];
    for (Atom::bond_iterator cat = AtB.bondbegin(); cat != AtB.bondend(); ++cat) {
      if (*cat != at) {
        Atom const& AtC = topIn[*cat];
        mprintf("Building torsion INTERNALs for: %s  around: %s - %s\n",
               topIn.LeapName(at).c_str(),
               topIn.LeapName(*bat).c_str(),
               topIn.LeapName(*cat).c_str());
        std::vector<InternalCoords> iTorsions = getExistingInternals(*bat, *cat);
        int iShouldBe = (AtB.Nbonds() - 1) * (AtC.Nbonds() - 1);
        mprintf("ISHOULDBE= %i ITORSIONS= %zu\n", iShouldBe, iTorsions.size());
        if (iShouldBe != (int)iTorsions.size()) {
          Zmatrix tmpz; // FIXME
          assignTorsionsAroundBond(tmpz, *bat, *cat, frameIn, topIn, hasPosition);
        }
      }
    }
  }

  return 0;
}
 
/** Generate internal coordinates around a bond linking two residues
  * in the same manner as LEaP.
  * \param zmatrix Hold output ICs
  * \param at0 Atom in residue we are linking to (i.e. the current residue).
  * \param at1 Atom in resiude we are linking from.
  */
int Builder::GenerateInternalsAroundLink(Zmatrix& zmatrix, int at0, int at1,
                                         Frame const& frameIn, Topology const& topIn,
                                         Barray const& hasPosition)
{
  mprintf("DEBUG: ----- Entering Builder::GenerateInternalsAroundLink. -----\n");
  // Sanity check
  Atom const& A0 = topIn[at0];
  Atom const& A1 = topIn[at1];
  if (A0.ResNum() == A1.ResNum()) {
    mprinterr("Internal Error: Builder::GenerateInternalsAroundLink(): Atoms are in the same residue.\n");
    return 1;
  }
  // In order to mimic the way LEaP does things, mark all atoms before
  // this residue as having position, the first known atom of this residue
  // has having position, and all other atoms as not having position.
  Residue const& R0 = topIn.Res(A0.ResNum());
  Barray tmpHasPosition( topIn.Natom(), false );
  for (int at = 0; at < R0.FirstAtom(); at++)
    tmpHasPosition[at] = true;
  for (int at = R0.FirstAtom(); at != R0.LastAtom(); at++) {
    if (hasPosition[at]) {
      tmpHasPosition[at] = true;
      break;
    }
  }
  // Create spanning tree across the link
  std::vector<int> span_atoms = GenerateSpanningTree(at0, at1, 4, topIn.Atoms());
  for (std::vector<int>::const_iterator it = span_atoms.begin(); 
                                        it != span_atoms.end(); ++it)
  {
    //mprintf("SPANNING TREE ATOM: %s\n", topIn.LeapName(*it).c_str());
    if (generateAtomInternals(*it, frameIn, topIn, tmpHasPosition)) {
      mprinterr("Error: Could not generate internals for atom %s\n", topIn.AtomMaskName(*it).c_str());
      return 1;
    }
  }
  // Create torsions around the link
  //if (assignTorsionsAroundBond( zmatrix, at0, at1, frameIn, topIn, hasPosition )) {
  //  mprinterr("Error Assign torsions around link %s - %s failed.\n",
  //            topIn.AtomMaskName(at0).c_str(),
  //            topIn.AtomMaskName(at1).c_str());
  //  return 1;
  //}
  mprintf("DEBUG: ----- Leaving Builder::GenerateInternalsAroundLink. -----\n");
  return 0;
}
