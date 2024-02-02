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
  currentZmatrix_(0)
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
/** Store info for modelling torsions around X-Y */
class Cpptraj::Structure::Builder::ModelTorsion {
  public:
    /// CONSTRUCTOR
    ModelTorsion() : ax_(-1), ay_(-1), dAbsolute_(0), Xorientation_(0), Yorientation_(0) {}
    /// Set up torsions around bonded atoms
    int SetupTorsion(int, int, AtomType::HybridizationType, AtomType::HybridizationType,
                     Frame const&, Topology const&, std::vector<bool> const&);

  private:
    static int LeapAtomWeight(Atom const&);
    static inline std::vector<int> SiftBondedAtomsLikeLeap(unsigned int&, Atom const&, std::vector<bool> const&);
    static inline std::vector<int> SortBondedAtomsLikeLeap(unsigned int&, Atom const&,
                                                           Topology const& topIn, int ignoreAtom,
                                                           std::vector<bool> const& hasPosition);

    int ax_;              ///< Atom X
    int ay_;              ///< Atom Y
    Iarray sorted_ax_;    ///< Hold the leap-sorted bonds for atom X
    Iarray sorted_ay_;    ///< Hold the leap-sorted bonds for atom Y
    double dAbsolute_;    ///< Hold the value of the A-X-Y-D torsion in radians
    double Xorientation_; ///< Orientation around the X atom TODO make an enum?
    double Yorientation_; ///< Orientation around the Y atoms
};
    
/** \return the LEaP 'weight' of an atom.
  * Originally used to force the 'heaviest' atoms around a torsion trans to
  * each other. The 'weight' of an atom is defined as its element number,
  * unless the atom is CARBON, then it is 1000, making it the 'heaviest' atom.
  */
int Cpptraj::Structure::Builder::ModelTorsion::LeapAtomWeight(Atom const& At)
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
std::vector<int> 
  Cpptraj::Structure::Builder::ModelTorsion::SiftBondedAtomsLikeLeap(unsigned int& firstUnknownIdx,
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
}

/** Order atoms bonded to the given atom in a manner similar to LEaP's
  * zModelOrderAtoms. In that routine, first atoms were sorted into
  * known position > unknown position. Then the heaviest atom in each
  * subgroup was swapped with the first element of that list. 
  * The ignore atom is the index of the atom this atom is bonded to that
  * forms the torsion we are interested in.
  */
std::vector<int>
  Cpptraj::Structure::Builder::ModelTorsion::SortBondedAtomsLikeLeap(unsigned int& firstUnknownIdx,
                                                                     Atom const& At, Topology const& topIn,
                                                                     int ignoreAtom,
                                                                     std::vector<bool> const& hasPosition)
{
  std::vector<int> out;
  out.reserve( At.Nbonds() );
  // Sift so that atoms with known position are at the front
  std::vector<int> bondedAtoms = SiftBondedAtomsLikeLeap(firstUnknownIdx, At, hasPosition);
  // Find the index of the heaviest atom
  int iHighest = 0;
  int iPos = 0;
  for (unsigned int idx = 0; idx < bondedAtoms.size(); idx++) {
    int bat = bondedAtoms[idx];
    if (bat != ignoreAtom) {
      out.push_back( bat );
      int iWeight = LeapAtomWeight( topIn[bat] );
      if ( iHighest < iWeight ) {
        iHighest = iWeight;
        iPos = (int)out.size()-1;
      }
    }
  }
  // If highest weight atom not already in front, swap it there.
  if (iPos != 0) std::swap( out[0], out[iPos] );

  return out;
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
static inline double calculateOrientation(int iX, int iA, int iY, int iB, Frame const& frameIn, std::vector<bool> const& hasPosition)
{
  double dOrientation = 1.0;
  if (hasPosition[iX] &&
      hasPosition[iA] &&
      hasPosition[iY] &&
      hasPosition[iB])
  {
    dOrientation = VectorAtomChirality( frameIn.XYZ(iX), frameIn.XYZ(iA), frameIn.XYZ(iY), frameIn.XYZ(iB) );
  } else {
    mprinterr("Internal Error: Builder::calculateOrientation not yet set up for unknown positions.\n");
  }
  return dOrientation;
}

/** Set up model torsion for bonded atoms. */
int Cpptraj::Structure::Builder::ModelTorsion::SetupTorsion(int ax, int ay,
                                                            AtomType::HybridizationType Hx,
                                                            AtomType::HybridizationType Hy,
                                                            Frame const& frameIn,
                                                            Topology const& topIn,
                                                            std::vector<bool> const& hasPosition)
{
  if (Hx != AtomType::UNKNOWN_HYBRIDIZATION && Hy != AtomType::UNKNOWN_HYBRIDIZATION) {
    if (Hy > Hx) {
      mprinterr("Internal Error: :ModelTorsion::SetupTorsion() called with AX hybrid > AY hybrid.\n");
      return 1;
    }
  }
  ax_ = ax;
  ay_ = ay;
  Atom const& AX = topIn[ax];
  Atom const& AY = topIn[ay];
  // Sort AX bonds
  unsigned int firstUnknownIdxX = 0;
  sorted_ax_ = SortBondedAtomsLikeLeap(firstUnknownIdxX, AX, topIn, ay, hasPosition);
  // Sort AY bonds
  unsigned int firstUnknownIdxY = 0;
  sorted_ay_ = SortBondedAtomsLikeLeap(firstUnknownIdxY, AY, topIn, ax, hasPosition);
  // Calculate the chirality around atom X
  Xorientation_ = 0;
  if (Hx == AtomType::SP3) {
    Xorientation_ = calculateOrientation( ax, sorted_ax_[0], ay, sorted_ax_[1], frameIn, hasPosition );
  }
  // Calculate the chirality around atom Y
  Yorientation_ = 0;
  if (Hy == AtomType::SP3) {
    Yorientation_ = calculateOrientation( ay, sorted_ay_[0], ax, sorted_ay_[1], frameIn, hasPosition );
  }
  // DEBUG
  mprintf("Orientation around: %s = %f\n", *(AX.Name()), Xorientation_);
  //for (Atom::bond_iterator bat = AX.bondbegin(); bat != AX.bondend(); ++bat) mprintf(" %s", *(topIn[*bat].Name()));
  //mprintf("}\n");
  for (Iarray::const_iterator it = sorted_ax_.begin(); it != sorted_ax_.end(); ++it)
      mprintf("Atom %li: %s\n", it - sorted_ax_.begin(), *(topIn[*it].Name()));
  mprintf("Orientation around: %s = %f\n", *(AY.Name()), Yorientation_);
  //for (Atom::bond_iterator bat = AY.bondbegin(); bat != AY.bondend(); ++bat) mprintf(" %s", *(topIn[*bat].Name()));
  //mprintf("}\n");
  for (Iarray::const_iterator it = sorted_ay_.begin(); it != sorted_ay_.end(); ++it)
      mprintf("Atom %li: %s\n", it - sorted_ay_.begin(), *(topIn[*it].Name()));
  // Calculate the actual torsion angle between A-X-Y-D
  if (hasPosition[sorted_ax_[0]] &&
      hasPosition[ax] &&
      hasPosition[ay] &&
      hasPosition[sorted_ay_[0]])
  {
    dAbsolute_ = Torsion( frameIn.XYZ(sorted_ax_[0]),
                          frameIn.XYZ(ax),
                          frameIn.XYZ(ay),
                          frameIn.XYZ(sorted_ay_[0]) );
  } else {
    dAbsolute_ = 180.0 * Constants::DEGRAD;
  }
  mprintf("DABSOLUTE= %g\n", dAbsolute_);
  return 0;
}

// -----------------------------------------------
/** Create torsions around SP3-SP3. */
void Builder::createSp3Sp3Torsions() {
  return;
}

void Builder::createSp3Sp2Torsions() {
  return;
}

void Builder::createSp2Sp2Torsions() {
  return;
}

/** Assign torsions around bonded atoms in manner similar to LEaP's ModelAssignTorsionsAround. */
int Builder::assignTorsionsAroundBond(int a1, int a2, Frame const& frameIn, Topology const& topIn, Barray const& hasPosition)
{
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
  // Check if there is at least one atom on either side of the ax-ay pair
  // whose position is known.
  Atom const& AX = topIn[ax];
  Atom const& AY = topIn[ay];
  bool axHasKnownAtoms = false;
  for (Atom::bond_iterator bat = AX.bondbegin(); bat != AX.bondend(); ++bat) {
    if (*bat != ay && hasPosition[*bat]) {
      axHasKnownAtoms = true;
      break;
    }
  }
  bool ayHasKnownAtoms = false;
  for (Atom::bond_iterator bat = AY.bondbegin(); bat != AY.bondend(); ++bat) {
    if (*bat != ax && hasPosition[*bat]) {
      ayHasKnownAtoms = true;
      break;
    }
  }

  if (!axHasKnownAtoms && !ayHasKnownAtoms) {
    mprinterr("Internal Error: assignTorsionsAroundBond both not known not yet implemented.\n");
    return 1;
  } else {
    mprintf("DEBUG: Using externals to fit new torsions around: %s - %s\n",
            topIn.LeapName(ax).c_str(),
            topIn.LeapName(ay).c_str());

    ModelTorsion mT;
    if (mT.SetupTorsion(ax, ay, Hx, Hy, frameIn, topIn, hasPosition)) {
      mprinterr("Error: Could not set up torsions around %s - %s\n",
                topIn.LeapName(ax).c_str(),
                topIn.LeapName(ay).c_str());
      return 1;
    } 

    // Build the new internals
    if (Hx == AtomType::SP3 && Hy == AtomType::SP3) {
      mprintf("SP3 SP3\n");
      createSp3Sp3Torsions();
    } else if (Hx == AtomType::SP3 && Hy == AtomType::SP2) {
      mprintf("SP3 SP2\n");
      createSp3Sp2Torsions();
    } else if (Hx == AtomType::SP2 && Hy == AtomType::SP2) {
      mprintf("SP2 SP2\n");
      createSp2Sp2Torsions();
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
    if (assignTorsionsAroundBond( bnd->A1(), bnd->A2(), frameIn, topIn, hasPosition )) {
      mprinterr("Error Assign torsions around bond %s - %s failed.\n",
                topIn.AtomMaskName(bnd->A1()).c_str(),
                topIn.AtomMaskName(bnd->A2()).c_str());
      return 1;
    }
/*
    Atom const& A2 = topIn[bnd->A1()];
    Atom const& A3 = topIn[bnd->A2()];
    if (A2.Nbonds() > 1 && A3.Nbonds() > 1) {
      //Residue const& R2 = topIn.Res(A2.ResNum());
      //Residue const& R3 = topIn.Res(A3.ResNum());
      mprintf("Building torsion INTERNALs around: %s - %s\n",
              topIn.LeapName(bnd->A1()).c_str(), topIn.LeapName(bnd->A2()).c_str());
      Iarray sorted_a2 = SortBondedAtomsLikeLeap(A2, topIn, bnd->A2());
      Iarray sorted_a3 = SortBondedAtomsLikeLeap(A3, topIn, bnd->A1());
      mprintf("Orientation around: %s {", *(A2.Name()));
      for (Atom::bond_iterator bat = A2.bondbegin(); bat != A2.bondend(); ++bat) mprintf(" %s", *(topIn[*bat].Name()));
      mprintf("}\n");
      for (Iarray::const_iterator it = sorted_a2.begin(); it != sorted_a2.end(); ++it)
        mprintf("Atom %li: %s\n", it - sorted_a2.begin(), *(topIn[*it].Name()));
      mprintf("Orientation around: %s {", *(A3.Name()));
      for (Atom::bond_iterator bat = A3.bondbegin(); bat != A3.bondend(); ++bat) mprintf(" %s", *(topIn[*bat].Name()));
      mprintf("}\n");
      for (Iarray::const_iterator it = sorted_a3.begin(); it != sorted_a3.end(); ++it)
        mprintf("Atom %li: %s\n", it - sorted_a3.begin(), *(topIn[*it].Name()));
      // Build the torsions
      int aj = bnd->A1();
      int ak = bnd->A2();
      for (Iarray::const_iterator ai = sorted_a2.begin(); ai != sorted_a2.end(); ++ai) {
        for (Iarray::const_iterator al = sorted_a3.begin(); al != sorted_a3.end(); ++al) {
          //double dval = Torsion(frameIn.XYZ(*ai), frameIn.XYZ(aj), frameIn.XYZ(ak), frameIn.XYZ(*al));
          zmatrix.AddIC(*ai, aj, ak, *al, frameIn);
          mprintf("++++Torsion INTERNAL: %f to %s - %s - %s - %s\n",
                  zmatrix.back().Phi(),
                  topIn.LeapName(*ai).c_str(),
                  topIn.LeapName(aj).c_str(),
                  topIn.LeapName(ak).c_str(),
                  topIn.LeapName(*al).c_str());
          zmatrix.AddIC(*al, ak, aj, *ai, frameIn);
        }
      }
    }*/
  }
  mprintf("DEBUG: ----- Leaving Builder::GenerateInternals. ------\n");
  return 0;
}
 
