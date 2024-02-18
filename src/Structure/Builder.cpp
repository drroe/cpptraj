#include "Builder.h"
#include "Chirality.h"
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
#include <cmath> // fabs

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Builder::Builder() :
  debug_(0),
  params_(0),
  currentTop_(0)
{}

/** Set optional parameter set. */
void Cpptraj::Structure::Builder::SetParameters(ParameterSet const* paramsIn) {
  if (paramsIn == 0) {
    mprinterr("Internal Error: Builder::SetParmaters called with null set.\n");
    return;
  }
  params_ = paramsIn;
}

/** Set an atoms chirality */
void Cpptraj::Structure::Builder::SetAtomChirality(int at, double chi) {
  int cidx = getExistingChiralityIdx(at);
  if (cidx == -1) {
    internalChirality_.push_back( InternalChirality(at, chi) );
  } else {
    mprintf("Warning: Overriding existing chirality %f for atom %i with %f\n",
            internalChirality_[cidx].ChiralVal(), at+1, chi);
    internalChirality_[cidx].SetChiralVal( chi );
  }
}

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

/** \return Index of atom with longest 'depth' around an atom, ignoring one bonded atom. */
int Builder::get_depths_around_atom(int at0, int at1, Topology const& topIn) {
  int largest_depth = 0;
  int largest_depth_idx = -1;
  Atom const& Atm0 = topIn[at0];
  mprintf("DEBUG: Depths around %s\n", topIn.AtomMaskName(at0).c_str());
  for (Atom::bond_iterator bat = Atm0.bondbegin(); bat != Atm0.bondend(); ++bat)
  {
    if (*bat != at1) {
      Barray visited(topIn.Natom(), false);
      visited[at0] = true;
      visited[at1] = true;
      // Get depth from here
      int currentDepth = 0;
      int depth = atom_depth(currentDepth, *bat, topIn, visited, 10);
      mprintf("DEBUG:\t\t%s = %i\n", topIn.AtomMaskName(*bat).c_str(), depth);
      if (depth > largest_depth) {
        largest_depth = depth;
        largest_depth_idx = *bat;
      }
    }
  }
  return largest_depth_idx;
}

/** Adjust internals around a link so that atoms with longest 'depth' are trans. */
int Builder::AdjustIcAroundLink(int at0, int at1, Frame const& frameIn, Topology const& topIn)
{
  int atA = get_depths_around_atom(at0, at1, topIn);
  int atD = get_depths_around_atom(at1, at0, topIn);
  if (atA < 0 || atD < 0) return 0;
  mprintf("DEBUG: Highest depth torsion is %s - %s - %s - %s\n",
          topIn.AtomMaskName(atA).c_str(),
          topIn.AtomMaskName(at0).c_str(),
          topIn.AtomMaskName(at1).c_str(),
          topIn.AtomMaskName(atD).c_str());
  // Get that torsion
  int tidx = getExistingTorsionIdx(atA, at0, at1, atD);
  if (tidx < 0) {
    mprintf("Warning: Could not adjust internal torsion %s - %s - %s - %s, not present.\n",
          topIn.AtomMaskName(atA).c_str(),
          topIn.AtomMaskName(at0).c_str(),
          topIn.AtomMaskName(at1).c_str(),
          topIn.AtomMaskName(atD).c_str());
    return 0;
  }
  // Figure out the delta from 180
  double dTorsion = 180.0 * Constants::DEGRAD;
  double dInternalValue = internalTorsions_[tidx].PhiVal();
  double tDiff = (dTorsion - dInternalValue);
  mprintf("\tdTorsion= %f  dInternalValue= %f  tDiff= %f\n",
          dTorsion*Constants::RADDEG, dInternalValue*Constants::RADDEG, tDiff*Constants::RADDEG);
  if (fabs(tDiff) > Constants::SMALL) {
    // Find all ICs that share atoms 1 and 2 (J and K)
    Iarray iTorsions = getExistingTorsionIdxs(at0, at1);
    mprintf("Twisting torsions centered on %s - %s by %f degrees\n",
            topIn.LeapName(at0).c_str(),
            topIn.LeapName(at1).c_str(),
            tDiff*Constants::RADDEG);
    for (Iarray::const_iterator idx = iTorsions.begin(); idx != iTorsions.end(); ++idx)
    {
      InternalTorsion& dih = internalTorsions_[*idx];
      double dNew = dih.PhiVal() + tDiff;
      mprintf("Twisting torsion for atoms: %s-%s-%s-%s\n",
              *(topIn[dih.AtI()].Name()),
              *(topIn[dih.AtJ()].Name()),
              *(topIn[dih.AtK()].Name()),
              *(topIn[dih.AtL()].Name()));
      mprintf("------- From %f to %f\n", dih.PhiVal()*Constants::RADDEG, dNew*Constants::RADDEG);
      dih.SetPhiVal( dNew );
    }
  }

  return 0;
}

/** For existing torsions, see if all coordinates in that torsion
  * exist. If so, update the torsion from the existing coordinates.
  */
int Builder::UpdateICsFromFrame(Frame const& frameIn, Topology const& topIn, Barray const& hasPosition)
{
  // Create a list of residues that have atoms with internals 
  std::vector<Residue> residues;
  std::vector<int> Rnums;
  for (Tarray::const_iterator dih = internalTorsions_.begin();
                              dih != internalTorsions_.end(); ++dih)
  {
    int rnum = topIn[dih->AtI()].ResNum();
    bool has_rnum = false;
    for (std::vector<int>::const_iterator it = Rnums.begin(); it != Rnums.end(); ++it) {
      if (*it == rnum) {
        has_rnum = true;
        break;
      }
    }
    if (!has_rnum) {
      mprintf("DEBUG: Residue %s has internals.\n", topIn.TruncResNameNum(rnum).c_str());
      residues.push_back( topIn.Res(rnum) );
      Rnums.push_back(rnum);
    }
  }
  // Get list of bonds.
  BondArray myBonds = GenerateBondArray( residues, topIn.Atoms() );
  for (BondArray::const_iterator bnd = myBonds.begin(); bnd != myBonds.end(); ++bnd) {
    //if (topIn[bnd->A1()].ResNum() == ires && topIn[bnd->A2()].ResNum() == ires) {
      mprintf("Looking at torsions around: %s - %s\n", topIn.LeapName(bnd->A1()).c_str(), topIn.LeapName(bnd->A2()).c_str());
      // Find all ICs that share atoms 1 and 2 (J and K)
      bool needsUpdate = false;
      double tDiff = 0;
      Iarray iTorsions = getExistingTorsionIdxs(bnd->A1(), bnd->A2());
      for (Iarray::const_iterator idx = iTorsions.begin(); idx != iTorsions.end(); ++idx)
      {
        InternalTorsion const& dih = internalTorsions_[*idx];
        if (hasPosition[dih.AtI()] &&
            hasPosition[dih.AtJ()] &&
            hasPosition[dih.AtK()] &&
            hasPosition[dih.AtL()])
        {
          mprintf("Measuring torsion of fixed atoms: %s - %s - %s - %s\n",
                  topIn.LeapName(dih.AtI()).c_str(),
                  topIn.LeapName(dih.AtJ()).c_str(),
                  topIn.LeapName(dih.AtK()).c_str(),
                  topIn.LeapName(dih.AtL()).c_str());
          double dTorsion = Torsion(frameIn.XYZ(dih.AtI()),
                                    frameIn.XYZ(dih.AtJ()),
                                    frameIn.XYZ(dih.AtK()),
                                    frameIn.XYZ(dih.AtL()));
          double dInternalValue = dih.PhiVal();
          tDiff = (dTorsion - dInternalValue);
          mprintf("\tdTorsion= %f  dInternalValue= %f\n", dTorsion, dInternalValue);
          needsUpdate = true;
        } // END all coords present
      } // END loop over torsions matching current bond
      // If any difference was found, shift all of the torsions
      if (needsUpdate) {
        mprintf("Twisting torsions centered on %s - %s by %f degrees\n",
                topIn.LeapName(bnd->A1()).c_str(),
                topIn.LeapName(bnd->A2()).c_str(),
                tDiff*Constants::RADDEG);
        for (Iarray::const_iterator idx = iTorsions.begin(); idx != iTorsions.end(); ++idx)
        {
          InternalTorsion& dih = internalTorsions_[*idx];
          double dNew = dih.PhiVal() + tDiff;
          mprintf("Twisting torsion for atoms: %s-%s-%s-%s\n",
                  *(topIn[dih.AtI()].Name()),
                  *(topIn[dih.AtJ()].Name()),
                  *(topIn[dih.AtK()].Name()),
                  *(topIn[dih.AtL()].Name()));
          mprintf("------- From %f to %f\n", dih.PhiVal()*Constants::RADDEG, dNew*Constants::RADDEG);
          dih.SetPhiVal( dNew );
        }
      } // END ICs need update
    //} // END both bond atoms belong to this residue
  } // END loop over bonds
  return 0;
}

// -----------------------------------------------------------------------------
/// Used to track atoms for mock externals
class MockAtom {
  public:
    /// CONSTRUCTOR
    MockAtom() : idx_(-1), pos_(0.0), known_(false), buildInternals_(false) {}
    /// CONSTRUCTOR - index, position
    MockAtom(int idx, Vec3 const& pos) : idx_(idx), pos_(pos), known_(true), buildInternals_(false) {}
    /// CONSTRUCTOR - index, unknown position
    MockAtom(int idx) : idx_(idx), pos_(0.0), known_(false), buildInternals_(false) {}
    /// Set position
    void SetPos(Vec3 const& p) { pos_ = p; known_ = true; }
    /// Set position status to unknown
    void SetUnknown() { pos_ = Vec3(0.0); known_ = false; }
    /// Set build internals status
    void SetBuildInternals(bool b) { buildInternals_ = b; }

    int Idx()         const { return idx_; }
    Vec3 const& Pos() const { return pos_; }
    bool Known()      const { return known_; }
    bool BuildInternals() const { return buildInternals_; }
  private:
    int idx_;    ///< Atom index
    Vec3 pos_;   ///< Atom position
    bool known_; ///< True if atom position is known
    bool buildInternals_; ///< True if internals should be built for this atom
};

// -----------------------------------------------------------------------------
/** Store info for modelling torsions around X-Y */
class Cpptraj::Structure::Builder::TorsionModel {
  public:
    typedef std::vector<MockAtom> Marray;
    /// CONSTRUCTOR
    TorsionModel() : dAbsolute_(0), Xorientation_(0), Yorientation_(0), axHasKnownAtoms_(false), ayHasKnownAtoms_(false) {}
    /// Initialize torsions around bonded atoms
    int InitTorsion(int, int, Frame const&, Topology const&, std::vector<bool> const&, int);
    /// Set up torsions around bonded atoms
    int SetupTorsion(AtomType::HybridizationType, AtomType::HybridizationType, Topology const& topIn, double, double);
    /// Build mock externals from given internals
    int BuildMockExternals(Iarray const&, Tarray const&, Topology const&);

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
    /// Sort an array of MockAtoms the way that LEaP does
    static Marray SortBondedAtomsLikeLeap(unsigned int&, Topology const& topIn, Marray const&);
  private:
    static int LeapAtomWeight(Atom const&);
    //static inline std::vector<int> SiftBondedAtomsLikeLeap(unsigned int&, Atom const&, std::vector<bool> const&);
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

/** Assuming atoms have been ordered with SortBondedAtomsLikeLeap,
  * calculate the orientation of the iB atom with respect to the
  * triangle (iA, iX, iY). This orientation will be used by the
  * CreateSpXSpX routines to determine which torsion values to use.
  */
static inline double calculateOrientation(MockAtom const& iX, double chiX, Atom const& AX, MockAtom const& iA, MockAtom const& iY, MockAtom const& iB)
{
  using namespace Cpptraj::Structure::Chirality;
  double dOrientation = 1.0;
  if (iX.Known() &&
      iA.Known() &&
      iY.Known() &&
      iB.Known())
  {
    dOrientation = VectorAtomChirality( iX.Pos(), iA.Pos(), iY.Pos(), iB.Pos() );
    mprintf("ORIENTATION: known = %f\n", dOrientation);
  } else {
    double dChi = chiX;
    if (fabs(dChi) < Constants::SMALL) {
      mprintf("default chirality on\n");
      dChi = 1.0;
    }
    mprintf("ORIENTATION: Chirality %f\n", dChi);
    dOrientation = chiralityToOrientation(dChi, AX, iA.Idx(), iY.Idx(), iB.Idx(), -1);
  }
  return dOrientation;
}

/** Initialize model torsion for bonded atoms. */
int Cpptraj::Structure::Builder::TorsionModel::InitTorsion(int ax, int ay,
                                                           Frame const& frameIn,
                                                           Topology const& topIn,
                                                           std::vector<bool> const& hasPosition,
                                                           int aAtomIdx)
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
      sorted_ax_.back().SetBuildInternals( (aAtomIdx == -1 || aAtomIdx == *bat) );
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
      sorted_ay_.back().SetBuildInternals( (aAtomIdx == -1 || aAtomIdx == *bat) );
    }
  }
  return 0;
}

/** Set up model torsion for bonded atoms. */
int Cpptraj::Structure::Builder::TorsionModel::SetupTorsion(AtomType::HybridizationType Hx,
                                                            AtomType::HybridizationType Hy,
                                                            Topology const& topIn,
                                                            double chiX, double chiY)
//                                                            double orientX, double orientY)
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
    if (sorted_ax_.size() < 2)
      Xorientation_ = 1.0;
    else
      Xorientation_ = calculateOrientation( atX_, chiX, topIn[atX_.Idx()], sorted_ax_[0], atY_, sorted_ax_[1] );
  }
  // Calculate the chirality around atom Y
  Yorientation_ = 0;
  if (Hy == AtomType::SP3) {
    if (sorted_ay_.size() < 2)
      Yorientation_ = 1.0;
    else
      Yorientation_ = calculateOrientation( atY_, chiY, topIn[atY_.Idx()], sorted_ay_[0], atX_, sorted_ay_[1] );
  }
  // DEBUG
  Atom const& AX = topIn[atX_.Idx()];
  Atom const& AY = topIn[atY_.Idx()];
  //mprintf("Orientation around: %s = %f (chiX= %f)\n", *(AX.Name()), Xorientation_, chiX);
  mprintf("Orientation around: %s = %f\n", *(AX.Name()), Xorientation_);
  //for (Atom::bond_iterator bat = AX.bondbegin(); bat != AX.bondend(); ++bat) mprintf(" %s", *(topIn[*bat].Name()));
  //mprintf("}\n");
  for (Marray::const_iterator it = sorted_ax_.begin(); it != sorted_ax_.end(); ++it)
      mprintf("Atom %li: %s (%i) (build=%i)\n", it - sorted_ax_.begin(), *(topIn[it->Idx()].Name()), (int)it->Known(), (int)it->BuildInternals());
  //mprintf("Orientation around: %s = %f (chiY= %f)\n", *(AY.Name()), Yorientation_, chiY);
  mprintf("Orientation around: %s = %f\n", *(AY.Name()), Yorientation_);
  //for (Atom::bond_iterator bat = AY.bondbegin(); bat != AY.bondend(); ++bat) mprintf(" %s", *(topIn[*bat].Name()));
  //mprintf("}\n");
  for (Marray::const_iterator it = sorted_ay_.begin(); it != sorted_ay_.end(); ++it)
      mprintf("Atom %li: %s (%i) (build=%i)\n", it - sorted_ay_.begin(), *(topIn[it->Idx()].Name()), (int)it->Known(), (int)it->BuildInternals());
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
int Cpptraj::Structure::Builder::TorsionModel::BuildMockExternals(Iarray const& iaTorsions,
                                                                  Tarray const& internalTorsionsIn,
                                                                  Topology const& topIn) // DEBUG topIn for debug only
{
  if (iaTorsions.empty()) {
    mprinterr("Internal Error: TorsionModel::BuildMockExternals() called with no internal torsions.\n");
    return 1;
  }
  mprintf("=======  Started mock coords from: %s\n", topIn.LeapName(internalTorsionsIn[iaTorsions.front()].AtI()).c_str());
  mprintf("========  %zu Torsions to build mock coords from:\n", iaTorsions.size());
  for (Iarray::const_iterator idx = iaTorsions.begin(); idx != iaTorsions.end(); ++idx) {
    InternalTorsion const& ic = internalTorsionsIn[*idx];
    mprintf("------- Known torsion: %s - %s - %s - %s  %f\n",
            topIn.LeapName(ic.AtI()).c_str(),
            topIn.LeapName(ic.AtJ()).c_str(),
            topIn.LeapName(ic.AtK()).c_str(),
            topIn.LeapName(ic.AtL()).c_str(),
            ic.PhiVal()*Constants::RADDEG);
  }

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
  for (Iarray::const_iterator idx = iaTorsions.begin(); idx != iaTorsions.end(); ++idx)
  {
    InternalTorsion const& ic = internalTorsionsIn[*idx];
    if (idx == iaTorsions.begin()) {
      // Define first outer atom as being in the XY plane
      outerAtoms.push_back( MockAtom(ic.AtI(), Vec3(1, 1, 0)) );
    } else {
      Marray::iterator mi = find_mock_atom( outerAtoms, ic.AtI() );
      if (mi == outerAtoms.end())
        outerAtoms.push_back( MockAtom(ic.AtI()) );
    }
    Marray::iterator ml = find_mock_atom( outerAtoms, ic.AtL() );
    if (ml == outerAtoms.end())
      outerAtoms.push_back( MockAtom(ic.AtL()) );
  }
# ifdef CPPTRAJ_DEBUG_BUILDER
  mprintf("DEBUG: Outer atoms:\n");
  for (Marray::const_iterator it = outerAtoms.begin(); it != outerAtoms.end(); ++it)
    mprintf("DEBUG:\t\t%i %4s (%i) {%f %f %f}\n", it->Idx()+1, topIn.AtomMaskName(it->Idx()).c_str(),
            (int)it->Known(), it->Pos()[0], it->Pos()[1], it->Pos()[2]);
# endif
  // Loop through the known torsions looking for those that
  // have one position defined, then build coords for the
  // other atom and mark the torsion as used.
  std::vector<bool> used( iaTorsions.size(), false );
  unsigned int nused = 0;
  for (unsigned int idx = 0; idx != iaTorsions.size(); idx++) {
    //bool gotOne = false;
    for (unsigned int jdx = 0; jdx != iaTorsions.size(); jdx++) {
      if (!used[jdx]) {
        InternalTorsion const& iInt = internalTorsionsIn[iaTorsions[jdx]];
        Marray::iterator tmpAt1 = find_mock_atom(outerAtoms, iInt.AtI());
        Marray::iterator tmpAt4 = find_mock_atom(outerAtoms, iInt.AtL());
        if (tmpAt1 == outerAtoms.end()) {
          mprinterr("Internal Error: TorsionModel::BuildMockExternals(): Outer atom I %i not found.\n", iInt.AtI()+1);
          return 1;
        }
        if (tmpAt4 == outerAtoms.end()) {
          mprinterr("Internal Error: TorsionModel::BuildMockExternals(): Outer atom L %i not found.\n", iInt.AtL()+1);
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
          tgt->SetPos( Zmatrix::AtomIposition(maPC1, maPC2, knownAt->Pos(), 1.0, 90.0, iInt.PhiVal()*Constants::RADDEG) );
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
# ifdef CPPTRAJ_DEBUG_BUILDER
  mprintf("DEBUG: Final outer atoms:\n");
  for (Marray::const_iterator it = outerAtoms.begin(); it != outerAtoms.end(); ++it) {
    mprintf("DEBUG:\t\t%i %4s (%i) {%f %f %f}\n", it->Idx()+1, topIn.AtomMaskName(it->Idx()).c_str(),
            (int)it->Known(), it->Pos()[0], it->Pos()[1], it->Pos()[2]);
  }
# endif
  // Update the outer atom positions for this torsion
  for (Marray::const_iterator it = outerAtoms.begin(); it != outerAtoms.end(); ++it) {
    Marray::iterator itx = find_mock_atom(sorted_ax_, it->Idx());
    if (itx != sorted_ax_.end()) {
      itx->SetPos( it->Pos() );
    } else {
      itx = find_mock_atom(sorted_ay_, it->Idx());
      if (itx != sorted_ay_.end()) {
        itx->SetPos( it->Pos() );
      } else {
        mprinterr("Internal Error: TorsionModel::BuildMockExternals(): Could not update mock atom.\n");
        return 1;
      }
    }
  }

  return 0;
}

// -----------------------------------------------------------------------------
/** Determine the chirality around a single atom. Assume all positions are known. */
double Builder::DetermineChiralityAroundAtom(int at, Frame const& frameIn, Topology const& topIn)
{
  double dChi = 0;
  determineChirality(dChi, at, frameIn, topIn, std::vector<bool>(topIn.Natom(), true));
  return dChi;
}

/** Update all indices in internals according to the given offset. */
void Builder::UpdateIndicesWithOffset(int atomOffset) {
  for (Tarray::iterator it = internalTorsions_.begin(); it != internalTorsions_.end(); ++it)
    it->OffsetIndices( atomOffset );
  for (Aarray::iterator it = internalAngles_.begin(); it != internalAngles_.end(); ++it)
    it->OffsetIndices( atomOffset );
  for (Larray::iterator it = internalBonds_.begin(); it != internalBonds_.end(); ++it)
    it->OffsetIndices( atomOffset );
  for (Carray::iterator it = internalChirality_.begin(); it != internalChirality_.end(); ++it)
    it->OffsetIndices( atomOffset );
}

/** Find any existing torsions around ax-ay. */
Builder::Iarray Builder::getExistingTorsionIdxs(int ax, int ay) const {
  Iarray iTorsions;
  for (Tarray::const_iterator it = internalTorsions_.begin(); it != internalTorsions_.end(); ++it)
  {
    if ((it->AtJ() == ax && it->AtK() == ay) ||
        (it->AtJ() == ay && it->AtK() == ax))
    {
      iTorsions.push_back( it - internalTorsions_.begin() );
    }
  }
  return iTorsions;
}

/** \return index of existing torsion matching the 4 given atoms, -1 for no match. */
int Builder::getExistingTorsionIdx(int ai, int aj, int ak, int al) const {
  //mprintf("SEARCHING FOR %i %i %i %i\n", ai, aj, ak, al);
  int idx = -1;
  for (Tarray::const_iterator it = internalTorsions_.begin(); it != internalTorsions_.end(); ++it)
  {
    //mprintf("\t\t%i %i %i %i\n", it->AtI(), it->AtJ(), it->AtK(), it->AtL());
    if ( (it->AtI() == ai && it->AtJ() == aj && it->AtK() == ak && it->AtL() == al) ||
         (it->AtI() == al && it->AtJ() == ak && it->AtK() == aj && it->AtL() == ai) )
    {
      idx = (int)(it - internalTorsions_.begin());
      //mprintf("FOUND at index %i\n", idx);
      break;
    }
  }
  return idx;
}

/** \return Index of existing angle matching the 3 given atoms, -1 for no match. */
int Builder::getExistingAngleIdx(int ai, int aj, int ak) const {
  int idx = -1;
  for (Aarray::const_iterator it = internalAngles_.begin(); it != internalAngles_.end(); ++it)
  {
    if ((it->AtI() == ai && it->AtJ() == aj && it->AtK() == ak) ||
        (it->AtI() == ak && it->AtJ() == aj && it->AtK() == ai))
    {
      idx = (int)(it - internalAngles_.begin());
      break;
    }
  }
  return idx;
}

/** \return Index of existing bond matching the 2 given atoms, -1 for no match. */
int Builder::getExistingBondIdx(int ai, int aj) const {
  int idx = -1;
  for (Larray::const_iterator it = internalBonds_.begin(); it != internalBonds_.end(); ++it)
  {
    if ((it->AtI() == ai && it->AtJ() == aj) ||
        (it->AtI() == aj && it->AtJ() == ai))
    {
      idx = (int)(it - internalBonds_.begin());
      break;
    }
  }
  return idx;
}

/** \return Index of existing chirality value matching given atom, 1 for no match. */
int Builder::getExistingChiralityIdx(int ai) const {
  int idx = -1;
  for (Carray::const_iterator it = internalChirality_.begin(); it != internalChirality_.end(); ++it)
  {
    if (it->AtI() == ai) {
      idx = (int)(it - internalChirality_.begin());
      break;
    }
  }
  return idx;
}

/** Determine hybridization in the same manner as leap */
AtomType::HybridizationType Builder::getAtomHybridization(Atom const& aAtom, std::vector<Atom> const& atoms) const {
  AtomType::HybridizationType H1 = AtomType::UNKNOWN_HYBRIDIZATION;
  // Check if hybridization is defined in the parameter set.
  if (params_ != 0) {
    ParmHolder<AtomType>::const_iterator it;
    if (aAtom.HasType()) {
      it = params_->AT().GetParam( TypeNameHolder(aAtom.Type()) );
      if (it != params_->AT().end())
        H1 = it->second.Hybridization();
    }
  }
  // Use the cpptraj guess.
  if (H1 == AtomType::UNKNOWN_HYBRIDIZATION) {
    H1 = GuessAtomHybridization(aAtom, atoms);
    // FIXME this is a bit of a hack. If we get an SP type here its likely we
    //       are generating internals for a fragment that has not been bonded
    //       yet. Do SP2 instead.
    if (H1 == AtomType::SP)
      H1 = AtomType::SP2;
  }
  // If still unknown, guess in the same way leap does.
  if (H1 == AtomType::UNKNOWN_HYBRIDIZATION) {
    // TODO bond orders?
    int iSingle = 0;
    int iDouble = 0;
    int iTriple = 0;
    int iAromatic = 0;
    for (int i = 0; i < aAtom.Nbonds(); i++)
      iSingle++;
    //printf("iAtomHybridization: %s isingle=%i\n", iSingle);
    if ( iAromatic != 0 )       H1 = AtomType::SP2;
                // one or more triple bonds makes the atom SP1
    else if ( iTriple != 0 )    H1 = AtomType::SP;
                // Two or more double bonds makes the atom linear, SP1
    else if ( iDouble >= 2 )    H1 = AtomType::SP;
                // One double bond makes the atom SP2
    else if ( iDouble != 0 )    H1 = AtomType::SP2;
                // Otherwise the atom is SP3
    else                        H1 = AtomType::SP3;
  }

  return H1;
}

/** Model bond */
double Builder::ModelBondLength(int ai, int aj, Topology const& topIn) const {
  // First look up parameter
  double dist = 0;
  if (getLengthParam(dist, ai, aj, topIn)) {
    return dist;
  }
  Atom const& AI = topIn[ai];
  Atom const& AJ = topIn[aj];
  AtomType::HybridizationType hybridI = getAtomHybridization( AI, topIn.Atoms() );
  AtomType::HybridizationType hybridJ = getAtomHybridization( AJ, topIn.Atoms() );

  if (hybridI == AtomType::UNKNOWN_HYBRIDIZATION ||
      hybridJ == AtomType::UNKNOWN_HYBRIDIZATION)
  {
    // Default to bond length based on elements
    dist = Atom::GetBondLength( AI.Element(), AJ.Element() );
  } else {
    // Use leap method based on atom hybridization
    AtomType::HybridizationType hybrid1, hybrid2;
    if (hybridI < hybridJ) {
      hybrid1 = hybridI;
      hybrid2 = hybridJ;
    } else {
      hybrid1 = hybridJ;
      hybrid2 = hybridI;
    }
    if (hybrid1 == AtomType::SP3 && hybrid2 == AtomType::SP3)
      dist = 1.5;
    else if (hybrid1 == AtomType::SP2 && hybrid2 == AtomType::SP3)
      dist = 1.4;
    else if (hybrid1 == AtomType::SP && hybrid2 == AtomType::SP3)
      dist = 1.3;
    else if (hybrid1 == AtomType::SP2 && hybrid2 == AtomType::SP2)
      dist = 1.35;
    else if (hybrid1 == AtomType::SP && hybrid2 == AtomType::SP2)
      dist = 1.3;
    else if (hybrid1 == AtomType::SP && hybrid2 == AtomType::SP)
      dist = 1.1;
    else
      dist =  Atom::GetBondLength( AI.Element(), AJ.Element() );
  }
  return dist;
}

/** Model angle */
double Builder::ModelBondAngle(int ai, int aj, int ak, Topology const& topIn) const {
  // First look up parameter
  double theta = 0;
  if (getAngleParam(theta, ai, aj, ak, topIn)) {
    return theta;
  }
  Atom const& AJ = topIn[aj];
  // Figure out hybridization and chirality of atom j.
  if (debug_ > 0) {
    mprintf("DEBUG:\t\tI %s Nbonds: %i\n", topIn[ai].ElementName(), topIn[ai].Nbonds());
    mprintf("DEBUG:\t\tJ %s Nbonds: %i\n", AJ.ElementName(), AJ.Nbonds());
    mprintf("DEBUG:\t\tK %s Nbonds: %i\n", topIn[ak].ElementName(), topIn[ak].Nbonds());
  }
  AtomType::HybridizationType hybrid = getAtomHybridization( AJ, topIn.Atoms() );

  // Guess hybrid if needed
  //if (hybrid == AtomType::UNKNOWN_HYBRIDIZATION)
  //  hybrid = GuessAtomHybridization(AJ, topIn.Atoms());
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
  return theta;
}

/** Model torsion */
void Builder::ModelTorsion(TorsionModel const& MT, unsigned int iBondX, unsigned int iBondY, double dvalIn)
{
  if (iBondX >= MT.SortedAx().size() ||
      iBondY >= MT.SortedAy().size())
    return;
  mprintf("CALLING ModelTorsion for iBondX=%u iBondY=%u dVal=%g\n",iBondX,iBondY,dvalIn*Constants::RADDEG);
  MockAtom const& AA = MT.SortedAx()[iBondX];
  MockAtom const& AX = MT.AtX();
  MockAtom const& AY = MT.AtY();
  MockAtom const& AD = MT.SortedAy()[iBondY];
  if ( !(AA.BuildInternals() || AD.BuildInternals()) ) {
    mprintf("%s does not need internals.\n", *((*currentTop_)[AA.Idx()].Name()));
    return;
  }
  int aa = AA.Idx();
  int ax = AX.Idx();
  int ay = AY.Idx();
  int ad = AD.Idx();
  // If the coordinates for the atoms are defined then
  // measure the torsion angle between them and use that for
  // the internal.
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
    mprintf(" %s replacing dval with %f\n", *((*currentTop_)[aa].Name()), phiVal*Constants::RADDEG);
  }
  // Look for an existing internal
  int icIdx = getExistingTorsionIdx( aa, ax, ay, ad );
  if (icIdx < 0) {
    mprintf("++++Torsion INTERNAL: %f to %s - %s - %s - %s\n", phiVal*Constants::RADDEG,
            currentTop_->LeapName(aa).c_str(),
            currentTop_->LeapName(ax).c_str(),
            currentTop_->LeapName(ay).c_str(),
            currentTop_->LeapName(ad).c_str());
    internalTorsions_.push_back( InternalTorsion(aa, ax, ay, ad, phiVal) );
  } else {
    mprintf( "Torsional INTERNAL already exists: %f\n", internalTorsions_[icIdx].PhiVal()*Constants::RADDEG );
  }
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
  mprintf("In ModelCreateSp2Sp2Torsions, dAbsolute= %g, dADOffset= %g, d180= %g, d0= %g\n",
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

/** Assign torsions around bonded atoms in manner similar to LEaP's ModelAssignTorsionsAround. */
int Builder::assignTorsionsAroundBond(int a1, int a2, Frame const& frameIn, Topology const& topIn, Barray const& hasPosition, int aAtomIdx)
{
  // Save address of current topology.
  // These are required for the ModelTorsion routine. TODO zero at the end?
  currentTop_ = &topIn;
  // No need to do this if either atom only has 1 bond.
  if (topIn[a1].Nbonds() < 2 || topIn[a2].Nbonds() < 2)
    return 0;
  // Get atom hybridizations
  AtomType::HybridizationType H1 = getAtomHybridization( topIn[a1], topIn.Atoms() );
  AtomType::HybridizationType H2 = getAtomHybridization( topIn[a2], topIn.Atoms() );//AtomType::UNKNOWN_HYBRIDIZATION;
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
  static const char* hstr[] = { "SP1", "SP2", "SP3", "Unknown" };
  mprintf("DEBUG: assignTorsionsAroundBond: AX= %s (%s)  AY= %s (%s), aAtomIdx= %i",
          *(topIn[ax].Name()), hstr[Hx],
          *(topIn[ay].Name()), hstr[Hy], aAtomIdx+1);
  if (aAtomIdx != -1) mprintf(" %s", topIn.AtomMaskName(aAtomIdx).c_str()); // DEBUG
  mprintf("\n"); // DEBUG
  TorsionModel mT;
  if (mT.InitTorsion( ax, ay, frameIn, topIn, hasPosition, aAtomIdx )) {
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
    Iarray iTorsions = getExistingTorsionIdxs(ax, ay);
    if (!iTorsions.empty()) {
      mprintf("Using INTERNALs to fit new torsions around: %s - %s\n",
              topIn.LeapName(ax).c_str(), topIn.LeapName(ay).c_str());
      if (mT.BuildMockExternals(iTorsions, internalTorsions_, topIn)) {
        mprinterr("Error: Building mock externals around %s - %s failed.\n",
                  topIn.AtomMaskName(ax).c_str(), topIn.AtomMaskName(ay).c_str());
        return 1;
      }
      //return 0; // FIXME DEBUG 
    } else {
      mprintf("Completely free in assigning new torsions for: %s - %s\n",
              topIn.LeapName(ax).c_str(), topIn.LeapName(ay).c_str());
      //return 0; // FIXME DEBUG 
    }
  } else {
    // Use existing atoms to determine torsions
    mprintf("Using externals to fit new torsions around: %s - %s\n",
            topIn.LeapName(ax).c_str(),
            topIn.LeapName(ay).c_str());
  }
  // Get chiralities around X and Y if they exist
  double chiX = 0;
  int Xcidx = getExistingChiralityIdx( ax );
  if (Xcidx != -1) chiX = internalChirality_[Xcidx].ChiralVal();
  double chiY = 0;
  int Ycidx = getExistingChiralityIdx( ay );
  if (Ycidx != -1) chiY = internalChirality_[Ycidx].ChiralVal();
  // Set up the torsion model
  if (mT.SetupTorsion(Hx, Hy, topIn, chiX, chiY))
  {
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

/** Build angle internal. */
void Builder::buildAngleInternal(int a1, int a2, int a3, Frame const& frameIn, Topology const& topIn,
                                 Barray const& hasPosition)
{
    double dValue = 0;
    if (hasPosition[a1] &&
        hasPosition[a2] &&
        hasPosition[a3])
    {
      dValue = CalcAngle( frameIn.XYZ(a1), frameIn.XYZ(a2), frameIn.XYZ(a3) );
    } else {
      dValue = ModelBondAngle( a1, a2, a3, topIn );
    }
    internalAngles_.push_back( InternalAngle(a1, a2, a3, dValue) );
    mprintf("++++Angle INTERNAL: %f  for %s - %s - %s\n", dValue*Constants::RADDEG,
            topIn.LeapName(a1).c_str(),
            topIn.LeapName(a2).c_str(),
            topIn.LeapName(a3).c_str());
}

/** Build bond internal. */
void Builder::buildBondInternal(int a1, int a2, Frame const& frameIn, Topology const& topIn,
                                Barray const& hasPosition)
{
    double dValue = 0;
    if (hasPosition[a1] &&
        hasPosition[a2])
    {
      dValue = sqrt(DIST2_NoImage( frameIn.XYZ(a1), frameIn.XYZ(a2) ) );
    } else {
      dValue = ModelBondLength( a1, a2, topIn );
    }
    internalBonds_.push_back( InternalBond(a1, a2, dValue) );
    mprintf("++++Bond INTERNAL: %f  for %s - %s\n", dValue,
            topIn.LeapName(a1).c_str(),
            topIn.LeapName(a2).c_str());
}


/** Determine chirality around a single atom.
  * \return 1 if chirality was determined, 0 if left undefined.
  */
int Builder::determineChirality(double& dChi, int at, Frame const& frameIn, Topology const& topIn, Barray const& hasPosition)
{
  using namespace Cpptraj::Structure::Chirality;
  dChi = 0.0;
  if (!hasPosition[at]) {
    return 0;
  }
  // Only check atoms with 3 or 4 bonds
  Atom const& A0 = topIn[at];
  mprintf("CHIRALITY CALCULATION FOR %s (nbonds= %i)\n", *(A0.Name()), A0.Nbonds());
  for (Atom::bond_iterator bat = A0.bondbegin(); bat != A0.bondend(); ++bat)
    mprintf("\tneighbor %s (id= %i) [%i]\n", *(topIn[*bat].Name()), *bat, (int)hasPosition[*bat]);
  if ( A0.Nbonds() == 3 ||
       A0.Nbonds() == 4 )
  {
    int aAtomA = -1;
    int aAtomB = -1;
    int aAtomC = -1;
    int aAtomD = -1;
    chiralityOrderNeighbors(A0, aAtomA, aAtomB, aAtomC, aAtomD);

    bool knowA = (aAtomA != -1 && hasPosition[aAtomA]);
    bool knowB = (aAtomB != -1 && hasPosition[aAtomB]);
    bool knowC = (aAtomC != -1 && hasPosition[aAtomC]);
    bool knowD = (aAtomD != -1 && hasPosition[aAtomD]);

    if (knowA) mprintf("Chirality order A %s %i\n", *(topIn[aAtomA].Name()), aAtomA);
    if (knowB) mprintf("Chirality order B %s %i\n", *(topIn[aAtomB].Name()), aAtomB);
    if (knowC) mprintf("Chirality order C %s %i\n", *(topIn[aAtomC].Name()), aAtomC);
    if (knowD) mprintf("Chirality order D %s %i\n", *(topIn[aAtomD].Name()), aAtomD);

    Vec3 vPA, vPB, vPC, vPD;
    if (knowA) vPA = Vec3(frameIn.XYZ(aAtomA));
    if (knowB) vPB = Vec3(frameIn.XYZ(aAtomB));
    if (knowC) vPC = Vec3(frameIn.XYZ(aAtomC));
    if (knowD) vPD = Vec3(frameIn.XYZ(aAtomD));

    dChi = VectorAtomNormalizedChirality( Vec3(frameIn.XYZ(at)),
                                          vPA, knowA,
                                          vPB, knowB,
                                          vPC, knowC,
                                          vPD, knowD );
    return 1;
  }
/*
  BuildAtom bldAt;
  if (topIn[at].Nbonds() > 2 && hasPosition[at]) {
    // All bonded atoms must have position
    bool bonded_atoms_have_position = true;
    for (Atom::bond_iterator bat = topIn[at].bondbegin(); bat != topIn[at].bondend(); ++bat)
    {
      if (!hasPosition[*bat]) {
        bonded_atoms_have_position = false;
        break;
      }
    }
    if (bonded_atoms_have_position) {
      if (bldAt.DetermineChirality(at, topIn, frameIn, 0)) {// FIXME debug level
        mprinterr("Error: Problem determining chirality of atom %s\n",
                  topIn.AtomMaskName(at).c_str());
        return 1;
      }
    }
    mprintf("Got chirality from external coordinates\n" );
    mprintf("++++Chirality INTERNAL: %f  for %s\n", bldAt.TorsionVal()*Constants::RADDEG,
            topIn.LeapName(at).c_str());
  } else {
    mprintf("Left chirality undefined for %s\n",topIn.LeapName(at).c_str() );
  }*/
  return 0;
}

/** Generate internal coordinates in the same manner as LEaP's
  * BuildInternalsForContainer.
  */
int Builder::GenerateInternals(Frame const& frameIn, Topology const& topIn, Barray const& hasPosition)
{
  mprintf("DEBUG: ----- Entering Builder::GenerateInternals. -----\n");
  // First generate the bond array for use in determining torsions.
  BondArray bonds = GenerateBondArray( topIn.Residues(), topIn.Atoms() );
  // Loop over bonds to determine torsions.
  for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
  {
    if (topIn[bnd->A1()].Nbonds() > 1 && topIn[bnd->A2()].Nbonds() > 1) mprintf("Building torsion INTERNALs around: %s - %s\n", topIn.LeapName(bnd->A1()).c_str(), topIn.LeapName(bnd->A2()).c_str()); // DEBUG
    if (assignTorsionsAroundBond( bnd->A1(), bnd->A2(), frameIn, topIn, hasPosition, -1 )) {
      mprinterr("Error Assign torsions around bond %s - %s failed.\n",
                topIn.AtomMaskName(bnd->A1()).c_str(),
                topIn.AtomMaskName(bnd->A2()).c_str());
      return 1;
    }
  }
  // Loop over angles.
  AngleArray angles = GenerateAngleArray( topIn.Residues(), topIn.Atoms() );
  for (AngleArray::const_iterator ang = angles.begin(); ang != angles.end(); ++ang)
  {
    buildAngleInternal(ang->A1(), ang->A2(), ang->A3(), frameIn, topIn, hasPosition);
  }
  // Loop over bonds
  for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
  {
    buildBondInternal(bnd->A1(), bnd->A2(), frameIn, topIn, hasPosition);
  }
  // Chirality
  Iarray atomIdxs = GenerateAtomArray(topIn.Residues(), topIn.Atoms());
  for (Iarray::const_iterator it = atomIdxs.begin(); it != atomIdxs.end(); ++it) {
    int at = *it; // FIXME just use *it
  //for (int at = 0; at != topIn.Natom(); ++at) {
    double dValue = 0;
    if (determineChirality(dValue, at, frameIn, topIn, hasPosition)) {
      mprintf("Got chirality from external coordinates\n" );
      mprintf("++++Chirality INTERNAL: %f  for %s\n", dValue,
              topIn.LeapName(at).c_str());
    } else {
      mprintf("Left chirality undefined for %s\n",topIn.LeapName(at).c_str() );
    }
    internalChirality_.push_back( InternalChirality(at, dValue) );
  }
  mprintf("DEBUG: ----- Leaving Builder::GenerateInternals. ------\n");
  return 0;
}

/** Build internal coordinates around an atom. */
int Builder::generateAtomInternals(int at, Frame const& frameIn, Topology const& topIn, Barray const& hasPosition)
{
  mprintf( "Building internals for: %s\n", topIn.LeapName(at).c_str());
  // Torsions
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
        Iarray iTorsions = getExistingTorsionIdxs(*bat, *cat);
        int iShouldBe = (AtB.Nbonds() - 1) * (AtC.Nbonds() - 1);
        mprintf("ISHOULDBE= %i ITORSIONS= %zu\n", iShouldBe, iTorsions.size());
        if (iShouldBe != (int)iTorsions.size()) {
          assignTorsionsAroundBond(*bat, *cat, frameIn, topIn, hasPosition, at);
        }
      }
    }
  }
  // Angles
  for (Atom::bond_iterator bat = AtA.bondbegin(); bat != AtA.bondend(); ++bat) {
    Atom const& AtB = topIn[*bat];
    for (Atom::bond_iterator cat = AtB.bondbegin(); cat != AtB.bondend(); ++cat) {
      if (*cat != at) {
        mprintf("Building angle INTERNAL for: %s - %s - %s\n",
                topIn.LeapName(at).c_str(),
                topIn.LeapName(*bat).c_str(),
                topIn.LeapName(*cat).c_str());
        int aidx = getExistingAngleIdx(at, *bat, *cat);
        if (aidx < 0) {
          double dValue = 0;
          if (hasPosition[at] &&
              hasPosition[*bat] &&
              hasPosition[*cat])
          {
            mprintf("Got bond angle from externals\n");
            dValue = CalcAngle(frameIn.XYZ(at), frameIn.XYZ(*bat), frameIn.XYZ(*cat));
          } else {
            mprintf("Got bond angle from model builder\n");
            dValue = ModelBondAngle(at, *bat, *cat, topIn);
          }
          mprintf("++++Angle INTERNAL: %f  for %s - %s - %s\n", dValue*Constants::RADDEG,
                  topIn.LeapName(at).c_str(),
                  topIn.LeapName(*bat).c_str(),
                  topIn.LeapName(*cat).c_str());
          internalAngles_.push_back(InternalAngle(at, *bat, *cat, dValue));
        } else {
          mprintf("Angle INTERNAL was already defined\n");
        }
      }
    } // END loop over atoms bonded to B
  } // END loop over atoms bonded to A
  // Bonds
  for (Atom::bond_iterator bat = AtA.bondbegin(); bat != AtA.bondend(); ++bat) {
    mprintf("Building bond INTERNAL for: %s - %s\n",
            topIn.LeapName(at).c_str(),
            topIn.LeapName(*bat).c_str());
    int bidx = getExistingBondIdx(at, *bat);
    if (bidx < 0) {
      double dValue = 0;
      if (hasPosition[at] &&
          hasPosition[*bat])
      {
        mprintf("Got bond length from externals\n");
        dValue = sqrt(DIST2_NoImage(frameIn.XYZ(at), frameIn.XYZ(*bat)));
      } else {
        mprintf("Got bond length from the model builder\n");
        dValue = ModelBondLength(at, *bat, topIn);
      }
      mprintf("++++Bond INTERNAL: %f  for %s - %s\n", dValue,
              topIn.LeapName(at).c_str(),
              topIn.LeapName(*bat).c_str());
      internalBonds_.push_back(InternalBond(at, *bat, dValue));
    } else {
      mprintf( "Bond length INTERNAL already defined\n" );
    }
  } // END loop over atoms bonded to A
  // Chirality
  double dChi = 0;
  int cidx = getExistingChiralityIdx(at);
  if (determineChirality(dChi, at, frameIn, topIn, hasPosition)) {
    mprintf("Got chirality from external coordinates\n" );
    mprintf("++++Chirality INTERNAL: %f  for %s\n", dChi,
              topIn.LeapName(at).c_str());
    if (cidx == -1)
      internalChirality_.push_back( InternalChirality(at, dChi) );
    else {
      // Check that this chirality matches previously determine chirality
      if (!internalChirality_[cidx].ChiralityMatches(dChi))
        mprintf("Warning: Atom %s chirality (%f) does not match previous chirality (%f)\n",
                topIn.AtomMaskName(at).c_str(), dChi, internalChirality_[cidx].ChiralVal());
    }
  } else {
    if (cidx == -1)
      mprintf("Left chirality undefined for %s\n",topIn.LeapName(at).c_str() );
    else
      mprintf("Using already-defined chirality (%f).\n", internalChirality_[cidx].ChiralVal());
  }

  return 0;
}
 
/** Generate internal coordinates around a bond linking two residues
  * in the same manner as LEaP.
  * \param zmatrix Hold output ICs
  * \param at0 Atom in residue we are linking to (i.e. the current residue).
  * \param at1 Atom in resiude we are linking from.
  */
int Builder::GenerateInternalsAroundLink(int at0, int at1,
                                         Frame const& frameIn, Topology const& topIn,
                                         Barray const& hasPosition,
                                         BuildType btype)
{
  mprintf("DEBUG: ----- Entering Builder::GenerateInternalsAroundLink. -----\n");
  mprintf("DEBUG: Link: %s to %s\n", topIn.AtomMaskName(at0).c_str(), topIn.AtomMaskName(at1).c_str());
  // Sanity check
  Atom const& A0 = topIn[at0];
  Atom const& A1 = topIn[at1];
  if (A0.ResNum() == A1.ResNum()) {
    mprintf("Warning: Builder::GenerateInternalsAroundLink(): Atoms are in the same residue.\n");
    //return 1;
  }
  // In order to mimic the way LEaP does things, mark all atoms before
  // this residue as having position, and all other atoms as not having
  // position.
  Residue const& R0 = topIn.Res(A0.ResNum());
  Barray tmpHasPosition( topIn.Natom(), false );
  for (int at = 0; at < R0.FirstAtom(); at++)
    tmpHasPosition[at] = true;
  // Create spanning tree across the link
  int actualAt1;
  if (btype == BUILD)
    // Want the full tree
    actualAt1 = -1;
  else
    // Only want the 'forward' tree
    actualAt1 = at1;
  Iarray span_atoms = GenerateSpanningTree(at0, actualAt1, 4, topIn.Atoms());
  for (Iarray::const_iterator it = span_atoms.begin(); it != span_atoms.end(); ++it)
  {
    mprintf("SPANNING TREE ATOM: %s\n", *(topIn[*it].Name()));
    if (generateAtomInternals(*it, frameIn, topIn, tmpHasPosition)) {
      mprinterr("Error: Could not generate internals for atom %s\n", topIn.AtomMaskName(*it).c_str());
      return 1;
    }
  }
  // FIXME this is a hack to make certain we have all the angle/bond terms we need
  mprintf("DEBUG: LOOKING FOR MISSING ANGLE/BOND PARAMS.\n");
  for (Tarray::const_iterator dih = internalTorsions_.begin(); dih != internalTorsions_.end(); ++dih)
  {
    int idx = getExistingAngleIdx(dih->AtI(), dih->AtJ(), dih->AtK());
    if (idx < 0) buildAngleInternal( dih->AtI(), dih->AtJ(), dih->AtK(), frameIn, topIn, hasPosition );
    idx = getExistingAngleIdx(dih->AtJ(), dih->AtK(), dih->AtL());
    if (idx < 0) buildAngleInternal( dih->AtJ(), dih->AtK(), dih->AtL(), frameIn, topIn, hasPosition );
    idx = getExistingBondIdx(dih->AtI(), dih->AtJ());
    if (idx < 0) buildBondInternal( dih->AtI(), dih->AtJ(), frameIn, topIn, hasPosition );
    idx = getExistingBondIdx(dih->AtK(), dih->AtL());
    if (idx < 0) buildBondInternal( dih->AtK(), dih->AtL(), frameIn, topIn, hasPosition );
  }

  mprintf("DEBUG: ----- Leaving Builder::GenerateInternalsAroundLink. -----\n");
  return 0;
}

/** For debugging, print all the complete internals associated with the given atom. */
void Builder::printAllInternalsForAtom(int at, Topology const& topIn, Barray const& hasPosition) const
{
  mprintf("DEBUG: All internals for atom %s\n", topIn.LeapName(at).c_str());
  for (Tarray::const_iterator dih = internalTorsions_.begin(); dih != internalTorsions_.end(); ++dih)
  {
    if (at == dih->AtI()) {
      if (hasPosition[dih->AtJ()] &&
          hasPosition[dih->AtK()] &&
          hasPosition[dih->AtL()])
      {
        mprintf("DEBUG:\t\t%s - %s - %s Phi= %f\n",
                topIn.LeapName(dih->AtJ()).c_str(),
                topIn.LeapName(dih->AtK()).c_str(),
                topIn.LeapName(dih->AtL()).c_str(),
                dih->PhiVal()*Constants::RADDEG);
      }
    } else if (at == dih->AtL()) {
      if (hasPosition[dih->AtK()] &&
          hasPosition[dih->AtJ()] &&
          hasPosition[dih->AtI()])
      {
        mprintf("DEBUG:\t\t%s - %s - %s Phi= %f\n",
                topIn.LeapName(dih->AtK()).c_str(),
                topIn.LeapName(dih->AtJ()).c_str(),
                topIn.LeapName(dih->AtI()).c_str(),
                dih->PhiVal()*Constants::RADDEG);
      }
    }
  }
}

/** Find internal coordinates for given atom.
  * Find torsion that contains this atom as one of the end atoms. The other
  * three atoms must have known position.
  * There must be a bond and angle that lie on the torsion and include
  * the given atom as a terminal atom.
  * \return 1 if complete internal coords were found, 0 if not.
  */
int Builder::getIcFromInternals(InternalCoords& icOut, int at, Barray const& hasPosition) const
{
  for (Tarray::const_iterator dih = internalTorsions_.begin(); dih != internalTorsions_.end(); ++dih)
  {
    if (at == dih->AtI()) {
      if (hasPosition[dih->AtJ()] &&
          hasPosition[dih->AtK()] &&
          hasPosition[dih->AtL()])
      {
        int bidx = getExistingBondIdx(dih->AtI(), dih->AtJ());
        if (bidx > -1) {
          int aidx = getExistingAngleIdx(dih->AtI(), dih->AtJ(), dih->AtK());
          if (aidx > -1) {
            icOut = InternalCoords(dih->AtI(), dih->AtJ(), dih->AtK(), dih->AtL(),
                                   internalBonds_[bidx].DistVal(),
                                   internalAngles_[aidx].ThetaVal()*Constants::RADDEG,
                                   dih->PhiVal()*Constants::RADDEG);
            return 1;
          }
        }
      }
    } else if (at == dih->AtL()) {
      if (hasPosition[dih->AtK()] &&
          hasPosition[dih->AtJ()] &&
          hasPosition[dih->AtI()])
      {
        int bidx = getExistingBondIdx(dih->AtL(), dih->AtK());
        if (bidx > -1) {
          int aidx = getExistingAngleIdx(dih->AtL(), dih->AtK(), dih->AtJ());
          if (aidx > -1) {
            icOut = InternalCoords(dih->AtL(), dih->AtK(), dih->AtJ(), dih->AtI(),
                                   internalBonds_[bidx].DistVal(),
                                   internalAngles_[aidx].ThetaVal()*Constants::RADDEG,
                                   dih->PhiVal()*Constants::RADDEG);
            return 1;
          }
        }
      }
    }
  } // END loop over internal torsions
  return 0;
}

/** \return Array of residues with atoms that need positions. */
std::vector<Residue> Builder::residuesThatNeedPositions(Topology const& topIn,
                                                        Barray const& hasPosition)
const
{
  std::vector<Residue> residues;
  std::vector<int> Rnums;
  for (Tarray::const_iterator dih = internalTorsions_.begin();
                              dih != internalTorsions_.end(); ++dih)
  {
    if (!hasPosition[dih->AtI()]) {
      int rnum = topIn[dih->AtI()].ResNum();
      bool has_rnum = false;
      for (std::vector<int>::const_iterator it = Rnums.begin(); it != Rnums.end(); ++it) {
        if (*it == rnum) {
          has_rnum = true;
          break;
        }
      }
      if (!has_rnum) {
        mprintf("DEBUG: Need to build for residue %s\n", topIn.TruncResNameNum(rnum).c_str());
        residues.push_back( topIn.Res(rnum) );
        Rnums.push_back(rnum);
      }
    }
  }
  return residues;
} 

/** Build coordinates for any atom with an internal that does
  * not have its position set. Use same atom order as leap sequence.
  */
int Builder::BuildSequenceFromInternals(Frame& frameOut, Topology const& topIn,
                                        Barray& hasPosition, int at0, int at1)
const
{
  // Create a list of residues that have atoms that need positions
  //std::vector<Residue> residues = residuesThatNeedPositions(topIn, hasPosition);
  Iarray atomIndices = GenerateSpanningTree(at0, at1, -1, topIn.Atoms());

  return buildExternalsForAtoms( atomIndices, frameOut, topIn, hasPosition ); 
  return 0;
}

/** Build coordinates for any atom with an internal that does
  * not have its position set.
  */ 
int Builder::BuildFromInternals(Frame& frameOut, Topology const& topIn, Barray& hasPosition)
const
{
  // Create a list of residues that have atoms that need positions
  std::vector<Residue> residues = residuesThatNeedPositions(topIn, hasPosition);
//  std::vector<int> Rnums;
//  for (Tarray::const_iterator dih = internalTorsions_.begin();
//                              dih != internalTorsions_.end(); ++dih)
//  {
//    if (!hasPosition[dih->AtI()]) {
//      int rnum = topIn[dih->AtI()].ResNum();
//      bool has_rnum = false;
//      for (std::vector<int>::const_iterator it = Rnums.begin(); it != Rnums.end(); ++it) {
//        if (*it == rnum) {
//          has_rnum = true;
//          break;
//        }
//      }
//      if (!has_rnum) {
//        mprintf("DEBUG: Need to build for residue %s\n", topIn.TruncResNameNum(rnum).c_str());
//        residues.push_back( topIn.Res(rnum) );
//        Rnums.push_back(rnum);
//      }
//    }
//  }
  // Generate array over residue in same order that leap would do
  Iarray atomIndices = GenerateAtomArray(residues, topIn.Atoms());
  residues.clear();
//  Rnums.clear();
  return buildExternalsForAtoms( atomIndices, frameOut, topIn, hasPosition ); 
}

/** Build externals for atoms in the given array using internals. */
int Builder::buildExternalsForAtoms(Iarray const& atomIndices,
                                    Frame& frameOut, Topology const& topIn,
                                    Barray& hasPosition)
const
{
  // Count how many atoms need their positions set
  unsigned int nAtomsThatNeedPositions = 0;
  for (std::vector<int>::const_iterator it = atomIndices.begin();
                                        it != atomIndices.end(); ++it)
    if (!hasPosition[*it])
      nAtomsThatNeedPositions++;
  mprintf("DEBUG: %u atoms need positions.\n", nAtomsThatNeedPositions);
  if (nAtomsThatNeedPositions == 0) return 0;

  // Loop over residue atoms
  while (nAtomsThatNeedPositions > 0) {
    unsigned int nAtomsBuilt = 0;
    for (std::vector<int>::const_iterator idx = atomIndices.begin();
                                          idx != atomIndices.end(); ++idx)
    {
      int at = *idx;
      int atToBuildAround = -1;
      if (!hasPosition[at]) {
        // Position of atom is not known.
        //mprintf("BUILD: Position of %s is not known.\n", *(topIn[at].Name()));
        // Is this bonded to an atom with known position?
        for (Atom::bond_iterator bat = topIn[at].bondbegin(); bat != topIn[at].bondend(); ++bat) {
          if (hasPosition[*bat]) {
            atToBuildAround = *bat;
            break;
          }
        }
      } else {
        // Position of atom is known.
        //mprintf("BUILD: Position of %s is known.\n", *(topIn[at].Name()));
        atToBuildAround = at;
      }
      // Build unknown positions around known atom
      if (atToBuildAround != -1) {
        mprintf("Building externals from %s\n", topIn.LeapName(atToBuildAround).c_str());
        Atom const& bAtom = topIn[atToBuildAround];
        for (Atom::bond_iterator bat = bAtom.bondbegin(); bat != bAtom.bondend(); ++bat)
        {
          if (!hasPosition[*bat]) {
            // Find an internal coordinate.
            InternalCoords ic;
            if (getIcFromInternals(ic, *bat, hasPosition)) {
              //printAllInternalsForAtom(*bat, topIn, hasPosition); // DEBUG

              mprintf("Building atom %s using torsion/angle/bond\n", topIn.LeapName(*bat).c_str());
              mprintf("Using - %s - %s - %s\n",
                      topIn.LeapName(ic.AtJ()).c_str(),
                      topIn.LeapName(ic.AtK()).c_str(),
                      topIn.LeapName(ic.AtL()).c_str());
              mprintf( "Torsion = %f\n", ic.Phi() );
              mprintf( "Angle   = %f\n", ic.Theta() );
              mprintf( "Bond    = %f\n", ic.Dist() );
              Vec3 posI = Zmatrix::AtomIposition(ic, frameOut);
              mprintf( "ZMatrixAll:  %f,%f,%f\n", posI[0], posI[1], posI[2]);
              frameOut.SetXYZ( ic.AtI(), posI );
              hasPosition[ ic.AtI() ] = true;
              nAtomsBuilt++;
              nAtomsThatNeedPositions--;
            }
          }
        } // END loop over atoms bonded to atom with known position
      }
    } // END loop over residue atoms
    // If we built no atoms this is a problem
    if (nAtomsBuilt < 1) {
      mprinterr("Error: No more atoms could be built for %s\n", topIn.c_str());
      return 1;
    }
  } // END loop while atoms need position
  return 0;
}

/** Generate a Zmatrix from the current internals. TODO only for atoms that need it? */
/*
int Builder::GetZmatrixFromInternals(Zmatrix& zmatrix, Topology const& topIn) const {
  mprintf("DEBUG: ----- Enter GetZmatrixFromInternals -----\n");
  zmatrix.clear();

  for (Tarray::const_iterator dih = internalTorsions_.begin(); dih != internalTorsions_.end(); ++dih)
  {
    // Get angles i-j-k and j-k-l
    int aidx0 = getExistingAngleIdx(dih->AtI(), dih->AtJ(), dih->AtK());
    int aidx1 = getExistingAngleIdx(dih->AtJ(), dih->AtK(), dih->AtL());
    if (aidx0 < 0) {
      mprinterr("Error: Missing angle0 internal for %s - %s - %s\n",
                topIn.AtomMaskName(dih->AtI()).c_str(),
                topIn.AtomMaskName(dih->AtJ()).c_str(),
                topIn.AtomMaskName(dih->AtK()).c_str());
      return 1;
    }
    if (aidx1 < 0) {
      mprinterr("Error: Missing angle1 internal for %s - %s - %s\n",
                topIn.AtomMaskName(dih->AtJ()).c_str(),
                topIn.AtomMaskName(dih->AtK()).c_str(),
                topIn.AtomMaskName(dih->AtL()).c_str());
      return 1;
    }
    // Get Bonds i-j and k-l
    int bidx0 = getExistingBondIdx(dih->AtI(), dih->AtJ());
    int bidx1 = getExistingBondIdx(dih->AtK(), dih->AtL());
    if (bidx0 < 0) {
      mprinterr("Error: Missing bond0 internal for %s - %s\n",
                topIn.LeapName(dih->AtI()).c_str(),
                topIn.LeapName(dih->AtJ()).c_str());
      //mprintf("DEBUG: Internal %s - %s - %s - %s\n",
      //        topIn.LeapName(dih->AtI()).c_str(),
      //        topIn.LeapName(dih->AtJ()).c_str(),
      //        topIn.LeapName(dih->AtK()).c_str(),
      //        topIn.LeapName(dih->AtL()).c_str());
    }
    if (bidx1 < 0) {
    mprinterr("Error: Missing bond1 internal for %s - %s\n",
                topIn.AtomMaskName(dih->AtK()).c_str(),
                topIn.AtomMaskName(dih->AtL()).c_str());
    }
    // Add internal coordinates
    zmatrix.AddIC( InternalCoords(dih->AtI(), dih->AtJ(), dih->AtK(), dih->AtL(),
                                  internalBonds_[bidx0].DistVal(),
                                  internalAngles_[aidx0].ThetaVal()*Constants::RADDEG,
                                  dih->PhiVal()*Constants::RADDEG) );
    zmatrix.AddIC( InternalCoords(dih->AtL(), dih->AtK(), dih->AtJ(), dih->AtI(),
                                  internalBonds_[bidx1].DistVal(),
                                  internalAngles_[aidx1].ThetaVal()*Constants::RADDEG,
                                  dih->PhiVal()*Constants::RADDEG) );
  } // END loop over internal torsions
  mprintf("DEBUG: ----- Exit GetZmatrixFromInternals -----\n");
  return 0;
}*/
