#include "HbCalc.h"
#include <cmath> //sqrt
//#include <algorithm>
//#include <utility>
#include "../ArgList.h"
#include "../Constants.h"
#include "../CpptrajStdio.h"
#include "../DistRoutines.h"
#include "../Topology.h"

using namespace Cpptraj::HB;

/** CONSTRUCTOR */
HbCalc::HbCalc() :
  dcut2_(0),
  acut_(0)
{}

const char* HbCalc::TypeStr_[] = {
  "Solute Donor",
  "Solute Acceptor",
  "Solute Both",
  "Solvent Donor",
  "Solvent Acceptor",
  "Solvent Both",
  "Unknown"
};

/** Initialize */
int HbCalc::InitHbCalc(ArgList& argIn, int debugIn) {
  double dcut = argIn.getKeyDouble("dist",3.0);
  dcut = argIn.getKeyDouble("distance", dcut); // for PTRAJ compat.
  dcut2_ = dcut * dcut;
  acut_ = argIn.getKeyDouble("angle", 135.0 * Constants::DEGRAD);

  generalMask_.SetMaskString( argIn.GetMaskNext() );

  pairList_.InitPairList( 8.0, 0.1, debugIn );

  return 0;
}

/** Print current options */
void HbCalc::PrintHbCalcOpts() const {
  mprintf("\tSearching for atoms in mask '%s'\n", generalMask_.MaskString());
  mprintf("\tHeavy atom distance cutoff= %g Ang.\n", sqrt(dcut2_));
  if (acut_ > -1)
    mprintf("\tAngle cutoff= %g deg.\n", acut_*Constants::RADDEG);
  else
    mprintf("\tNo angle cutoff.\n");
}

/** Set up calculation */
int HbCalc::SetupHbCalc(Topology const& topIn, Box const& boxIn) {
  if (setupPairlistAtomMask( topIn )) return 1;

  if (pairList_.SetupPairList( boxIn )) return 1;

  return 0;
}

/// \return True if given atom is F, O, or N
bool HbCalc::IsFON( Atom const& at ) {
  return ( at.Element() == Atom::OXYGEN ||
           at.Element() == Atom::NITROGEN ||
           at.Element() == Atom::FLUORINE );
}

/** Set up mask for pair list. */
int HbCalc::setupPairlistAtomMask(Topology const& topIn) {
  if (topIn.SetupIntegerMask( generalMask_ )) {
    mprinterr("Error: Could not set up mask '%s'\n", generalMask_.MaskString());
    return 1;
  }
  if (generalMask_.None()) {
    mprintf("Warning: No atoms selected by mask '%s' \n", generalMask_.MaskString());
    return 1;
  }
  generalMask_.MaskInfo();
  // Decide what each atom is
  //typedef std::pair<int, Type> Ptype;
  //typedef std::vector<Ptype> Parray;
  //Parray IdxTypes;
  //IdxTypes.reserve( generalMask_.Nselected() ); // TODO not reserve?

  plMask_ = AtomMask( std::vector<int>(), topIn.Natom() );
  plTypes_.clear();
  plNames_.clear();
  plHatoms_.clear();

  int nsites = 0;
  unsigned int NacceptorOnly = 0;
  unsigned int Nboth = 0;
  unsigned int NdonorOnly = 0;
  unsigned int NumH = 0;
  for (AtomMask::const_iterator at = generalMask_.begin(); at != generalMask_.end(); ++at) {
    Atom const& currentAtom = topIn[*at];
    if (IsFON( currentAtom )) {
      // Check if there are any hydrogens bonded to this atom
      Iarray h_atoms;
      for (Atom::bond_iterator bat = currentAtom.bondbegin(); bat != currentAtom.bondend(); ++bat) {
        if (topIn[*bat].Element() == Atom::HYDROGEN) {
          h_atoms.push_back( *bat );
        }
      }
      int molnum = currentAtom.MolNum();
      Type currentType;
      if ( topIn.Mol(molnum).IsSolvent()) {
        // Solvent atom
        if (h_atoms.empty())
          currentType = VACCEPTOR;
        else
          currentType = VBOTH;
      } else {
        // Solute atom
        if (h_atoms.empty())
          currentType = ACCEPTOR;
         else {
          NumH += h_atoms.size();
          currentType = BOTH;
        }
        // Solute Count
        if (currentType == ACCEPTOR)
          NacceptorOnly++;
        else if (currentType == BOTH)
          Nboth++;
        else if (currentType == DONOR)
          NdonorOnly++;
      }
      nsites++;

      //IdxTypes.push_back( Ptype(*at, currentType) );
      plMask_.AddSelectedAtom( *at );
      plTypes_.push_back( currentType );
      plNames_.push_back( topIn.TruncResAtomName( *at ) );
      plHatoms_.push_back( h_atoms );
    }
  }
  mprintf("\tNumber of heavy atom sites: %i\n", nsites);
  mprintf("\tAcceptor-only atoms: %u\n", NacceptorOnly);
  mprintf("\tDonor/acceptor sites: %u\n", Nboth);
  mprintf("\tDonor-only sites: %u\n", NdonorOnly);
  mprintf("\t%u solute hydrogens.\n", NumH);

//  std::sort( IdxTypes.begin(), IdxTypes.end() );

//  plTypes_.reserve( IdxTypes.size() );
//  for (Parray::const_iterator it = IdxTypes.begin(); it != IdxTypes.end(); ++it) {
//    plMask_.AddSelectedAtom( it->first );
//    plTypes_.push_back( it->second );
//    plNames_.push_back( topIn.TruncResAtomName(it->first) );
//  }

  for (int idx = 0; idx != plMask_.Nselected(); idx++) {
    //mprintf("\t%8i %4s %s\n", plMask_[idx]+1, *(topIn[plMask_[idx]].Name()), TypeStr_[plTypes_[idx]]);
    mprintf("\t%8i", plMask_[idx]+1);
    mprintf(" %4s", *(topIn[plMask_[idx]].Name()));
    mprintf(" %s\n", TypeStr_[plTypes_[idx]]);
  }

  return 0;
}

/** Determine if the interaction is valid. */
bool HbCalc::validInteraction(Type t0, Type t1) {
  if (t0 == BOTH || t0 == VBOTH || t1 == BOTH || t1 == VBOTH) return true;
  // If we are here, t0/t1 must be either a DONOR or ACCEPTOR.
  if ((t0 == DONOR    || t0 == VDONOR)    && (t1 == ACCEPTOR || t1 == VACCEPTOR)) return true;
  if ((t0 == ACCEPTOR || t0 == VACCEPTOR) && (t1 == DONOR    || t1 == VDONOR)   ) return true;
  return false;
}

/** Calculate angle in radians between 3 atoms with imaging. */
double HbCalc::Angle(const double* XA, const double* XH, const double* XD, Box const& boxIn) const
{ 
  //if (imageOpt_.ImagingType() == ImageOption::NO_IMAGE)
  //  return (CalcAngle(XA, XH, XD));
  //else {
    double angle;
    Vec3 VH = Vec3(XH);
    Vec3 H_A = MinImagedVec(VH, Vec3(XA), boxIn.UnitCell(), boxIn.FracCell());
    Vec3 H_D = Vec3(XD) - VH;
    double rha = H_A.Magnitude2();
    double rhd = H_D.Magnitude2();
    if (rha > Constants::SMALL && rhd > Constants::SMALL) {
      angle = (H_A * H_D) / sqrt(rha * rhd);
      if      (angle >  1.0) angle =  1.0;
      else if (angle < -1.0) angle = -1.0;
      angle = acos(angle);
    } else
      angle = 0.0;
    return angle;
  //}
}

/** Calculate hydrogen bonds between given solute donor site and 
  * solute acceptor atom.
  * The distance cutoff should already be satisfied between donor and
  * acceptor heavy atoms.
  */
void HbCalc::CalcSiteHbonds(int frameNum, double dist2,
                            int d_idx, Iarray const& Hatoms,
                            int a_idx,
                            Frame const& frmIn, int& numHB,
                            int trajoutNum)
{
  int d_atom = plMask_[d_idx];
  int a_atom = plMask_[a_idx]; 
  // Determine if angle cutoff is satisfied
  for (Iarray::const_iterator h_atom = Hatoms.begin(); h_atom != Hatoms.end(); ++h_atom)
  {
    //double angle = 0;
    //if (acut_ > -1)
    //double angle = Angle(XYZA.Dptr(), frmIn.XYZ(*h_atom), XYZD.Dptr(), frmIn.BoxCrd());
    double angle = Angle(frmIn.XYZ(a_atom), frmIn.XYZ(*h_atom), frmIn.XYZ(d_atom), frmIn.BoxCrd());
    if ( !(angle < acut_) )
    {
      mprintf("DBG: %12s %12i %12s %12.4f %12.4f\n", plNames_[a_idx].c_str(), *h_atom + 1, plNames_[d_idx].c_str(), sqrt(dist2), angle*Constants::RADDEG);
/*#     ifdef _OPENMP
      // numHB holds thread number, will be counted later on.
      thread_HBs_[numHB].push_back( Hbond(sqrt(dist2), angle, a_atom, *h_atom, d_atom) );
#     else
      ++numHB;
      AddUU(sqrt(dist2), angle, frameNum, a_atom, *h_atom, d_atom, trajoutNum);
#     endif*/
    }
  }
}

/** Calculate hbonds between two atoms. */
void HbCalc::CalcHbonds(int frameNum, double dist2,
                        int a0idx,
                        int a1idx,
                        Frame const& frmIn, int& numHB,
                        int trajoutNum)
{
  // BOTH ACCEPTOR
  // DONOR ACCEPTOR
  // BOTH DONOR
  // ACCEPTOR DONOR
  // BOTH BOTH
  // DONOR BOTH
  // ACCEPTOR BOTH
  if ((plTypes_[a0idx] == BOTH || plTypes_[a0idx] == DONOR)    && (plTypes_[a1idx] == BOTH || plTypes_[a1idx] == ACCEPTOR)) {
    CalcSiteHbonds(frameNum, dist2, a0idx, plHatoms_[a0idx], a1idx, frmIn, numHB, trajoutNum);
  } 
  if ((plTypes_[a0idx] == BOTH || plTypes_[a0idx] == ACCEPTOR) && (plTypes_[a1idx] == BOTH || plTypes_[a1idx] == DONOR)) {
    CalcSiteHbonds(frameNum, dist2, a1idx, plHatoms_[a1idx], a0idx, frmIn, numHB, trajoutNum);
  }// else if (plTypes_[a0idx] == BOTH && plTypes_[a1idx] == BOTH) {
  //  CalcSiteHbonds(frameNum, dist2, a0idx, plHatoms_[a0idx], a0xyz, a1idx, a1xyz, frmIn, numHB, trajoutNum);
  //  CalcSiteHbonds(frameNum, dist2, a1idx, plHatoms_[a1idx], a1xyz, a0idx, a0xyz, frmIn, numHB, trajoutNum);
  //}
}

/** HB calc loop with a pairlist */
int HbCalc::RunCalc_PL(Frame const& currentFrame)
{
  int frameNum = 0; // FIXME
  int trajoutNum = 0; // FIXME
  int retVal = pairList_.CreatePairList(currentFrame,
                                        currentFrame.BoxCrd().UnitCell(),
                                        currentFrame.BoxCrd().FracCell(), plMask_);
  if (retVal < 0) {
    mprinterr("Error: Grid setup failed.\n");
    return 1;
  } else if (retVal > 0) {
    mprintf("Warning: %i atoms are off the grid.\n", retVal);
  }
  //problemAtoms_.clear();

  int Ninteractions = 0; // DEBUG
  int numHB = 0;
  int cidx;
# ifdef _OPENMP
  int mythread;
# pragma omp parallel private(cidx,mythread) 
  {
  mythread = omp_get_thread_num();
//  thread_problemAtoms_[mythread].clear();
# pragma omp for
# endif 
  for (cidx = 0; cidx < pairList_.NGridMax(); cidx++)
  {
    PairList::CellType const& thisCell = pairList_.Cell( cidx );
    if (thisCell.NatomsInGrid() > 0)
    {
      // cellList contains this cell index and all neighbors.
      PairList::Iarray const& cellList = thisCell.CellList();
      // transList contains index to translation for the neighbor.
      PairList::Iarray const& transList = thisCell.TransList();
      // Loop over all atoms of thisCell.
      for (PairList::CellType::const_iterator it0 = thisCell.begin();
                                              it0 != thisCell.end(); ++it0)
      {
        Vec3 const& xyz0 = it0->ImageCoords();
        // Calc interaction of atom to all other atoms in thisCell.
        for (PairList::CellType::const_iterator it1 = it0 + 1;
                                                it1 != thisCell.end(); ++it1)
        {
          if (validInteraction(plTypes_[it0->Idx()], plTypes_[it1->Idx()]))
          {
            Vec3 const& xyz1 = it1->ImageCoords();
            Vec3 dxyz = xyz1 - xyz0;
            double D2 = dxyz.Magnitude2();
            if (D2 < dcut2_) {
              Ninteractions++; // DEBUG
              CalcHbonds(frameNum, D2, it0->Idx(), it1->Idx(), currentFrame, numHB, trajoutNum);
              //mprintf("DBG: %12s %12s %12.4f\n", plNames_[it0->Idx()].c_str(), plNames_[it1->Idx()].c_str(), sqrt(D2));
              //mprintf("DBG: %i %s to %i %s %g\n", plMask_[it0->Idx()]+1, TypeStr_[plTypes_[it0->Idx()]],
              //                                  plMask_[it1->Idx()]+1, TypeStr_[plTypes_[it1->Idx()]], sqrt(D2));

            }
          }
        } // END loop over all other atoms in thisCell
        // Loop over all neighbor cells
        for (unsigned int nidx = 1; nidx != cellList.size(); nidx++)
        {
          PairList::CellType const& nbrCell = pairList_.Cell( cellList[nidx] );
          // Translate vector for neighbor cell
          Vec3 const& tVec = pairList_.TransVec( transList[nidx] );
          // Loop over every atom in nbrCell
          for (PairList::CellType::const_iterator it1 = nbrCell.begin();
                                                  it1 != nbrCell.end(); ++it1)
          {
            if (validInteraction(plTypes_[it0->Idx()], plTypes_[it1->Idx()]))
            {
              Vec3 const& xyz1 = it1->ImageCoords();
              Vec3 dxyz = xyz1 + tVec - xyz0;
              double D2 = dxyz.Magnitude2();
              if (D2 < dcut2_) {
                Ninteractions++; // DEBUG
                CalcHbonds(frameNum, D2, it0->Idx(), it1->Idx(), currentFrame, numHB, trajoutNum);
                //mprintf("DBG: %12s %12s %12.4f\n", plNames_[it0->Idx()].c_str(), plNames_[it1->Idx()].c_str(), sqrt(D2));
                //mprintf("DBG: %i %s to %i %s %g\n", plMask_[it0->Idx()]+1, TypeStr_[plTypes_[it0->Idx()]],
                //                                    plMask_[it1->Idx()]+1, TypeStr_[plTypes_[it1->Idx()]], sqrt(D2));
              }
            }
          } // END loop over atoms in neighbor cell
        } // END loop over neighbor cells
      } // END loop over atoms in thisCell
    } // END cell not empty
  } // END loop over cells
# ifdef _OPENMP
  } // END omp parallel
# endif
  //ConsolidateProblems();
  mprintf("DEBUG: %i interactions.\n", Ninteractions);
  return 0;
}

