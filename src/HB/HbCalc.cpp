#include "HbCalc.h"
#include <cmath> //sqrt
#include "../ArgList.h"
#include "../Constants.h"
#include "../CpptrajStdio.h"
#include "../DataFileList.h"
#include "../DistRoutines.h"
#include "../Topology.h"
#ifdef _OPENMP
# include <omp.h>
#endif

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

/** Set debug level. */
void HbCalc::SetDebug(int debugIn) {
  hbdata_.SetDebug( debugIn );
}

/** Initialize */
int HbCalc::InitHbCalc(ArgList& argIn, DataSetList* masterDslPtr, DataFileList& DFL, int debugIn) {
  double dcut = argIn.getKeyDouble("dist",3.0);
  dcut = argIn.getKeyDouble("distance", dcut); // for PTRAJ compat.
  dcut2_ = dcut * dcut;
  acut_ = argIn.getKeyDouble("angle", 135.0 * Constants::DEGRAD);

  if (hbdata_.ProcessArgs(argIn, DFL)) {
    mprinterr("Error: Could not process hydrogen bond data args.\n");
    return 1;
  } 

  generalMask_.SetMaskString( argIn.GetMaskNext() );

  // Setup datasets
  std::string hbsetname = argIn.GetStringNext();

  // ----- All arguments should be processed now. ----------
  pairList_.InitPairList( 8.0, 0.1, debugIn );

  if (hbdata_.InitHbData( masterDslPtr, hbsetname )) {
    mprinterr("Error: Could not initialize hydrogen bond data.\n");
    return 1;
  }
# ifdef _OPENMP
  // Each thread needs temp. space to store found hbonds every frame
  // to avoid memory clashes when adding/updating in map.
# pragma omp parallel
  {
# pragma omp master
  {
  thread_HBs_.resize( omp_get_num_threads() );
  }
  }
# endif

  return 0;
}

/** Print current options */
void HbCalc::PrintHbCalcOpts() const {
# ifdef _OPENMP
  if (thread_HBs_.size() > 1)
    mprintf("\tParallelizing calculation with %zu threads.\n", thread_HBs_.size());
# endif
  mprintf("\tSearching for atoms in mask '%s'\n", generalMask_.MaskString());
  mprintf("\tHeavy atom distance cutoff= %g Ang.\n", sqrt(dcut2_));
  if (acut_ > -1)
    mprintf("\tAngle cutoff= %g deg.\n", acut_*Constants::RADDEG);
  else
    mprintf("\tNo angle cutoff.\n");
  hbdata_.PrintHbDataOpts();
}

/** Set up calculation */
int HbCalc::SetupHbCalc(Topology const& topIn, Box const& boxIn) {
  if (setupPairlistAtomMask( topIn )) return 1;

  if (pairList_.SetupPairList( boxIn )) return 1;

  // For backwards compatibility, if saving time series we need to store
  // the donor hydrogen indices and acceptor heavy atom indices.
  // We also need to do this if we are saving an interaction matrix.
  Iarray acceptorOnly_indices, donor_h_indices, both_indices, donorOnly_indices;
  if (hbdata_.Series() || hbdata_.InteractionMatrix()) {
    for (int idx = 0; idx != plMask_.Nselected(); idx++) {
      if (plTypes_[idx] == BOTH) {
        both_indices.push_back( plMask_[idx] );
        //acceptor_indices.push_back( plMask_[idx] );
        for (Iarray::const_iterator ht = plHatoms_[idx].begin(); ht != plHatoms_[idx].end(); ++ht)
          donor_h_indices.push_back( *ht );
      } else if (plTypes_[idx] == ACCEPTOR) {
        acceptorOnly_indices.push_back( plMask_[idx] );
      } else if (plTypes_[idx] == DONOR) {
        donorOnly_indices.push_back( plMask_[idx] );
        for (Iarray::const_iterator ht = plHatoms_[idx].begin(); ht != plHatoms_[idx].end(); ++ht)
          donor_h_indices.push_back( *ht );
      }
    }
  }
  if (hbdata_.SetCurrentParm( &topIn, both_indices,
                              donorOnly_indices, donor_h_indices,
                              acceptorOnly_indices ))
  {
    mprinterr("Error: Could not set up hbond Topology data.\n");
    return 1;
  }

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
  plId_.clear();
//  plNames_.clear();
  plHatoms_.clear();

  int nsites = 0;
  unsigned int NacceptorOnly = 0;
  unsigned int Nboth = 0;
  unsigned int NdonorOnly = 0;
  unsigned int NumH = 0;
  unsigned int NV_acceptorOnly = 0;
  unsigned int NV_both = 0;
  unsigned int NV_donorOnly = 0;
  unsigned int NV_H = 0;
  unsigned int NIons = 0;
  for (AtomMask::const_iterator at = generalMask_.begin(); at != generalMask_.end(); ++at) {
    Atom const& currentAtom = topIn[*at];
    int molnum = currentAtom.MolNum();
    bool isSolventAtom = topIn.Mol(molnum).IsSolvent();
    // If we do not care about solvent hydrogen bonds and this is solvent, skip it.
    if (isSolventAtom && !hbdata_.CalcSolvent()) continue;
    // Determine atom ID
    int atid;
    if (hbdata_.NoIntramol())
      atid = molnum;
    else
      atid = *at;
    if (IsFON( currentAtom )) {
      // Check if there are any hydrogens bonded to this atom
      Iarray h_atoms;
      for (Atom::bond_iterator bat = currentAtom.bondbegin(); bat != currentAtom.bondend(); ++bat) {
        if (topIn[*bat].Element() == Atom::HYDROGEN) {
          h_atoms.push_back( *bat );
        }
      }
      Type currentType;
      if ( isSolventAtom ) {
        // Solvent atom
        if (h_atoms.empty())
          currentType = VACCEPTOR;
        else {
          NV_H += h_atoms.size();
          currentType = VBOTH;
        }
        // Solvent count
        if (currentType == VACCEPTOR)
          NV_acceptorOnly++;
        else if (currentType == VBOTH)
          NV_both++;
        else if (currentType == VDONOR)
          NV_donorOnly++;
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
      plId_.push_back( atid );
//      plNames_.push_back( topIn.TruncResAtomName( *at ) );
      plHatoms_.push_back( h_atoms );
    } else if (hbdata_.CalcSolvent()) {
      if (currentAtom.Nbonds() == 0) {
        // If no bonds to this atom assume it is an ion. Set the H atom
        // to be the same as D atom; this will skip the angle calc.
        plMask_.AddSelectedAtom( *at );
        plTypes_.push_back( VDONOR );
        plId_.push_back( atid );
        // TODO check charge to see if it can be an acceptor?
        plHatoms_.push_back( Iarray(1, *at) );
        NIons++;
        nsites++;
      } // END atom has no bonds
    } 
  }
  mprintf("\tTotal Number of heavy atom sites: %i\n", nsites);
  mprintf("\t  Solute acceptor-only atoms: %u\n", NacceptorOnly);
  mprintf("\t  Solute donor/acceptor sites: %u\n", Nboth);
  mprintf("\t  Solute donor-only sites: %u\n", NdonorOnly);
  mprintf("\t  %u solute hydrogens.\n", NumH);
  if (hbdata_.CalcSolvent()) {
    mprintf("\t  Solvent acceptor-only atoms: %u\n", NV_acceptorOnly);
    mprintf("\t  Solvent donor/acceptor sites: %u\n", NV_both);
    mprintf("\t  Solvent donor-only sites: %u\n", NV_donorOnly);
    mprintf("\t  %u solvent hydrogens, %u ions.\n", NV_H, NIons);
  }

  unsigned int uuSize = NumH * (Nboth + NacceptorOnly);
  unsigned int uvSize = NumH * (NV_both + NV_acceptorOnly + NIons);
  uvSize += (Nboth + NacceptorOnly);
  mprintf("\tEstimated max potential memory usage: %s\n",
          hbdata_.MemoryUsage( uuSize,
                               uvSize,
                               0 ).c_str());
//  std::sort( IdxTypes.begin(), IdxTypes.end() );

//  plTypes_.reserve( IdxTypes.size() );
//  for (Parray::const_iterator it = IdxTypes.begin(); it != IdxTypes.end(); ++it) {
//    plMask_.AddSelectedAtom( it->first );
//    plTypes_.push_back( it->second );
//    plNames_.push_back( topIn.TruncResAtomName(it->first) );
//  }

  if (hbdata_.Debug() > 1) {
    for (int idx = 0; idx != plMask_.Nselected(); idx++) {
      //mprintf("\t%8i %4s %s\n", plMask_[idx]+1, *(topIn[plMask_[idx]].Name()), TypeStr_[plTypes_[idx]]);
      mprintf("\t%8i", plMask_[idx]+1);
      mprintf(" %4s", *(topIn[plMask_[idx]].Name()));
      mprintf(" %s\n", TypeStr_[plTypes_[idx]]);
    }
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
void HbCalc::calc_UU_Hbonds(int frameNum, double dist2,
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
//      mprintf("DBG: %12s %12i %12s %12.4f %12.4f\n", plNames_[a_idx].c_str(), *h_atom + 1, plNames_[d_idx].c_str(), sqrt(dist2), angle*Constants::RADDEG);
#     ifdef _OPENMP
      // numHB holds thread number, will be counted later on. -1 indicates UU hbond.
      thread_HBs_[numHB].push_back( Hbond(sqrt(dist2), angle, a_atom, *h_atom, d_atom, -1) );
#     else
      ++numHB;
      hbdata_.AddUU(sqrt(dist2), angle, frameNum, a_atom, *h_atom, d_atom, trajoutNum);
#     endif
    }
  }
}

/** Calculate hydrogen bonds between solute/solvent pair.
  * The distance cutoff should already be satisfied between donor and
  * acceptor heavy atoms.
  */
void HbCalc::calc_UV_Hbonds(int frameNum, double dist2,
                            int d_idx, Iarray const& Hatoms,
                            int a_idx,
                            Frame const& frmIn, int& numHB, bool soluteDonor,
                            int trajoutNum)
{
  int d_atom = plMask_[d_idx];
  int a_atom = plMask_[a_idx]; 
  // Determine if angle cutoff is satisfied
  for (Iarray::const_iterator h_atom = Hatoms.begin(); h_atom != Hatoms.end(); ++h_atom)
  {
    double angle = 0.0;
    bool angleSatisfied = true;
    // For ions, donor atom will be same as h atom so no angle needed.
    if (d_atom != *h_atom) {
      //angle = Angle(XYZA, frmIn.XYZ(*h_atom), XYZD, frmIn.BoxCrd());
      angle = Angle(frmIn.XYZ(a_atom), frmIn.XYZ(*h_atom), frmIn.XYZ(d_atom), frmIn.BoxCrd());
      angleSatisfied = !(angle < acut_);
    }
    if (angleSatisfied)
    {
#     ifdef _OPENMP
      // numHB holds thread number, will be counted later on.
      thread_HBs_[numHB].push_back( Hbond(sqrt(dist2), angle, a_atom, *h_atom, d_atom, (int)soluteDonor) );
#     else
      ++numHB;
      hbdata_.AddUV(sqrt(dist2), angle, frameNum, a_atom, *h_atom, d_atom, soluteDonor, trajoutNum);
#     endif
    }
  }
}

/** Calculate hbonds between two heavy atoms. */
void HbCalc::CalcHbonds(int frameNum, double dist2,
                        int a0idx,
                        int a1idx,
                        Frame const& frmIn, int& numHB,
                        int trajoutNum)
{
  Type a0type = plTypes_[a0idx];
  Type a1type = plTypes_[a1idx];
  bool a0IsSolvent = (a0type == VBOTH || a0type == VDONOR || a0type == VACCEPTOR);
  bool a1IsSolvent = (a1type == VBOTH || a1type == VDONOR || a1type == VACCEPTOR);
  if (!a0IsSolvent && !a1IsSolvent) {
    // Solute-solute
    // BOTH ACCEPTOR
    // DONOR ACCEPTOR
    // BOTH DONOR
    // ACCEPTOR DONOR
    // BOTH BOTH
    // DONOR BOTH
    // ACCEPTOR BOTH
    if ((a0type == BOTH || a0type == DONOR)    && (a1type == BOTH || a1type == ACCEPTOR)) {
      calc_UU_Hbonds(frameNum, dist2, a0idx, plHatoms_[a0idx], a1idx, frmIn, numHB, trajoutNum);
    } 
    if ((a0type == BOTH || a0type == ACCEPTOR) && (a1type == BOTH || a1type == DONOR)) {
      calc_UU_Hbonds(frameNum, dist2, a1idx, plHatoms_[a1idx], a0idx, frmIn, numHB, trajoutNum);
    }
  } else {
    if (!hbdata_.CalcSolvent()) return;
    // Check solvent-solvent and solvent-solute
    if (a0IsSolvent && a1IsSolvent) {
      return; // FIXME
    } else {
      // Solvent - solute.  By convention, put the solvent atom first.
      int i0, i1;
      Type t0, t1;
      if (a0IsSolvent) {
        i0 = a0idx;
        t0 = a0type;
        i1 = a1idx;
        t1 = a1type;
      } else {
        i0 = a1idx;
        t0 = a1type;
        i1 = a0idx;
        t1 = a0type;
      }
      // VBOTH ACCEPTOR
      // VBOTH DONOR
      // VBOTH BOTH
      // VDONOR ACCEPTOR
      // VDONOR BOTH
      // VACCEPTOR DONOR
      // VACCEPTOR BOTH
      if ((t0 == VBOTH || t0 == VDONOR)  && (t1 == BOTH || t1 == ACCEPTOR)) {
        calc_UV_Hbonds(frameNum, dist2, i0, plHatoms_[i0], i1, frmIn, numHB, false, trajoutNum);
      }
      if ((t0 == VBOTH || t0 == VACCEPTOR) && (t1 == BOTH || t1 == DONOR)) {
        calc_UV_Hbonds(frameNum, dist2, i1, plHatoms_[i1], i0, frmIn, numHB, true, trajoutNum);
      }
    }
  }
}

/** HB calc loop with a pairlist */
int HbCalc::RunCalc_PL(Frame const& currentFrame, int frameNum, int trajoutNum)
{
# ifdef TIMER
  t_action_.Start();
# endif
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
# ifdef TIMER
  t_hbcalc_.Start();
# endif
//  int Ninteractions = 0; // DEBUG
  int numHB = 0;
  int cidx;
# ifdef _OPENMP
//  int mythread;
# pragma omp parallel private(cidx, numHB) 
  {
  // In OpenMP, numHB is used to track thread number
  numHB = omp_get_thread_num();
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
          if (plId_[it0->Idx()] != plId_[it1->Idx()] &&
              validInteraction(plTypes_[it0->Idx()], plTypes_[it1->Idx()]))
          {
            Vec3 const& xyz1 = it1->ImageCoords();
            Vec3 dxyz = xyz1 - xyz0;
            double D2 = dxyz.Magnitude2();
            if (D2 < dcut2_) {
//              Ninteractions++; // DEBUG
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
            if (plId_[it0->Idx()] != plId_[it1->Idx()] && 
                validInteraction(plTypes_[it0->Idx()], plTypes_[it1->Idx()]))
            {
              Vec3 const& xyz1 = it1->ImageCoords();
              Vec3 dxyz = xyz1 + tVec - xyz0;
              double D2 = dxyz.Magnitude2();
              if (D2 < dcut2_) {
//                Ninteractions++; // DEBUG
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
  // Add all found hydrogen bonds
  for (std::vector<Harray>::iterator it = thread_HBs_.begin();
                                     it != thread_HBs_.end(); ++it)
  {
    for (Harray::const_iterator hb = it->begin(); hb != it->end(); ++hb)
    {
      if (hb->Frames() < 0)
        hbdata_.AddUU(hb->Dist(), hb->Angle(), frameNum, hb->A(), hb->H(), hb->D(), trajoutNum);
      else
        hbdata_.AddUV(hb->Dist(), hb->Angle(), frameNum, hb->A(), hb->H(), hb->D(), (bool)hb->Frames(), trajoutNum);
    }
    it->clear();
  }
# endif /* _OPENMP */
//  mprintf("DEBUG: %i interactions.\n", Ninteractions);
# ifdef TIMER
  t_hbcalc_.Stop();
# endif
  hbdata_.IncrementNframes(frameNum, trajoutNum);
# ifdef TIMER
  t_action_.Stop();
# endif
  return 0;
}

/** Finish HB calc and do output. */
void HbCalc::FinishHbCalc() {
  hbdata_.PrintHbData();
# ifdef TIMER
  t_hbcalc_.WriteTiming(2, "Hydrogen Bond Calc. :", t_action_.Total());
  t_action_.WriteTiming(1, "Total :");
# endif
}
