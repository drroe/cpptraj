#include "HbCalc.h"
#include <cmath>
#include <algorithm>
#include <utility>
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"

using namespace Cpptraj::HB;

HbCalc::HbCalc() :
  dcut2_(0)
{}

bool HbCalc::IsFON( Atom const& at ) {
  return ( at.Element() == Atom::OXYGEN ||
           at.Element() == Atom::NITROGEN ||
           at.Element() == Atom::FLUORINE );
}

const char* HbCalc::TypeStr_[] = {
  "Hydrogen",
  "Solute Donor",
  "Solute Acceptor",
  "Solute Both",
  "Solvent Donor",
  "Solvent Acceptor",
  "Solvent Both"
};

/** Initialize */
int HbCalc::InitHbCalc(ArgList& argIn, int debugIn) {
  double dcut = argIn.getKeyDouble("dist",3.0);
  dcut = argIn.getKeyDouble("distance", dcut); // for PTRAJ compat.
  dcut2_ = dcut * dcut;

  generalMask_.SetMaskString( argIn.GetMaskNext() );

  pairList_.InitPairList( 8.0, 0.1, debugIn );

  return 0;
}

/** Print current options */
void HbCalc::PrintHbCalcOpts() const {
  mprintf("\tSearching for atoms in mask '%s'\n", generalMask_.MaskString());
  mprintf("\tHeavy atom distance cutoff= %g Ang.\n", sqrt(dcut2_));
}

/** Set up calculation */
int HbCalc::SetupHbCalc(Topology const& topIn, Box const& boxIn) {
  if (setupPairlistAtomMask( topIn )) return 1;

  if (pairList_.SetupPairList( boxIn )) return 1;

  return 0;
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
  typedef std::pair<int, Type> Ptype;
  typedef std::vector<Ptype> Parray;
  Parray IdxTypes;
  IdxTypes.reserve( generalMask_.Nselected() ); // TODO not reserve?

  plMask_ = AtomMask( std::vector<int>(), topIn.Natom() );
  plTypes_.clear();

  int nsites = 0;
  Both_.clear();
  Acceptor_.clear();
  //Iarray Donor;
  for (AtomMask::const_iterator at = generalMask_.begin(); at != generalMask_.end(); ++at) {
    Atom const& currentAtom = topIn[*at];
    if (IsFON( currentAtom )) {
      Iarray h_atoms;
      for (Atom::bond_iterator bat = currentAtom.bondbegin(); bat != currentAtom.bondend(); ++bat) {
        if (topIn[*bat].Element() == Atom::HYDROGEN) {
          h_atoms.push_back( *bat );
          IdxTypes.push_back( Ptype(*bat, HYDROGEN) );
        }
      }
      int molnum = currentAtom.MolNum();
      Type currentType;
      if ( topIn.Mol(molnum).IsSolvent()) {
        if (h_atoms.empty()) {
          currentType = VACCEPTOR;
          Acceptor_.push_back( *at );
        } else {
          currentType = VBOTH;
          Both_.push_back( Site(*at, h_atoms) );
        }
      } else {
        if (h_atoms.empty()) {
          currentType = ACCEPTOR;
          Acceptor_.push_back( *at );
        } else {
          currentType = BOTH;
          Both_.push_back( Site(*at, h_atoms) );
        }
      }
      nsites++;
      IdxTypes.push_back( Ptype(*at, currentType) );
    }
  }
  mprintf("\tNumber of heavy atom sites: %i\n", nsites);

  std::sort( IdxTypes.begin(), IdxTypes.end() );

  plTypes_.reserve( IdxTypes.size() );
  for (Parray::const_iterator it = IdxTypes.begin(); it != IdxTypes.end(); ++it) {
    plMask_.AddSelectedAtom( it->first );
    plTypes_.push_back( it->second );
  }

  for (int idx = 0; idx != plMask_.Nselected(); idx++) {
    //mprintf("\t%8i %4s %s\n", plMask_[idx]+1, *(topIn[plMask_[idx]].Name()), TypeStr_[plTypes_[idx]]);
    mprintf("\t%8i", plMask_[idx]+1);
    mprintf(" %4s", *(topIn[plMask_[idx]].Name()));
    mprintf(" %s\n", TypeStr_[plTypes_[idx]]);
  }

  return 0;
}

/** Place atom into cell TODO put into pairlist */
int HbCalc::GridAtom(HbCellArray& cells, int atomIdx, Vec3 const& frac, Vec3 const& cart) const {
  int i1 = (int)((frac[0]) * (double)pairList_.NX());
  int i2 = (int)((frac[1]) * (double)pairList_.NY());
  int i3 = (int)((frac[2]) * (double)pairList_.NZ());
  int idx = (i3*pairList_.NX()*pairList_.NY())+(i2*pairList_.NX())+i1;
# ifdef DEBUG_PAIRLIST
  mprintf("DBGPL: Atom %6i assigned to cell %6i %10.5f%10.5f%10.5f\n", atomIdx+1, idx, frac[0], frac[1], frac[2]);
# endif
  if (idx < 0 || idx >= (int)cells.size()) {
    // This can happen for e.g. NaN coords
    //mprinterr("Internal Error: Grid %i is out of range (>= %zu || < 0)\n",
    //          idx, cells_.size());
    return 1;
  }
  cells[idx].push_back( PairList::AtmType(atomIdx, cart) );

  return 0;
}

/** Orthogonal imaging TODO put into PairList */
int HbCalc::grid_orthogonal(HbCellArray& cells, int idx, const double* XYZ, Matrix_3x3 const& ucell,  Matrix_3x3 const& recip) const {
      Vec3 fc( XYZ[0]*recip[0],    XYZ[1]*recip[4],    XYZ[2]*recip[8]   );
      Vec3 fcw(fc[0]-floor(fc[0]), fc[1]-floor(fc[1]), fc[2]-floor(fc[2]));
      Vec3 ccw(fcw[0]*ucell[0],    fcw[1]*ucell[4],    fcw[2]*ucell[8]   );
#     ifdef DEBUG_PAIRLIST
      mprintf("DBG: o %6i fc=%7.3f%7.3f%7.3f  fcw=%7.3f%7.3f%7.3f  ccw=%7.3f%7.3f%7.3f\n",
              *atom+1, fc[0], fc[1], fc[2], fcw[0], fcw[1], fcw[2], ccw[0], ccw[1], ccw[2]);
#     endif
      return GridAtom( cells, idx, fcw, ccw );
}

/** Non-orthogonal imaging TODO put into PairList */
int HbCalc::grid_nonOrthogonal(HbCellArray& cells, int idx, const double* XYZ, Matrix_3x3 const& ucell,  Matrix_3x3 const& recip) const {
      Vec3 fc = recip * Vec3(XYZ);
      Vec3 fcw(fc[0]-floor(fc[0]), fc[1]-floor(fc[1]), fc[2]-floor(fc[2]));
      Vec3 ccw = ucell.TransposeMult( fcw );
#     ifdef DEBUG_PAIRLIST
      mprintf("DBG: n %6i fc=%7.3f%7.3f%7.3f  fcw=%7.3f%7.3f%7.3f  ccw=%7.3f%7.3f%7.3f\n",
              *atom+1, fc[0], fc[1], fc[2], fcw[0], fcw[1], fcw[2], ccw[0], ccw[1], ccw[2]);
#     endif
      return GridAtom( cells, idx, fcw, ccw );
}

/** Place donor/acceptor sites and acceptor sites into their respective grids. */
int HbCalc::PlaceSitesOnGrid(Frame const& frmIn, Matrix_3x3 const& ucell,
                             Matrix_3x3 const& recip)
{
  for (HbCellArray::iterator cell = bothCells_.begin(); cell != bothCells_.end(); ++cell)
    cell->clear();
  for (HbCellArray::iterator cell = acceptorCells_.begin(); cell != acceptorCells_.end(); ++cell)
    cell->clear();
  int nOffGrid = 0;
  if (frmIn.BoxCrd().Is_X_Aligned_Ortho()) {
    // Orthogonal imaging
    for (Sarray::const_iterator site = Both_.begin(); site != Both_.end(); ++site)
    {
      const double* XYZ = frmIn.XYZ(site->Idx());
      nOffGrid += grid_orthogonal( bothCells_, site->Idx(), XYZ, ucell, recip );
    }
    for (Iarray::const_iterator it = Acceptor_.begin(); it != Acceptor_.end(); ++it)
    {
      const double* XYZ = frmIn.XYZ( *it );
      nOffGrid += grid_orthogonal( acceptorCells_, *it, XYZ, ucell, recip );
    }
  } else {
    // Non-orthogonal imaging
    for (Sarray::const_iterator site = Both_.begin(); site != Both_.end(); ++site)
    {
      const double* XYZ = frmIn.XYZ(site->Idx());
      nOffGrid += grid_nonOrthogonal( bothCells_, site->Idx(), XYZ, ucell, recip );
    }
    for (Iarray::const_iterator it = Acceptor_.begin(); it != Acceptor_.end(); ++it)
    {
      const double* XYZ = frmIn.XYZ( *it );
      nOffGrid += grid_nonOrthogonal( acceptorCells_, *it, XYZ, ucell, recip );
    }
  }
  return nOffGrid;
}

/** HB calc loop with a pairlist */
int HbCalc::RunCalc_PL(Frame const& currentFrame)
{
  //int retVal = pairList_.CreatePairList(currentFrame,
  //                                      currentFrame.BoxCrd().UnitCell(),
  //                                      currentFrame.BoxCrd().FracCell(), plMask_);
  int retVal = pairList_.PreparePairList(currentFrame.BoxCrd().UnitCell(), currentFrame.BoxCrd().RecipLengths());
  if (retVal < 0) {
    mprinterr("Error: Grid setup failed.\n");
    return 1;
  }
  retVal = PlaceSitesOnGrid(currentFrame, currentFrame.BoxCrd().UnitCell(), currentFrame.BoxCrd().FracCell());
  if (retVal < 0) {
    mprinterr("Error: Placing sites on grids failed.\n");
  } else if (retVal > 0) {
   mprintf("Warning: %i atoms are off the grid.\n", retVal);
  }
  //problemAtoms_.clear();

  int Ninteractions = 0; // DEBUG

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
    // cellList contains this cell index and all neighbors.
    PairList::Iarray const& cellList = thisCell.CellList();
    // transList contains index to translation for the neighbor.
    PairList::Iarray const& transList = thisCell.TransList();

    PairList::Aarray const& bothList = bothCells_[ cidx ];
    if (!bothList.empty()) {
      PairList::Aarray const& accList = acceptorCells_[ cidx ];

      // Loop over donor/acceptor atoms of thisCell
      for (PairList::Aarray::const_iterator it0 = bothList.begin();
                                            it0 != bothList.end(); ++it0)
      {
        Vec3 const& xyz0 = it0->ImageCoords();
        // Loop over other donor/acceptor atoms of thisCell
        for (PairList::Aarray::const_iterator it1 = it0 + 1; // FIXME needs to be changed if donor only atoms added
                                              it1 != bothList.end(); ++it1)
        {
          Vec3 const& xyz1 = it1->ImageCoords();
          Vec3 dxyz = xyz1 - xyz0;
          double D2 = dxyz.Magnitude2();
          if (D2 < dcut2_) {
            Ninteractions++; // DEBUG
            mprintf("DBG: %i to %i %g\n", it0->Idx()+1, it1->Idx()+1, sqrt(D2));
          }
        } // END loop over other donor/acceptor atoms of thisCell
        // Loop over acceptor only atoms of thisCell
        for (PairList::Aarray::const_iterator it1 = accList.begin();
                                              it1 != accList.end(); ++it1)
        {
          Vec3 const& xyz1 = it1->ImageCoords();
          Vec3 dxyz = xyz1 - xyz0;
          double D2 = dxyz.Magnitude2();
          if (D2 < dcut2_) {
            Ninteractions++; // DEBUG
            mprintf("DBG: %i to %i %g\n", it0->Idx()+1, it1->Idx()+1, sqrt(D2));
          }
        } // END loop over acceptor only atoms of thisCell

        // Loop over all neighbor cells
        for (unsigned int nidx = 1; nidx != cellList.size(); nidx++)
        {
          PairList::Aarray const& nbrBothList = bothCells_[ cellList[nidx] ];
          PairList::Aarray const& nbrAccList = acceptorCells_[ cellList[nidx] ];
          // Translate vector for neighbor cell
          Vec3 const& tVec = pairList_.TransVec( transList[nidx] );
          // Loop over donor/acceptor atoms of neighbor cell
          for (PairList::Aarray::const_iterator it1 = nbrBothList.begin();
                                                it1 != nbrBothList.end(); ++it1)
          {
            Vec3 const& xyz1 = it1->ImageCoords();
            Vec3 dxyz = xyz1 + tVec - xyz0;
            double D2 = dxyz.Magnitude2();
            if (D2 < dcut2_) {
              Ninteractions++; // DEBUG
              mprintf("DBG: %i to %i %g\n", it0->Idx()+1, it1->Idx()+1, sqrt(D2));
            }
          } // END loop over donor/acceptor atoms in neighbor cell
          // Loop over acceptor only atoms of neighbor cell
          for (PairList::Aarray::const_iterator it1 = nbrAccList.begin();
                                                it1 != nbrAccList.end(); ++it1)
          {
            Vec3 const& xyz1 = it1->ImageCoords();
            Vec3 dxyz = xyz1 + tVec - xyz0;
            double D2 = dxyz.Magnitude2();
            if (D2 < dcut2_) {
              Ninteractions++; // DEBUG
              mprintf("DBG: %i to %i %g\n", it0->Idx()+1, it1->Idx()+1, sqrt(D2));
            }
          } // END loop over acceptor only atoms of neighbor cell
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

