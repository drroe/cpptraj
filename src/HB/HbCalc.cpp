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

  pairList_.InitPairList( dcut, 1.0, debugIn );

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

/** HB calc loop with a pairlist */
int HbCalc::RunCalc_PL(Frame const& currentFrame)
{
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

  int cidx;
# ifdef _OPENMP
  int mythread;
# pragma omp parallel private(cidx,mythread) 
  {
  mythread = omp_get_thread_num();
  thread_problemAtoms_[mythread].clear();
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
        if (plTypes_[it0->Idx()] == HYDROGEN) continue;
        Vec3 const& xyz0 = it0->ImageCoords();
        // Exclusion list for this atom
        //ExclusionArray::ExListType const& excluded = Excluded_[it0->Idx()];
        // Calc interaction of atom to all other atoms in thisCell.
        for (PairList::CellType::const_iterator it1 = it0 + 1;
                                                it1 != thisCell.end(); ++it1)
        {
          if (plTypes_[it1->Idx()] == HYDROGEN) continue;
          Vec3 const& xyz1 = it1->ImageCoords();
          Vec3 dxyz = xyz1 - xyz0;
          double D2 = dxyz.Magnitude2();
          if (D2 < dcut2_) {
            Ninteractions++; // DEBUG
            mprintf("DBG: %i %s to %i %s %g\n", plMask_[it0->Idx()]+1, TypeStr_[plTypes_[it0->Idx()]],
                                                plMask_[it1->Idx()]+1, TypeStr_[plTypes_[it1->Idx()]], sqrt(D2));
/*
#           ifdef _OPENMP
            thread_problemAtoms_[mythread]
#           else
            problemAtoms_
#           endif
              .push_back(Problem(Mask1_[it0->Idx()], Mask1_[it1->Idx()], sqrt(D2)));
*/
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
            if (plTypes_[it1->Idx()] == HYDROGEN) continue;
            Vec3 const& xyz1 = it1->ImageCoords();
            Vec3 dxyz = xyz1 + tVec - xyz0;
            double D2 = dxyz.Magnitude2();
            if (D2 < dcut2_) {
              Ninteractions++; // DEBUG
              mprintf("DBG: %i %s to %i %s %g\n", plMask_[it0->Idx()]+1, TypeStr_[plTypes_[it0->Idx()]],
                                                  plMask_[it1->Idx()]+1, TypeStr_[plTypes_[it1->Idx()]], sqrt(D2));
/*
#             ifdef _OPENMP
              thread_problemAtoms_[mythread]
#             else
              problemAtoms_
#             endif
                .push_back(Problem(Mask1_[it0->Idx()], Mask1_[it1->Idx()], sqrt(D2)));
*/
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

