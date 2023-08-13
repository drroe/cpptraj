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
  "Donor",
  "Acceptor",
  "Both"
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

  for (AtomMask::const_iterator at = generalMask_.begin(); at != generalMask_.end(); ++at) {
    Atom const& currentAtom = topIn[*at];
    if (IsFON( currentAtom )) {
      int nh = 0;
      for (Atom::bond_iterator bat = currentAtom.bondbegin(); bat != currentAtom.bondend(); ++bat) {
        if (topIn[*bat].Element() == Atom::HYDROGEN) {
          nh++;
          IdxTypes.push_back( Ptype(*bat, HYDROGEN) );
        }
      }
      Type currentType;
      if (nh == 0)
        currentType = ACCEPTOR;
       else
        currentType = BOTH;
      IdxTypes.push_back( Ptype(*at, currentType) );
    }
  }

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
