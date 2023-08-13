#include "HbCalc.h"
#include <algorithm>
#include <utility>
#include "../CpptrajStdio.h"
#include "../Topology.h"

using namespace Cpptraj::HB;

HbCalc::HbCalc() {}

bool HbCalc::IsFON( Atom const& at ) {
  return ( at.Element() == Atom::OXYGEN ||
           at.Element() == Atom::NITROGEN ||
           at.Element() == Atom::FLUORINE );
}

/** Set up mask for pair list. */
int HbCalc::SetupPairlistAtomMask(Topology const& topIn) {
  if (topIn.SetupIntegerMask( generalMask_ )) {
    mprinterr("Error: Could not set up mask '%s'\n", generalMask_.MaskString());
    return 1;
  }
  if (generalMask_.None()) {
    mprintf("Warning: No atoms selected by mask '%s' \n", generalMask_.MaskString());
    return 1;
  }
  // Decide what each atom is
  typedef std::pair<int, Type> Ptype;
  typedef std::vector<Ptype> Parray;
  Parray IdxTypes;
  IdxTypes.reserve( generalMask_.Nselected() ); // TODO not reserve?

  plMask_.ClearSelected();
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

  return 0;
}
