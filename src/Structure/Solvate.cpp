#include "Solvate.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../Frame.h"
#include "../Topology.h"
#include "../Parm/ParameterSet.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Solvate::Solvate() :
  debug_(0)
{
}

/** Initialize arguments. */
int Solvate::InitSolvate(ArgList& argIn, int debugIn) {
  debug_ = debugIn;

  return 0;
}

/** Create box, fill with solvent */
int Solvate::SolvateBox(Topology& topOut, Frame& frameOut, Cpptraj::Parm::ParameterSet const& set0) {
  using namespace Cpptraj::Parm;
  // Sanity check
  if (topOut.Natom() != frameOut.Natom()) {
    mprinterr("Internal Error: Solvate::SolvateBox(): Topology %s #atoms %i != frame #atoms %i\n",
              topOut.c_str(), topOut.Natom(), frameOut.Natom());
    return 1;
  }
  // TODO Remove any existing box info?
  //if (frameOut.BoxCrd.HasBox())

  // Set vdw bounding box
  double Xmin, Ymin, Zmin;
  double Xmax, Ymax, Zmax;

  for (int at = 0; at < topOut.Natom(); at++)
  {
    // Get radius
    double atom_radius = 0.0;
    bool has_vdw = false;
    if (topOut[at].HasType()) {
      ParmHolder<AtomType>::const_iterator it = set0.AT().GetParam( TypeNameHolder(topOut[at].Type()) );
      if (it != set0.AT().end() && it->second.HasLJ()) {
        atom_radius = it->second.LJ().Radius();
        has_vdw = true;
      }
    }
    mprintf("DEBUG: Atom %s has_vdw= %i VDW=%f\n", topOut.AtomMaskName(at).c_str(), (int)has_vdw, atom_radius);
  }

  return 0;
}

