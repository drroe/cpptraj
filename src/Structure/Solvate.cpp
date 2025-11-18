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

/** Atom default radius in Angstroms from LEaP */
const double Solvate::ATOM_DEFAULT_RADIUS_ = 1.5;

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
  // TODO principal align

  // Set vdw bounding box
  double Xmin = 0;
  double Ymin = 0;
  double Zmin = 0;
  double Xmax = 0;
  double Ymax = 0;
  double Zmax = 0;

  for (int at = 0; at < topOut.Natom(); at++)
  {
    // Get radius
    double atom_radius = 0.0;
    //bool has_vdw = false;
    if (topOut[at].HasType()) {
      ParmHolder<AtomType>::const_iterator it = set0.AT().GetParam( TypeNameHolder(topOut[at].Type()) );
      if (it != set0.AT().end() && it->second.HasLJ()) {
        atom_radius = it->second.LJ().Radius();
        //has_vdw = true;
      }
    }
    if (atom_radius < 0.1) {
      if (topOut[at].Element() == Atom::HYDROGEN)
        atom_radius = 1.0;
      else
        atom_radius = ATOM_DEFAULT_RADIUS_;
    }
    //mprintf("DEBUG: Atom %s has_vdw= %i VDW=%f\n", topOut.AtomMaskName(at).c_str(), (int)has_vdw, atom_radius);

    const double* XYZ = frameOut.XYZ(at);
    double dXp = XYZ[0] + atom_radius;
    double dYp = XYZ[1] + atom_radius;
    double dZp = XYZ[2] + atom_radius;
    double dXm = XYZ[0] - atom_radius;
    double dYm = XYZ[1] - atom_radius;
    double dZm = XYZ[2] - atom_radius;
    if (at == 0) {
      Xmin = dXm;
      Ymin = dYm;
      Zmin = dZm;
      Xmax = dXp;
      Ymax = dYp;
      Zmax = dZp;
    } else {
      if (dXm < Xmin) Xmin = dXm;
      if (dYm < Ymin) Ymin = dYm;
      if (dZm < Zmin) Zmin = dZm;
      if (dXp > Xmax) Xmax = dXp;
      if (dYp > Ymax) Ymax = dYp;
      if (dZp > Zmax) Zmax = dZp;
    }
  }

  // Define box
  double boxX = Xmax - Xmin;
  double boxY = Ymax - Ymin;
  double boxZ = Zmax - Zmin;
  mprintf("  Solute vdw bounding box:              %-5.3lf %-5.3lf %-5.3lf\n", boxX, boxY, boxZ);

  // Define center
  Vec3 toCenter( -(Xmin + 0.5 * boxX),
                 -(Ymin + 0.5 * boxY),
                 -(Zmin + 0.5 * boxZ) );
  // Translate to origin
  frameOut.Translate(toCenter);

  

  return 0;
}

