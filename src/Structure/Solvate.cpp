#include "Solvate.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DataSet_Coords.h"
#include "../DataSetList.h"
#include "../Frame.h"
#include "../Topology.h"
#include "../Parm/ParameterSet.h"
#include <algorithm> //std::max

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Solvate::Solvate() :
  debug_(0),
  bufferX_(0),
  bufferY_(0),
  bufferZ_(0),
  isotropic_(false),
  clip_(true)
{
}

/** Initialize arguments. */
int Solvate::InitSolvate(ArgList& argIn, int debugIn) {
  debug_ = debugIn;

  if (argIn.Contains("buffer")) {
    bufferX_ = argIn.getKeyDouble("buffer", -1.0);
    bufferY_ = bufferX_;
    bufferZ_ = bufferX_;
  } else {
    bufferX_ = argIn.getKeyDouble("bufx", -1.0);
    bufferY_ = argIn.getKeyDouble("bufy", -1.0);
    bufferZ_ = argIn.getKeyDouble("bufz", -1.0);
  }
  if (bufferX_ < 0 || bufferY_ < 0 || bufferZ_ < 0) {
    mprinterr("Error: Either 'buffer' or 'bufx/bufy/bufx' must be specified and >= 0\n");
    return 1;
  }

  isotropic_ = argIn.hasKey("iso");

  solventBoxName_ = argIn.GetStringKey("solventbox");
  if (solventBoxName_.empty()) {
    mprinterr("Error: Specify solvent box unit name with 'solventbox'\n");
    return 1;
  }

  return 0;
}

/** Get solvent unit box from DataSetList */
DataSet_Coords* Solvate::GetSolventUnit(DataSetList const& DSL) const {
  if (solventBoxName_.empty()) {
    mprinterr("Internal Error: Solvate::GetSolventUnit() called before solventBoxName_ set.\n");
    return 0;
  }
  DataSetList sets = DSL.SelectGroupSets( "*", DataSet::COORDINATES ); // TODO specific set type for units?
  // First try to match aspect, then match name
  DataSet_Coords* solventUnit = 0;
  // Aspect
  for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
  {
    DataSet_Coords* ds = static_cast<DataSet_Coords*>( *it );
    if (!ds->Meta().Aspect().empty()) {
      if (solventBoxName_ == ds->Meta().Aspect()) {
        solventUnit = ds;
        break;
      }
    }
  }
  // Name
  if (solventUnit == 0) {
    for (DataSetList::const_iterator it = sets.begin(); it != sets.end(); ++it)
    {
      DataSet_Coords* ds = static_cast<DataSet_Coords*>( *it );
      if (solventBoxName_ == ds->Meta().Name()) {
        solventUnit = ds;
        break;
      }
    }
  }
  if (solventUnit != 0)
    mprintf("\tSolvent unit: %s\n", solventUnit->legend());
  else
    mprinterr("Error: Could not get solvent unit named %s\n", solventBoxName_.c_str());

  return solventUnit;
}

/** Atom default radius in Angstroms from LEaP */
const double Solvate::ATOM_DEFAULT_RADIUS_ = 1.5;

/** Set VDW bounding box. */
int Solvate::setVdwBoundingBox(double& boxX, double& boxY, double& boxZ,
                               Topology const& topOut, Frame& frameOut,
                               Cpptraj::Parm::ParameterSet const& set0)
const
{
  using namespace Cpptraj::Parm;
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
    bool has_vdw = false;
    if (topOut[at].HasType()) {
      ParmHolder<AtomType>::const_iterator it = set0.AT().GetParam( TypeNameHolder(topOut[at].Type()) );
      if (it != set0.AT().end() && it->second.HasLJ()) {
        atom_radius = it->second.LJ().Radius();
        has_vdw = true;
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
    //mprintf("DBG: %12.4f %12.4f %12.4f %12.4f\n", XYZ[0], XYZ[1], XYZ[2], atom_radius);
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
  mprintf("Min= %f %f %f  Max= %f %f %f\n", Xmin, Ymin, Zmin, Xmax, Ymax, Zmax);
  // Define box
  boxX = Xmax - Xmin;
  boxY = Ymax - Ymin;
  boxZ = Zmax - Zmin;

  // Define center
  Vec3 toCenter( -(Xmin + 0.5 * boxX),
                 -(Ymin + 0.5 * boxY),
                 -(Zmin + 0.5 * boxZ) );
  // Translate to origin
  frameOut.Translate(toCenter);

  return 0;
}

/** Create box, fill with solvent */
int Solvate::SolvateBox(Topology& topOut, Frame& frameOut, Cpptraj::Parm::ParameterSet const& set0,
                        DataSet_Coords& SOLVENTBOX)
const
{
  // Sanity check
  if (topOut.Natom() != frameOut.Natom()) {
    mprinterr("Internal Error: Solvate::SolvateBox(): Topology %s #atoms %i != frame #atoms %i\n",
              topOut.c_str(), topOut.Natom(), frameOut.Natom());
    return 1;
  }
  // TODO Remove any existing box info?
  //if (frameOut.BoxCrd.HasBox())
  // TODO principal align

  // Set vdw box
  double boxX, boxY, boxZ;
  if (setVdwBoundingBox(boxX, boxY, boxZ, topOut, frameOut, set0)) {
    mprinterr("Error: Setting vdw bounding box for %s failed.\n", topOut.c_str());
    return 1;
  }
  mprintf("  Solute vdw bounding box:              %-5.3f %-5.3f %-5.3f\n", boxX, boxY, boxZ);

  double dXWidth = boxX + bufferX_ * 2;
  double dYWidth = boxY + bufferY_ * 2;
  double dZWidth = boxZ + bufferZ_ * 2;

  if (isotropic_) {
    double dTemp = dXWidth * dYWidth * dZWidth;

    double dMax = std::max(dXWidth, dYWidth);
    dMax = std::max(dMax, dZWidth);
    dXWidth = dYWidth = dZWidth = dMax;

    dTemp = (dMax * dMax * dMax - dTemp ) / dTemp;

    mprintf("  Total bounding box for atom centers:  %5.3f %5.3f %5.3f\n", 
            dXWidth, dYWidth, dZWidth );
    mprintf("      (box expansion for 'iso' is %5.1lf%%)\n", dTemp * 100.0 );

     // To make the actual clip right, 'iso' the solute box
    dTemp = std::max(boxX, boxY);
    dTemp = std::max(dTemp, boxZ);
    //dXBox = dYBox = dZBox = dTemp;
    boxX = boxY = boxZ = dTemp;
  } else
    mprintf("  Total bounding box for atom centers:  %5.3f %5.3f %5.3f\n", 
            dXWidth, dYWidth, dZWidth );

  /*if ( bClip ) {
        //  If the solvated system should be clipped to the exact
        //      size the user specified then note the criterion
        //      & dimensions (for 0,0,0-centered system)
        iCriteria |= TOOLOUTSIDEOFBOX;
        cCriteria.dX = 0.5 * dXBox + dXW;
        cCriteria.dY = 0.5 * dYBox + dYW;
        cCriteria.dZ = 0.5 * dZBox + dZW;
    }
    if ( bOct ) {
        // maybe allow oct clip on integer boxes someday.. but for now:
        if ( !bClip )
                DFATAL(( "oct but no clip\n" ));
        iCriteria |= TOOLOUTSIDEOFOCTBOX;
    }
*/
  // TODO check COORDS size
  Frame solventFrame = SOLVENTBOX.AllocateFrame();
  SOLVENTBOX.GetFrame(0, solventFrame);

  // Set vdw box for solvent
  double solventX, solventY, solventZ;
  if (solventFrame.BoxCrd().HasBox()) {
    // TODO check ortho?
    solventX = solventFrame.BoxCrd().Param(Box::X);
    solventY = solventFrame.BoxCrd().Param(Box::Y);
    solventZ = solventFrame.BoxCrd().Param(Box::Z);
  } else {
    if (setVdwBoundingBox(solventX, solventY, solventZ, SOLVENTBOX.Top(), solventFrame, set0)) {
      mprinterr("Error: Setting vdw bounding box for %s failed.\n", topOut.c_str());
      return 1;
    }
  }
  mprintf("  Solvent unit box:                     %5.3f %5.3f %5.3f\n", solventX, solventY, solventZ);

  // See how many solvent boxes required in each dimension

  int iX = (int)( dXWidth / solventX ) + 1;
  int iY = (int)( dYWidth / solventY ) + 1;
  int iZ = (int)( dZWidth / solventZ ) + 1;

  //  Calculate the center of the first solvent box 
  //  (the one that goes in the max XYZ corner), given
  //  that the solute is centered at 0,0,0
  double dXStart = 0.5 * solventX * (double) (iX-1);
  double dYStart = 0.5 * solventY * (double) (iY-1);
  double dZStart = 0.5 * solventZ * (double) (iZ-1);

             /* If the caller wants a solvent shell then */
             /* make sure that the box used to find interesting solute */
             /* spheres takes into account the dFarness parameter */
             /* so that there are at least some solute spheres in */
             /* the interesting list to check against solvent */

 //if ( bShell ) 
 //    dBuffer = dFarness;
 //else 
 //    dBuffer = 0.0;

  addSolventUnits(iX, iY, iZ, dXStart, dYStart, dZStart, solventX, solventY, solventZ,
                  solventFrame, SOLVENTBOX.Top(), frameOut, topOut);

  return 0;
}

int Solvate::addSolventUnits(int numX, int numY, int numZ,
                             double dXStart, double dYStart, double dZStart,
                             double dXSolvent, double dYSolvent, double dZSolvent,
                             Frame& solventFrame, Topology const& solventTop,
                             Frame& frameOut, Topology& topOut)
const
{

  mprintf( "The number of boxes:  x=%2d  y=%2d  z=%2d\n", numX, numY, numZ );
  int NboxesToAdd = numX * numY * numZ;
  int NatomsToAdd = NboxesToAdd * solventTop.Natom();
  mprintf("Will add %i boxes, %i atoms.\n", NboxesToAdd, NatomsToAdd);
  frameOut.IncreaseX( NatomsToAdd );

  // Current solvent unit center
  Vec3 currentSolventCenter = solventFrame.VGeometricCenter();

  int firstSolventAtom = topOut.Natom();
  mprintf("DEBUG: First solvent atom is %i\n", firstSolventAtom+1);

  std::vector<int> bondedAtoms;
  bondedAtoms.reserve(12); // Reserve for 6 bonds

  double dX = dXStart;
  for ( int ix=0; ix < numX; ix++, dX -= dXSolvent ) {
    double dY = dYStart;
    for ( int iy=0; iy < numY; iy++, dY -= dYSolvent ) {
      double dZ = dZStart;
      for ( int iz=0; iz < numZ; iz++, dZ -= dZSolvent ) {
        mprintf( "Adding box at: x=%d  y=%d  z=%d\n", ix, iy, iz);

        mprintf( "Center of solvent box is: %lf, %lf, %lf\n",
                                dX, dY, dZ );
        Vec3 trans( dX - currentSolventCenter[0],
                    dY - currentSolventCenter[1],
                    dZ - currentSolventCenter[2] );
        solventFrame.Translate(trans);
        Vec3 debugVec = solventFrame.VGeometricCenter();
        debugVec.Print("DEBUG: check solvent center");
        currentSolventCenter[0] = dX;
        currentSolventCenter[1] = dY;
        currentSolventCenter[2] = dZ;
        // Add solvent unit to output topology for this cube
        int atomOffset = topOut.Natom();
        int currentResNum = topOut.Nres();
        for (int ires = 0; ires != solventTop.Nres(); ires++) {
          Residue solventRes = solventTop.Res(ires);
          solventRes.SetOriginalNum( currentResNum + ires );
          bondedAtoms.clear();
          for (int iat = solventRes.FirstAtom(); iat != solventRes.LastAtom(); iat++)
          {
            Atom solventAtom = solventTop[iat];
            // Save solvent bonds
            for (Atom::bond_iterator bat = solventAtom.bondbegin(); bat != solventAtom.bondend(); ++bat) {
              if (*bat > iat) {
                bondedAtoms.push_back(  iat + atomOffset );
                bondedAtoms.push_back( *bat + atomOffset );
              }
            }
            solventAtom.ClearBonds(); // FIXME AddTopAtom should clear
            topOut.AddTopAtom( solventAtom, solventRes );
          } // END loop over solvent atoms
          // Add bonds
          for (std::vector<int>::const_iterator it = bondedAtoms.begin(); it != bondedAtoms.end(); ++it) {
            int at0 = *it;
            ++it;
            topOut.AddBond( at0, *it );
          }
        } // END loop over solvent unit residues
        // Append solvent frame
        frameOut.AppendFrame( solventFrame );
      } // END loop over Z
    } // END loop over Y
  } // END loop over X

  return 0;
}
