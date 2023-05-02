#include "MeadGrid.h"
#include "MeadError.h"
#include "../CpptrajStdio.h"
#include "../Frame.h"
#include "../DistRoutines.h"
// MEAD includes
#include "../../mead/FinDiffMethod.h"
#include "../../mead/MEADexcept.h"
// For debug print
#include <iostream>
#include <cmath>

using namespace Cpptraj::Mead;
/** CONSTRUCTOR */
MeadGrid::MeadGrid() :
  fdm_(0)
{
  levelOpts_.reserve(3);
}

/** DESTRUCTOR */
MeadGrid::~MeadGrid() {
  if (fdm_ != 0) delete fdm_;
}

/** COPY CONSTRUCTOR */
MeadGrid::MeadGrid( MeadGrid const& rhs ) :
  fdm_(0)
{
  fdm_ = new FinDiffMethod();
  for (std::vector<GridOpt>::const_iterator it = rhs.levelOpts_.begin();
                                            it != rhs.levelOpts_.end(); ++it)
    addGrid( *it );
}

/** ASSIGNMENT OPERATOR */
MeadGrid& MeadGrid::operator=(MeadGrid const& rhs) {
  if (this == &rhs) return *this;
  if (fdm_ != 0) delete fdm_;
  levelOpts_.clear();
  fdm_ = new FinDiffMethod();
  for (std::vector<GridOpt>::const_iterator it = rhs.levelOpts_.begin();
                                            it != rhs.levelOpts_.end(); ++it)
    addGrid( *it );
  return *this;
}

/** Corresponds to enum Center_Mode */
const char* MeadGrid::Center_ModeStr_[] = {
  "ON_ORIGIN", "ON_CENT_OF_INTR", "ON_GEOM_CENT", "SPECIFIED"
};

/** \return Char string corresponding to Center_Mode */
const char* MeadGrid::Center_ModeStr(Center_Mode gc) {
  return Center_ModeStr_[gc];
}

/** Add a grid to the FDM object. */
int MeadGrid::addGrid(GridOpt const& go) {
  if (go.mode_ == C_SPECIFIED)
    return AddGrid(go.npoints_, go.spacing_, go.coord_);
  else
    return AddGrid(go.npoints_, go.spacing_, go.mode_);
}

/** Add a grid to the finite difference method object with explicit centering. */
int MeadGrid::AddGrid(int ngrd, float spc, Vec3 const& cntr)
{
  if (fdm_ == 0) {
    fdm_ = new FinDiffMethod();
  }

  try { 
    fdm_->add_level( ngrd, spc, Coord(cntr[0], cntr[1], cntr[2]) );
  }
  catch (MEADexcept& e) {
    return ERR("AddGrid(coord)", e);
  }
  levelOpts_.push_back(GridOpt(ngrd, spc, cntr));
  return 0;
}

/** Add a grid to the finite difference method object with centering string. */
int MeadGrid::AddGrid(int ngrd, float spc, Center_Mode ctrmode)
{
  if (fdm_ == 0) {
    fdm_ = new FinDiffMethod();
  }

  CenteringStyle censtl = ON_ORIGIN;
  switch (ctrmode) {
    case C_ON_ORIGIN       : censtl = ON_ORIGIN; break;
    case C_ON_CENT_OF_INTR : censtl = ON_CENT_OF_INTR; break;
    case C_ON_GEOM_CENT    : censtl = ON_GEOM_CENT; break;
    case C_SPECIFIED:
      mprinterr("Internal Error: MeadGrid::AddGrid called with C_SPECIFIED and no coords.\n");
      return 1;
  }
  try {
    fdm_->add_level( ngrd, spc, censtl );
  }
  catch (MEADexcept& e) {
    return ERR("AddGrid(style)", e);
  }
  levelOpts_.push_back(GridOpt(ngrd, spc, ctrmode));
  return 0;
}

/** Print to stdout. */
void MeadGrid::Print() const {
  mprintf("Grid Options:\n");
  if (fdm_ != 0) {
    std::cout << *fdm_;
  }
}

/** Resolve on given coords. */
int MeadGrid::Resolve(Coord const& geom_center, Coord const& site_of_interest)
const
{
  //try {
    fdm_->resolve( geom_center, site_of_interest );
  //}
  //catch (MEADexcept& e) {
  //  return ERR("Resolve()");
  //}
  return 0;
}

/** Set up MEAD grid based on furthest atom-atom distance. */
int MeadGrid::SetupGridFromCoords(Frame const& frameIn, int maxGridIn) {
  int MAXGRID;
  if (maxGridIn < 1)
    MAXGRID = 41;
  else
    MAXGRID = maxGridIn;

  double maxdist2 = 0;
  int maxat1 = -1;
  int maxat2 = -1;
  for (int at1 = 0; at1 < frameIn.Natom(); at1++) {
    const double* xyz1 = frameIn.XYZ(at1);
    for (int at2 = at1 + 1; at2 < frameIn.Natom(); at2++) {
      const double* xyz2 = frameIn.XYZ(at2);
      double dist2 = DIST2_NoImage(xyz1, xyz2);
      if (dist2 > maxdist2) {
        maxdist2 = dist2;
        maxat1 = at1;
        maxat2 = at2;
      }
    }
  }
  double max = sqrt(maxdist2);
  mprintf("\tMax distance is %g Ang. between atoms %i and %i.\n", max, maxat1+1, maxat2+1);
  max = int(max + 5);
  int ival = (int)max;
  ival = ival % 2;
  if (ival == 0)
    max = max + 1;

  mprintf("\tMAX= %g\n", max);

  if (max <= MAXGRID) {
    AddGrid(max, 8.0, MeadGrid::C_ON_GEOM_CENT);
    AddGrid(max, 2.0, MeadGrid::C_ON_CENT_OF_INTR);
    AddGrid(max, 0.5, MeadGrid::C_ON_CENT_OF_INTR);
  } else if (max > MAXGRID && max <= MAXGRID*4 ) {
    // For larger compounds use  grid with MAXGRID  grid points.
    AddGrid(MAXGRID, 8.0, MeadGrid::C_ON_GEOM_CENT);
    AddGrid(MAXGRID, 2.0, MeadGrid::C_ON_CENT_OF_INTR);
    AddGrid(MAXGRID, 0.5, MeadGrid::C_ON_CENT_OF_INTR);
  } else if (max > MAXGRID*4 && max <= MAXGRID*16 ) {
    // For even larger ones, add one more level of focusing.
    AddGrid(MAXGRID, 32.0, MeadGrid::C_ON_GEOM_CENT);
    AddGrid(MAXGRID,  8.0, MeadGrid::C_ON_CENT_OF_INTR);
    AddGrid(MAXGRID,  2.0, MeadGrid::C_ON_CENT_OF_INTR);
    AddGrid(MAXGRID,  0.5, MeadGrid::C_ON_CENT_OF_INTR);
  } else if (max > MAXGRID*16 && max <= MAXGRID*64 ) {
    // For even larger ones, add yet one more level of focusing.
    AddGrid(MAXGRID, 128.0, MeadGrid::C_ON_GEOM_CENT);
    AddGrid(MAXGRID,  32.0, MeadGrid::C_ON_CENT_OF_INTR);
    AddGrid(MAXGRID,   8.0, MeadGrid::C_ON_CENT_OF_INTR);
    AddGrid(MAXGRID,   2.0, MeadGrid::C_ON_CENT_OF_INTR);
    AddGrid(MAXGRID,   0.5, MeadGrid::C_ON_CENT_OF_INTR);
  } else {
    // For gargantuan ones, re-scale grids, but it is pretty bad.
    mprintf("\tWarning: Grid is extremely large, calculations may not be as accurate.\n");
    double fac = (double)(max) / (double)(MAXGRID);
    AddGrid(MAXGRID, 2.0       * fac, MeadGrid::C_ON_GEOM_CENT);
    AddGrid(MAXGRID, 0.5       * fac, MeadGrid::C_ON_CENT_OF_INTR);
    AddGrid(MAXGRID, 0.125     * fac, MeadGrid::C_ON_CENT_OF_INTR);
    AddGrid(MAXGRID, 0.03125   * fac, MeadGrid::C_ON_CENT_OF_INTR);
    AddGrid(MAXGRID, 0.0078125 * fac, MeadGrid::C_ON_CENT_OF_INTR);
  }
  return 0;
}

/** Set up MEAD grid based on furthest atom-atom distance. */
int MeadGrid::SetupGridFromCoords(Frame const& frameIn) {
  return SetupGridFromCoords(frameIn, -1);
}
