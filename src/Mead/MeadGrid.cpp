#include "MeadGrid.h"
#include "MeadInterface.h"
#include "../CpptrajStdio.h"
// MEAD includes
#include "../../mead/FinDiffMethod.h"
#include "../../mead/MEADexcept.h"
// For debug print
#include <iostream>

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

/** Corresponds to enum Center_Mode */
const char* MeadGrid::Center_ModeStr_[] = {
  "ON_ORIGIN", "ON_CENT_OF_INTR", "ON_GEOM_CENT", "SPECIFIED"
};

/** \return Char string corresponding to Center_Mode */
const char* MeadGrid::Center_ModeStr(Center_Mode gc) {
  return Center_ModeStr_[gc];
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

  CenteringStyle censtl;
  switch (ctrmode) {
    case C_ON_ORIGIN       : censtl = ON_ORIGIN; break;
    case C_ON_CENT_OF_INTR : censtl = ON_CENT_OF_INTR; break;
    case C_ON_GEOM_CENT    : censtl = ON_GEOM_CENT; break;
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
