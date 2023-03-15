#include "MeadInterface.h"
#include "Vec3.h"
#include "CpptrajStdio.h"
// MEAD includes
#include "../mead/FinDiffMethod.h"
#include "../mead/MEADexcept.h"

using namespace Cpptraj;

/** CONSTRUCTOR */
MeadInterface::MeadInterface() :
  fdm_(0)
{
  fdm_ = new FinDiffMethod();
}

/** DESTRUCTOR */
MeadInterface::~MeadInterface() {
  if (fdm_ != 0) delete fdm_;
}

/** Add a grid to the finite difference method object. */
int MeadInterface::AddGrid(int ngrd, float spc, Vec3 const& cntr)
{
  try { 
    fdm_->add_level( ngrd, spc, Coord(cntr[0], cntr[1], cntr[2]) );
  }
  catch (MEADexcept& e) {
    mprinterr("Error: MEAD error in AddGrid(): '%s' '%s' '%s'\n",
              e.get_error1().c_str(),
              e.get_error2().c_str(),
              e.get_error3().c_str());
    return 1;
  }
  return 0;
}
