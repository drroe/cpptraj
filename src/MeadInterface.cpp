#include "MeadInterface.h"
// MEAD includes
#include "../mead/FinDiffMethod.h"

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
