#ifndef INC_HB_HBPARALLEL_H
#define INC_HB_HBPARALLEL_H
#ifdef MPI
#include "../Parallel.h"
#endif
namespace Cpptraj {
namespace HB {
/// Hold routines for syncing up HB data in parallel (MPI)
class HbParallel {
  public:
    HbParallel();
  private:
#   ifdef MPI
    Parallel::Comm trajComm_; ///< Across-trajectory communicator
#   endif
};
}
}
#endif
