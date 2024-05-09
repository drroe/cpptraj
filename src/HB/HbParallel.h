#ifndef INC_HB_HBPARALLEL_H
#define INC_HB_HBPARALLEL_H
#ifdef MPI
#include "../Parallel.h"
#include <vector>
#endif
namespace Cpptraj {
namespace HB {
class Hbond;
/// Hold routines for syncing up HB data in parallel (MPI)
class HbParallel {
  public:
    HbParallel();
  private:
#   ifdef MPI
    /// Determine number of hydrogen bonds on each rank
    static std::vector<int> GetRankNhbonds( int, Parallel::Comm const& );
    /// Flatten given Hbond into arrays
    static void HbondToArray(std::vector<double>&, std::vector<int>&, Hbond const&);

    Parallel::Comm trajComm_; ///< Across-trajectory communicator
#   endif
};
}
}
#endif
