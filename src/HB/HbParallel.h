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
    typedef std::vector<int> Iarray;
  public:
    HbParallel();
    /// Sync all data to the master process
    int SyncToMaster(int&, Iarray const&) const;
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
