#ifndef INC_REPLICATE_H
#define INC_REPLICATE_H
#include <string>
#include "Unit.h"
//namespace Cpptraj {
/// Define atoms that belong to the same "replicate"
class Replicate {
  public:
    Replicate();
    /// CONSTRUCTOR - replicate is atom0 up to atom1
    Replicate(int, int);
  private:
    Unit replicate_; ///< All atoms in the unit belong to the same replicate
    std::string id_; ///< Unique identifier for atoms in the replicate
};
//}
#endif
