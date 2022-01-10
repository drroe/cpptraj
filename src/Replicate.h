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

    Unit& ModifyRepUnit()         { return replicate_; }
    Unit const& RepUnit()   const { return replicate_; }
    std::string const& ID() const { return id_;   }
  private:
    Unit replicate_; ///< All atoms in the unit belong to the same replicate
    std::string id_; ///< Unique identifier for atoms in the replicate
};
//}
#endif
