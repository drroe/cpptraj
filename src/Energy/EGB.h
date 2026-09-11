#ifndef INC_ENERGY_EGB_H
#define INC_ENERGY_EGB_H
#include <vector>
namespace Cpptraj {
namespace Energy {
/// Implement generalized Born implicit solvent models
class EGB {
  public:
    typedef std::vector<double> Darray;

    enum GBtype { HCT = 0, ///< igb = 1
                  OBC2,    ///< igb = 2
                  OBC5,    ///< igb = 5
                  NECK,    ///< igb = 7
                  NECK2,   ///< igb = 6
                  GB66,    ///< igb = 66
                  NGBTYPES };
    /// CONSTRUCTOR
    EGB();
    /// DESTRUCTOR - virtual since inherited
    virtual ~EGB() {}
  private:
    GBtype gbtype_;
};
}
}
#endif
