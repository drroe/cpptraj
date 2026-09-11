#ifndef INC_ENERGY_EGB_H
#define INC_ENERGY_EGB_H
namespace Cpptraj {
namespace Energy {
/// Implement generalized Born implicit solvent models
class EGB {
  public:
    /// CONSTRUCTOR
    EGB();
    /// DESTRUCTOR - virtual since inherited
    virtual ~EGB() {}
  private:
};
}
}
#endif
