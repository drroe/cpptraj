#ifndef INC_ENERGY_EGB_OBC_H
#define INC_ENERGY_EGB_OBC_H
#include "EGB.h"
namespace Cpptraj {
namespace Energy {
/// Implement Onufriev-Bashfor-Case variants of GB
class EGB_OBC : public EGB {
  public:
    /// CONSTRUCTOR
    EGB_OBC();

    int InitABG(GBtype, int);
  private:
    Darray gbalpha_;
    Darray gbbeta_;
    Darray gbgamma_;
};
}
}
#endif
