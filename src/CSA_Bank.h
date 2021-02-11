#ifndef INC_CSA_BANK_H
#define INC_CSA_BANK_H
#include <vector>
namespace Cpptraj {
namespace CSA {

// Forward declare
class Trial;

/** A CSA bank holds 1 or more Trials to be optimized. */
class Bank {
  public:
    Bank() {}
  private:
    typedef std::vector<Trial*> Btype;

    Btype trials_;
};

}
}
#endif
