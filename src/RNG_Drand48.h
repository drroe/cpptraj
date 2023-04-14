#ifndef INC_RNG_DRAND48_H
#define INC_RNG_DRAND48_H
#include "RNG.h"
namespace Cpptraj {
/// Drand48 RNG from stdlib
class RNG_Drand48 : public RNG {
  public:
    RNG_Drand48();

    unsigned int Number();
    double Generate();
  private:
    int SetupRng();
};
}
#endif
