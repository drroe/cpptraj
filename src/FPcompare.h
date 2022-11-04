#ifndef INC_FPCOMPARE_H
#define INC_FPCOMPARE_H
#include "Constants.h"
namespace Cpptraj {
/// Floating point comparison routines
namespace FPcompare {
/// \return true if 2 floating point numbers are equal within a tolerance
template <typename T> bool FEQ(T v1, T v2) {
  T delta = v1 - v2;
  if (delta < 0.0) delta = -delta;
  return (delta < Constants::SMALL);
}
/// \return true if 2 floating point numbers are not equal outside a tolerance
template <typename T> bool FNE(T v1, T v2) {
  T delta = v1 - v2;
  if (delta < 0.0) delta = -delta;
  return (delta > Constants::SMALL);
}

}
}
#endif
