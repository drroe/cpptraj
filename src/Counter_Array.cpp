#include "Counter_Array.h"

using namespace Cpptraj;
/** CONSTRUCTOR */
Counter_Array::Counter_Array() {}

/** CONSTRUCTOR. Take array. */
Counter_Array::Counter_Array(Iarray const& arr) : numbers_(arr) {}

/** \return string containin counter info. */
std::string Counter_Array::CounterInfo() const {
  // FIXME fill this in
  return std::string("Array");
}
