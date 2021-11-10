#include "Counter_Array.h"
#include "StringRoutines.h"

using namespace Cpptraj;
/** CONSTRUCTOR */
Counter_Array::Counter_Array() {}

/** CONSTRUCTOR. Take array. */
Counter_Array::Counter_Array(Iarray const& arr) : numbers_(arr) {}

/** \return string containing counter info. */
std::string Counter_Array::CounterInfo() const {
  // FIXME fill this in
  return std::string("Array");
}

/** \return string with more verbose info. */
std::string Counter_Array::Verbose(int total_frames) const{
  std::string out;
  out.assign("reading " + integerToString(numbers_.size()) + " of " +
             integerToString(total_frames));
  return out;
}
