#include "Counter_Regular.h"
#include "StringRoutines.h"

using namespace Cpptraj;
/** CONSTRUCTOR */
Counter_Regular::Counter_Regular() :
  start_(0),
  stop_(-1),
  offset_(1),
  current_(0),
  total_read_frames_(0)
{}

/** Determine total number of elements represented based on start, stop, and offset. */
int Counter_Regular::determineTotal() const {
  int nTotal = -1;
  if (stop_ != -1) {
    int Nframes = stop_ - start_;
    nTotal = Nframes / offset_;
    // Round up
    if ( (Nframes % offset_) > 0 )
      ++nTotal;
  }
  return nTotal;
}

/** CONSTRUCTOR - start stop offset */
Counter_Regular::Counter_Regular(int start, int stop, int offset) :
  start_(start),
  stop_(stop),
  offset_(offset),
  current_(0)
{
  total_read_frames_ = determineTotal();
}

/** \return String with counter info. */
std::string Counter_Regular::CounterInfo() const {
  std::string sstop;
  if (stop_ != -1)
    sstop = integerToString(stop_);
  else
    sstop.assign("EOF");
  return integerToString(start_) + "-" + sstop + ", " + integerToString(offset_);
}
