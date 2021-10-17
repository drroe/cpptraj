#include "Counter_Regular.h"

using namespace Cpptraj;
/** CONSTRUCTOR */
Counter_Regular::Counter_Regular() :
  start_(0),
  stop_(-1),
  offset_(1),
  current_(0)
{}

/** CONSTRUCTOR - start stop offset */
Counter_Regular::Counter_Regular(int start, int stop, int offset) :
  start_(start),
  stop_(stop),
  offset_(offset),
  current_(0)
{}
