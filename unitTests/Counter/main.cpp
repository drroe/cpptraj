#include <cstdio>
#include <cmath>
#include "Counter_Regular.h"
#include "Counter_Array.h"

using namespace Cpptraj;

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

//const double SMALL        = 0.00000000000001;

//bool f_not_equals(double d1, double d2) {
//  return (fabs(d1 - d2) > SMALL);
//}

int regular(int start, int stop, int offset) {
  Counter_Regular regular(start, stop, offset);
  regular.StartCounter();
  int frame = start;
  printf("Regular Counter: (%i-%i, %i)\n", start, stop, offset);
  while (!regular.IsFinished()) {
    printf("\t%8i %8i\n", frame, regular.CurrentNumber());
    if (frame != regular.CurrentNumber()) {
      fprintf(stderr,"Error: %i != %i\n", frame, regular.CurrentNumber());
      return Err("Regular counter value does not match.");
    }
    regular.UpdateCounter();
    frame += offset;
    if (frame > stop*2) {
      return Err("Regular counter IsFinished() failed.");
    }
  }
  return 0;
}

int main() {
  if (regular(0,10,1)) return 1;
  if (regular(5,23,3)) return 1;
  return 0;
}
