#include <cstdio>
#include <cmath>
#include <vector>
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
    printf("\t[%6i] %8i %8i\n", regular.CurrentIdx(), frame, regular.CurrentNumber());
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

int array(int start, int stop, int offset) {
  std::vector<int> numbers;
  int frame = start;
  while (frame < stop) {
    numbers.push_back(frame);
    frame += offset;
  }

  Counter_Array regular(numbers);
  regular.StartCounter();
  frame = start;
  printf("Array Counter: (%i-%i, %i)\n", start, stop, offset);
  while (!regular.IsFinished()) {
    printf("\t[%6i] %8i %8i\n", regular.CurrentIdx(), frame, regular.CurrentNumber());
    if (frame != regular.CurrentNumber()) {
      fprintf(stderr,"Error: %i != %i\n", frame, regular.CurrentNumber());
      return Err("Array counter value does not match.");
    }
    regular.UpdateCounter();
    frame += offset;
    if (frame > stop*2) {
      return Err("Array counter IsFinished() failed.");
    }
  }
  return 0;
}

int main() {
  if (regular(0,10,1)) return 1;
  if (regular(5,23,3)) return 1;
  if (array(0,10,1)) return 1;
  if (array(5,23,3)) return 1;

  // Test NumberAtIdx()
  Counter_Regular c0(5,23,3);
  if (c0.NumberAtIdx(2) != 11)
    return Err("NumberAtIdx() failed for Counter_Regular");
  std::vector<int> numbers;
  numbers.push_back(5);
  numbers.push_back(2);
  numbers.push_back(8);
  numbers.push_back(7);
  Counter_Array c1(numbers);
  if (c1.NumberAtIdx(2) != 8)
    return Err("NumberAtIdx() failed for Counter_Array");

  return 0;
}
