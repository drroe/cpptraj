// Unit test for Range class
#include <cstdio>
#include "Range.h"

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

int main() {
  Range range1("8,1-4,15-20,12");
  Range::const_iterator it = range1.begin();
  if (*(it++) != 1) return Err("Expected 1");
  if (*(it++) != 2) return Err("Expected 2");
  if (*(it++) != 3) return Err("Expected 3");
  if (*(it++) != 4) return Err("Expected 4");
  if (*(it++) != 8) return Err("Expected 8");
  if (*(it++) != 12) return Err("Expected 12");
  if (*(it++) != 15) return Err("Expected 15");
  if (*(it++) != 16) return Err("Expected 16");
  if (*(it++) != 17) return Err("Expected 17");
  if (*(it++) != 18) return Err("Expected 18");
  if (*(it++) != 19) return Err("Expected 19");
  if (*(it++) != 20) return Err("Expected 20");

  return 0;
}
