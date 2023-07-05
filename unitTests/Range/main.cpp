// Unit test for Range class
#include <cstdio>
#include "Range.h"

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

int BasicRangeTest() {
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

  Range range2;
  range2.SetRange("22,10-12,8,27", Range::UNSORTED);
  it = range2.begin();
  if (*(it++) != 22) return Err("Expected 22");
  if (*(it++) != 10) return Err("Expected 10");
  if (*(it++) != 11) return Err("Expected 11");
  if (*(it++) != 12) return Err("Expected 12");
  if (*(it++) !=  8) return Err("Expected 8");
  if (*(it++) != 27) return Err("Expected 27");
  return 0;
}

int TestRange(Range const& myRange) {
  Range::const_iterator it = myRange.begin();
  for (int i = 3; i < 7; i++, ++it) {
    if (*it != i) return Err("Iterator range 3-6 does not match.");
  }
  if (*it != 9) return Err("Iterator 9 does not match.");
  ++it;
  for (int i = 12; i < 15; i++, ++it) {
    if (*it != i) return Err("Iterator range 12-14 does not match.");
  }
  return 0;
}

int main() {
  if (BasicRangeTest()) return 1;

  Range myRange;
  if (myRange.SetRange("3-6,9,12-14")) return Err("SetRange failed.");
  myRange.PrintToStdout();
  if (myRange.Size() != 8) return Err("myRange size is not 8.");
  printf("myRange (+ 1) %s\n", myRange.PrintRange(1).c_str());
  printf("Range arg: '%s'\n", myRange.RangeArg());
  if (TestRange(myRange)) return 1;

  // Test range with duplicates
  myRange.Clear();
  if (myRange.SetRange("5,12,9,3-6,9,12-14")) return Err("SetRange duplicate failed.");
  myRange.PrintToStdout();
  if (myRange.Size() != 8) return Err("myRange duplicate size is not 8.");
  if (TestRange(myRange)) return 1;

  // Test copying and assignment
  Range rangeCopy( myRange );
  if (TestRange(rangeCopy)) return 1;
  Range rangeAssign;
  rangeAssign = myRange;
  if (TestRange(rangeAssign)) return 1;

  // Test basic range from numbers
  Range range2;
  if (range2.SetRange(10, 15)) return Err("SetRange range2 failed.");
  printf("range2 (should be 10-14 minus 1) %s\n", range2.PrintRange(-1).c_str());
  printf("Range2 arg: '%s'\n", range2.RangeArg());
  range2.ShiftBy(-1);
  if (range2.Size() != 5) return Err("range2 size is not 5.");
  Range::const_iterator it = range2.begin();
  for (int i = 9; i < 14; i++, ++it) {
    if (*it != i) return Err("Iterator range 9-14 does not match.");
  }

  // Test AddToRange
  Range range3;
  range3.SetRange("6,9,12");
  range3.AddToRange(4);
  range3.AddToRange(3);
  range3.AddToRange(14);
  range3.AddToRange(5);
  range3.AddToRange(9);
  range3.AddToRange(13);
  printf("range3 %s\n", range3.PrintRange(0).c_str());
  if (TestRange(myRange)) return Err("AddToRange test failed.");

  return 0;
}
