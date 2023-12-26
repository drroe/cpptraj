// Unit test for NameType class
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>
#include "TypeNameHolder.h"

static const int Err(const char* msg) {
  fprintf(stderr, "Error: %s\n", msg);
  return 1;
}

/*struct reverseSort {
  bool operator() (const NameType& lhs, const NameType& rhs) const {
    return ( lhs > rhs );
  }
} reverseSortObj;*/

int main() {
  TypeNameHolder type1;
  type1.AddName("CT");
  TypeNameHolder type2("O");
  bool ismatch = type1.Match_NoWC( type2 );
  if (ismatch) return Err("TypeNameHolder single type match failed.");
  return 0;
}
